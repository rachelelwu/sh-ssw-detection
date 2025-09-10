# src/sh_ssw_methods/u_wind_decel.py
from __future__ import annotations
import numpy as np
import pandas as pd
import xarray as xr

from .utils.event_iden_funcs import (
    build_events_df,
)

def detect_u_decel_events(
    u: xr.DataArray,
    *,
    # Algorithm controls
    wdw: int = 10,                 # event window length (days)
    drop_per_day: float | None = None,
    mode: str = "decel",           # "decel" or "accel"
    min_gap_days: int = 20,        # next event must be at least this far after the previous
    ratio_thres: float = 1.2,      # max(|window range| / |window end-start|)
    time_dim: str = "time",
    # Metadata (for build_events_df)
    method: str = "u_wind_decel",
    data_source: str | None = None,
    level_hpa: float | None = None,
    latitude: float | None = None,
    lat_band: str | None = None,
) -> tuple[pd.DataFrame, pd.DatetimeIndex, np.ndarray]:
    """
    Identify deceleration/acceleration events in a zonal-mean wind time series.

    Reimplements the logic from Wu et al. (2022)-style detection:
      - event magnitude = u[t+wdw-1] - u[t] (over a wdw-day window)
      - for decel, require diff <= drop_per_day * wdw
        for accel, require diff >= drop_per_day * wdw
      - ratio = (max(window)-min(window)) / |diff|  <= ratio_thres
      - enforce separation by min_gap_days AND, within that gap window,
        keep only the stronger event (more negative for decel, more positive for accel)

    Returns
    -------
    events_df : pd.DataFrame
        tidy table via build_events_df (with only `date` = onset)
    event_dates : pd.DatetimeIndex
        onset dates
    drop_vals : np.ndarray
        event magnitudes (end - start over the window)
    """

    # make sure input is xr.DataArray
    if not isinstance(u, xr.DataArray):
        raise TypeError(
            "detect_u_decel_events expects an xarray.DataArray. "
            "If you have a Dataset, select a variable first, e.g. ds['u']."
        )
    
    if time_dim not in u.dims:
        raise ValueError(f"Expected `{time_dim}` dimension in input DataArray.")
    if mode not in {"decel", "accel"}:
        raise ValueError("mode must be 'decel' or 'accel'.")

    # Choose a sensible default threshold if not provided (matches your original)
    if drop_per_day is None:
        # If level metadata is present, use it; else default to -2 m/s/day (10 hPa-like)
        if level_hpa in (10, 10.0):
            drop_per_day = -2.0 if mode == "decel" else 2.0
        elif level_hpa in (50, 50.0):
            drop_per_day = -1.0 if mode == "decel" else 1.0
        else:
            drop_per_day = -2.0 if mode == "decel" else 2.0

    # 1) Prepare a clean 1-D time series
    t = pd.to_datetime(u[time_dim].to_index())
    s = pd.Series(u.values, index=t)

    # 2) Compute window start→end diff and window range (vectorized)
    # Align "diff" with window START index
    s_end = s.shift(-(wdw - 1))         # value at t+wdw-1
    diff = s_end - s                    # window difference

    # Rolling range over forward-looking window:
    # We can compute a centered rolling on a reversed series to align at window start,
    # or use expanding on slices. Simpler: compute rolling max/min on forward windows via numpy stride.
    # For clarity (and still efficient), do this with pandas rolling on the *reversed* series, then flip back.
    sr = s[::-1]
    rmax = sr.rolling(window=wdw, min_periods=wdw).max()[::-1]
    rmin = sr.rolling(window=wdw, min_periods=wdw).min()[::-1]
    window_range = rmax - rmin

    # Keep only indices where a full window exists
    valid = ~diff.isna()
    diff = diff[valid]
    window_range = window_range[valid]

    # 3) Apply threshold and ratio filters
    drop_thres_total = drop_per_day * wdw
    if mode == "decel":
        cand_mask = (diff <= drop_thres_total)
        stronger_is = np.argmin   # more negative is stronger
        better = lambda a, b: a < b
    else:
        cand_mask = (diff >= drop_thres_total)
        stronger_is = np.argmax   # more positive is stronger
        better = lambda a, b: a > b

    ratio = (window_range / diff.abs()).replace([np.inf, -np.inf], np.nan)
    ratio_mask = ratio <= ratio_thres
    cand_idx = diff.index[cand_mask & ratio_mask & diff.notna() & ratio.notna()]

    if len(cand_idx) == 0:
        events_df = build_events_df(
            dates=[],
            method=method,
            definition=definition or f"{mode}_wdw{wdw}_ratio{ratio_thres:g}_thr{drop_per_day:g}perday",
            data_source=data_source or "",
            threshold=drop_per_day,
            level_hpa=level_hpa,
            latitude=latitude,
            lat_band=lat_band,
            notes=notes or f"min_gap={min_gap_days}d",
        )
        return events_df, pd.DatetimeIndex([]), np.array([])

    # 4) Enforce min-gap and “stronger-in-window” replacement
    cand = pd.DataFrame({
        "onset": cand_idx,
        "diff": diff.loc[cand_idx].values,            # event magnitude over window
    }).sort_values("onset")

    selected_onsets: list[pd.Timestamp] = []
    selected_diffs: list[float] = []

    # We walk candidates in time; if a candidate is within min_gap_days of the last kept event,
    # we keep only the stronger one (compare magnitudes according to mode).
    for onset, dval in zip(cand["onset"].values, cand["diff"].values):
        if not selected_onsets:
            selected_onsets.append(onset)
            selected_diffs.append(dval)
            continue

        gap_days = int((pd.Timestamp(onset) - pd.Timestamp(selected_onsets[-1])) / pd.Timedelta(days=1))

        if gap_days >= min_gap_days:
            # far enough — accept as new event
            selected_onsets.append(onset)
            selected_diffs.append(dval)
        else:
            # within gap — replace if stronger
            if better(dval, selected_diffs[-1]):
                selected_onsets[-1] = onset
                selected_diffs[-1] = dval
            # else keep existing

    event_dates = pd.DatetimeIndex(selected_onsets)
    drop_vals = np.asarray(selected_diffs)

    # 5) Build unified table
    events_df = build_events_df(
        dates=event_dates,
        method=method,
        definition= f"{mode}_wdw{wdw}_ratio{ratio_thres:g}_thr{drop_per_day:g}perday",
        data_source=data_source or "",
        threshold=f"{drop_per_day}ms-1/day",
        level_hpa=level_hpa,
        latitude=latitude,
        lat_band=lat_band,
        notes=f"min_gap={min_gap_days}d; window={wdw}d",
        extra_cols={"drop_val": drop_vals}, 
    )

    return events_df, event_dates
