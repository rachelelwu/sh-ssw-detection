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
    # Window definition (in DAYS; code infers samples from the time step)
    wdw: int = 10,
    # Choose ONE of the two thresholding modes:
    drop_per_day: float | None = None,   # fixed threshold (m/s per day)
    percentile: float | None = None,     # e.g., 60 -> keep strongest 40% (by magnitude)
    # Algorithm controls
    mode: str = "decel",                 # "decel" or "accel"
    min_gap_days: int = 20,
    ratio_thres: float = 1.2,            # (max-min)/|diff| must be <= ratio_thres
    time_dim: str = "time",
    # Metadata
    method: str = "u_wind_decel",
    definition: str | None = None,
    data_source: str | None = None,
    level_hpa: float | None = None,
    latitude: float | None = None,
    lat_band: str | None = None,
    notes: str | None = None,
) -> tuple[pd.DataFrame, pd.DatetimeIndex, np.ndarray]:
    """
    Detect acceleration/deceleration events with either a fixed per-day threshold
    or a data-driven percentile threshold on event magnitudes.

    Window diff: diff(t) = u[t+window_end] - u[t] over a forward window of length `wdw` days.
    For decel, negative diffs are candidates; for accel, positive diffs.
    """
    # --- validation
    if not isinstance(u, xr.DataArray):
        raise TypeError("detect_u_decel_events expects an xarray.DataArray.")
    if time_dim not in u.dims:
        raise ValueError(f"Expected `{time_dim}` in input DataArray.")
    if mode not in {"decel", "accel"}:
        raise ValueError("mode must be 'decel' or 'accel'.")
    if (drop_per_day is None) and (percentile is None):
        # sensible defaults if nothing specified
        if level_hpa in (10, 10.0):
            drop_per_day = -2.0 if mode == "decel" else 2.0
        elif level_hpa in (50, 50.0):
            drop_per_day = -1.0 if mode == "decel" else 1.0
        else:
            drop_per_day = -2.0 if mode == "decel" else 2.0

    # --- 1) 1-D series + infer sampling; map wdw days -> n samples
    t = pd.to_datetime(u[time_dim].to_index())
    s = pd.Series(u.values, index=t)

    if len(s) < 2:
        # too short for any window
        empty = build_events_df(
            dates=[],
            method=method,
            definition=definition or f"{mode}_wdw{wdw}",
            data_source=data_source or "",
            level_hpa=level_hpa, latitude=latitude, lat_band=lat_band, notes=notes,
        )
        return empty, pd.DatetimeIndex([]), np.array([])

    dt = pd.Series(t).diff().dropna().median()
    dt_days = float(dt / pd.Timedelta(days=1))
    n = max(1, int(round(wdw / dt_days)))        # samples in window
    window_days = n * dt_days                     # actual days spanned by window

    # --- 2) forward-window diff and range aligned to window START
    s_end = s.shift(-(n - 1))                    # value at t + (n-1 steps)
    diff = s_end - s

    sr = s[::-1]
    rmax = sr.rolling(window=n, min_periods=n).max()[::-1]
    rmin = sr.rolling(window=n, min_periods=n).min()[::-1]
    window_range = (rmax - rmin)

    valid = ~diff.isna()
    diff = diff[valid]
    window_range = window_range[valid]

    # --- 3) candidates by sign + coherence ratio (no magnitude threshold yet)
    if mode == "decel":
        sign_mask = (diff <= 0)
        better = lambda a, b: a < b  # more negative is stronger
    else:
        sign_mask = (diff >= 0)
        better = lambda a, b: a > b  # more positive is stronger

    ratio = (window_range / diff.abs()).replace([np.inf, -np.inf], np.nan)
    ratio_mask = ratio <= ratio_thres
    cand_idx = diff.index[sign_mask & ratio_mask & diff.notna() & ratio.notna()]

    if len(cand_idx) == 0:
        empty = build_events_df(
            dates=[],
            method=method,
            definition=definition or f"{mode}_wdw{wdw}_ratio{ratio_thres:g}",
            data_source=data_source or "",
            level_hpa=level_hpa, latitude=latitude, lat_band=lat_band,
            notes=(notes or "") + f"; window={wdw}d",
        )
        return empty, pd.DatetimeIndex([]), np.array([])

    # --- 4) min-gap with "keep stronger within the gap" selection
    cand = pd.DataFrame({"onset": cand_idx, "diff": diff.loc[cand_idx].values}) \
            .sort_values("onset")

    selected_onsets: list[pd.Timestamp] = []
    selected_diffs: list[float] = []

    for onset, dval in zip(cand["onset"].values, cand["diff"].values):
        if not selected_onsets:
            selected_onsets.append(onset); selected_diffs.append(dval); continue
        gap_days = int((pd.Timestamp(onset) - pd.Timestamp(selected_onsets[-1])) / pd.Timedelta(days=1))
        if gap_days >= min_gap_days:
            selected_onsets.append(onset); selected_diffs.append(dval)
        else:
            if better(dval, selected_diffs[-1]):
                selected_onsets[-1] = onset
                selected_diffs[-1] = dval

    event_dates = pd.DatetimeIndex(selected_onsets)
    drop_vals = np.asarray(selected_diffs)  # signed: negative (decel) / positive (accel)

    # --- 5) Percentile mode: compute magnitude cutoff and filter
    used_threshold_per_day: float | None = None
    used_threshold_total: float | None = None
    used_mode = "fixed"
    if percentile is not None:
        if len(drop_vals) == 0:
            # nothing to percentile over
            pass
        else:
            mag = np.abs(drop_vals)                 # magnitude in m/s over the window
            cutoff_mag = np.percentile(mag, percentile)
            keep = mag >= cutoff_mag                # keep strongest (100 - p)% events
            event_dates = event_dates[keep]
            drop_vals = drop_vals[keep]

            used_threshold_total = float(cutoff_mag)                 # m/s over window
            # convert to per-day, set signed according to mode
            signed = -1.0 if mode == "decel" else 1.0
            used_threshold_per_day = signed * (cutoff_mag / window_days)
            used_mode = f"percentile_{percentile:g}"

    # --- 6) Fixed-threshold mode (if requested in addition or instead)
    # If percentile is None, apply fixed threshold here:
    if percentile is None and drop_per_day is not None:
        total = drop_per_day * window_days          # signed total over window
        if mode == "decel":
            keep = drop_vals <= total               # total is negative
        else:
            keep = drop_vals >= total
        event_dates = event_dates[keep]
        drop_vals = drop_vals[keep]
        used_threshold_per_day = float(drop_per_day)
        used_threshold_total = float(abs(total))    # report magnitude for clarity
        used_mode = "fixed"

    # --- 7) Build unified table
    definition_auto = (
        f"{mode}_{used_mode}_wdw{wdw}_ratio{ratio_thres:g}"
        + ("" if used_threshold_per_day is None else f"_thr{used_threshold_per_day:.3g}/day")
    )

    events_df = build_events_df(
        dates=event_dates,
        method=method,
        definition=definition or definition_auto,
        data_source=data_source or "",
        threshold=used_threshold_per_day,          # per-day signed threshold used
        level_hpa=level_hpa,
        latitude=latitude,
        lat_band=lat_band,
        notes=((notes or "") + f"; min_gap={min_gap_days}d; window≈{window_days:.2f}d"),
        extra_cols={
            "drop_val": drop_vals,                 # signed total change over the window (m/s)
            "window_days": np.repeat(window_days, len(event_dates)),
            # "threshold_total_mag": (
            #     [] if used_threshold_total is None else np.repeat(used_threshold_total, len(event_dates))
            #),
        },
    )

    return events_df, event_dates#, drop_vals
