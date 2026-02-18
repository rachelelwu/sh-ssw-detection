import pandas as pd

import numpy as np
import xarray as xr

import os
from pathlib import Path

from .utils.event_iden_funcs import (
    build_events_df
)


def detect_pc1_events(lv=10, basepath="sh_ssw_methods/processed", thresh=1.0, persistence=None,
    data_source: str | None = "era5", 
    latitude: str | None = "-20_to_-90",
    minsep=60):
    """
    Detect daily EOF1-based events from standardized PC1 time series.

    Parameters
    ----------
    lv : int
        Pressure level (default 50 hPa).
    basepath : str
        Directory containing daily.pc1.z{lv}.{YYYY}{MM}.txt files.
    thresh : float
        Threshold for event detection (default 1.0).
    persistence : int or None
        Number of consecutive days PC1 must stay above threshold.
        If None: defaults are 10 (for 50 hPa) or 14 (for 10/1 hPa).
    minsep : int
        Minimum days between independent events (default 60).
    outfile : str or None
        Path to save events. If None, auto-generate.

    Returns
    -------
    events : list of (year, month, day)
        Detected independent events.
    counts : dict
        Number of events per year.
    """

    # --- persistence defaults by level ---
    if persistence is None:
        if lv == 50:
            persistence = 10
        elif lv in (10, 1):
            persistence = 14
        else:
            raise ValueError(f"No default persistence set for level {lv}. Please provide persistence manually.")

    years = np.arange(1979, 2024)
    eday = [30,31,31,30,31,30,31]  # days in June–Dec
    cumdays = np.cumsum([0]+eday)
    tot_day = cumdays[-1]

    # --- load daily PC1 into array (year, day) ---
    data = np.full((len(years), tot_day), np.nan)
    for yi, yr in enumerate(years):
        for mi, nd in enumerate(eday):
            mon = mi+6
            basepath = Path(basepath)
            zeof_dir = basepath / "zeof" ####
            fname = os.path.join(zeof_dir, f"daily.pc1.z{lv}.{yr}{mon:02d}.txt")
            arr = np.loadtxt(fname, ndmin=2)
            vals = arr[0,1:nd+1] if arr.ndim==2 else arr[1:nd+1]
            data[yi, cumdays[mi]:cumdays[mi+1]] = vals

    events = []
    counts = {}

    # --- detection loop ---
    for yi, yr in enumerate(years):
        candidates = []
        for d in range(tot_day - persistence):
            # require day d and next `persistence` days > thresh
            if np.all(data[yi, d:d+persistence+1] > thresh):
                candidates.append(d)

        kept = []
        if candidates:
            lastd = candidates[0]
            kept.append(lastd)
            for dd in candidates[1:]:
                if dd - lastd >= minsep:
                    kept.append(dd)
                    lastd = dd

        # convert to (year, month, day)
        for dd in kept:
            mon_idx = np.searchsorted(cumdays, dd, side='right')-1
            mon = mon_idx+6
            day_in_mon = dd - cumdays[mon_idx] + 1
            events.append(pd.Timestamp(year=yr, month=mon, day=day_in_mon))

        counts[yr] = len(kept)

    event_dates = pd.to_datetime(events)

    # --- build events DataFrame ---
    min_gap_days = minsep
    
    # --- build standardized events DataFrame ---
    events_df = build_events_df(
        dates=event_dates,
        method="zeof1", ####
        definition=f"EOF1_PC1_gt{thresh:.1f}_{int(lv)}hPa_{latitude}",
        data_source=data_source or "",
        threshold=f">{thresh:.1f}σ",
        level_hpa=str(lv),
        latitude=str(latitude),
        notes=(
            f"PC1 > {thresh:.1f}σ for at least {persistence} consecutive days; "
            f"events separated by ≥{min_gap_days} days"
        ),
        # extra_cols={"n_events": len(event_dates)}  # optional if build_events_df supports it
    )

    return events_df, event_dates


    
