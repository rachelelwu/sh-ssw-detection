import pandas as pd

import numpy as np
import os

def select_pc1_events(lv=50, basepath=".", thresh=1.0, persistence=None, minsep=60, output_csv=True):
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
            fname = os.path.join(basepath, f"daily.pc1.z{lv}.{yr}{mon:02d}.txt")
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

    event_dates = np.array(events, dtype="datetime64[D]")

    # --- optional CSV output ---
    if output_csv:
        outfile = f"EOF1.Z{lv}.gt{thresh:.1f}.csv"
        with open(outfile, "w") as fout:
            fout.write("# Daily data projected onto EOF1 of GPHA calculated for each month\n")
            fout.write(f"# (June–December) at {lv} hPa over the domain of 20–90S\n")
            fout.write(f"# standardized daily PC1 > {thresh}, which persists for {persistence+1} days\n")
            fout.write(f"# two events should be separated by {minsep} days\n")
            pd.Series(event_dates).to_csv(fout, index=False, header=["Date"])
        print(f"Saved {len(event_dates)} events to {outfile}")

    return event_dates
