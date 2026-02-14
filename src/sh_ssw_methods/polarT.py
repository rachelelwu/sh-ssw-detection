import numpy as np
import pandas as pd
import glob

import xarray as xr

from pathlib import Path

from .utils import (
    build_events_df
)

from.preprocessing import (
prep_input_polarT
)


# iden algorithm
def polarT1_lim(years, T_augnov, fw_idx, delT, cdelT, stddv, polarT):
    
    xtime1 = T_augnov['date']
    
    count = np.full((len(years), 122), -999, dtype=int)
    sel_dates = []
    
    for i, yr in enumerate(years):
        if fw_idx[i] < 0:
            continue
    
        # Loop j = 0 .. eidx-20  (avoid last 20 days before FW)
        jmax = int(fw_idx[i]) - 20
        jmax = min(jmax, delT.shape[1] - 5)   # prevent partial windows
        if jmax < 0:
            continue
    
        for j in range(0, jmax + 1):
            # 1) sign reversal (positive) & 5-day persistence
            if delT[i, j].item() > 0:
                win = delT[i, j:j+5].values  # length 5
                if win.shape[0] != 5:
                    continue                  # safety, but jmax should prevent this
    
                if np.all(win >= 0):
                    # print('cond1: %s' % (win))
                    maxidx = int(np.argmax(win))   # index in the 5-day window
                    midx   = j + maxidx            # absolute time index within season
    
                    # 2) anomaly vs climatology stddev
                    anom = (delT[i, midx] - cdelT[midx]).item()
                    
                    if anom >= float(stddv[midx]):
                        # print(xtime1[i,j].values)
                        # print('cond2: anom=%s, std=%s' % (anom, stddv[midx].values))
                        # print('polarT=%s' % polarT[i, j].item())
                        # 3) polar cap temperature positive at the window start j
                        if polarT[i, j].item() > 0:
                            
                            # print('cond3: polarT=%s' % polarT[i, j].item())
                            count[i, j] = 1
                            sel_dates.append(xtime1[i,j].values)

    # make sure events are 60 days apart
    starts = np.r_[True, np.diff(sel_dates) >= np.timedelta64(60, 'D')]
    event_dates = np.array(sel_dates)[starts]
    
    return event_dates


def polarT2_lim(years, T_augnov, fw_idx, delT, cdelT, stddv, polarT, pers=5):
    
    xtime1 = T_augnov['date']
    
    count = np.full((len(years), 122), -999, dtype=int)
    sel_dates = []
    
    for i, yr in enumerate(years):
        if fw_idx[i] < 0:
            continue
    
        # Loop j = 0 .. eidx-20  (avoid last 20 days before FW)
        jmax = int(fw_idx[i]) - 20
        jmax = min(jmax, delT.shape[1] - 5)   # prevent partial windows
        if jmax < 0:
            continue
    
        for j in range(0, jmax + 1):
            # 1) sign reversal (positive) & 5-day persistence
            if delT[i, j].item() > 0:
                w_delT = delT[i, j:j+pers].values  # length 5
                if w_delT.shape[0] != 5:
                    continue                  # safety, but jmax should prevent this
    
                if np.all(w_delT >= 0):
                    print('cond1: %s' % (w_delT))

                    # 2) anomaly ≥ 1σ for all 5 days
                    w_c   = cdelT[j:j+pers].values
                    w_std = stddv[j:j+pers].values
                    
                    anom = (w_delT - w_c) / w_std
                    
                    if not np.all(anom >= 1):
                        continue

                    print(xtime1[i,j].values)
                    print('cond2: anom=%s, std=%s' % (anom, w_std))
                    print('polarT=%s' % polarT[i, j].item())
                    
                    # 3) polar cap temperature positive at the window start j
                    w_polar = polarT[i, j + pers].values
                    if not (np.all(np.isfinite(w_polar)) and np.all(w_polar >= 0)):
                        continue    
                    print('cond3: polarT=%s' % polarT[i, j].item())
                    count[i, j] = 1
                    sel_dates.append(xtime1[i,j].values)

    # make sure events are 60 days apart
    starts = np.r_[True, np.diff(sel_dates) >= np.timedelta64(60, 'D')]
    event_dates = np.array(sel_dates)[starts]
    
    return event_dates


def detect_ssw_polarT(
    Tband_file, 
    polarT_file, 
    fsw_file,
    istart_date,
    iend_date,
    *,
    algo: str = "T1",               # choose between "T1" or "T2"
    persist_days: int = 5,
    min_gap_days: int = 60,
    thres_sigma: float = 1.0,
    level_hpa: float | None = None,
    latitude: str | None = "-80_to_-90_-60_to_-70",
    data_source: str | None = None,
) -> tuple[pd.DataFrame, pd.DatetimeIndex]:
    """
    Wrapper for Shen et al. (2022) temperature-gradient–based SSW detection.

    Calls either `polarT1_lim` or `polarT2_lim` depending on `algo`.
    Returns standardized output (events_df, event_dates) compatible with `detect_ssw_u_anom`.

    Parameters
    ----------
     Tband_file : str or Path
        Path to the NetCDF file containing temperature-band data (e.g., T80–90S and T60–70S).Used by `output_Tvars()` to compute ΔT, climatological mean, and standard deviation.
    polarT_file : str or Path
        Path to the NetCDF file containing polar-cap (60–90S) mean temperature.
    fsw_file : str or Path
        Path to the NetCDF or .npy file containing final stratospheric  warming (FSW) dates. Used by `output_fwidx()` to produce per-year index positions for each event.
    algo : str, optional
        "T1" or "T2" algorithm choice.

    persist_days : int, optional
        Persistence requirement (used by T2).
    min_gap_days : int, optional
        Minimum separation between detected events.
    thres_sigma : float, optional
        Standard deviation threshold (used by T2).
    level_hpa : float, optional
        Pressure level label.
    latitude : str, optional
        Latitude range label.
    data_source : str, optional
        Metadata label for data source.

    Returns
    -------
    events_df : pd.DataFrame
        Standardized event metadata table.
    event_dates : pd.DatetimeIndex
        Detected event onset dates.
    """

    iyear = int(istart_date[:4])
    iiyear = int(iend_date[:4])
    years = np.arange(iyear, iiyear+1, 1)
    
    T_augnov, cdelT, delT, stddv, polarT, fw_idx = prep_input_polarT(Tband_file, polarT_file, fsw_file, istart_date, iend_date)
    
    # --- Dispatch to chosen algorithm ---
    if algo.upper() == "T1":
        event_dates = polarT1_lim(years, T_augnov, fw_idx, delT, cdelT, stddv, polarT)
        method_name = "polar_T1"
    elif algo.upper() == "T2":
        event_dates = polarT2_lim(years, T_augnov, fw_idx, delT, cdelT, stddv, polarT, pers=persist_days)
        method_name = "polar_T2"
    else:
        raise ValueError("Invalid algo: choose 'T1' or 'T2'")

    # Ensure DatetimeIndex output
    event_dates = pd.to_datetime(event_dates)

    # --- Build standardized event dataframe ---
    events_df = build_events_df(
        dates=event_dates,
        method=method_name,
        definition=f"{method_name}_{thres_sigma}σ_{persist_days}d_persist",
        data_source=data_source or "",
        threshold=f"{thres_sigma}σ",
        level_hpa=str(level_hpa),
        latitude=str(latitude),
        notes=f"{method_name}: ≥{persist_days}d persistence; min_gap={min_gap_days}d",
    )

    return events_df, event_dates

