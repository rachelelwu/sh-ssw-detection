# src/sh_ssw_methods/u_tend.py
from __future__ import annotations
import numpy as np
import pandas as pd
import xarray as xr

from .utils.event_iden_funcs import (
    build_events_df
)

from .utils.xr_operators import (
    remove_doy_climatology,
)


def detect_ssw_u_tend(
    u_daily: xr.DataArray,
    *,
    time_dim: str = "time",
    thres: float = -35.0,          # e.g., -35 (10 hPa), -19 (50 hPa)
    season_month_min: int = 5,     # months condition: (m > 4) & (m < 11) -> 5..10
    season_month_max: int = 10,
    min_gap_days: int = 20,
    persist_days: int = 10,
    rng: int = 7, # +/- 7 days around current day to define u tendency
    level_hpa: float | None = None,
    data_source: str | None = None,
    latitude: str | None = "-60",
) -> tuple[pd.DataFrame, pd.DatetimeIndex]:


    u_anom = remove_doy_climatology(u_daily)

    # To pandas for simple window logic
    t = pd.to_datetime(u_daily[time_dim].to_index())
    s_full = pd.Series(u_daily.values, index=t)
    s_anom = pd.Series(u_anom.values,  index=t)

    mns = t.month
    ndd = len(s_full)

     #Calculate delta
    delta=np.zeros(ndd)
    for i in range(rng,ndd-1-rng):
        delta[i]=s_full.iloc[i+rng]-s_full.iloc[i-rng] 
    
    dates_list = []  # to store datetime objects

    for i in range(rng, ndd - 1 - rng):
        if (delta[i] < thres) & (delta[i-1] >= thres) & (mns[i+1] > 4) & (mns[i+1] < 11):
            min1 = np.amin(delta[i-min_gap_days:i])
            min2 = s_full.iloc[i:i+persist_days].min()
    
            if (min1 > thres) & (min2 > 0):
                date_obj = t[i+1]
                dates_list.append(date_obj)

    event_dates = pd.to_datetime(dates_list)
    
    # output into df
    events_df = build_events_df(
        dates=event_dates,
        method="u_tend", #### change here the name later
        definition=f"u_tend_{int(thres)}m/s_{int(level_hpa)}hPa_{int(np.abs(latitude))}S",
        data_source=data_source or "",
        threshold=f"{int(thres)}m/s",
        level_hpa=str(level_hpa),
        latitude=str(latitude),
        notes=f"u_tend < thres at_least_{persist_days}d; u_tend defined as u_full difference between +/-{rng} days from current day; min_gap={min_gap_days}d",
        # extra_cols={
        # "u_anom": u_anom_event, }, 
    )

    return events_df, event_dates