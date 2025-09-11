# src/sh_ssw_methods/u_anom.py

from __future__ import annotations
import xarray as xr
import pandas as pd
import numpy as np
from .utils import (
    remove_doy_climatology,
    build_events_df
)


def detect_ssw_u_anom(
    u_daily: xr.DataArray,
    *,
    time_dim: str = "time",
    thres: float = -20.0,          # e.g., -20 (10 hPa), -11 (50 hPa)
    season_month_min: int = 5,     # months condition: (m > 4) & (m < 11) -> 5..10
    season_month_max: int = 10,
    min_gap_days: int = 20,
    persist_days: int = 10,
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

    
    dates_list = []  # to store datetime objects
    u_min1 = []

    for i in range(ndd - 1):
        if (s_anom.iloc[i+1] < thres) & (s_anom.iloc[i] >= thres) & (mns[i+1] > 4) & (mns[i+1] < 11):
            min1 = s_anom.iloc[i-min_gap_days:i].min()
            min2 = s_full.iloc[i:i+persist_days].min()
    
            if (min1 > thres) & (min2 > 0):
                date_obj = t[i+1]
                dates_list.append(date_obj)
                u_min1.append(min1)

    event_dates = pd.to_datetime(dates_list)
    u_anom_event = np.array(u_min1)
    
    # output into df
    # Ozone threshold events (dates and ozone value at onset)
    events_df = build_events_df(
        dates=event_dates,
        method="u_anom", #### change here the name later
        definition=f"u_anom_{int(thres)}m/s_{int(level_hpa)}hPa_{int(abs(latitude))}S",
        data_source=data_source or "",
        threshold=f"{int(thres)}m/s",
        level_hpa=str(level_hpa),
        latitude=str(latitude),
        notes=f"at_least_{persist_days}d_positive_after_ssw; min_gap={min_gap_days}d",
        # extra_cols={
        # "u_anom": u_anom_event, }, 
    )

    return events_df, event_dates