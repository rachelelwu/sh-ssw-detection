# src/sh_ssw_methods/methods/ozone_threshold.py
from __future__ import annotations
import xarray as xr
import pandas as pd
from .utils import (
    cos_weighted_mean,
    remove_doy_climatology,
    remove_yearly_mean, 
    build_events_df, 
    group_consecutive_dates,
    enforce_min_gap,
    apply_persistence
)


def detect_ozone_threshold(
    tco3: xr.DataArray,
    *,
    years=("1980","2020"),
    lat_band=(-90.0,-65.0),
    thresh_du=40.0,
    min_persist_days=3,
    min_gap_days=20,
    time_dim="time",
    lat_dim="lat", 
    data_source: str | None = None
):
    # 1. select years
    da = tco3.sel({time_dim: slice(*years)})

    # 2. polar cap weighted
    cap = cos_weighted_mean(da.sel({lat_dim: slice(*lat_band)}), dim=lat_dim)

    # 3. compute anom from daily climatology
    anom = remove_doy_climatology(cap, time_dim)

    # 4. remove trend by removing yearly mean
    fin = remove_yearly_mean(anom)
    
    # 5. filter out dates that satisfy threshold
    ind = fin.where(fin > thresh_du, drop=True)

    # 6. group consecutive days
    counts = group_consecutive_dates(ind['time'].to_index())

    # 7. test persistence criterion, need to last at least for 3 days
    persist = apply_persistence(counts, min_persist_days = 3)

    # 8. enforce minimum gap, events between need to separate at least 20 days
    tpersist = enforce_min_gap(persist, min_gap_days=20)

    event_dates = tpersist["first"]
    
    # 9. output into df
    # Ozone threshold events (dates only)
    events_df = build_events_df(
        dates=event_dates,
        method="ozone_threshold", #### change here the name later
        definition=f"ozone_du{int(thresh_du)}_{int(lat_band[0])}to{int(lat_band[1])}S",
        data_source=data_source or "",
        threshold=thresh_du,
        lat_band=str(list(lat_band)),
        notes=f"min_persist={min_persist_days}d; min_gap={min_gap_days}d",
    )
    
    return events_df, event_dates.values