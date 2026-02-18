# src/sh_ssw_methods/ozone_threshold.py
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
    lat_band=(-90.0,-65.0),
    thresh_du=40.0,
    min_persist_days=3,
    min_gap_days=20,
    time_dim="time",
    lat_dim="lat", 
    data_source: str | None = None
):

    # raise error if input is not xr.DataArray
    if not isinstance(tco3, xr.DataArray):
        raise TypeError(
            "detect_ozone_threshold expects an xarray.DataArray. "
            "If you have a Dataset, select a variable first, e.g. ds['tcO3']."
        )

    # polar cap weighted
    cap = cos_weighted_mean(tco3.sel({lat_dim: slice(*lat_band)}), dim=lat_dim)

    # compute anom from daily climatology
    anom = remove_doy_climatology(cap, time_dim)

    # remove trend by removing yearly mean
    fin = remove_yearly_mean(anom)
    
    # filter out dates that satisfy threshold
    ind = fin.where(fin > thresh_du, drop=True)

    # group consecutive days
    counts = group_consecutive_dates(ind['time'].to_index())

    # test persistence criterion, need to last at least for 3 days
    persist = apply_persistence(counts, min_persist_days = 3)

    # enforce minimum gap, events between need to separate at least 20 days
    tpersist = enforce_min_gap(persist, min_gap_days=20)

    event_dates = pd.to_datetime(tpersist["first"])

    # also look up corresponding ozone polar cap anomalies
    series_for_value = fin          
    
    # robust alignment: convert to pandas Series and reindex on event dates
    s_val = pd.Series(
        series_for_value.values,
        index=pd.to_datetime(series_for_value[time_dim].to_index())
    )
    ozone_onset = s_val.reindex(event_dates).to_numpy()  # value at event 'first' day
    
        
    # output into df
    # Ozone threshold events (dates and ozone value at onset)
    events_df = build_events_df(
        dates=event_dates,
        method="ozone_threshold", #### change here the name later
        definition=f"ozone_du{int(thresh_du)}_{int(lat_band[0])}to{int(lat_band[1])}S",
        data_source=data_source or "",
        threshold=f"{thresh_du}DU",
        lat_band=str(list(lat_band)),
        notes=f"min_persist={min_persist_days}d; min_gap={min_gap_days}d",
        extra_cols={
        "ozone_onset_DU": ozone_onset, }, 
    )
    
    return events_df, event_dates.values