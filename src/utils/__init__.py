from .xr_operators import cos_weighted_mean, remove_doy_climatology, remove_yearly_mean
from .event_iden_funcs import (
    group_consecutive_dates, 
    apply_persistence, 
    enforce_min_gap, 
    build_events_df
)



__all__ = [
    "cos_weighted_mean", "remove_doy_climatology",
    "remove_yearly_mean",
    "group_consecutive_dates", "apply_persistence", "enforce_min_gap", "build_events_df",
]

