# /src/sh_ssw_methods/utils/xr_operators.py


from __future__ import annotations
import numpy as np


import xarray as xr


# Data processing 

def cos_weighted_mean(da: xr.DataArray, dim: str | None = None) -> xr.DataArray:
    """
    Cosine-weighted mean along latitude.
    If dim is None, tries common names ["lat", "latitude"].

    Input:
    da (xr.DataArray) - with pre-selected time period and latitude
    """
    if dim is None:
        for candidate in ["lat", "latitude"]:
            if candidate in da.dims:
                dim = candidate
                break
        else:
            raise ValueError("Latitude dimension not found. Expected one of ['lat', 'latitude'].")

    w = np.cos(np.deg2rad(da[dim]))
    w.name = "weights"
    return da.weighted(w).mean(dim)

def remove_doy_climatology(
    da: xr.DataArray, time_dim: str = "time"
) -> xr.DataArray:
    """Anomaly = value - long-term mean for that day-of-year."""
    return da.groupby(f"{time_dim}.dayofyear") - da.groupby(f"{time_dim}.dayofyear").mean(time_dim)


def remove_yearly_mean(
    da: xr.DataArray, time_dim: str = "time"
) -> xr.DataArray:
    """Remove the mean of each calendar year to mitigate trends."""
    return da.groupby(f"{time_dim}.year") - da.groupby(f"{time_dim}.year").mean(time_dim)

