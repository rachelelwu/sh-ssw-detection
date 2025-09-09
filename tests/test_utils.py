from __future__ import annotations
import numpy as np
import pandas as pd
import xarray as xr

def polar_cap_cos_weighted_mean(
    da: xr.DataArray,
    lat_min: float = -90,
    lat_max: float = -65,
    lat_dim: str = "lat",
) -> xr.DataArray:
    """Cosine-weighted polar-cap mean over [lat_min, lat_max]."""
    sub = da.sel({lat_dim: slice(lat_min, lat_max)})
    w = np.cos(np.deg2rad(sub[lat_dim]))
    w.name = "weights"
    return sub.weighted(w).mean(lat_dim)