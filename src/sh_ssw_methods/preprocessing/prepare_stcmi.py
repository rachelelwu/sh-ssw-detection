import xarray as xr
import numpy as np
from eofs.xarray import Eof

import pandas as pd
from pathlib import Path


import os



def zonal_area_mean(
    da: xr.DataArray | xr.Dataset,
    lat_band: tuple[float, float] = (-65, -55),
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """
    Compute zonal mean and area-weighted mean over a latitude band.
    Automatically sorts latitude to ascending order.

    Parameters
    ----------
    da : xr.DataArray or xr.Dataset
        Input data with dimensions including `lat` and `lon`.
        e.g. uwind(time, plev, lat, lon)
    lat_band : tuple(float, float)
        Latitude range (south, north), e.g. (-65, -55) for 55–65°S.
    lat_name, lon_name : str
        Coordinate names (defaults are 'lat' and 'lon').

    Returns
    -------
    da_mean : xr.DataArray
        Zonal + area-weighted mean, dimensions (time, plev).
    """

    # --- Step 0: ensure latitude ascending (S -> N) ---
    if da[lat_name][0] > da[lat_name][-1]:
        da = da.sortby(lat_name)

    # --- Step 1: take zonal mean ---
    if lon_name in da.dims:
        da_zm = da.mean(lon_name)
    else:
        da_zm = da

    # --- Step 2: select latitude band ---
    latmin, latmax = sorted(lat_band)
    da_sel = da_zm.sel({lat_name: slice(latmin, latmax)})

    # --- Step 3: cosine latitude weights ---
    lat_radians = np.deg2rad(da_sel[lat_name])
    weights = np.cos(lat_radians)
    weights /= weights.sum()  # normalize so they sum to 1

    # --- Step 4: weighted mean over latitude ---
    da_wt = (da_sel * weights).sum(lat_name)

    # add year and month coordinates
    da_wt = da_wt.assign_coords(
        year=da_wt.time.dt.year,
        month=da_wt.time.dt.month
    )

    
    # reshape: create MultiIndex and unstack to get (year, month, plev)
    da_year_month = (
        da_wt
        .set_index(time=('year', 'month'))
        .unstack('time')
        .transpose('year', 'month', 'plev')
    )
    

    return da_year_month


def output_stcmI(ds, outfname, save_txt=True):

    anom = ds

    # Rename 'year' → 'time' so Eof sees it as the sample axis
    anom_year_space = anom.rename(year='time')
    
    solver = Eof(anom_year_space)
    
    eofs = solver.eofs(neofs=3)        # (mode, space)
    pcs  = solver.pcs(npcs=3, pcscaling=1)  # (time=44, mode)
    varfrac = solver.varianceFraction(neigs=3)

    # -------------------------
    # PC1 time series (standardized & sign-adjusted)
    # -------------------------
    pc1 = pcs[:,0] * -1.0
    pc1 = (pc1 - pc1.mean()) / pc1.std()

    if save_txt:
    # -------------------------
    # Save PC1 time series
    # -------------------------
        np.savetxt(outfname, pc1.values)

    return pc1



def compute_stcmI(xu_era5, outdir, save_txt=True):
    """
    Compute the STC Mode Index (STCMI) from ERA5 zonal wind data.

    Steps:
      1. Compute zonal and latitude-band means (creates uwind anomalies).
      2. Compute EOFs of (plev × month) anomalies.
      3. Extract standardized PC1 as STCMI.
      4. Optionally save PC1 as text file named by the actual year range.
      5. Skip computation if file already exists.

    Parameters
    ----------
    xu_era5 : str or xarray.Dataset
        Input ERA5 uwind dataset (time, plev, lat, lon) or path to NetCDF.
    outdir : str
        Directory to save output file (e.g., "results/").
    save_txt : bool, default=True
        If True, saves standardized PC1 time series as "stcmI.<start>-<end>.txt".

    Returns
    -------
    dict
        {
            'pc1'     : xarray.DataArray, standardized leading PC
            'pcs'     : xarray.DataArray, first 3 PCs
            'eofs'    : xarray.DataArray, first 3 EOFs
            'varfrac' : xarray.DataArray, variance fractions
        }
    """

    # -----------------------------
    # 0. Load dataset
    # -----------------------------

    ds = xu_era5
    # Check that 'time' exists
    if 'time' not in ds.coords:
        raise ValueError("Input dataset must have a 'time' coordinate.")

    # -----------------------------
    # 1. Get year range dynamically
    # -----------------------------
    years = ds['time'].dt.year
    start_year = int(years.min())
    end_year = int(years.max())

    # -----------------------------
    # 2. Prepare output directory & file name
    # -----------------------------
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"stcmI.{start_year}-{end_year}.txt")

    if os.path.exists(outfile):
        print(f" File already exists: {outfile}")
        print("→ Skipping computation. To recompute, delete the file first.")
        return None

    # -----------------------------
    # 3. Compute zonal & lat-band means
    # -----------------------------
    anom = zonal_area_mean(ds)  # returns uwind anomalies (time, plev)

    # -----------------------------
    # 4. Compute EOFs and PC1
    # -----------------------------
    result = output_stcmI(anom, outfile, save_txt=save_txt)

    print(f" STCMI successfully computed for {start_year}–{end_year}")
    print(f"   Saved to: {outfile}")

    return result



