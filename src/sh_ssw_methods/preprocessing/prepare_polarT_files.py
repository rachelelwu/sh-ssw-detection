"""
prepare_polarT.py
-----------------
Utilities to prepare polar-cap temperature diagnostics for SSW detection.

Includes:
    - return_temp_bands()
    - return_temp_polarCap()
    - daily_anom()
    - build_multiyear_temp()
"""

import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path


# --- 1. Latitude-band temperature time series ---

def return_temp_bands(
    temp_sh_era5: xr.DataArray,
    *,
    level: float = 10.0,
    save_file: bool = False,
    out_dir: str | Path = ".",
    data_source: str = "era5",
) -> xr.Dataset:
    """
    Compute daily, area-weighted mean temperature for 80–90S and 60–70S bands
    from an ERA5 temperature DataArray (e.g. var130). Optionally save to NetCDF.

    Parameters
    ----------
    temp_sh_era5 : xr.DataArray
        ERA5 temperature field with dims (time, lat, lon).
    level : float, optional
        Pressure level label (for filename).
    save_file : bool, optional
        If True, saves output as NetCDF.
    out_dir : str or Path, optional
        Directory to save file into (default current directory).
    data_source : str, optional
        Data source label for filename (default "era5").

    Returns
    -------
    xr.Dataset
        Dataset containing two daily mean time series:
          - T8090S : temperature averaged over 80–90S
          - T6070S : temperature averaged over 60–70S
    """

    T = temp_sh_era5.copy()
    
    # --- Ensure latitude is sorted from -90 to 0 (south to north) ---
    if np.any(np.diff(T.lat) > 0):
        # already increasing (−90 → 0)
        pass
    else:
        T = T.sortby("lat")

    # --- 80–90S band ---
    band1 = T.sel(lat=slice(-90, -80))
    w1 = np.cos(np.deg2rad(np.abs(band1["lat"])))
    ts_80_90S = band1.weighted(w1).mean(dim=("lat", "lon")).rename("T8090S")

    # --- 60–70S band ---
    band2 = T.sel(lat=slice(-70, -60))
    w2 = np.cos(np.deg2rad(np.abs(band2["lat"])))
    ts_60_70S = band2.weighted(w2).mean(dim=("lat", "lon")).rename("T6070S")

    # --- Combine ---
    ds_all = xr.Dataset({"T8090S": ts_80_90S, "T6070S": ts_60_70S})

    # --- Optional: save to file ---
    if save_file:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Find start and end year from time coordinate
        start_year = pd.to_datetime(ds_all.time.values[0]).year
        end_year = pd.to_datetime(ds_all.time.values[-1]).year

        fname = f"Tmid_polar_{int(level)}hPa_daily_{start_year}_{end_year}_{data_source}.nc"
        fpath = out_dir / fname

        ds_all.to_netcdf(fpath)
        print(f" Saved to {fpath}")

    return ds_all



# --- 2. Polar-cap (60–90S) mean temperature ---

def return_temp_polarCap(
    temp_sh_era5: xr.DataArray,
    *,
    level_hpa: float = 10.0,
    save_file: bool = False,
    out_dir: str | Path = ".",
    data_source: str = "era5",
) -> xr.Dataset:
    """
    Compute daily, area-weighted mean temperature over the 60–90°S polar cap
    from an already-loaded ERA5 temperature DataArray (e.g. var130).

    Parameters
    ----------
    temp_sh_era5 : xr.DataArray
        ERA5 temperature field with dims (time, lat, lon).
    level_hpa : float, optional
        Pressure level label (for filename).
    save_file : bool, optional
        If True, saves the output to NetCDF.
    out_dir : str or Path, optional
        Directory to save file into (default current directory).
    data_source : str, optional
        Data source label for filename (default "era5").

    Returns
    -------
    xr.Dataset
        Dataset containing one daily mean time series:
          - T6090S : temperature averaged over 60–90S
    """
    T = temp_sh_era5.copy()
    
    # --- Ensure latitude is sorted from -90 to 0 (south to north) ---
    if np.any(np.diff(T.lat) > 0):
        # already increasing (−90 → 0)
        pass
    else:
        T = T.sortby("lat")

    # --- Polar cap average (60–90S) ---
    band = T.sel(lat=slice(-90, -60))
    w = np.cos(np.deg2rad(np.abs(band["lat"])))  # area weights ~ cos(lat)
    ts_60_90S = band.weighted(w).mean(dim=("lat", "lon")).rename("T6090S")

    # --- Combine into a Dataset ---
    ds_all = xr.Dataset({"T6090S": ts_60_90S})

    # --- Optional save ---
    if save_file:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        start_year = pd.to_datetime(ds_all.time.values[0]).year
        end_year = pd.to_datetime(ds_all.time.values[-1]).year

        fname = f"polarT_{int(level_hpa)}hPa_daily_{start_year}_{end_year}_{data_source}.nc"
        fpath = out_dir / fname

        ds_all.to_netcdf(fpath)
        print(f" Saved to {fpath}")

    return ds_all



# --- 3. Daily anomaly computation ---

def daily_anom(
    ds: xr.Dataset,
    *,
    baseline=('1979-01-01', '2023-12-31'),
    drop_feb29=True,
    level_hpa: float = 10.0,
    save_file: bool = False,
    out_dir: str | Path = ".",
    data_source: str = "era5",
) -> xr.DataArray:
    """
    Compute daily anomalies relative to a baseline climatology for the polar-cap
    temperature (T6090S), optionally saving the output to a NetCDF file.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset containing variable 'T6090S'.
    baseline : tuple of str, optional
        (start_date, end_date) defining baseline period.
    drop_feb29 : bool, optional
        If True, remove Feb 29 to keep consistent DOY climatology.
    level_hpa : float, optional
        Pressure level label for filename.
    save_file : bool, optional
        If True, saves anomalies to NetCDF.
    out_dir : str or Path, optional
        Directory to save NetCDF file.
    data_source : str, optional
        Data source label for filename.

    Returns
    -------
    xr.DataArray
        Daily anomalies of T6090S relative to baseline climatology.
    """

    ds = ds.sortby("time")

    # --- Optionally drop Feb 29 ---
    if drop_feb29:
        is_feb29 = (ds.time.dt.month == 2) & (ds.time.dt.day == 29)
        ds = ds.sel(time=~is_feb29)

    # --- Baseline slice ---
    ds_base = ds.sel(time=slice(*baseline))

    # --- Daily climatology and stddev ---
    clim = ds_base["T6090S"].groupby("time.dayofyear").mean("time")
    stdv = ds_base["T6090S"].groupby("time.dayofyear").std("time")

    # --- Compute anomalies ---
    anom = ds["T6090S"].groupby("time.dayofyear") - clim
    anom = anom.rename("T6090S_anom")
    anom.attrs["description"] = f"Daily anomaly relative to {baseline[0]}–{baseline[1]}"
    anom.attrs["baseline"] = f"{baseline[0]} to {baseline[1]}"

    # --- Optional: save to file ---
    if save_file:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        start_year = pd.to_datetime(anom.time.values[0]).year
        end_year = pd.to_datetime(anom.time.values[-1]).year

        fname = f"polarT_anom_{int(level_hpa)}hPa_daily_{start_year}_{end_year}_{data_source}.nc"
        fpath = out_dir / fname

        anom.to_netcdf(fpath)
        print(f"Saved to {fpath}")

    return anom




# --- 4. Optional builder to merge multi-year temperature data ---

def output_polarT_prefiles(
    temp_sh_era5,
    start_year=1979,
    end_year=2023,
    ilev=10,
    out_dir="./sh_ssw_methods/processed",
    save_file=True,
):
    """
    Build and (optionally) save all key preprocessed polar-temperature datasets:
        - Temperature bands (T8090S, T6070S)
        - Polar-cap mean (T6090S)
        - Polar-cap anomaly relative to baseline

    Automatically checks whether output files already exist and skips recomputation
    when possible.

    Parameters
    ----------
    temp_sh_era5 : xr.DataArray
        ERA5 temperature field (time, lat, lon).
    start_year, end_year : int
        Start and end years for file naming and baseline period.
    ilev : int or float
        Pressure level in hPa.
    out_dir : str or Path, optional
        Directory where processed files are saved.
    save_file : bool, optional
        If True, saves results to NetCDF files.

    Returns
    -------
    tuple of xr.Dataset/xr.DataArray
        (ds_bands, polarT, polarT_anom)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- expected file paths ---
    bands_path = out_dir / f"Tmid_polar_{int(ilev)}hPa_daily_{start_year}_{end_year}_era5.nc"
    polar_path = out_dir / f"polarT_{int(ilev)}hPa_daily_{start_year}_{end_year}_era5.nc"
    anom_path  = out_dir / f"polarT_anom_{int(ilev)}hPa_daily_{start_year}_{end_year}_era5.nc"

    # --- Temperature bands (T8090S, T6070S) ---
    if bands_path.exists():
        print(f" Found existing file: {bands_path}")
        ds_bands = xr.open_dataset(bands_path)
    else:
        print(f"Creating: {bands_path.name}")
        ds_bands = return_temp_bands(
            temp_sh_era5,
            level_hpa=ilev,
            save_file=save_file,
            out_dir=out_dir,
        )

    # --- Polar-cap mean temperature (60–90S) ---
    if polar_path.exists():
        print(f"Found existing file: {polar_path}")
        polarT = xr.open_dataset(polar_path)
    else:
        print(f"Creating: {polar_path.name}")
        polarT = return_temp_polarCap(
            temp_sh_era5,
            level_hpa=ilev,
            save_file=save_file,
            out_dir=out_dir,
        )

    # --- Polar-cap temperature anomalies ---
    if anom_path.exists():
        print(f"Found existing file: {anom_path}")
        polarT_anom = xr.open_dataarray(anom_path)
    else:
        print(f"Creating: {anom_path.name}")
        polarT_anom = daily_anom(
            ds=polarT,
            baseline=(f"{start_year}-01-01", f"{end_year}-12-31"),
            level_hpa=ilev,
            save_file=save_file,
            out_dir=out_dir,
        )

    return ds_bands, polarT, polarT_anom

