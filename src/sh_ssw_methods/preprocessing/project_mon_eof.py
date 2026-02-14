import numpy as np
import xarray as xr

from pathlib import Path



def compute_daily_clim(
    ds: xr.Dataset | xr.DataArray,
    mon: int,
    *,
    drop_feb29: bool = False,
    lat_slice: tuple[float, float] = (-90, -20),
    ilev: int | None = None,
    varname: str = "var129",
    save_file: bool = False,
    out_dir: str | Path = "./processed",
) -> xr.DataArray:
    """
    Compute daily (day-of-month) climatology for a given month
    from an already loaded multi-year geopotential dataset.

    Parameters
    ----------
    ds : xr.Dataset or xr.DataArray
        Multi-year geopotential height dataset or dataarray.
    mon : int
        Calendar month (1–12).
    drop_feb29 : bool, optional
        Whether to drop Feb 29 (default False).
    lat_slice : tuple, optional
        Latitude range (default: (-90, -20)).
    ilev : int or None, optional
        Pressure level in hPa or Pa (used for subsetting and filename).
    varname : str, optional
        Variable name to use if `ds` is a Dataset (default: 'var129').
    save_file : bool, optional
        Save the climatology to NetCDF (default False).
    out_dir : str or Path, optional
        Directory for output NetCDF (default './processed').

    Returns
    -------
    xr.DataArray
        Daily climatology for that month: (day, lat, lon)
    """

    # --- Select data variable ---
    if isinstance(ds, xr.Dataset):
        if varname not in ds:
            raise KeyError(f"Variable '{varname}' not found in dataset.")
        z = ds[varname]
    elif isinstance(ds, xr.DataArray):
        z = ds
    else:
        raise TypeError("Input must be an xarray Dataset or DataArray.")

    # Convert geopotential to meters if needed
    if "units" in z.attrs and "m**2 s**-2" in z.attrs["units"]:
        z = z / 9.81
    z.attrs["units"] = "m"

    # --- Sort and subset latitude before processing ---
    if "lat" not in z.coords:
        raise ValueError("Latitude coordinate 'lat' not found in dataset.")
    if not np.all(np.diff(z["lat"]) > 0):
        z = z.sortby("lat")
    z = z.sel(lat=slice(*lat_slice))

    # --- Handle level selection if present ---
    if ilev is not None and "plev" in z.dims:
        lev_val = ilev 
        if lev_val not in z["plev"]:
            raise ValueError(
                f"Requested level {lev_val} not found in dataset. Available: {z['plev'].values}"
            )
        z = z.sel(plev=lev_val)
        # Keep plev dimension explicitly
        z = z.expand_dims("plev") if "plev" not in z.dims else z
        z = z.assign_coords(plev=[lev_val])
        # Ensure consistent dimension order
        if set(["time", "plev", "lat", "lon"]).issubset(z.dims):
            z = z.transpose("time", "plev", "lat", "lon")

    z.attrs["level_hPa"] = ilev
    z = z.sortby("time")

    # --- Select month ---
    z_mon = z.where(z.time.dt.month == mon, drop=True)
    if z_mon.time.size == 0:
        raise ValueError(f"No data found for month {mon}.")

    # --- Handle leap year day ---
    if drop_feb29 and mon == 2:
        z_mon = z_mon.sel(time=~((z_mon.time.dt.month == 2) & (z_mon.time.dt.day == 29)))

    # --- Compute daily climatology across years ---
    clim = z_mon.groupby("time.day").mean("time")
    clim.name = "z_daily_clim"
    clim.attrs.update(
        dict(
            description=f"Daily climatology (DoM) for month={mon}, 1979–2023",
            units="m",
        )
    )

    # --- Optional save ---
    if save_file:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        start_year = pd.to_datetime(z.time.values[0]).year
        end_year = pd.to_datetime(z.time.values[-1]).year
        ilev_label = f"_{int(ilev)}hPa" if ilev is not None else ""
        out_path = out_dir / f"daily_clim_z{ilev_label}_Mon{mon:02d}_{start_year}_{end_year}_era5.nc"

        if not out_path.exists():
            clim.to_netcdf(out_path)
            print(f"Saved climatology to {out_path}")
        else:
            print(f"File already exists: {out_path}")

    return clim


def project_daily_pc1(
    monthly_eof: xr.DataArray,
    gh_data: xr.Dataset | xr.DataArray,
    *,
    month: int,
    lat_hi: float = -20,
    lat_lo: float = -90,
    varname: str = "var129",
    ilev: int | None = None,
) -> xr.DataArray:
    """
    Project daily geopotential height anomalies onto a monthly EOF pattern
    to obtain daily PC1 time series for a given level and month.

    Parameters
    ----------
    monthly_eof : xr.DataArray
        EOF1 spatial pattern (lat, lon) for a specific month and level.
    gh_data : xr.Dataset or xr.DataArray
        Multi-year geopotential height dataset containing 'time', 'lat', and 'lon'.
    gh_clim : xr.DataArray
        Daily climatology for that same month (day, lat, lon).
    month : int
        Calendar month (1–12) corresponding to the EOF pattern and climatology.
    lat_hi, lat_lo : float, optional
        Latitude limits for projection (default: -20 to -90).
    varname : str, optional
        Variable name if `gh_data` is a Dataset (default: 'var129').
    ilev : int or None, optional
        Level in hPa or Pa to select if multi-level data are provided.

    Returns
    -------
    xr.DataArray
        Daily standardized PC1 time series for that month (time,).
    """

    # --- Select variable ---
    if isinstance(gh_data, xr.Dataset):
        if varname not in gh_data:
            raise KeyError(f"Variable '{varname}' not found in dataset.")
        z = gh_data[varname]
    elif isinstance(gh_data, xr.DataArray):
        z = gh_data
    else:
        raise TypeError("gh_data must be an xarray Dataset or DataArray.")

    # --- Convert geopotential units to meters ---
    z = z / 9.81
    z.attrs["units"] = "m"

    # --- Select pressure level if available ---
    if ilev is not None and "plev" in z.dims:
        lev_val = ilev
        if lev_val not in z["plev"]:
            raise ValueError(
                f"Requested level {lev_val} not found. Available: {z['plev'].values}"
            )
        z = z.sel(plev=lev_val)

    # --- Sort latitude and restrict to region ---
    if not np.all(np.diff(z["lat"]) > 0):
        z = z.sortby("lat")
    z = z.sel(lat=slice(lat_lo, lat_hi))

    # --- Select only the target month ---
    z_mon = z.sel(time=z.time.dt.month == month)
    if z_mon.time.size == 0:
        raise ValueError(f"No data found for month {month} in provided dataset.")
        
    gh_clim = compute_daily_clim(z_mon,
            mon=month,
            ilev=ilev, # Pa
            varname=varname,
            drop_feb29=True,
            save_file=False)
        
    # --- Compute anomalies relative to given month’s climatology ---
    # Align day-of-month to climatology day coordinate
    gh_anom = z_mon.groupby("time.day") - gh_clim

    # --- Apply sqrt(cos(lat)) weighting ---
    w = np.sqrt(np.cos(np.deg2rad(gh_anom["lat"])))
    w2d = w.broadcast_like(monthly_eof)
    wgt_gh_anom = gh_anom * w2d

    # --- Project anomalies onto EOF pattern ---
    # pc(t) = ∑_lat ∑_lon [ wgt_gh_anom(t,lat,lon) * EOF1(lat,lon) ]
    pc1 = (wgt_gh_anom * monthly_eof).sum(dim=("lat", "lon"))

    return pc1


def compute_pc1_all_months(
    lv: int,
    *,
    base_dir: str | Path = "./processed",
    months: list[int] | None = None,
    lat_range: tuple[float, float] = (-90, -20),
    varname: str = "var129",
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """
    Project daily geopotential height anomalies onto monthly EOF1 patterns
    for all years, returning (year, day) arrays for each month.

    Parameters
    ----------
    lv : int
        Pressure level in Pa.
    base_dir : str or Path, optional
        Directory containing processed EOFs and geopotential data.
    months : list[int], optional
        List of months (default: [1..12]).
    lat_range : tuple, optional
        Latitude limits for projection (default: (-90, -20)).
    varname : str, optional
        Variable name in geopotential dataset (default: 'var129').

    Returns
    -------
    dict[int, tuple[np.ndarray, np.ndarray]]
        {month: (pc1_all, years)} where
        - pc1_all: ndarray (n_years, n_days)
        - years: array of corresponding year labels
    """

    base_dir = Path(base_dir)
    results = {}
    months = months or range(1, 13)

    # --- Load geopotential height dataset once ---
    gh_path = base_dir / f"era5_geopot_{int(lv/100)}hPa_1979_2023.nc"
    if not gh_path.exists():
        raise FileNotFoundError(f"Missing geopotential file: {gh_path}")
    gh_data = xr.open_dataset(gh_path)

    # --- Loop through requested months ---
    for mon in months:
        print(f"\n=== Computing PC1 for z{int(lv/100)}hPa Month {mon:02d} ===")

        eof_path = base_dir / f"zeof/ERA5.EOF1.z{int(lv/100)}.Mon{mon:02d}.nc"
        if not eof_path.exists():
            print(f" Missing EOF file → {eof_path.name}, skipping month {mon:02d}")
            continue

        eof = xr.open_dataarray(eof_path)

        try:
            pc1 = project_daily_pc1(
                monthly_eof=eof,
                gh_data=gh_data,
                month=mon,
                ilev=lv,  
                lat_hi=lat_range[1],
                lat_lo=lat_range[0],
                varname=varname,
            )
        except Exception as e:
            print(f" Failed to compute PC1 for month {mon:02d}: {e}")
            continue

        # --- Convert to (year, day) array ---
        grouped = list(pc1.groupby("time.year"))
        years = np.array([int(y) for y, _ in grouped])
        n_days = max(len(g) for _, g in grouped)
        pc1_all = np.full((len(years), n_days), np.nan)

        for i, (_, g) in enumerate(grouped):
            vals = g.values
            pc1_all[i, :len(vals)] = vals

        
        # standardise
        pc1_all_std = (pc1_all - np.mean(pc1_all, axis=0))/ np.std(pc1_all, axis=0)

        results[mon] = (pc1_all_std, years)

    
    print("\n Finished computing PC1 for all requested months.")
    return results


def write_daily_pc1(pc1_std_all, years, mon, lv, out_dir="."):
    """
    Write standardized daily PC1 for each year to text files like:
    daily.pc1.z10.197906.txt

    Parameters
    ----------
    pc1_std_all : (year, day) array-like
        Standardized PC1 values for a given month (e.g., shape (45, 30))
    years : list or array of int
        Year labels corresponding to axis 0
    mon : int
        Month number (1–12)
    lv : int
        Pressure level (e.g., 10 for 10 hPa)
    out_dir : str or Path
        Directory to save text files
    """
    pc1_std_all = np.asarray(pc1_std_all)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for i, yr in enumerate(years):
        vals = pc1_std_all[i, :]
        # format: year + 30 space-padded floats (width=9, precision=3)

        base = Path(out_dir)
        zeof_dir = base / "zeof"
        zeof_dir.mkdir(parents=True, exist_ok=True)

        line = f"{yr}" + "".join([f"{v:9.3f}" for v in vals])
        fname = zeof_dir / f"daily.pc1.z{lv}.{yr}{mon:02d}.txt"
        with open(fname, "w") as f:
            f.write(line + "\n")



def save_pc1_ascii_all_months(pc1_dict, lv: int, *, out_dir: str | Path = "./processed"):
    """
    Save standardized daily PC1 values for each month to ASCII files.

    Parameters
    ----------
    pc1_dict : dict[int, tuple[np.ndarray, np.ndarray]]
        Output from `compute_pc1_all_months`.
        Keys are months; values are (pc1_all, years).
    lv : int
        Pressure level in hPa (e.g., 10).
    out_dir : str or Path, optional
        Output base directory (default: './processed').
    """

    out_dir = Path(out_dir)
    for mon, (pc1_all, years) in pc1_dict.items():
        print(f"\n--- Writing standardized PC1 files for Month {mon:02d} ---")
        write_daily_pc1(pc1_all, years, mon, lv, out_dir=out_dir)
    print("\n All PC1 text files written successfully.")



def output_daily_pc1(ilev, months, base_dir):

    # Step 1: Compute PC1 arrays for all months
    pc1_all_std = compute_pc1_all_months(
        lv=ilev,
        base_dir=base_dir,
        months=months,  # June–December
    )
    
    # Step 2: Write to ASCII text files
    save_pc1_ascii_all_months(pc1_all_std, lv=ilev, out_dir=base_dir)


