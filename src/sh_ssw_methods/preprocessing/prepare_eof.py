import numpy as np
import xarray as xr

from eofs.xarray import Eof

import os
import glob
from pathlib import Path


def eof(X,n=-1,detrend='constant',eof_in=None):
	"""Principal Component Analysis / Empirical Orthogonal Functions / SVD

		Uses Singular Value Decomposition to find the dominant modes of variability.
		The field X can be reconstructed with Y = dot(EOF,PC) + X.mean(axis=time)

		INPUTS:
			X	-- Field, shape (time x space).
			n	-- Number of modes to extract. All modes if n < 0
			detrend -- detrend with global mean ('constant')
						  or linear trend ('linear')
		    eof_in  -- If not None, compute PC by projecting eof onto X.
		OUTPUTS:
			EOF - Spatial modes of variability
			PC  - Temporal evolution of EOFs - only output if eof_in is not None
			E   - Explained value of variability
			u   - spatial modes
			s   - variances
			v   - temporal modes

    from aostools from Martin Jucker
	"""
	is_xr = False
	try:
		import xarray as xr
		if isinstance(X,xr.DataArray):
			is_xr = True
			dims = []
			for dim in X.dims[1:]:
				dims.append(X[dim])
			tdim = [X[X.dims[0]]]
			X = X.values
	except:
		pass
	import scipy.signal as sg
	# make sure we have a matrix time x space
	shpe = X.shape
	if len(shpe) > 2:
		X = X.reshape([shpe[0],np.prod(shpe[1:])])
		if eof_in is not None:
			if len(eof_in.shape) > 2:
				eof_in = eof_in.reshape([np.prod(eof_in.shape[:-1]),eof_in.shape[-1]])
			else:
				eof_in = eof_in.reshape([np.prod(eof_in.shape),1])
	# take out the time mean or trend
	X = sg.detrend(X.transpose(),type=detrend)
	if eof_in is not None:
		if eof_in.shape[-1] == X.shape[0]:
			PC =  np.matmul(eof_in, X)
			eof_norm = np.dot(eof_in.transpose(),eof_in)
			return np.dot(PC,np.linalg.inv(eof_norm))
		else:
			PC = np.matmul(eof_in.transpose(), X)
			eof_norm = np.dot(eof_in.transpose(),eof_in)
			return np.dot(PC.transpose(),np.linalg.inv(eof_norm)).transpose()
		# return sg.detrend(PC,type='constant')
	# perform SVD - v is actually V.H in X = U*S*V.H
	u,s,v = np.linalg.svd(X, full_matrices=False)
	# now, u contains the spatial, and v the temporal structures
	# s contains the variances, with the same units as the input X
	# u.shape = (space, modes(space)), v.shape = (modes(space), time)

	# get the first n modes, in physical units
	#  we can either project the data onto the principal component, X*V
	#  or multiply u*s. This is the same, as U*S*V.H*V = U*S
	if n < 0:
		n = s.shape[0]
	EOF = np.dot(u[:,:n],np.diag(s)[:n,:n])
	# time evolution is in v
	PC  = v[:n,:]
	# EOF wants \lambda = the squares of the eigenvalues,
	#  but SVD yields \gamma = \sqrt{\lambda}
	s2 = s*s
	E   = s2[:n]/sum(s2)
	# now we need to make sure we get everything into the correct shape again
	u = u[:,:n]
	s = s[:n]
	v = v.transpose()[:,:n]
	if len(shpe) > 2:
		# replace time dimension with modes at the end of the array
		newshape = list(shpe[1:])+[n]
		EOF = EOF.reshape(newshape)
		u   = u	 .reshape(newshape)
	if is_xr: # return xarray dataarrays
		mode = [('n',np.arange(1,n+1))]
		EOF = xr.DataArray(EOF,coords=dims+mode,name='EOF')
		PC  = xr.DataArray(PC,coords=tdim+mode,name='PC')
		E   = xr.DataArray(E,coords=mode,name='E')
	return EOF,PC,E,u,s,v


def eof_analysis_eof1(
    pwtanom: xr.DataArray,
    *,
    lv: int,
    mon: int,
    neofs: int = 2,
    lat_name: str | None = None,
    lon_name: str | None = None,
    sign_lat: float = -60.0,
    sign_lon: float = 90.0,
    out_dir: str | Path = ".",
    save: bool = True,
    overwrite: bool = False,     # NEW: optional overwrite control
):
    """
    EOF analysis on weighted anomalies (pwtanom):
      - EOFs over (lat, lon), PCs over time (unit variance)
      - sign convention: ensure EOF1(lat≈-60, lon≈90E) > 0
      - regression pattern = <pwtanom * PC1_std> over time
      - saves EOF1 to ERA5.EOF1.z{lv}.Mon{mm}.nc if save=True

    Parameters
    ----------
    pwtanom : xr.DataArray
        Weighted anomalies with dims (time, lat, lon). Already demeaned across years for a month.
    lv, mon : int
        Pressure level (for filename) and calendar month.
    neofs : int
        Number of EOF modes to compute (default 2).
    lat_name, lon_name : str or None
        Coordinate names; auto-detected if None.
    sign_lat, sign_lon : float
        Reference point for sign convention (nearest grid point is used).
    out_dir : path-like
        Directory to save output NetCDF.
    save : bool
        Whether to write EOF1 to file.
    overwrite : bool
        Whether to overwrite an existing file (default False).

    Returns
    -------
    eof1 : xr.DataArray       (lat, lon)
    pc1_std : xr.DataArray    (time,)
    varfrac1 : float          variance fraction of EOF1 (0..1)
    """
    # --- coord names ---
    if lat_name is None:
        lat_name = "lat" if "lat" in pwtanom.dims else "latitude"
    if lon_name is None:
        lon_name = "lon" if "lon" in pwtanom.dims else "longitude"

    # ensure correct order (time, lat, lon)
    da = pwtanom.transpose("time", lat_name, lon_name)

    # --- output path ---
    out_dir = Path(out_dir)
    zeof_dir = out_dir / "zeof"
    zeof_dir.mkdir(parents=True, exist_ok=True)
    lv_label = int(lv / 100)  # hPa
    out_path = zeof_dir / f"ERA5.EOF1.z{lv_label}.Mon{mon:02d}.nc"

    # --- if file exists ---
    if save and out_path.exists() and not overwrite:
        print(f"Existing EOF1 file found: {out_path}")
        eof1 = xr.open_dataarray(out_path)
        return eof1

    # --- compute EOFs ---
    solver = Eof(da, center=False)
    eofs = solver.eofs(neofs=neofs, eofscaling=0)
    pcs = solver.pcs(npcs=neofs, pcscaling=1)
    varf = solver.varianceFraction()

    print(f"Explained variance: {varf.values[:neofs]}")

    eof1 = eofs.isel(mode=0)
    pc1_std = pcs.isel(mode=0)
    varfrac1 = float(varf.isel(mode=0).values)

    # --- sign convention ---
    lat_vals = eof1[lat_name].values
    lon_vals = eof1[lon_name].values
    ilat = int(np.argmin(np.abs(lat_vals - sign_lat)))
    lon360 = (lon_vals % 360 + 360) % 360
    sign_lon360 = sign_lon % 360
    ilon = int(np.argmin(np.abs(lon360 - sign_lon360)))

    flip = -1.0 if float(eof1.values[ilat, ilon]) < 0 else 1.0
    eof1 *= flip
    pc1_std *= flip

    # --- save to file ---
    if save:
        eof1.rename("eof1").to_dataset().to_netcdf(out_path)
        print(f"Saved EOF1 to {out_path}")

    return eof1


# --- function to loop over months ---

def compute_monthly_eof_series(
    da: xr.DataArray,
    *,
    lv: int,
    mons: list[int] | np.ndarray = range(1, 13),
    lat_min: float = -90,
    lat_max: float = -20,
    out_dir: str | Path = "./sh_ssw_methods/processed",
    save: bool = True,
    overwrite: bool = False,
):
    """
    Compute monthly EOF1 patterns from a multi-year geopotential height file.

    For each month:
      - Selects the corresponding subset from the full dataset
      - Checks if requested level matches data
      - Sorts latitudes (ascending)
      - Applies √cos(lat) weighting
      - Computes anomalies across years
      - Runs EOF1 analysis via `eof_analysis_eof1`

    Parameters
    ----------
    da : xr.DataArray
        Multi-year geopotential height (time, [plev,] lat, lon) [m].
    lv : int
        Pressure level in Pa (used for file naming and consistency check).
    mons : list[int] or np.ndarray, optional
        List of months (1–12) to process.
    lat_min, lat_max : float, optional
        Latitude range to select (default: -90 to -20 for SH).
    out_dir : str or Path, optional
        Directory for output NetCDFs.
    save : bool, optional
        Whether to save EOF1 NetCDF files.
    overwrite : bool, optional
        Whether to overwrite existing files.

    Returns
    -------
    dict[int, xr.DataArray]
        Dictionary of EOF1 patterns for each processed month.
    """

    results = {}
    out_dir = Path(out_dir)

    # --- Level consistency check ---
    if "plev" in da.dims:
        da_levels = da["plev"].values
        lv_pa = lv 

        if lv_pa not in da_levels:
            raise ValueError(
                f"Requested level {lv_pa} Pa not found in DataArray. "
                f"Available levels: {np.round(da_levels, 1)}"
            )
        else:
            da = da.sel(plev=lv_pa)
            print(f"Using level {lv_pa} Pa from data.")
    else:
        print("No 'plev' dimension in DataArray. Assuming input is single-level data. ")

    # --- Sort latitude if needed ---
    if not np.all(np.diff(da["lat"]) > 0):
        da = da.sortby("lat")

    # --- Restrict latitude range ---
    da = da.sel(lat=slice(lat_min, lat_max))

    # --- Loop over requested months ---
    for mon in mons:
        print(f"\nProcessing month: {mon:02d}")

        # Subset by month
        da_mon = da.sel(time=da.time.dt.month == mon)
        if da_mon.time.size == 0:
            print(f"No data for month {mon:02d}; skipping.")
            continue

        # Compute yearly mean for this month
        da_m = da_mon.groupby("time.year").mean("time")

        # Apply √cos(lat) weighting
        w = np.sqrt(np.cos(np.deg2rad(da_m["lat"])))
        da_weighted = da_m * w

        # Compute anomalies across years
        clim = da_weighted.mean("year")
        anom = (da_weighted - clim).rename(year="time").astype("float32")

        # Run EOF analysis
        eof1 = eof_analysis_eof1(
            anom,
            lv=lv,
            mon=mon,
            out_dir=out_dir,
            save=save,
            overwrite=overwrite,
        )

        results[mon] = eof1

    print("\nAll requested months processed.")
    return results


