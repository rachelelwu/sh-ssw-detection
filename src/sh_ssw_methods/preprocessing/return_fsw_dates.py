import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path


# --- identify FSW dates ---

# ---- helpers ----
def _cummax_time_bool(da_bool: xr.DataArray) -> xr.DataArray:
    """Cumulative max along 'time' for booleans, compatible with older xarray."""
    out = xr.apply_ufunc(
        lambda a: np.maximum.accumulate(a, axis=-1),
        da_bool.astype(np.int8),
        input_core_dims=[["time"]],
        output_core_dims=[["time"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.int8],
    )
    return out.astype(bool)

def _future_any_true(da_bool: xr.DataArray) -> xr.DataArray:
    """For each t, return True if any True occurs at/after t (inclusive)."""
    rev = da_bool.isel(time=slice(None, None, -1))
    rev_cum = _cummax_time_bool(rev)
    return rev_cum.isel(time=slice(None, None, -1))


# ---- main detector ----
def detect_fsw_SH(
    da: xr.DataArray,
    *,
    max_westerly_return: int = 10,   # no ≥(max+1)-day westerly return allowed
    smooth: int | None = None,       # optional rolling smoothing
    level_hpa: float = 10.0,
    save_file: bool = False,
    out_dir: str | Path = ".",
    data_source: str = "era5",
    overwrite: bool = False,         # NEW: optionally overwrite existing files
) -> xr.DataArray:
    """
    Detect final stratospheric warming (FSW) dates in the SH following Butler & Domeisen (2020).

    If save_file=True, checks whether the corresponding NetCDF file already exists.
    If it does and overwrite=False, the existing file is loaded and returned.

    Parameters
    ----------
    da : xr.DataArray
        Daily zonal-mean zonal wind (time) at ~60°S and specified level (e.g., 10 hPa).
    max_westerly_return : int, optional
        Maximum allowed consecutive westerly days after reversal (default 10).
    smooth : int or None, optional
        Length of optional running mean to smooth the time series.
    level_hpa : float, optional
        Pressure level label for filename and metadata.
    save_file : bool, optional
        If True, saves the detected FSW dates to NetCDF.
    out_dir : str or Path, optional
        Directory to save file into.
    data_source : str, optional
        Data source label for filename (default: "era5").
    overwrite : bool, optional
        If True, recomputes and overwrites existing file.

    Returns
    -------
    xr.DataArray
        FSW dates (datetime64) per season-year.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Construct expected filename ---
    start_year = pd.to_datetime(da.time.values[0]).year
    end_year = pd.to_datetime(da.time.values[-1]).year
    fname = f"FSW_dates_SH_{int(level_hpa)}hPa_60S_{start_year}_{end_year}_{data_source}.nc"
    fpath = out_dir / fname

    # --- Check if file already exists ---
    if save_file and fpath.exists() and not overwrite:
        print(f"Found existing file: {fpath}")
        return xr.open_dataarray(fpath)

    # --- Ensure DataArray and sort time ---
    da = da if isinstance(da, xr.DataArray) else da.to_array().squeeze()
    da = da.sortby("time")

    # Optional smoothing to reduce day-to-day noise
    if smooth not in (None, 0, 1):
        da = da.rolling(time=smooth, center=True, min_periods=max(1, smooth // 2)).mean()

    # Define "season year" (May–Apr)
    month = da["time"].dt.month
    sy = xr.where(month <= 4, da["time"].dt.year - 1, da["time"].dt.year).astype("int32").values
    da = da.assign_coords(season_year=("time", sy))

    fsw_dates, seasons = [], []
    L = int(max_westerly_return) + 1  # forbid ≥L-day westerly runs after reversal

    for sy_val, g in da.groupby("season_year"):
        # Candidate window: May 1 – Dec 31
        t0 = np.datetime64(f"{int(sy_val)}-05-01")
        t1 = np.datetime64(f"{int(sy_val)}-12-31")
        win = g.sel(time=slice(t0, t1))
        if win.time.size == 0:
            continue

        westerly = (win >= 0).fillna(False)
        easterly = (win < 0).fillna(False)

        had_westerly_before = _cummax_time_bool(westerly).shift(time=1, fill_value=False)
        long_west_end = (
            westerly.rolling(time=L, min_periods=L)
            .construct("win")
            .all("win")  # True at END of ≥L-day westerly runs
        )
        any_long_west_ahead = _future_any_true(long_west_end)
        candidate_mask = easterly & had_westerly_before & (~any_long_west_ahead)

        if bool(candidate_mask.any()):
            fsw_time = win.time.where(candidate_mask, drop=True)[0].item()
            fsw_dates.append(fsw_time)
            seasons.append(int(sy_val))

    fsw_da = xr.DataArray(
        np.array(fsw_dates, dtype="datetime64[ns]"),
        coords={"season_year": ("season_year", np.array(seasons, dtype="int32"))},
        dims=["season_year"],
        name=f"FSW_date_SH_{int(level_hpa)}hPa_60S_ret≤{max_westerly_return}d",
    )

    fsw_da.attrs.update({
        "title": "Final Stratospheric Warming (SH)",
        "definition": "Following Butler & Domeisen (2020)",
        "latitude": "-60°",
        "level_hpa": level_hpa,
        "max_westerly_return_days": max_westerly_return,
        "data_source": data_source,
    })

    # --- Save output if requested ---
    if save_file:
        fsw_da.to_netcdf(fpath)
        print(f"💾 Saved FSW dates to {fpath}")

    return fsw_da


