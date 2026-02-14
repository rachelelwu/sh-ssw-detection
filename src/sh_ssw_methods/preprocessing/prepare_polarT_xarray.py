import numpy as np
import pandas as pd
import glob

import xarray as xr

from pathlib import Path

from pathlib import Path


# T6080, T8090
def aug_nov_by_year(ds):
    """
    ds: xarray.Dataset with daily 'time' and variables 'T8090S', 'T6070S'
    returns: Dataset with dims ('year','day') and a 'date' variable for reference.
    """
    ds = ds.sortby('time')

    # Keep Aug–Nov only
    sub = ds.sel(time=ds.time.dt.month.isin([8, 9, 10, 11]))

    # Index within each year: 0..121 (Aug–Nov has 122 days in non-leap years)
    ti = sub.indexes['time']  # pandas.DatetimeIndex
    day_idx = pd.Series(np.arange(ti.size), index=ti).groupby(ti.year).cumcount().to_numpy()

    # Add coords for pivoting
    sub = sub.assign_coords(
        year=('time', sub.time.dt.year.values),
        day=('time', day_idx),
        date=('time', sub.time.values)  # keep the actual datestamp as a variable
    )

    # Pivot time -> (year, day)
    wide = sub[['T8090S', 'T6070S', 'date']].set_index(time=['year', 'day']).unstack('time')

    # Tidy output
    out = xr.Dataset(
        {
            'T8090S': wide['T8090S'],
            'T6070S': wide['T6070S'],
            'date':   wide['date'],  # same shape (year, day), useful for debugging / alignment
        }
    )
    # Optional: name dims explicitly (xarray may already name them 'year' and 'day')
    out = out.transpose('year', 'day')

    # (Optional) Sanity check: every year should have 122 days
    # If you have missing days, this will vary; handle reindexing per-year if needed.
    # counts = out['date'].count(dim='day').to_pandas()

    return out


def output_Tvars(fpath, start_date="1979-01-01", end_date="2021-12-31"):
    xT = xr.open_dataset(fpath)

    # select time period
    xT = xT.sel(time=slice(start_date, end_date))
    
    # extract only aug to nov
    T_augnov = aug_nov_by_year(xT)
    
    # compute climatologies
    cT8090 = T_augnov['T8090S'].isel(year=slice(0, 43)).mean('year')
    cT6070 = T_augnov['T6070S'].isel(year=slice(0, 43)).mean('year')
    
    # delta climatology (80–90S minus 60–70S)
    cdelT = cT8090 - cT6070
    
    # Anomalies (broadcast day climatology across years)
    T8090anom = T_augnov['T8090S'] - cT8090
    T6070anom = T_augnov['T6070S'] - cT6070

    # Indices where climatological delta ≤ 0
    cidx = np.where(cdelT.values <= 0)[0]
    ncidx = cidx.size
    
    if ncidx > 0:
        first_idx, last_idx = cidx[0], cidx[-1]
        # print(first_idx, last_idx)
    
    # Year-by-year daily difference
    delT = T_augnov['T8090S'] - T_augnov['T6070S']          # (year, day)
    
    # Std dev across years for each day
    stddv = delT.std(dim='year')                          # (day,)
    
    return T_augnov, cdelT, delT, stddv


# final warming dates
def output_fwidx(fwdate_path, T_augnov):
    fwdate = xr.load_dataarray(fwdate_path)
    
    fwmask = T_augnov['date'] == fwdate.rename({'season_year':'year'}).broadcast_like(T_augnov['date'])
    fw_idx = fwmask.argmax(dim='day') 

    return fw_idx


# polarT
def return_polarT(polarT_path, start_date, end_date):
    da = xr.load_dataarray(polarT_path)

    # Assume da is your DataArray (name 'T6090S'), with a daily 'time' coord
    da_aug_nov = da.sel(time=da.time.dt.month.isin([8, 9, 10, 11]))
    da_aug_nov = da_aug_nov.sel(time=slice(start_date, end_date)) 

    ti = da_aug_nov.indexes['time']  # DatetimeIndex
    day_idx = pd.Series(range(ti.size), index=ti).groupby(ti.year).cumcount().to_numpy()
    
    polarT = da_aug_nov.assign_coords(
        year=('time', da_aug_nov.time.dt.year.values),
        day =('time', day_idx),
    ).set_index(time=['year','day']).unstack('time')  # dims now ('year','day')

    return polarT


# --- main to return arrays needed for polarT detect function ---

def prep_input_polarT(Tband_file, polarT_file, fsw_file, istart_date, iend_date):
    """
    Prepare all input variables required for the SSW detection algorithms based on
    polar-cap temperature gradients (e.g., for detect_ssw_tgrad).

    This function loads preprocessed datasets from disk — including the temperature
    band data (for computing ΔT and climatology), the polar-cap temperature, and the
    final stratospheric warming (FSW) date indices — and returns them as ready-to-use
    xarray and numpy objects within a specified time range.

    Parameters
    ----------
    Tband_file : str or Path
        Path to the NetCDF file containing temperature-band data (e.g., T80–90S and T60–70S).
        Used by `output_Tvars()` to compute ΔT, climatological mean, and standard deviation.
    polarT_file : str or Path
        Path to the NetCDF file containing polar-cap (60–90S) mean temperature.
    fsw_file : str or Path
        Path to the NetCDF or .npy file containing final stratospheric warming (FSW) dates.
        Used by `output_fwidx()` to produce per-year index positions for each event.
    istart_date : str
        Start date (inclusive) of the analysis period (e.g., '1979-08-01').
    iend_date : str
        End date (inclusive) of the analysis period (e.g., '2023-11-30').

    Returns
    -------
    tuple
        (T_augnov, cdelT, delT, stddv, polarT, fw_idx)
        where:
            T_augnov : xr.Dataset
                Temperature dataset for the selected time window.
            cdelT : xr.DataArray
                Day-of-year climatological mean of ΔT.
            delT : xr.DataArray
                Daily temperature difference between 80–90S and 60–70S.
            stddv : xr.DataArray
                Day-of-year standard deviation of ΔT.
            polarT : xr.DataArray
                Daily polar-cap (60–90S) mean temperature.
            fw_idx : np.ndarray
                Array of FSW index positions per year for event alignment.

    Notes
    -----
    This function is typically called before SSW detection algorithms such as
    `detect_ssw_tgrad()` or `polarT1_lim()` to assemble all necessary temperature
    and reference data in memory.
    """

    T_augnov, cdelT, delT, stddv = output_Tvars(Tband_file, start_date=istart_date, end_date=iend_date)
    polarT = return_polarT(polarT_file, istart_date, iend_date)
    fw_idx = output_fwidx(fsw_file, T_augnov)

    return T_augnov, cdelT, delT, stddv, polarT, fw_idx
    