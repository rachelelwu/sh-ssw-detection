import numpy as np
import pandas as pd
import glob

import xarray as xr

from pathlib import Path


# read inputs
years = np.arange(1979, 2021+1,1)

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
    fw_idx = fwmask.argmax(dim='day')  # same as j in NCL

    return fw_idx


# polarT
def return_polarT(polarT_path):
    da = xr.load_dataarray(polarT_path)

    # Assume da is your DataArray (name 'T6090S'), with a daily 'time' coord
    da_aug_nov = da.sel(time=da.time.dt.month.isin([8, 9, 10, 11]))
    da_aug_nov = da_aug_nov.sel(time=slice("1979-01-01", "2021-12-31")) 

    ti = da_aug_nov.indexes['time']  # DatetimeIndex
    day_idx = pd.Series(range(ti.size), index=ti).groupby(ti.year).cumcount().to_numpy()
    
    polarT = da_aug_nov.assign_coords(
        year=('time', da_aug_nov.time.dt.year.values),
        day =('time', day_idx),
    ).set_index(time=['year','day']).unstack('time')  # dims now ('year','day')

    return polarT

# iden algorithm
def polarT1_lim(T_augnov, fw_idx, delT, cdelT, stddv, polarT):
    
    xtime1 = T_augnov['date']
    
    count = np.full((len(years), 122), -999, dtype=int)
    sel_dates = []
    
    for i, yr in enumerate(years):
        if fw_idx[i] < 0:
            continue
    
        # Loop j = 0 .. eidx-20  (avoid last 20 days before FW)
        jmax = int(fw_idx[i]) - 20
        jmax = min(jmax, delT.shape[1] - 5)   # prevent partial windows
        if jmax < 0:
            continue
    
        for j in range(0, jmax + 1):
            # 1) sign reversal (positive) & 5-day persistence
            if delT[i, j].item() > 0:
                win = delT[i, j:j+5].values  # length 5
                if win.shape[0] != 5:
                    continue                  # safety, but jmax should prevent this
    
                if np.all(win >= 0):
                    # print('cond1: %s' % (win))
                    maxidx = int(np.argmax(win))   # index in the 5-day window
                    midx   = j + maxidx            # absolute time index within season
    
                    # 2) anomaly vs climatology stddev
                    anom = (delT[i, midx] - cdelT[midx]).item()
                    
                    if anom >= float(stddv[midx]):
                        # print(xtime1[i,j].values)
                        # print('cond2: anom=%s, std=%s' % (anom, stddv[midx].values))
                        # print('polarT=%s' % polarT[i, j].item())
                        # 3) polar cap temperature positive at the window start j
                        if polarT[i, j].item() > 0:
                            
                            # print('cond3: polarT=%s' % polarT[i, j].item())
                            count[i, j] = 1
                            sel_dates.append(xtime1[i,j].values)

    # make sure events are 60 days apart
    starts = np.r_[True, np.diff(sel_dates) >= np.timedelta64(60, 'D')]
    event_dates = np.array(sel_dates)[starts]
    
    return event_dates
                        


def polarT2_lim(T_augnov, fw_idx, delT, cdelT, stddv, polarT, pers=5):
    
    xtime1 = T_augnov['date']
    
    count = np.full((len(years), 122), -999, dtype=int)
    sel_dates = []
    
    for i, yr in enumerate(years):
        if fw_idx[i] < 0:
            continue
    
        # Loop j = 0 .. eidx-20  (avoid last 20 days before FW)
        jmax = int(fw_idx[i]) - 20
        jmax = min(jmax, delT.shape[1] - 5)   # prevent partial windows
        if jmax < 0:
            continue
    
        for j in range(0, jmax + 1):
            # 1) sign reversal (positive) & 5-day persistence
            if delT[i, j].item() > 0:
                w_delT = delT[i, j:j+pers].values  # length 5
                if w_delT.shape[0] != 5:
                    continue                  # safety, but jmax should prevent this
    
                if np.all(w_delT >= 0):
                    print('cond1: %s' % (w_delT))

                    # 2) anomaly ≥ 1σ for all 5 days
                    w_c   = cdelT[j:j+pers].values
                    w_std = stddv[j:j+pers].values
                    
                    anom = (w_delT - w_c) / w_std
                    
                    if not np.all(anom >= 1):
                        continue

                    print(xtime1[i,j].values)
                    print('cond2: anom=%s, std=%s' % (anom, w_std))
                    print('polarT=%s' % polarT[i, j].item())
                    
                    # 3) polar cap temperature positive at the window start j
                    w_polar = polarT[i, j + pers].values
                    if not (np.all(np.isfinite(w_polar)) and np.all(w_polar >= 0)):
                        continue    
                    print('cond3: polarT=%s' % polarT[i, j].item())
                    count[i, j] = 1
                    sel_dates.append(xtime1[i,j].values)

    # make sure events are 60 days apart
    starts = np.r_[True, np.diff(sel_dates) >= np.timedelta64(60, 'D')]
    event_dates = np.array(sel_dates)[starts]
    
    return event_dates


def main_polarT1_lim(Tmid_file, fw_dates_file, polarT_file, write_to_file=True):

    T_augnov, cdelT, delT, stddv = output_Tvars(Tmid_file)
    fw_idx = output_fwidx(fw_dates_file, T_augnov)

    polarT = return_polarT(polarT_file)

    event_dates = polarT1_lim(T_augnov, fw_idx, delT, cdelT, stddv, polarT)
    
    if write_to_file:
        iso_strs = np.datetime_as_string(event_dates, unit='D')
        
        # Prepare the header text
        header = """# Shen et al. (2022)
        # if T80-90s >= T6070 for 5 days
        #  if the maximum of Tdiff anomalies >= 1 stddev of the corresponding date
        #   if the corresponding polar cap T of the maximum Tdiff >= 0
        #    if the detected event occurs 20 days earlier than the final warming date
        # of the year and it is not detected within 60 days of an earlier event of the
        # year
        #"""
        
        # Write to file
        with open("ShenEtAl2022_org_Tgrad_1979_2021_era5.txt", "w") as f:
            f.write(header + "\n")
            for d in iso_strs:
                f.write(f"{d}\n")

            print("ShenEtAl2022_org_Tgrad_1979_2021_era5.txt saved.")

    return event_dates
    