import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
# ================= Ozone_Butler =====================

def sh_vortex_ozone_dates(tcO3_filepath, thresh, output_csv=None):
    """
    Define SH SSW based on ozone concentration. Original script: Butler.
    (rachel wu: only put together the lines of scripts into a single function)

    tcO3_filepath (str) - filepath for the tco3 file, from MERRA2
    """
    o = xr.open_dataset(tcO3_filepath)
    
    # The SSWC has odd time units (Julian days) which Xarray doesn't recognize, so here I replace the time dimension.
    # Need to make sure that these are the correct dates for a given dataset.
    
    oldtime = o.time
    newtime_beg = pd.to_datetime("1979-07-01")
    newtime = pd.date_range(newtime_beg, periods=np.shape(oldtime)[0], freq='D')
    o = o['tcO3'].assign_coords(time=newtime)


    # take polar-cap average, cosine-weighted 
    otot = o.sel(time=slice('1980','2020'),lat=slice(-90,-65))#.mean(dim="lon")
    
    weights = np.cos(np.deg2rad(otot.lat))
    weights.name = "weights"
    obs_weighted = otot.weighted(weights)
    sh_o3 = obs_weighted.mean("lat")

    #remove seasonal cycle (could get more fancy with this, right now it's just removing long-term mean of each day)
    anom= sh_o3.groupby("time.dayofyear") - sh_o3.groupby("time.dayofyear").mean("time")

    # Because of ozone depletion/recovery, there is a trend particularly from 1980-2000.
    # This means choosing anomalies above some positive threshold may select more events in the 1980s, for example.
    # To try to account for this, I remove the mean of each year from each year of data
    
    # old = xr.DataArray(dims=["time"],coords={"time": [np.datetime64('1979-12-31')]})
    old = xr.DataArray(dims=["time"], coords={"time": [np.datetime64("1979-12-31", "ns")]})

    for y in range(1980,2021,1):
        newd = anom.sel(time=str(y)).drop('dayofyear')
        anompol = newd - newd.mean("time") #mean of each year removed
        fin = xr.concat([old,anompol],dim="time")
        old = fin.copy()
    
    fin = fin[1:]

    # Find anomalies exceeding some positive threshold (weak vortex events)
    ind = fin.where(fin > thresh, drop=True)

    # Convert to pandas to deal with dates
    indp = ind['time'].to_dataframe()

    # Take difference between dates to determine which events are consecutive
    # Dates that are consecutive make up an "event" 
    indp['event']=indp.time.diff().dt.days.ne(1).cumsum()

    # Find first and last date of each "event". "count" shows how many days each event lasts
    data = indp.groupby(['event'])['time'].agg(['first','last','count'])

    # Only retain events that last at least 3 days
    persist = data.loc[data['count'].ge(3)].copy()

    # Now determine which of remaining dates are still too close together to be considered "independent" events
    # cevent is shifted so that the difference between the last date of event #1 and the first date of event #2 can be compared ("consec")
    cevent = persist.copy()
    cevent['first'] = cevent['first'].shift(periods=-1)
    cevent['consec'] = cevent['first'].sub(cevent['last'], axis = 0)

    
    # Now put "consec" info back into non-shifted array and search for at least 20 day separation
    persist.loc[:,('consec')] = cevent['consec'].dt.days.gt(20) # will return "True" where "consec" is greater than 20 days
    
    persist["consec"] = persist["consec"].astype("boolean")
    persist.loc[:,('consec')] = persist.loc[:,('consec')].shift()
    
    persist = persist.fillna(value=True) #First value will be NaN but we want to retain that as an event so set to True

    # Finally pull out the events that are "independent", as marked by "True"
    fevent = persist.loc[persist['consec']==True]

    final = fevent.rename(columns={"first": "date"}).drop(columns=['last','count','consec'])

    final['Year'] = final.date.dt.year
    final['Month'] = final.date.dt.month
    final['Day'] = final.date.dt.day
    final = final.drop(columns="date")
    
    # after you did: final = final.drop(columns="date")
    event_dates = pd.to_datetime(
        dict(year=final['Year'], month=final['Month'], day=final['Day']),
        errors='coerce'
    ).values.astype('datetime64[D]')


    if output_csv:

        # if output_csv is True, write in current folder
        if output_csv is True:
            out_dir = Path(".")
        else:
            out_dir = Path(output_csv)
            out_dir.mkdir(parents=True, exist_ok=True)
            
        filepath = out_dir / f"sh_ssw_totcol_o3_{thresh}DU.csv"
 
        with open(filepath, 'w') as fout:
              fout.write('# Description\n')
              fout.write('# total column ozone 65-90S anomalies > '+str(thresh)+' DU\n')
              fout.write('# Must persist for at least 3 days to be event\n')
              fout.write('# Events must be separated by >20 days\n')
              fout.write('# MERRA2\n')
              final.to_csv(fout,index=False)

              print("%s saved." % f"sh_ssw_totcol_o3_{thresh}DU.csv")

    return final, event_dates

def main_ozone_butler(tco3_file, thresh=None, output_csv=None):
    """
    Detect Southern Hemisphere weak vortex (ozone-based SSW) events.

    Parameters
    ----------
    tcO3_filepath : str
        Path to the NetCDF file containing total column ozone (tco3).
    thresh : float, optional
        Positive anomaly threshold in Dobson Units (DU) used to define
        weak vortex events. Default is 40 DU.
    output_dates : bool, optional
        If True, return only event dates. If False, return the full 
        processed dataset along with event dates. Default is True.

    Returns
    -------
    np.ndarray of np.datetime64[D]
        Array of detected event dates, representing independent 
        Southern Hemisphere weak vortex events.

    Notes
    -----
    The detection algorithm proceeds as follows:

    - Polar cap ozone anomalies (65°S–90°S) are used to detect weak vortex events.
    - Events are identified when ozone anomalies exceed a positive threshold (default 40 DU).
    - Events need to persist for at least 3 days 


    References/ Author
    ----------
    author: Butler 

    """
    if thresh is None:
        thresh = 40

    # print(thresh)
    final, event_dates = sh_vortex_ozone_dates(tco3_file, thresh, output_csv)

    return event_dates



