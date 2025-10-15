#!/usr/bin/env python
# coding: utf-8
#
# This code computes 2D polar vortex moments following Mitchell et al (2011) and Seviour et al (2013)

import xarray as xr
import numpy as np


def ComputeVortexMoments(da,edge=30.5,time='time',verbose=0):
    '''Computes Aspect Ration and Centroid Latitude from geopotential height.
        It is important that latitude is sorted from south to north, not north to south as in ERA5.

       INPUTS:
          da  : xr.DataArray containing 2D geopotential height at a given pressure level
          edge: Height of vortex edge in km.
          time: name of time dimension in da
          verbose: print progress if >0
     
       OUTPUTS:
          AR,CL: aspect ration and centroid latitude
    '''
    from vortex_moments import vor
    from aostools import climate as ac

    # make sure lat is from south to north, and lon from 0 to 360
    #  this also renames latitude and longitude from whatever name to lat and lon
    da = ac.StandardGrid(da,rename=True)
    da = da.transpose(time,'lat','lon')
    z = da.values
    lons = da['lon'].values
    lats = da['lat'].values
    times= da[time].values
    nt = len(times)
    aspects = np.zeros(nt,)
    latc = np.zeros_like(aspects)
    # compute moments for each timestep
    for t in range(nt):
        if verbose > 0:
            ac.update_progress(t/nt)
        moms = vor.calc_moments(z[t,:],lats,lons,hemisphere='SH',field_type='GPH',edge=edge*1000)
        aspects[t] = moms['aspect_ratio']
        latc[t] = moms['centroid_latitude']
    if verbose > 0:
        ac.update_progress(1)

    # put it all together into an xr.Dataset
    aspx = xr.DataArray(aspects,coords=[da[time]],name='aspect_ratio')
    latx = xr.DataArray(latc,coords=[da[time]],name='centroid_latitude')

    return xr.merge([aspx,latx])

def DetectEvents(da,ar_cl,thresh,period,time='time',edge=30.5,distinct=20,is_mom=False,verbose=0):
    '''Detect splits or displacement events from geopotential height. We assume daily data.

      INPUTS:
        da :  daily geopotential height as a function of time,lon,lat, i.e. at a specific level.
              OR, daily vortex moment if is_mom=True.
        ar_cl: which method to use: 
                 'ar' for Aspect Ratio
                 'cl' for Centroid Latitude. 
        thresh:  detection threshold
        period:  how many timesteps (==days) beyond threshold?
        time:    name of time dimension
        edge:    polar vortex edge in km
        distinct:number of days between events to count as independent event
        is_mom:  whether da is actually the wanted vortex moment rather than geopotential height.
                  this is mostly useful for debugging or checking different thresholds
        verbose: add a progress bar for moment calculations
      OUTPUTS:
    '''
    if not is_mom:
        # get timeseries of aspect ratio and centroid latitude
        asp = ComputeVortexMoments(da,edge=edge,time='time',verbose=verbose)
        # get the right field
        if ar_cl.lower() == 'ar':
            da = asp['aspect_ratio']
        elif ar_cl.lower() == 'cl':
            # Note that centroid latitude is positive even if it's for the Southern Hemisphere
            #  we want it to be negative so we can search for maxima
            da = -asp['centroid_latitude']
        else:
            raise ValueError('ar_cl has to be "ar" or "cl" but got '+ar_cl)
    
    # here comes the detection:
    #  THIS IS THE SAME FOR ALL METHODS, SO SHOULD PROBABLY BE DONE IN A SEPARATE, COMMON FUNCTION

    #  if we want N days *above* a value,
    #  that means we want the *minimum* of a rolling N days to be *above* that value
    events = da.rolling(time=period).min() > thresh
    # da.isel(time=events) now gives all instances when this is true.
    #  but, a single event going for longer than period will result in a sequence of True.
    #  we only want the first day when this happens
    #  we find this by looking for the first time the index goes from 0 (False) to 1 (True)
    starts = np.insert(np.diff(events.astype(int)) > 0,0,False)
    # similarly, to find the end of each events, we do the inverse as it goes from 1 to 0
    ends   = np.insert(np.diff(events.astype(int)) < 0,0,False)
    # the times when events happen are
    events = da.isel({time:starts})[time]
    event_ends = da.isel({time:ends})[time]
    # but this is after the rolling mean, so the start of the event is the rolling period before that
    events = events - np.timedelta64(period,'D')
    # now we also need to remove events which happen within `distinct` days of an earlier event
    #  this is the same as retaining the ones which have a larger delta-t to the next event
    #    and then shifting one to the right
    #  note we are comparing the end of the previous with the start of the subsequent event
    # this only works if the data does not end in the middle of an event
    day_diffs = np.insert((events.values[1:]-event_ends.values[:-1])/np.timedelta64(1,'D') > distinct, 0,True)
    #day_diffs = np.insert(np.diff(events)/np.timedelta64(1,'D') > distinct,0,True)
    # return the final list of onset dates
    return events.isel({time:day_diffs})

def WriteCSV(onset_dates,init_text,filename):
    '''
    Write CSV files with all onset dates. Will write each date on one line following the format year,month,day
    INPUTS:
       onset_dates:  xarrat.DataArray of type time. Defines the onset dates to be written as lines
       init_text:    any file header to be written before the onset dates
       filename:     name of file to be created
    '''
    import os
    with open(filename,'w') as csvfile:
        csvfile.write(init_text)
        for event in onset_dates:
            csvfile.writelines(str(event.dt.strftime('%Y,%m,%d').values)+os.linesep)
    print(filename)



if __name__ == "__main__":
    # this is specific to MJ's machine at UNSW
    datadir = '/srv/ccrc/AtmMJ/shared/ERA5/'
    ds = xr.open_mfdataset(datadir+'ERA5_dm.*.z.nc')
    # from here it's generic - note ERA5 z is geopotential, not geopotential height
    z10 = ds.sel(pressure_leve=10).z/9.81
    moms = ComputeVortexMoments(z10,edge=30.5,time='time',verbose=1)
    events = DetectEvents(moms.aspect_ratio,'AR',1.62,7)
    init_text = "Vortex moment definition: geopotential height at 10 hPa, vortex edge = 30.5 km\n # individual 7-day periods of aspect_ratio above 1.62. Events are considered the same if spaced by less than 20 days\n"
    WriteCSV(events,init_text,'Split.csv')
    events = DetectEvents(moms.centroid_latitude,'CL',-72.3,7)
    init_text = "# Vortex moment definition: geopotential height at 10 hPa, vortex edge = 30.5 km\n# individual 7-day periods of centroid_latitude below -72.3. Events are considered the same if spaced by less than 20 days\n"
    WriteCSV(events,init_text,'Displace.csv')
