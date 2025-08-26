#!/usr/bin/env python
# coding: utf-8
#
# Module for functions related to deceleration event analysis
# - Algorithm used in Wu et al. (2022) to identify wind acceleration and deceleraiton events
#
# Rachel Wu, rachel.wu@env.ethz.ch

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as dates

import pandas as pd

import xarray as xr

from datetime import datetime, timedelta

# Identify deceleration events

def find_decel(U_data, U_time, drop, nxtevnt, wdw, mode):
    """
    Identify deceleration event given time series of 
    zonal mean zonal wind for a specific level, and
    given criteria
    
    Identification criteria:
        - over the e.g. 10-window there's a deceleration 
        - the deceleration is over the specified value (drop)
        - after an identified event, another event is only identified e.g. 20 days later (nxtevnt)
        - in the 10 days moving window after an event, if the next identified event is stronger 
          than the first identified event, pick the stronger event in the window
        - the ratio between the maximum wind change in the event window
          to the identified wind change for the given event window
          has to be less 1.2, to avoid identifying with variations that
          is much shorter than the event window, e.g. 10 days
    
    Inputs:
        U_data (array)         - time series of u1060
        U_time (datenum array) - datenum of the period of the data
        drop (float)           - wind change threshold over the specified event window (wdw) to exceed such that an event is considered (units: m/s per day) 
        nxtevnt (int)          - number of days after the last identified event to identify the next
        wdw (int)              - window width of events
        mode (str)             - 'decel'/ 'accel' to identify wind deceleration/ acceleration events
    
    Outputs:
        onset_nums (array) - numbers specifying the onset date of the event (convert to date using num2date) 
        drop_vals (array)  - values of the drop in wind, (first - last day) 
    """

    # parameter settings
    step = wdw
    drop_thres = drop*(wdw)
    ratio_thres = 1.2
    
    # output arrays
    onset_nums = np.ones([200])* np.nan
    drop_vals = np.ones([200])* np.nan

    diff = [] 
    ratio = []
    
    data = U_data
    data = data[~np.isnan(data)] # filter out data if nan values present
    
    time_num = U_time
    time_num = time_num[~np.isnan(data)]
    
    # == compute deceleration over 10 day time step - rolling ===
    for i in range(len(U_data)-step): # !!! not rolling the last 10, not sure if this could cause a problem in the future
        
        window = np.roll(data,-i)[:step] # roll data, investigation window: first 10 elements of the rolled data
        window_diff = window[-1] - window[0] # compute difference over 10 time steps
        max_diff = np.abs(np.max(window) - np.min(window))
        signal_ratio = max_diff / np.abs(window_diff)

        diff.append(window_diff)
        ratio.append(signal_ratio)

    # == identify events ==
    
    # add criteria: decrease/ increase cannot be greater than a value
    # check ratio of max diff/ event magnitude
    
    diff = np.array(diff)
    ratio = np.array(ratio)
       
    if mode == 'decel':

        dropind = np.where(diff<=drop_thres)[0] # find where diff less than threshold, i.e. fulfilling conditions
        dropval = diff[diff<=drop_thres]
        dropratio = ratio[diff<=drop_thres]
        
    elif mode == 'accel':
        
        dropind = np.where(diff>=drop_thres)[0] # find where diff less than threshold, i.e. fulfilling conditions
        dropval = diff[diff>=drop_thres]
        dropratio = ratio[diff>=drop_thres]
        
    # only identify event/ save data when the event is e.g. 10 days after the first event
    # so, diff btw dropind is more than 10, ind here is per day
    ind_ref = -nxtevnt 
    lastdrop = 0
    
    count = 0 # number of events count
    for j in range(len(dropind)):
        ind = dropind[j]  # index in data array
        check = ind-ind_ref # check against reference index, i.e. the index of last identified event
        ratio_check = dropratio[j]
        
        if check >= nxtevnt: # if difference in time between this event and the last identified event is more than or equal to 10 days
            
            if ratio_check <= ratio_thres: # criteria 2: ratio less than 1
                
                # save data
                onset_nums[count] = time_num[ind]
                drop_vals[count]  = dropval[j]

                # reset index to identified event index
                ind_ref = ind
                lastdrop = dropval[j]

                count+=1
            
        elif check < nxtevnt: # if event is at less than e.g. 5 days between the last identified event 
            # if deceleration of next scan window is stronger, replace the last identified event with this event
            
            if ratio_check <= ratio_thres: # criteria 2: ratio less than 1
  
                if mode == 'decel':
                    if dropval[j] < lastdrop:
                        lastind = np.where(drop_vals == lastdrop)[0] # find the index where the last event was saved

                        # replace event
                        onset_nums[lastind] = time_num[ind] 
                        drop_vals[lastind] = dropval[j]

                        # reset index
                        ind_ref = ind
                        lastdrop = dropval[j]

                elif mode == 'accel':
                    if dropval[j] > lastdrop:
                        lastind = np.where(drop_vals == lastdrop)[0] # find the index where the last event was saved

                        # replace event
                        onset_nums[lastind] = time_num[ind] 
                        drop_vals[lastind] = dropval[j]

                        # reset index
                        ind_ref = ind
                        lastdrop = dropval[j]

    onset_nums = onset_nums[~np.isnan(onset_nums)]
    drop_vals = drop_vals[~np.isnan(drop_vals)]

    return onset_nums, drop_vals

# def run_find_decel(u1060f, wp0, drop=None, mode='decel', varname=None):
#     """
#     Wrapper to apply find_decel directly on u1060f ERA5 file.
    
#     Parameters
#     ----------
#     u1060f : str
#         Path to ERA5 zonal mean zonal wind NetCDF file.
#     wp0 : str or int
#         Pressure level (10 or 50 hPa).
#     drop, nxtevnt, wdw, mode : see find_decel
#     varname : str, optional
#         Variable name inside NetCDF. If None, assumed 'u{wp0}60S'.
#     """
#     # Set default drop if not given
#     if drop is None:
#         if str(wp0) == '10':
#             drop = -2
#         elif str(wp0) == '50':
#             drop = -1
#         else:
#             raise ValueError(f"Unsupported wp0: {wp0}. Expected '10' or '50'.")

#     nxtevnt = 20 # next event needs to be at least 20 days apart
#     wdw = 10 # window size to detect events
    
#     if varname is None:
#         varname = f"u{wp0}60S"

#     # Load data
#     ds = xr.open_dataset(u1060f)
#     da = ds[varname]

#     # Ensure sorted by time
#     da = da.sortby("time")

#     # Apply time subsetting if requested
#     if clim_range is not None:
#         start, end = pd.to_datetime(clim_range[0]), pd.to_datetime(clim_range[1])
#         da = da.sel(time=slice(start, end))


#     # Extract values
#     U_data = da.values
#     U_time = dates.date2num(da["time"].values)  # convert datetime64 → matplotlib datenums

#     # Call your existing function
#     onset_nums, drop_vals = find_decel(U_data, U_time, drop, nxtevnt, wdw, mode)

#     # Convert onset_nums back to datetime
#     onset_dates = dates.num2date(onset_nums)

#     return onset_dates, drop_vals


def detect_u_decel(u1060f, wp0, drop=None, clim_range=None, nxtevnt=20, wdw=10, mode='decel', varname=None):
    """
    Detect deceleration/acceleration events from zonal mean zonal wind.

    Parameters
    ----------
    u1060f : str
        Path to ERA5 zonal mean zonal wind NetCDF file.
    wp0 : str or int
        Pressure level (10 or 50 hPa).
    drop : float, optional
        Threshold for wind change per day (m/s/day). 
        If None, defaults are:
        -2 for 10 hPa
        -1 for 50 hPa
    nxtevnt : int, optional
        Minimum separation (days) between identified events.
    wdw : int, optional
        Window size in days.
    mode : str, optional
        'decel' or 'accel' (default 'decel').
    varname : str, optional
        Variable name inside NetCDF. If None, assumed 'u{wp0}60S'.

    Returns
    -------
    event_dates : list of datetime
        Onset dates of identified events.
    drop_vals : np.ndarray
        Magnitude of wind changes for events.
    """
    # Set default drop if not given
    if drop is None:
        if str(wp0) == '10':
            drop = -2
        elif str(wp0) == '50':
            drop = -1
        else:
            raise ValueError(f"Unsupported wp0: {wp0}. Expected '10' or '50'.")

    if varname is None:
        varname = f"u{wp0}60S"

    # Load data
    ds = xr.open_dataset(u1060f)
    da = ds[varname].sortby("time")

    # Apply time subsetting if requested
    if clim_range is not None:
        start, end = pd.to_datetime(clim_range[0]), pd.to_datetime(clim_range[1])
        da = da.sel(time=slice(start, end))

    U_data = da.values
    U_time = dates.date2num(da["time"].values)

    # Use your original core finder
    onset_nums, drop_vals = find_decel(U_data, U_time, drop, nxtevnt, wdw, mode)
    event_dates = dates.num2date(onset_nums)

    return event_dates, drop_vals
