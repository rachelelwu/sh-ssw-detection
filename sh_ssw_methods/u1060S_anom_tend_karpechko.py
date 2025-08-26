import pandas as pd
import numpy as np

from datetime import date, timedelta, datetime
import xarray as xr
import csv



def load_u1060S_with_anoms(nc_path, var='u1060S', base_start='1979-01-01', base_end='2023-12-31'):
    """
    take in .nc 6hourly era5 file and return dataframes of u_full and u_anom
    """
    # Open dataset
    ds = xr.open_dataset(nc_path)

    if var not in ds:
        raise ValueError(f"Variable '{var}' not found in {nc_path}")

    # Drop Feb 29
    is_feb29 = (ds.time.dt.month == 2) & (ds.time.dt.day == 29)
    ds = ds.sel(time=~is_feb29)

    # for now only take 1979-2021 to be consistent
    ds = ds.sel(time=slice(base_start, base_end))

    # Convert to DataFrame for easier month/day grouping
    df = ds[var].to_dataframe().reset_index()

    # Convert from 6-hourly to daily mean
    ds_daily = ds[var].resample(time='1D').mean()

    # Convert to DataFrame for detection
    df = ds_daily.to_dataframe().reset_index()
    df['year'] = df['time'].dt.year
    df['month'] = df['time'].dt.month
    df['day'] = df['time'].dt.day

    # Compute climatology over base period only (month/day groups)
    clim = (df.groupby(['month', 'day'], as_index=False)[var]
              .mean()
              .rename(columns={var: 'clim'}))

    # Merge climatology back to all years
    anom_df = df.merge(clim, on=['month', 'day'], how='left')

    # Compute anomaly
    anom_df[var] = anom_df[var] - anom_df['clim']
    anom_df.drop(columns='clim', inplace=True)

    # Full field DataFrame (no anomalies subtracted)
    data_df = df.copy()

    return data_df, anom_df

def iden_ssw_uanom(data_df, anom_df, wp0, thres, clim_start, clim_end, output_csv):
    dataset='era5'
    ndd=len(anom_df)
    
    u10a=anom_df['u%s60S' % wp0].to_numpy()
    u10=data_df['u%s60S' % wp0].to_numpy()
    yrs=anom_df['year'].to_numpy()
    mns=anom_df['month'].to_numpy()
    dds=anom_df['day'].to_numpy()
    
    csv_header = ['year','month','day','u1060S']
    
    outfile=f'u{wp0}60S_{dataset}_anom.csv'
    with open(outfile, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(csv_header)
        for i in range(ndd):
            row=[yrs[i],mns[i],dds[i],u10a[i]]
            writer.writerow(row)
    

    # -20 m/s for 10 hPa
    # -11 m/s for 50 hPa
    # if (wp0=='10'):
    #     thres=-20
    # if (wp0=='50'):
    #     thres=-11
    wthres=str(thres)
    
    
    # if (dataset=='jra55'):
    #     wperiod='1958-2021'
    # if (dataset=='era5'):
    #     wperiod='1979-2021'

        
    nssw=0
    dates_list = []  # to store datetime objects
        

    for i in range(ndd - 1):
        if (u10a[i+1] < thres) & (u10a[i] >= thres) & (mns[i+1] > 4) & (mns[i+1] < 11):
            min1 = np.amin(u10a[i-20:i])
            min2 = np.amin(u10[i:i+10])

            if (min1 > thres) & (min2 > 0):
                date_obj = datetime(yrs[i+1], mns[i+1], dds[i+1])
                dates_list.append(date_obj)
                nssw=nssw+1
                
    # convert to numpy datetime64 array
    dates_array = np.array(dates_list, dtype='datetime64[D]')
    print(dates_array)   
    
    wperiod = '%s-%s' % (clim_start[:4], clim_end[:4])

    if output_csv: 
        outfile=f'ssw_u{wp0}60S_{dataset}_anom.csv'
        
        with open(outfile, "w", newline="") as file:
            file.write(f"# data: {dataset} {wperiod}, daily means from 00,06,12,18 UTC instantaneous.\n")
            file.write(f"# definition: U_anomaly < {wthres} at {wp0} hPa and -60S. At least 10 days "
                       "of positive U_full required after SSW date to separate from the final warming.\n")
            file.write("# two events are the same if occurring within 20 days.\n")
            file.write(f"# U_anomaly defined as departure from the climatology over {wperiod} for the calendar day.\n")
    
            writer = csv.writer(file)
            writer.writerow(["year", "month", "day"])
    
            for date in dates_array:
                y, m, d = str(date).split("-")
                writer.writerow([y, m, d])                                       
                    
    # if (dataset=='jra55'):
    #     nyrs=2021.-1958.+1.
    if (dataset=='era5'):
        nyrs=int(clim_end[:4])-int(clim_start[:4])+1.
                    
    freq=nssw/nyrs
    print(nssw,10.*freq,nyrs)

    return dates_array


def iden_ssw_utend(data_df, anom_df, wp0, thres, clim_start, clim_end, output_csv):

    dataset='era5'
    
    ndd=len(anom_df)
    
    u10a=anom_df['u%s60S' % wp0].to_numpy()
    u10=data_df['u%s60S' % wp0].to_numpy()
    yrs=anom_df['year'].to_numpy()
    mns=anom_df['month'].to_numpy()
    dds=anom_df['day'].to_numpy()
    
    csv_header = ['year','month','day','u1060S']
        
    nssw=0
    # -35m/s for 10hPa
    # -19.m/s for 50hPa
    # if (wp0=='10'):
    #     thres=-35
    # if (wp0=='50'):
    #     thres=-19
    wthres=str(thres)
    rng=7
        
    #Calculate delta
    delta=np.zeros(ndd)
    for i in range(rng,ndd-1-rng):
        delta[i]=u10[i+rng]-u10[i-rng] 

    ndd = len(delta)
    dates_list = []

    for i in range(rng, ndd - 1 - rng):
        if (delta[i] < thres) & (delta[i-1] >= thres) & (4 < mns[i] < 11):
            min1 = np.amin(delta[i-20:i])
            min2 = np.amin(u10[i:i+10])

            if (min1 > thres) & (min2 > 0):
                date_obj = datetime(yrs[i+1], mns[i+1], dds[i+1])
                dates_list.append(date_obj)
                nssw=nssw+1
    
    # convert to numpy datetime64 array
    dates_array = np.array(dates_list, dtype='datetime64[D]')
    print(dates_array)

    

    # if (dataset=='jra55'):
    #     wperiod='1958-2021'
    # if (dataset=='era5'):
    #     wperiod='1979-2021'

    wperiod = '%s-%s' % (clim_start[:4], clim_end[:4])
    if output_csv: 
        outfile=f'ssw_u{wp0}60S_{dataset}_tend.csv'
        with open(outfile, "w", newline="") as file:
            file.write(f"# data: {dataset} {wperiod}, daily means from 00,06,12,18 UTC instantaneous.\n")
            file.write(f"# definition: U_tendency < {wthres} at {wp0} hPa and -60S. At least 10 days "
                       "of positive U_full required after SSW date to separate from the final warming.\n")
            file.write("# two events are the same if occurring within 20 days.\n")
            file.write("# U_tendency defined as the U_full difference between +/-7 days from the current day.\n")
    
            writer = csv.writer(file)
            writer.writerow(["year", "month", "day"])
    
            for date in dates_array:
                y, m, d = str(date).split("-")
                writer.writerow([y, m, d])



    # if (dataset=='jra55'):
    #     nyrs=2021.-1958.+1.
    if (dataset=='era5'):
        nyrs=int(clim_end[:4])-int(clim_start[:4])+1.
            
    freq=nssw/nyrs
    print(nssw,10.*freq,nyrs)

    return dates_array

def main_uanom(u1060f, wp0, thres=None, clim_range=['1979-01-01','2023-12-31'], output_csv=True):
    """
    Detect Southern Hemisphere SSW events based on zonal mean zonal wind at 10hPa/50hPa using anomalies.

    Parameters
    ----------
    u1060f : str
        Path to the zonal mean zonal wind file (now for ERA5).
    wp0 : str
        Pressure level to detect events, either '10' or '50' hPa.
    thres : float, optional
        Threshold value. If None, defaults are used:
        -20 for 10hPa, -11 for 50hPa.
    clim_range : list, optional
        Range of dates to compute climatologies and identify events.
    output_csv : bool, optional
        Whether to output results as CSV.

    Returns
    -------
    pandas.DatetimeIndex or xarray.Dataset
        Event dates if `output_csv=True`, otherwise full processed dataset.
    """

    # Set default threshold if not provided
    if thres is None:
        if str(wp0) == '10':
            thres = -20
        elif str(wp0) == '50':
            thres = -11
        else:
            raise ValueError(f"Unsupported wp0 level: {wp0}. Expected '10' or '50'.")

    clim_start, clim_end = clim_range

    data_df, anom_df = load_u1060S_with_anoms(
        u1060f, var=f'u{wp0}60S', base_start=clim_start, base_end=clim_end
    )

    event_dates_anom = iden_ssw_uanom(
        data_df, anom_df, wp0, thres, clim_start, clim_end, output_csv
    )

    return event_dates_anom



def main_utend(u1060f, wp0, thres=None, clim_range=['1979-01-01','2023-12-31'], output_csv=True):
    """
    Detect Southern Hemisphere SSW events based on zonal mean zonal wind at 10hPa/50hPa using tendencies. 

    Parameters
    ----------
    u1060f : str
        Path to the zonal mean zonal wind file (ERA5).
    wp0 : str or int
        Pressure level to detect events, either 10 or 50 hPa. 
    thres : float, optional
        Threshold value. If None, defaults are used:
        -35 for 10 hPa
        -19 for 50 hPa
    clim_range : list, optional
        Range of dates to compute climatologies and identify events.
    output_csv : bool, optional
        Whether to output results as CSV.

    Returns
    -------
    pandas.DatetimeIndex or xarray.Dataset
        Event dates if `output_csv=True`, otherwise full processed dataset.

    Notes
    -----
    The detection algorithm proceeds as follows:
    - tendency in zonal mean zonal wind < -35 m/s for 10 hPa  
    - tendency in zonal mean zonal wind < -19 m/s for 50 hPa  
    """

    # Handle default thresholds
    if thres is None:
        if str(wp0) == '10':
            thres = -35
        elif str(wp0) == '50':
            thres = -19
        else:
            raise ValueError(f"Unsupported wp0 level: {wp0}. Expected '10' or '50'.")

    clim_start, clim_end = clim_range

    data_df, anom_df = load_u1060S_with_anoms(
        u1060f, var=f'u{wp0}60S', base_start=clim_start, base_end=clim_end
    )

    event_dates_tend = iden_ssw_utend(
        data_df, anom_df, wp0, thres, clim_start, clim_end, output_csv
    )

    return event_dates_tend





# def main_uanom_tend_karpechko(u1060f, wp0, thres, clim_start='1979-01-01', clim_end='2023-12-31'):
#     """
#     Detect Southern Hemisphere SSW events based on zonal mean zonal wind at 10hPa/ 50hPa using anomalies and tendencies

#     """

#     data_df, anom_df = load_u1060S_with_anoms(u1060f, var='u%s60S'%wp0, base_start=clim_start, base_end=clim_end)

#     event_dates_anom = iden_ssw_uanom(data_df, anom_df, wp0, thres, clim_start, clim_end)

#     event_dates_tend = iden_ssw_utend(data_df, anom_df, wp0, thres, clim_start, clim_end)


#     return event_dates_anom, event_dates_tend
