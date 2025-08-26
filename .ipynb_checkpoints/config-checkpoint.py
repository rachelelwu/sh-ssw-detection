# config.py
data_paths = {
    "zeof": "data/zeof/",
    "tco3": "data/tco3_MERRA2.nc",
    "FN": "data/f2bym_100_20_90_SH_ERA5_1979_2023.bin",
    "Tmid": "data/Tmid_polar_10hPa_daily_1959_2023_era5.nc",
    "fw_dates": "data/FSW_dates_SH_60S_10hPa_era5_1959_2023.nc",
    "polarT": "data/polarT_anom_10hPa_daily_1979_2021_era5.nc",
    "pc1file": "data/stcmI.1979-2023.txt",
    # fully resolved u-files:
    "u10": "/net/cfc/s2s_nobackup/rachwu/ERA5/u1060s_era5_1959_2023.nc",
    "u50": "/net/cfc/s2s_nobackup/rachwu/ERA5/u5060s_era5_1959_2023.nc",
}

params = {
    # zeof1 for multiple levels
    "zeof1_lim_lv50": {"lv": 50, "basepath": data_paths["zeof"]},
    "zeof1_lim_lv10": {"lv": 10, "basepath": data_paths["zeof"]},
    "zeof1_lim_lv1":  {"lv": 1,  "basepath": data_paths["zeof"]},

    # ozone
    "ozone_butler": {"tco3_file": data_paths["tco3"], "thresh": None, "output_csv": True},

    # u_anom (10 & 50)
    "u_anom_karpechko_10": {
        "u1060f": data_paths["u10"], "wp0": "10",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'], "output_csv": True,
    },
    "u_anom_karpechko_50": {
        "u1060f": data_paths["u50"], "wp0": "50",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'], "output_csv": True,
    },

    # u_tend (10 & 50)
    "u_tend_karpechko_10": {
        "u1060f": data_paths["u10"], "wp0": "10",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'], "output_csv": True,
    },
    "u_tend_karpechko_50": {
        "u1060f": data_paths["u50"], "wp0": "50",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'], "output_csv": True,
    },

    # polarT (Shen)
    "polarT_shen": {
        "Tmid_file": data_paths["Tmid"], "fw_dates_file": data_paths["fw_dates"],
        "polarT_file": data_paths["polarT"], "write_to_file": True,
    },

    # EP flux (Reichler)
    "epflux_reichler": {"FN": data_paths["FN"], "DT": None, "DU": None, "output_csv": True},

    # u decel (Wu) 10 hPa example
    "u_decel_wu_10": {"u1060f": data_paths["u10"], "wp0": "10", "drop": None, # if None, drop = -2 m/s per 10 days
    "clim_range":['1979-01-01','2020-12-31']},

    # STMI
    "stmi_lim": {"inFile": data_paths["pc1file"], 
    "thresh": None, # if None, default = 0.8
    "save_txt":True,
    "save_mask":True},
}
