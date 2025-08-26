# config.py
from pathlib import Path

# Base directory = location of this config file
BASE_DIR = Path(__file__).resolve().parent

# Input data paths
data_paths = {
    "zeof": BASE_DIR / "data" / "zeof",
    "tco3": BASE_DIR / "data" / "tco3_MERRA2.nc",
    "FN": BASE_DIR / "data" / "f2bym_100_20_90_SH_ERA5_1979_2023.bin",
    "Tmid": BASE_DIR / "data" / "Tmid_polar_10hPa_daily_1959_2023_era5.nc",
    "fw_dates": BASE_DIR / "data" / "FSW_dates_SH_60S_10hPa_era5_1959_2023.nc",
    "polarT": BASE_DIR / "data" / "polarT_anom_10hPa_daily_1979_2021_era5.nc",
    "pc1file": BASE_DIR / "data" / "stcmI.1979-2023.txt",

    "u10": BASE_DIR / "data" / "u1060s_era5_1959_2023.nc",
    "u50": BASE_DIR / "data" / "u5060s_era5_1959_2023.nc",
}

# Output paths
output_paths = {
    "event_dates_csv": BASE_DIR / "outputs" / "event_dates_csv"
}

# Make sure output directories exist
for p in output_paths.values():
    p.mkdir(parents=True, exist_ok=True)

# Parameters
params = {
    # zeof1 for multiple levels
    "zeof1_lim_lv50": {"lv": 50, 
            "basepath": data_paths["zeof"],
            "output_csv":output_paths["event_dates_csv"]},
    "zeof1_lim_lv10": {"lv": 10, 
            "basepath": data_paths["zeof"],
            "output_csv":output_paths["event_dates_csv"]},
    "zeof1_lim_lv1":  {"lv": 1,  
            "basepath": data_paths["zeof"],
             "output_csv":output_paths["event_dates_csv"]},

    # ozone
    "ozone_butler": {
        "tco3_file": data_paths["tco3"],
        "thresh": None,
        "output_csv": output_paths["event_dates_csv"]
    },

    # u_anom (10 & 50)
    "u_anom_karpechko_10": {
        "u1060f": data_paths["u10"], "wp0": "10",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'],
        "output_csv": output_paths["event_dates_csv"]
    },
    "u_anom_karpechko_50": {
        "u1060f": data_paths["u50"], "wp0": "50",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'],
        "output_csv": output_paths["event_dates_csv"]
    },

    # u_tend (10 & 50)
    "u_tend_karpechko_10": {
        "u1060f": data_paths["u10"], "wp0": "10",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'],
        "output_csv": output_paths["event_dates_csv"]
    },
    "u_tend_karpechko_50": {
        "u1060f": data_paths["u50"], "wp0": "50",
        "thres": None, "clim_range": ['1979-01-01','2020-12-31'],
        "output_csv": output_paths["event_dates_csv"]
    },

    # polarT (Shen)
    "polarT_shen": {
        "Tmid_file": data_paths["Tmid"],
        "fw_dates_file": data_paths["fw_dates"],
        "polarT_file": data_paths["polarT"],
        "output_csv": output_paths["event_dates_csv"]
    },

    # EP flux (Reichler)
    "epflux_reichler": {
        "FN": data_paths["FN"],
        "DT": None, "DU": None,
        "output_csv": output_paths["event_dates_csv"]
    },

    # u decel (Wu) 10 hPa example
    "u_decel_wu_10": {
        "u1060f": data_paths["u10"], "wp0": "10",
        "drop": None,  # if None, drop = -2 m/s per 10 days
        "clim_range": ['1979-01-01','2020-12-31'],
        "output_csv": output_paths["event_dates_csv"]
    },

    # STMI
    "stmi_lim": {
        "inFile": data_paths["pc1file"],
        "thresh": None,  # if None, default = 0.8
        "save_txt": output_paths["event_dates_csv"],
        "save_mask": output_paths["event_dates_csv"]
    },
}
