#### Functions translation tracker: `sh_ssw_methods`

The python functions are stored in `src/sh_ssw_methods`, while `src/test_main.ipynb` is the script now to call the python functions.

```text
src/
├── test_main.ipynb               # Entry point for testing and validation
└── sh_ssw_methods/               # Main package containing ssw identification functions
    ├── __init__.py               # Exposes functions for easy import
    ├── zeof.py                   # (WIP) Zonal EOF analysis
    ├── ozone_threshold.py        # Ozone threshold detection
    ├── [other_modules].py        # (Add specific filenames here as created)
    └── preprocessing/            # Data preparation and formatting utilities
        ├── __init__.py
        └── prepare_polarT_files.py

```
#### 2. Module Status & Details

⚙️ **Preprocessing functions**


`preprocessing/prepare_polarT_files.py`


#### 3. Functions Status

🔴 **Not Started**

* FZ100
* Split
* Displace


🟠 **In Progress (WIP)**

**`zeof.py`**  
* **Method**: Identify evnets based on PC1 computed from daily Z anomalies
* **Current State**: ⚠️ Translation issues
    * Current implemenation (`detect_pc1_events(basepath="sh_ssw_methods/processed/")`) relies on pre-copmuted `daily.pc1.*.txt` files from eun-pa
    * _Need_: Modify to accept Z anomalies directly, compute EOF1 internally, project to PC1, and then detect events
* **Legacy Source**:
    * `ZEOF1_Lim/event_selection.ncl` (event identification from pc1 (translated in python as `detect_pc1_events`))
    * `ZEOF1_Lim/Z.month.EOF.ncl` (monthly EOF1 identification)
    * `ZEOF1_Lim/daily.pc1.ncl` (project monthly EOF1 to daily z anomalies to pc1)
* **Reference Notebook**: `src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb`
    (Contains current translation attempts and error logs)

🟢 **Completed**

**`ozone_threshold.py`** ✅
* **Function**: `detect_ozone_threshold(tcO3:xr.DataArray, thres_du)`
* **Method**: Detects events based on criterion on ozone threshold (`thres_du`)
* **Input Data**: `data/tco3_MERRA2.nc`
* **Legacy Source**: `Ozone_Butler/sh_vortex_ozone_date.ipynb`


**`u_wind_decel.py`** ✅
* **Function**: `detect_u_decel_events(u: xr.DataArray, percentile)`
* **Method**: Detects events based on percentiles (`percentile`) of wind deceleration of 10-day event window
* **Input Data**: `data/u1060s_era5_1959_2023.nc`
* **Legacy Source**: ``


**`u_anom.py`** ✅
* **Function**: `detect_ssw_u_anom(u_daily: xr.DataArray, thres)`
* **Method**: Detects events based on threshold (`thres`) (`thres`=-11 for 50hPa, `thres`=-20 for 10hPa)
* **Input Data**: 50hPa: `data/u5060s_era5_1959_2023.nc`; 10hPa: `data/u1060s_era5_1959_2023.nc`
* **Legacy Source**: `u1060s_Karpechko/u1060S_jra55_analyze.ipynb`


**`u_tend.py`** ✅
* **Function**: `detect_ssw_u_tend(u_daily: xr.DataArray, thres)`
* **Method**: Detects events based on threshold (`thres`) (`thres`=-19 for 50hPa, `thres`=-35 for 10hPa)
* **Input Data**: 50hPa: `data/u5060s_era5_1959_2023.nc`; 10hPa: `data/u1060s_era5_1959_2023.nc`
* **Legacy Source**: `u1060s_Karpechko/u1060S_jra55_analyze.ipynb`


**`polarT.py`** ✅
* **Function**: `detect_ssw_polarT(Tband_file, polarT_file, fsw_file)`
* **Method**: Detect events based on temperature files and final warming dates, temperature and final warming dates can be computed from temp and u files
* **Input Data**: temp: `processed/era5_temp_sh_10hPa_1979_2023.nc`; u: `data/u1060s_era5_1959_2023.nc`
* **Pre-processing**: Uses `output_polarT_prefiles`, `detect_fsw_SH` and output required files for event detection
* **Legacy Source**: `polarT_Lim/ShenEtAl2022.org.Tgrad.ncl` ??? check which function you translated 


**`STmode.py`** ✅
* **Function**: `detect_stmi_events(stcmI_file)`
* **Method**: 
* **Input Data**: monthly u: `processed/era5_uwind_all_lv_monthly_1979_2023.nc`
* **Pre-processing**: Uses `compute_stcmI` and output stcmI indices for event detection
* **Legacy Source**: `STCoupledMode_Lim/S-TcoupledMode.ncl`







| | Original Function | Author | Python Function | Status | Comments |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 1.| FZ100 | thomas |  | 🛑 Todo |  |
| 2.| STmode | eun-pa |  | ✅ Done |  |
| 3.| polarT1 | eun-pa | | ✅ Done |  |
| 4.| zeof1 | eun-pa | zeof.py | ⚠️ WIP | function translation attempt and problems described in src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb|
| 5.| ozone | amy| | ✅ Done  |  |
| 6.| Uanom |  | | ✅ Done |  |
| 7.| Utend |  | | ✅ Done |  |
| 8.| Udrop | rachel | | ✅ Done |  |
| 9.| Split | martin | | 🛑 Todo  |  |
| 10.| Displace | martin | |🛑 Todo |  |
