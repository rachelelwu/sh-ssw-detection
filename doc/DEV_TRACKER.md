#### Project Structure

The python functions are stored in `src/sh_ssw_methods`, while `src/test_main.ipynb` is the script now to call the python functions.

```text
src/
├── test_main.ipynb          # Entry point for testing and validation
└── sh_ssw_methods/          # Main package containing ssw identification functions
    ├── __init__.py          # Exposes functions for easy import
    ├── .py         
    └── .py


```
#### Available Functions

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
* **Function**: `detect_ozone_threshold(tcO3:xr.DataArray)`
* **Method**: Detects events based on criterion on ozone threshold (`thres_du`)
* **Data**: `data/tco3_MERRA2.nc`
* **Legacy Source**: `Ozone_Butler/sh_vortex_ozone_date.ipynb`



| Original Function | Author | Python Function | Status | Comments |
| :--- | :--- | :--- | :--- | :--- |
| FZ100 | thomas |  | 🛑 Todo |  |
| STmode | eun-pa |  | ✅ Done |  |
| polarT1 | eun-pa | | ✅ Done |  |
| zeof1 | eun-pa | zeof.py | ⚠️ WIP | function translation attempt and problems described in src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb|
| ozone | amy| | ✅ Done  |  |
| Uanom |  | | ✅ Done |  |
| Utend |  | | ✅ Done |  |
| Udrop | rachel | | ✅ Done |  |
| Split | martin | | 🛑 Todo  |  |
| Displace | martin | |🛑 Todo |  |
