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

**`zeof.py`**
* `detect_pc1_events(basepath="sh_ssw_methods/processed/")`: identify events from pc1 that are already computed and outputed as daily.pc1.*.txt
    * *original script:* `ZEOF1_Lim/event_selection.ncl`
* ⚠️ to modify: to directly take z anomalies to output pc1, and only then uses `detect_pc1_events` to identify events
    * *original script:*
    * `ZEOF1_Lim/Z.month.EOF.ncl` - to identify EOF1 from monthly z data
    * `ZEOF1_Lim/daily.pc1.ncl` - to project EOF1 identified from monthly data to daily z anomalies to get daily pc1  
    * working .ipynb on translating into Python: `src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb`
      
**`ozone_threshold.py`**
* `detect_ozone_threshold(tcO3:xr.DataArray)`: detect events based on criterion on ozone threshold (thres_du)
    * *original script:* `Ozone_Butler/sh_vortex_ozone_date.ipynb`

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
