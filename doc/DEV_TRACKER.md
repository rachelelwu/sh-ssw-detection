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
* `detect_pc1_events`: identify events from pc1 that are already computed and outputed as daily.pc1.*.txt
    * *original function:* `ZEOF1_Lim/event_selection.ncl`
* ⚠️ to modify: to directly take z anomalies to output pc1, and only then uses `detect_pc1_events` to identify events
    * *original functions:*
    * `ZEOF1_Lim/Z.month.EOF.ncl` - to identify EOF1 from monthly z data
    * `ZEOF1_Lim/daily.pc1.ncl` - to project EOF1 identified from monthly data to daily z anomalies to get daily pc1  
    * working .ipynb: `src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb`
      
**`thermodynamics.py`**
* `calc_potential_temp(ds)`: Converts temperature to theta.
    * *NCL original:* `pot_temp`

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
