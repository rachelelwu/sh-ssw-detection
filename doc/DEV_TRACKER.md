## Project Structure

The python functions are stored in `src/sh_ssw_methods`, while `src/test_main.ipynb` is the script now to call the python functions.

```text
src/
├── test_main.ipynb          # Entry point for testing and validation
└── sh_ssw_methods/          # Main package containing ssw identification functions
    ├── __init__.py          # Exposes functions for easy import
    ├── .py         
    └── .py


## Available Functions

### `vorticity.py`
* **`compute_vorticity(ds)`**: Calculates relative vorticity. Returns a Pandas DataFrame.
    * *NCL original:* `uv2vr_cfd`
* **`compute_circulation(ds)`**: Computes circulation over a specific region.
    * *NCL original:* Custom script

### `thermodynamics.py`
* **`calc_potential_temp(ds)`**: Converts temperature to theta.
    * *NCL original:* `pot_temp`

```
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
