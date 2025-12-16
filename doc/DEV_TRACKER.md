## Project Structure

The python functions are stored in `src/sh_ssw_methods`, while `src/test_main.ipynb` serves as the primary validation sandbox.

```text
src/
├── test_main.ipynb          # Entry point for testing and validation
└── sh_ssw_methods/          # Main package containing NCL-ported logic
    ├── __init__.py          # Exposes functions for easy import
    ├── vorticity.py         # Specific calculation modules
    └── thermodynamics.py

```
| Original Function | Author | Python Function | Status | Comments |
| :--- | :--- | :--- | :--- | :--- |
| FZ100 | thomas |  | 🛑 Todo |  |
| STmode | eun-pa |  | ✅ Done |  |
| polarT1 | eun-pa | | ✅ Done |  |
| zeof1 | eun-pa | | ⚠️ WIP | function translation attempt and problems described in src/sh_ssw_methods/zeof_dev_attempt_dec2025.ipynb|
| ozone | amy| | ✅ Done  |  |
| Uanom |  | | ✅ Done |  |
| Utend |  | | ✅ Done |  |
| Udrop | rachel | | ✅ Done |  |
| Split | martin | | 🛑 Todo  |  |
| Displace | martin | |🛑 Todo |  |
