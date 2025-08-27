### sh-ssw-detection

Module that contain different defintions to detect sudden stratospheric warmings in the Southern Hemisphere. 

Currently implemented functions include:

- **ZEOF index (Lim)**  
  - `zeof1_lim_lv1`  
  - `zeof1_lim_lv10`  
  - `zeof1_lim_lv50`  
- **EP flux (Reichler)**  
  - `epflux_reichler`  
- **Zonal wind anomalies (Karpechko)**  
  - `u_anom_karpechko_10`  
  - `u_anom_karpechko_50`  
- **Zonal wind tendencies (Karpechko)**  
  - `u_tend_karpechko_10`  
  - `u_tend_karpechko_50`  
- **Polar temperature gradient (Shen)**  
  - `polarT_shen`  
- **Ozone threshold (Butler)**  
  - `ozone_butler`  
- **Wind deceleration (Wu)**  
  - `u_decel_wu_10`  
- … .....

---

#### 📂 Repository structure

```
├── config.py # Central configuration: input/output paths, parameters
├── main.ipynb # Example notebook demonstrating detection methods
├── functions/ # Python modules for event detection
├── data/ # Input data (NetCDF, txt, bin) [not version controlled]
├── outputs/
│ └── event_dates_csv/ # CSV/ txt files with detected event dates
└── README.md
  
```

---

#### 🚀 Quickstart

Clone the repository:

```bash
git clone https://github.com/your-username/sh-ssw-detection.git
cd sh-ssw-detection
```

Edit config.py to point to your local data directories.
Detected events will be saved in outputs/event_dates_csv/.

Example usage (EP flux):
```
import sh_ssw_methods as ssw
from config import params

# zeof1_lim
ze_lv50 = ssw.zeof1_lim(**params["zeof1_lim_lv50"])
```
---

#### 📌 Required packages

`xarray`, `pandas`, `numpy`

`eofs` is included within the package.

---

#### 📌 Notes

For now, events are identified for the period 1979-2020.

---

#### 📜 References
....
