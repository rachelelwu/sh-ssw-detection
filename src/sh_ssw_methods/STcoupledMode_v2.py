import xarray as xr
import numpy as np
from .eofs.xarray import Eof

import pandas as pd
import re

from pathlib import Path
# ------------------------------
# Settings
# ------------------------------
hemi = "sh"
hlat, llat = -65, -55   # latitude band
yrs = np.arange(1979, 2023+1)
n_yrs = len(yrs)

# ufile = "u.lowr_era5.197904-202303.nc"

def compute_eof(ufile, save_txt=True):
    ds = xr.open_dataset(ufile)
    anom = ds['anom']  # [year=44, month=12, level=19]
    
    # Ensure dims order: (year, level, month)
    anom = anom.transpose('year','level','month')

    # Reshape into (year, space)
    anom_year_space = anom.stack(space=('level','month'))
    
    # Rename 'year' → 'time' so Eof sees it as the sample axis
    anom_year_space = anom_year_space.rename(year='time')
    
    solver = Eof(anom_year_space)
    
    eofs = solver.eofs(neofs=3)        # (mode, space)
    pcs  = solver.pcs(npcs=3, pcscaling=1)  # (time=44, mode)
    varfrac = solver.varianceFraction(neigs=3)

    eofs_unstacked = eofs.unstack('space').transpose('mode','level','month')

    # -------------------------
    # PC1 time series (standardized & sign-adjusted)
    # -------------------------
    pc1 = pcs[:,0] * -1.0
    pc1 = (pc1 - pc1.mean()) / pc1.std()

    if save_txt:
    # -------------------------
    # Save PC1 time series
    # -------------------------
        np.savetxt("stcmI.1979-2023.txt", pc1.values)

    return pc1


def classify_stmi_events(inFile, thresh=None, save_txt=True, save_mask=True):
    """
    Classify polar vortex events based on STCMI time series following Lim et al. (2018, 2019).
    
    Parameters
    ----------
    inFile : str
        Input file containing STCMI time series (one value per year).
        File name must contain year span, e.g. 'stcmI.1979-2023.txt'.
    thresh : float, optional
        Threshold (in std dev units) for defining events. Default is 0.8.
    save_txt : bool, optional
        If True, save a formatted text file with weak/strong event years. Default True.
    save_mask : bool, optional
        If True, save a CSV mask with year, STMI, and classification. Default True.
    
    Returns
    -------
    stmi : pandas.DataFrame
        DataFrame with columns: ['year','STMI','mask']
        mask is one of {'weak','strong','neutral'}.
    outFile, maskFile : str
        Paths to the saved files (or None if not saved).
    """

    # default
    if thresh is None:
        thresh = 0.8
        
    # --- Parse years from filename ---
    match = re.search(r'(\d{4})-(\d{4})', inFile)
    if not match:
        raise ValueError("Filename must contain year span like '1979-2023'")
    startYear, endYear = map(int, match.groups())

    # --- Read STMI time series ---
    stmi = pd.read_csv(inFile, header=None, names=['STMI'])
    stmi['year'] = np.arange(startYear, startYear+len(stmi))

    # --- Apply classification ---
    stmi['mask'] = 'neutral'
    stmi.loc[stmi['STMI'] >  thresh, 'mask'] = 'weak'
    stmi.loc[stmi['STMI'] < -thresh, 'mask'] = 'strong'


    # --- Save formatted text file ---
    if save_txt:
        if save_txt is True:
            out_dir = Path(".")
        else:
            out_dir = Path(save_txt)
            out_dir.mkdir(parents=True, exist_ok=True)

        # --- Prepare file names ---
        outFile  = out_dir / f"STmode_{startYear}-{endYear}_thr±{thresh:.1f}.txt" 
        

        weak_years   = stmi.loc[stmi['mask'] == 'weak', 'year'].tolist()
        strong_years = stmi.loc[stmi['mask'] == 'strong', 'year'].tolist()
        with open(outFile,'w') as fle:
            fle.write("# Polar Vortex weakening and strengthening events based on the multiple EOF (Lim et al. 2018)\n")
            fle.write(f"# Event threshold ±{thresh} stddev following Lim et al. (2019)\n")
            fle.write("# Weakening years\n")
            for y in weak_years:
                fle.write(f"{y}\n")
            fle.write("#--------------------------\n")
            fle.write("# Strengthening years\n")
            for y in strong_years:
                fle.write(f"{y}\n")
    
    # --- Save mask CSV ---
    if save_mask:
        if save_txt is True:
            out_dir = Path(".")
        else:
            out_dir = Path(save_txt)
            out_dir.mkdir(parents=True, exist_ok=True)
            
        maskFile = out_dir / f"STmode_mask_{startYear}-{endYear}_thr±{thresh:.1f}.csv" 
        stmi[['year','STMI','mask']].to_csv(maskFile, index=False)

    return stmi

    
