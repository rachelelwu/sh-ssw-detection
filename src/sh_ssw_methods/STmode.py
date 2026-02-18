import re
import numpy as np
import pandas as pd
from pathlib import Path

from .utils.event_iden_funcs import (
    build_events_df
)

def detect_stmi_events(inFile, thresh=None):
    """
    Classify STCMI events (strong/weak tropical circulation modes)
    following Lim et al. (2018, 2019).

    Output is standardized to match `detect_ssw_polarT`:
    returns (events_df, event_years).

    Parameters
    ----------
    inFile : str
        Input text file containing STCMI time series (one value per year).
        File name must contain year span, e.g. 'stcmI.1979-2023.txt'.
    thresh : float, optional
        Threshold (in std dev units) for defining events. Default = 0.8.
    save_txt : bool or str, optional
        If True, save formatted event-year text to current directory.
        If str (directory path), save there.
    save_mask : bool or str, optional
        If True, save CSV mask (year, STMI, mask) to current directory.
        If str (directory path), save there.

    Returns
    -------
    events_df : pandas.DataFrame
        Standardized event metadata table (built via `build_events_df`).
    event_years : pandas.Index
        Years of weak/strong events.
    """

    # --- Defaults ---
    if thresh is None:
        thresh = 0.8

    # --- Parse year span from filename ---
    match = re.search(r'(\d{4})-(\d{4})', inFile)
    if not match:
        raise ValueError("Filename must contain year span like 'stcmI.1979-2023.txt'")
    startYear, endYear = map(int, match.groups())

    # --- Load STCMI time series ---
    stmi = pd.read_csv(inFile, header=None, names=['STMI'])
    stmi['year'] = np.arange(startYear, startYear + len(stmi))

    # --- Classify ---
    stmi['mask'] = 'neutral'
    stmi.loc[stmi['STMI'] >  thresh, 'mask'] = 'weak'
    stmi.loc[stmi['STMI'] < -thresh, 'mask'] = 'strong'

    # -------------------------------
    #  Build standardized event table
    # -------------------------------
    # Build list of event years (like event_dates)
    event_years = pd.Index(stmi.loc[stmi['mask'] != 'neutral', 'year'])

    # Convert event years to datetime-like index (Jan 1st of each year)
    event_dates = pd.to_datetime(event_years.astype(str) + "-01-01")

    # Map to standardized format
    events_df = build_events_df(
        dates=event_dates,
        method="STCMI",
        definition=f"STCMI_±{thresh:.1f}σ",
        data_source="ERA5",
        threshold=f"±{thresh:.1f}σ",
        level_hpa=None,
        notes="STCMI weak/strong events following Lim et al. (2018, 2019)",
    )

    # Attach event type (weak/strong) to DataFrame
    events_df["mask"] = stmi.set_index("year").loc[event_years, "mask"].values

    return events_df, event_years
