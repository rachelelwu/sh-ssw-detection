# /src/sh_ssw_methods/utils/event_iden_funcs.py

from __future__ import annotations
import numpy as np
import pandas as pd
import xarray as xr

# event detection

def group_consecutive_dates(dates: pd.Index | pd.Series) -> pd.DataFrame:
    """
    Group consecutive dates and return events with a count of consecutive days (how long the event lasts) later for persistence.
    Returns a DataFrame with columns: first, last, count.
    Assumes daily frequency (i.e., 'consecutive' = 1 day apart).
    """
    d = pd.to_datetime(pd.Index(dates)).sort_values().unique()
    df = pd.DataFrame({"date": d})
    gap_days = df["date"].diff().dt.days
    df["event_id"] = (gap_days != 1).cumsum()
    events = (
        df.groupby("event_id")["date"]
          .agg(first="min", last="max", count="size")
          .reset_index(drop=True)
    )
    return events # columns: first, last, count

def apply_persistence(events: pd.DataFrame, min_persist_days: int = 3) -> pd.DataFrame:
    """
    Keep only events that last at least min_persist_days.
    """
    return events.loc[events["count"] >= min_persist_days].reset_index(drop=True)


def enforce_min_gap(runs: pd.DataFrame, min_gap_days: int = 20) -> pd.DataFrame:
    """
    Enforce independence: keep an event if the gap between its FIRST day
    and the PREVIOUS event's LAST day is > min_gap_days.
    Keeps the first event by definition.
    """
    r = runs.sort_values("first").reset_index(drop=True).copy()
    prev_last = r["last"].shift()
    gap = (r["first"] - prev_last).dt.days
    independent = prev_last.isna() | (gap > min_gap_days)
    return r.loc[independent].reset_index(drop=True)



# outputting event dates

def build_events_df(
    dates,
    *,
    method: str,
    definition: str,
    data_source: str,
    threshold: float | None = None,
    level_hpa: float | None = None,
    latitude: float | None = None,
    lat_band: str | None = None,
    notes: str | None = None,
    extra_cols: dict | None = None,   # NEW: per-event columns (e.g., {"drop_val": drop_vals})
) -> pd.DataFrame:
    """
    Create a unified event table from arrays/iterables.

    Parameters
    ----------
    dates : array-like of datetime-like
        Event 'date' (e.g., onset). Required.
    extra_cols : dict, optional
        Mapping of column_name -> scalar or array-like of length len(dates).
        Use this to add per-event metrics, e.g., {"drop_val": drop_vals, "window_days": 10}.
        Scalars will be broadcast to all rows.
    """

    df = pd.DataFrame({"date": pd.to_datetime(dates)})
    n = len(df)

    # Identity & provenance (constant per table)
    df["method"] = method
    df["definition"] = definition
    df["data_source"] = data_source

    # Optional metadata (constant per table unless you pass via extra_cols)
    df["threshold"] = threshold
    df["level_hpa"] = level_hpa
    df["latitude"] = latitude
    df["lat_band"] = lat_band
    df["notes"] = notes

    # Add per-event columns
    if extra_cols:
        for key, val in extra_cols.items():
            # Broadcast scalars
            if np.ndim(val) == 0:
                df[key] = val
            else:
                s = pd.Series(val)
                if len(s) != n:
                    raise ValueError(
                        f"extra_cols['{key}'] length {len(s)} != number of dates {n}"
                    )
                df[key] = s.values

    # Enforce dtypes for selected string columns
    for col in ["method", "definition", "data_source", "lat_band"]:
        df[col] = df[col].astype("string")

    return df.sort_values("date").reset_index(drop=True)

