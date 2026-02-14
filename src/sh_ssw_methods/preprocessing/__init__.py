# src/sh_ssw_methods/preprocessing/__init__.py

"""
sh_ssw_methods.preprocessing
----------------------------
Provides data-preparation utilities for SSW detection.
"""

from .prepare_polarT_files import (
    output_polarT_prefiles
)
from .prepare_polarT_xarray import prep_input_polarT

from .return_fsw_dates import detect_fsw_SH

from .prepare_eof import compute_monthly_eof_series

from .project_mon_eof import output_daily_pc1

from .prepare_stcmi import compute_stcmI


__all__ = [
    "output_polarT_prefiles", "prep_input_polarT",
    "detect_fsw_SH", "compute_monthly_eof_series",
    "output_daily_pc1", "compute_stcmI",
]

