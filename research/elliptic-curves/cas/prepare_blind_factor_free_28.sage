#!/usr/bin/env sage-python
"""Regenerate the complete control geometry from its original27 seed."""
from pathlib import Path
from importlib.machinery import SourceFileLoader
import sys
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import blind_factor_free_28_control as control
from research_runtime.store import checkpoint
geometry=SourceFileLoader('blind_geometry',str(CAS/'prepare_retained_native19_trial_v3.sage')).load_module()
geometry.control=control
geometry.main()
checkpoint(control.D/'geometry-data-access.json',sorted(control.READS))
