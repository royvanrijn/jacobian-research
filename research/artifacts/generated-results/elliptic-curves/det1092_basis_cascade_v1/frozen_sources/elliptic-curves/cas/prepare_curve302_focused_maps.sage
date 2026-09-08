#!/usr/bin/env sage-python
"""Only the configured seed supplies points to the unchanged geometry engine."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import curve302_focused_point_exposure as trial
geometry=SourceFileLoader('focused_geometry',str(CAS/'prepare_retained_native19_trial_v3.sage')).load_module()
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,required=True);a=ap.parse_args()
    trial.configure(a.index);geometry.control=trial;geometry.main()
