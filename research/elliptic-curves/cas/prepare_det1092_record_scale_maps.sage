#!/usr/bin/env sage-python
"""The calibrated own-subgroup geometry, on one frozen prospective wave."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import det1092_record_scale_points as trial
geometry=SourceFileLoader('record_scale_geometry',str(CAS/'prepare_retained_native19_trial_v3.sage')).load_module()
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,required=True);ap.add_argument('--wave',type=int,required=True);a=ap.parse_args()
    trial.configure(a.index,a.wave);geometry.control=trial;geometry.main()
