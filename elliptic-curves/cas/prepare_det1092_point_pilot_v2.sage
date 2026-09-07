#!/usr/bin/env sage-python
"""Reuse the calibrated factor-free own-subgroup geometry on the frozen pair."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import det1092_point_pilot_v2 as trial
geometry=SourceFileLoader('retained26_geometry',str(CAS/'prepare_retained_native19_trial_v3.sage')).load_module()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);args=p.parse_args();trial.configure(args.index);geometry.control=trial;geometry.main()
