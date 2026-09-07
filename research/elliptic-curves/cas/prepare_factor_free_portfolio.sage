#!/usr/bin/env sage-python
"""Reuse the calibrated factor-free own-subgroup geometry on the frozen pair."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
import prospective_factor_free_portfolio as trial
geometry=SourceFileLoader('retained26_geometry',str(CAS/'prepare_retained_native19_trial_v3.sage')).load_module()
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);args=p.parse_args();trial.configure(args.index);geometry.control=trial;
 if 'previous_maps' in trial.ROW:
  from research_runtime.store import checkpoint,digest
  import certify_compact_r17_candidates as cert
  maps=cert.read(trial.ROOT/trial.ROW['previous_maps']);fresh=cert.read(trial.ROOT/trial.ROW['factor_free_maps'])
  assert len(fresh['rows'])==49
  maps['rows']=[r['mapping'] for r in fresh['rows']]
  assert [m['centre'] for m in maps['rows']]==maps['centres']
  assert [r['parity'] for r in maps['sample']]==trial.masks(trial.protocol())
  maps['protocol_hash']=digest(trial.protocol());maps['status']='COMPLETE_DECLARED_MAPS'
  assert not (trial.D/'maps.json').exists();checkpoint(trial.D/'maps.json',maps);print('FROZEN49 audited factor-free boxes',flush=True)
 else:geometry.main()
