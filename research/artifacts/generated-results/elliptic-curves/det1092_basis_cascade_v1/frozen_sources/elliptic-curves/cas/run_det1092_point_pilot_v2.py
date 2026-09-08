#!/usr/bin/env python3
"""Finish the frozen portfolio, then run the independently motivated new-parent pilot."""
import sys,time,shutil
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';L=ROOT/'artifacts/local/elliptic-curves';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def main():
 import det1092_point_pilot_v2 as pilot
 d=pilot.BATCH;assert not (d/'run-ledger.json').exists();ledger={'status':'RUNNING','stages':[]};checkpoint(d/'run-ledger.json',ledger)
 checkpoint(d/'run-protocol.json',dict(driver_sha256=cert.hashed(Path(__file__)),pilot_sha256=cert.hashed(CAS/'det1092_point_pilot_v2.py'),maximum_workers=1,rss_bytes=2147483648,freeze_seconds=60,intake_seconds=180,seed_replay_seconds=120))
 def stage(name,cmd,seconds,cwd=ROOT):
  s=run(cmd,limits=Limits(seconds,2147483648),log_path=d/(name+'-driver.log'),checkpoint_path=d/(name+'-driver.supervisor.json'),cwd=cwd)
  ok=s['outcome']=='completed' and s['returncode']==0;ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(d/'run-ledger.json',ledger);print('DET1092',name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(d/'run-ledger.json',ledger);raise ArithmeticError('preserve failed new-parent stage')
 for name,seconds in [('freeze',60),('intake',180)]:stage(name,[sys.executable,str(CAS/'det1092_point_pilot_v2.py'),name],seconds)
 fresh=d/'seed-replay';fresh.mkdir(exist_ok=False)
 files=[CAS/'verify_det1092_point_seeds.sage',d/'parent-sections.json']
 for path in files:shutil.copy2(path,fresh/path.name)
 for ident in ['unit','outer']:shutil.copy2(d/ident/'seed.json',fresh/(ident+'-seed.json'))
 checkpoint(fresh/'protocol.json',dict(files={p.name:cert.hashed(p) for p in fresh.iterdir() if p.is_file()},seconds=120,rss_bytes=2147483648))
 stage('seed-replay',[SAGE,str(fresh/'verify_det1092_point_seeds.sage'),'--parent',str(fresh/'parent-sections.json'),'--seeds',str(fresh/'unit-seed.json'),str(fresh/'outer-seed.json'),'--output',str(fresh/'result.json')],120,fresh)
 pilot.launch()
 import complete_det1092_point_pilot_v2 as proof
 proof.main()
 ledger['status']='PASS';checkpoint(d/'run-ledger.json',ledger);print('COMPLETE all point exposure and certificates',flush=True)
if __name__=='__main__':main()
