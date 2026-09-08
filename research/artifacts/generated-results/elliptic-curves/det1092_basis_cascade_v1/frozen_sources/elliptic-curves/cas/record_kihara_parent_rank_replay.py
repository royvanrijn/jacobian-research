#!/usr/bin/env python3
"""Bind the completed independent first-Kihara-parent proof and stage costs."""
import argparse,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'kihara-rank-standalone-v1';OUT=ART/'kihara_first_parent_rank_replay_v1.json'
def compute():
 stages={};paths=[Path(__file__).resolve()]
 for name,folder,log in [('heights','kihara-section-involution-v1','worker.log'),('counts','kihara-picard-count-v1','worker.log'),('saturation','kihara-first-parent-rank-v1','worker.log'),('independent','kihara-rank-standalone-v1','replay.log')]:
  d=LOCAL/folder;s=json.loads((d/'supervisor.json').read_text());assert s['outcome']=='completed' and s['returncode']==0 and s['failure_reason'] is None
  assert hashlib.sha256((d/log).read_bytes()).hexdigest()==s['log_sha256'];stages[name]=s;paths.extend([d/'supervisor.json',d/'protocol.json',d/log])
 log=(D/'replay.log').read_text();assert log.splitlines()==[
 'PASS78 independent generic heights and exact6+6 symmetry; geometric NS>=18',
 'PASS complete finite quotient rank12 modulo 2',
 'PASS complete finite quotient rank12 modulo 3',
 'PASS p 53 complete Fp/Fp2 counts; rational/geometric Picard upper17/18',
 'PASS p 83 complete Fp/Fp2 counts; rational/geometric Picard upper17/18',
 'PASS full rational NSdet756; MW12/13; every Q-Jacobian fibration MW<=15']
 for p in [ROOT/'elliptic-curves/cas/verify_kihara_parent_rank.sage',ART/'kihara_first_parent_rank_bundle_v1.json']:
  assert p.read_bytes()==(D/p.name).read_bytes();paths.append(p)
 for n in ['kihara_section_involution_v1.json','kihara_first_parent_rank_v1.json']:paths.append(ART/n)
 return {'schema':'elliptic-curves.kihara-first-parent-rank-replay.v1','status':'PASS','stages':stages,'total_supervised_seconds':sum(s['wall_seconds'] for s in stages.values()),'transcript':log,
 'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'scope':'Completed standalone exact replay of the first new Kihara parent: Picard17/18, genericMW12/13, full rational NSdet756 and all-Q-fibration MW15 upper bound. Finite field counts and generic heights use independent arithmetic methods. No specialized-rank upper bound, new fibration, near-record curve or following search.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();r=compute()
 if args.check:assert r==json.loads(OUT.read_text())
 else:
  with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print('PASS independent Kihara parent proof bindings;',r['total_supervised_seconds'],'supervised seconds')
