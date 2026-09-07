#!/usr/bin/env python3
"""Report the original27-only, full-roster factor-free regression."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/blind-factor-free-28-control-v1';OUT=ROOT/'artifacts/generated-results/elliptic-curves/blind_factor_free_28_control_v1.json'
def expected():
 p=cert.read(D/'protocol.json');ledger=cert.read(D/'ledger.json');proof=cert.read(D/'certification-ledger.json')
 assert ledger['status']==proof['status']=='PASS' and proof['rank_lower_bound']==28 and proof['initial_rank']==27 and proof['completed_boxes']==49 and proof['odd_modulus_ranks']=={'3':28,'5':28}
 assert all(cert.hashed(ROOT/n)==h for n,h in p['sources'].items())
 original=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1/seed.json'
 assert cert.hashed(original)==p['original_seed_sha256'] and cert.read(original)==cert.read(D/'seed.json')
 inputs=[D/n for n in ['protocol.json','ledger.json','certification-protocol.json','certification-ledger.json','maps.json','result.json','seed.json']]+[original];access={}
 for stage in ['freeze','geometry','worker','replay']:
  path=D/(stage+'-data-access.json');inputs.append(path);reads=cert.read(path)
  assert all((ROOT/n).is_relative_to(D) or ROOT/n==original for n in reads)
  access[stage]=len(reads)
 cp=cert.read(D/'certification-protocol.json')
 assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items())
 return dict(schema='elliptic-curves.blind-factor-free-28-control.v1',status='PASS',initial_rank=27,rank_lower_bound=28,discovered_rank_gain=1,completed_boxes=49,retained_points=proof['point_count'],finite_ranks={'2':28,**proof['odd_modulus_ranks']},input_access_counts=access,inputs={str(x.relative_to(ROOT)):cert.hashed(x) for x in inputs},certificates=proof['certificates'],sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},total_stage_seconds=sum(r['supervision']['wall_seconds'] for r in ledger['stages']+proof['stages']),boundary='Known public28 curve recovered from original27 only. All49 charts regenerated before points; no public exceptional point, previous maps or successful chart enters geometry/execution. Python artifact-open audit is a process-level audit, not an OS sandbox. Prior one-chart success was known and disclosed. This is a method regression, not a new curve or new record.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS blind factor-free full49:27 ->28; mod2/3/5 and standalone Sage')
