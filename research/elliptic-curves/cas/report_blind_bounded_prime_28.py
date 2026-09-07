#!/usr/bin/env python3
"""Report the original27-only, full-roster factor-free regression."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/blind-bounded-prime-28-control-v2';OUT=ROOT/'artifacts/generated-results/elliptic-curves/blind_bounded_prime_28_control_v2.json'
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
 maps=cert.read(D/'maps.json');from collections import Counter
 from research_runtime.projective_box_change import classify
 assert all(m['minimization_primes']==p['minimization_primes'] and classify(m['first_matrix'],m['matrix'],p['height'])==m['comparison_to_factor_free'] for m in maps['rows'])
 failed=D.with_name('blind-bounded-prime-28-control-v1');fl=cert.read(failed/'ledger.json');assert fl['status']=='FAILED_OR_CENSORED' and not (failed/'result.json').exists();failure=sum(a['supervision']['wall_seconds'] for a in fl['stages']);inputs += [failed/'protocol.json',failed/'ledger.json']
 return dict(minimization_primes=p['minimization_primes'],box_comparison_counts=dict(Counter(m['comparison_to_factor_free']['status'] for m in maps['rows'])),failed_preparation_seconds=failure,schema='elliptic-curves.blind-bounded-prime-28-control.v2',status='PASS',initial_rank=27,rank_lower_bound=28,discovered_rank_gain=1,completed_boxes=49,retained_points=proof['point_count'],finite_ranks={'2':28,**proof['odd_modulus_ranks']},input_access_counts=access,inputs={str(x.relative_to(ROOT)):cert.hashed(x) for x in inputs},certificates=proof['certificates'],sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},total_stage_seconds=failure+sum(r['supervision']['wall_seconds'] for r in ledger['stages']+proof['stages']),boundary='Known public28 curve recovered from original27 only with explicit25-prime minimisation. No unrestricted minimisation or complete discriminant factorization requested. All49 charts regenerated before points; no public exceptional point, previous maps or successful chart enters geometry/execution. Python artifact-open audit is a process-level audit, not an OS sandbox. Prior one-chart success was known and disclosed. This is a method regression, not a new curve or new record.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS blind factor-free full49:27 ->28; mod2/3/5 and standalone Sage')
