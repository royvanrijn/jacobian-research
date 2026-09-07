#!/usr/bin/env python3
"""Exact completed exposure, ranks and cost accounting for the frozen11-curve trial."""
import argparse
from pathlib import Path
import prospective_factor_free_portfolio_v2 as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=batch.ROOT;D=batch.BATCH;ART=batch.ART;OUT=ART/'prospective_factor_free_portfolio_v2.json'
def supervision_seconds(path):
 s=cert.read(path);assert s['returncode'] is not None;return s['wall_seconds']
def expected():
 p=batch.protocol();search=cert.read(D/'ledger.json');post=cert.read(D/'post-verification-ledger.json')
 assert search['status']==post['status']=='PASS' and len(search['rows'])==len(post['rows'])==11 and len(search['maps'])==11
 inputs={str((D/n).relative_to(ROOT)):cert.hashed(D/n) for n in ['protocol.json','ledger.json','post-verification-ledger.json','post-verification-protocol.json','freeze-protocol.json','freeze.supervisor.json']};rows=[];models=[];total=supervision_seconds(D/'freeze.supervisor.json')
 for row in p['rows']:
  ident=row['id'];folder=D/ident;proof=cert.read(folder/'certification-ledger.json');data=cert.read(folder/'result.json');seed=cert.read(folder/'seed.json');maps=cert.read(folder/'maps.json')
  assert proof['status']=='PASS' and proof['initial_rank']==row['initial_rank'] and proof['completed_boxes']==49
  assert len(data['charts'])==49 and all(c['search']['status']=='bounded_search_complete' for c in data['charts'])
  assert data['rank_lower_bound']<=proof['rank_lower_bound']
  initial=len(seed['points']);assert data['initial_dimension']==initial
  # Distinct j proves distinct Q-isomorphism classes. Catalogue status stays the pinned intake status.
  model=tuple(map(cert.F,data['curve']));v=cert.weierstrass_invariants(model);j=v['c4']**3/v['discriminant'];assert j not in models;models.append(j)
  cp=cert.read(folder/'certification-protocol.json')
  assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items())
  inputs.update(proof['certificates'])
  for n in ['seed.json','maps.json','result.json','certification-protocol.json','certification-ledger.json']:
   f=folder/n;inputs[str(f.relative_to(ROOT))]=cert.hashed(f)
  entry=next(r for r in search['rows'] if r['id']==ident);m=next(r for r in search['maps'] if r['id']==ident)
  assert entry['status']==m['status']=='PASS'
  stages=[m['supervision']]+[s['supervision'] for s in entry['stages']]+[s['supervision'] for s in proof['stages']]
  assert all(s['outcome']=='completed' and s['returncode']==0 for s in stages)
  seconds=sum(s['wall_seconds'] for s in stages);total+=seconds
  rows.append(dict(id=ident,parent=row['parent'],family=row['family'],parameter=row['parameter'],initial_rank=initial,inherited_seed_dimension=row['generic_dimension'],rank_lower_bound=proof['rank_lower_bound'],discovered_rank_gain=proof['rank_lower_bound']-initial,completed_boxes=49,point_count=proof['point_count'],finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']},geometry_seconds=m['supervision']['wall_seconds'],point_worker_seconds=next(s['supervision']['wall_seconds'] for s in entry['stages'] if s['name']=='worker'),point_search_cpu_ms=sum(c['search']['search_cpu_ms'] or 0 for c in data['charts']),total_stage_seconds=seconds))
 old=D.with_name('prospective-factor-free-portfolio-v1');failed=cert.read(old/'ledger.json');assert failed['status']=='FAILED_OR_CENSORED' and not failed['rows']
 failure_cost=supervision_seconds(old/'freeze.supervisor.json')+sum(s['supervision']['wall_seconds'] for s in failed['maps'])
 for n in ['protocol.json','ledger.json','freeze.supervisor.json']:
  f=old/n;inputs[str(f.relative_to(ROOT))]=cert.hashed(f)
 total+=failure_cost
 return dict(schema='elliptic-curves.prospective-factor-free-portfolio.v2',status='PASS',rows=rows,distinct_curves=11,distinct_parent_surfaces=len({r['parent'] for r in rows}),completed_boxes=sum(r['completed_boxes'] for r in rows),retained_point_witnesses=sum(r['point_count'] for r in rows),discovered_rank_gains=sum(r['discovered_rank_gain'] for r in rows),curves_with_gain=sum(r['discovered_rank_gain']>0 for r in rows),retained26_27_rank_gains=sum(r['discovered_rank_gain'] for r in rows if r['initial_rank']>=26),failed_preparation_seconds=failure_cost,total_stage_seconds=total,inputs=inputs,sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},boundary='Fixed11-curve point exposure, no new parameter scan. Exact discovered directions relative to each prior certified subgroup. Finite quotient ranks are lower bounds, not upper bounds or proofs of dependence. No known public exceptional point or validation score enters selection/execution. Three parent surfaces; X948 fibration labels are not distinct parents. Prior intake novelty is relative to the pinned620 catalogue and201 inventory only. Failed V1 geometry is charged, with zero point boxes; no automatic following wave.')
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS prospective539 boxes;',r['discovered_rank_gains'],'directions;',r['retained26_27_rank_gains'],'on retained26/27;',r['total_stage_seconds'],'seconds')
