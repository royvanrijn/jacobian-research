#!/usr/bin/env python3
"""Certified point exposure and pinned-catalogue accounting for the new MW17 parent."""
import argparse
from pathlib import Path
import det1092_point_pilot as batch
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
from research_runtime.store import checkpoint
ROOT=batch.ROOT;D=batch.BATCH;ART=batch.ART;OUT=ART/'det1092_point_pilot_v1.json'
def expected():
 p=batch.protocol();search=cert.read(D/'ledger.json');post=cert.read(D/'post-verification-ledger.json');driver=cert.read(D/'run-ledger.json');transport=cert.read(D/'seed-replay/result.json')
 assert search['status']==post['status']==driver['status']==transport['status']=='PASS'
 parent=cert.read(D/'parent-sections.json');source=ART/'curve302_recovered_mw17_parent_v1.json';original=cert.read(source)
 assert set(parent)=={'source_sha256','a_invariants','basis_weierstrass_coordinates'} and parent['source_sha256']==cert.hashed(source)
 assert all(parent[k]==original[k] for k in ['a_invariants','basis_weierstrass_coordinates'])
 paths=[source]+[D/n for n in ['protocol.json','parent-sections.json','ledger.json','post-verification-ledger.json','post-verification-protocol.json','run-ledger.json','run-protocol.json','seed-replay/protocol.json','seed-replay/result.json']]
 dbpath=catalogue.D/'database.json';meta=catalogue.D/'metadata.json';indexpath=ART/'new_high_rank_curve_index_v22.json';paths += [dbpath,meta,indexpath];assert cert.hashed(dbpath)==cert.read(meta)['sha256'];db=cert.read(dbpath);index=cert.read(indexpath);assert len(db['curves'])==620 and len(index['curves'])==201
 byj={};local={}
 for row in db['curves']:byj.setdefault(catalogue.j(row['ainvs']),[]).append(row)
 for row in index['curves']:local.setdefault(catalogue.j(row['curve']),[]).append(row)
 fresh=D/'seed-replay';frozen=cert.read(fresh/'protocol.json')
 assert all(cert.hashed(fresh/n)==h for n,h in frozen['files'].items())
 assert cert.hashed(fresh/'parent-sections.json')==cert.hashed(D/'parent-sections.json')
 assert cert.hashed(fresh/'verify_det1092_point_seeds.sage')==cert.hashed(batch.CAS/'verify_det1092_point_seeds.sage')
 for ident in ['unit','outer']:assert cert.hashed(fresh/(ident+'-seed.json'))==cert.hashed(D/ident/'seed.json')
 rows=[];inputs={};total=sum(r['supervision']['wall_seconds'] for r in driver['stages'])
 for row in p['rows']:
  ident=row['id'];folder=D/ident;seed=cert.read(folder/'seed.json');data=cert.read(folder/'result.json');proof=cert.read(folder/'certification-ledger.json');cp=cert.read(folder/'certification-protocol.json')
  assert proof['status']=='PASS' and proof['initial_rank']==17 and proof['completed_boxes']==49 and len(data['charts'])==49
  assert all(c['search']['status']=='bounded_search_complete' for c in data['charts']) and data['rank_lower_bound']<=proof['rank_lower_bound']
  assert cert.hashed(folder/'seed.json')==search['seed_hashes'][ident] and seed['points']==seed['generic_points']
  assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items())
  inputs.update(proof['certificates'])
  paths += [folder/n for n in ['seed.json','maps.json','result.json','certification-protocol.json','certification-ledger.json']]
  entry=next(r for r in search['rows'] if r['id']==ident);m=next(r for r in search['maps'] if r['id']==ident)
  stages=[m['supervision']]+[s['supervision'] for s in entry['stages']]+[s['supervision'] for s in proof['stages']]
  assert all(s['outcome']=='completed' and s['returncode']==0 for s in stages);seconds=sum(s['wall_seconds'] for s in stages);total+=seconds
  j=catalogue.j(seed['curve'])
  rows.append(dict(id=ident,parameter=row['parameter'],curve=seed['curve'],initial_rank=17,rank_lower_bound=proof['rank_lower_bound'],discovered_rank_gain=proof['rank_lower_bound']-17,completed_boxes=49,point_count=proof['point_count'],model_coefficient_bits=seed['model_coefficient_bits'],finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']},geometry_seconds=m['supervision']['wall_seconds'],point_worker_seconds=next(s['supervision']['wall_seconds'] for s in entry['stages'] if s['name']=='worker'),point_search_cpu_ms=sum(c['search']['search_cpu_ms'] or 0 for c in data['charts']),total_stage_seconds=seconds,catalogue_matches=[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])],inventory_matches=[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])]))
 assert len(rows)==2 and len({catalogue.j(r['curve']) for r in rows})==2
 inputs.update({str(x.relative_to(ROOT)):cert.hashed(x) for x in paths})
 return dict(schema='elliptic-curves.det1092-point-pilot.v1',status='PASS',rows=rows,distinct_curves=2,distinct_parent_surfaces=1,completed_boxes=98,retained_point_witnesses=sum(r['point_count'] for r in rows),discovered_rank_gains=sum(r['discovered_rank_gain'] for r in rows),total_stage_seconds=total,inputs=inputs,sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},boundary='Two fixed non-anchor fibres of the new determinant1092 arithmetic MW17 parent, chosen before scores or point search. Execution input contains only the equation and17 generic sections; no public exceptional point or record-group embedding. Exact seed transport and full-cloud finite proofs precede any rank claim. Catalogue comparison is post-search against pinned620 equations and201 inventory only. No exact-rank, literature-wide novelty or parent-superiority claim; no automatic following wave.')
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS new determinant1092 parent98 boxes:',[(a['parameter'],a['rank_lower_bound'],a['point_count']) for a in r['rows']])
