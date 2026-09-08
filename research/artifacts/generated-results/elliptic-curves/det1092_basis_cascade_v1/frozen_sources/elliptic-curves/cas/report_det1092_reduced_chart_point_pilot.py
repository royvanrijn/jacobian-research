#!/usr/bin/env python3
"""Certified point exposure and pinned-catalogue accounting for the new MW17 parent."""
import argparse,gzip,json
from pathlib import Path
import det1092_reduced_chart_point_pilot_v2 as batch
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
from research_runtime.store import checkpoint
ROOT=batch.ROOT;D=batch.BATCH;ART=batch.ART;OUT=ART/'det1092_reduced_chart_point_pilot_v2.json'
def expected():
 p=batch.protocol();search=cert.read(D/'ledger.json');post=cert.read(D/'post-verification-ledger.json');driver=cert.read(D/'run-ledger.json');transport=cert.read(D/'seed-replay/result.json')
 assert search['status']==post['status']==driver['status']==transport['status']=='PASS'
 parent=cert.read(D/'parent-sections.json');source=ROOT/'artifacts/local/elliptic-curves/det1092-reduced-chart-proof-v2/parent-sections.json';original=cert.read(source)
 assert set(parent)=={'source_sha256','a_invariants','basis_weierstrass_coordinates'} and parent['source_sha256']==cert.hashed(source)
 assert all(parent[k]==original[k] for k in ['a_invariants','basis_weierstrass_coordinates'])
 paths=[source]+[D/n for n in ['protocol.json','parent-sections.json','ledger.json','post-verification-ledger.json','post-verification-protocol.json','run-ledger.json','run-protocol.json','seed-replay/protocol.json','seed-replay/result.json']]
 dbpath=ART/'icarm_catalogue_626_6a7ebb045f37.json.gz';indexpath=ART/'new_high_rank_curve_index_v22.json';paths += [dbpath,indexpath];db=json.loads(gzip.decompress(dbpath.read_bytes()));index=cert.read(indexpath);assert len(db['curves'])==626 and len(index['curves'])==201
 byj={};local={}
 for row in db['curves']:byj.setdefault(catalogue.j(row['ainvs']),[]).append(row)
 for row in index['curves']:local.setdefault(catalogue.j(row['curve']),[]).append(row)
 fresh=D/'seed-replay';frozen=cert.read(fresh/'protocol.json')
 assert all(cert.hashed(fresh/n)==h for n,h in frozen['files'].items())
 assert cert.hashed(fresh/'parent-sections.json')==cert.hashed(D/'parent-sections.json')
 assert cert.hashed(fresh/'verify_det1092_reduced_point_seeds.sage')==cert.hashed(batch.CAS/'verify_det1092_reduced_point_seeds.sage')
 for ident in ['unit','outer']:assert cert.hashed(fresh/(ident+'-seed.json'))==cert.hashed(D/ident/'seed.json')
 rows=[];inputs={};total=sum(r['supervision']['wall_seconds'] for r in driver['stages'])
 for row in p['rows']:
  ident=row['id'];folder=D/ident;seed=cert.read(folder/'seed.json');data=cert.read(folder/'result.json');proof=cert.read(folder/'certification-ledger.json');cp=cert.read(folder/'certification-protocol.json')
  assert proof['status']=='PASS' and proof['initial_rank']==row['initial_rank'] and proof['completed_boxes']==49 and len(data['charts'])==49
  assert all(c['search']['status']=='bounded_search_complete' for c in data['charts']) and data['rank_lower_bound']<=proof['rank_lower_bound']
  assert cert.hashed(folder/'seed.json')==search['seed_hashes'][ident] and seed['points']==seed['generic_points']
  assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items())
  inputs.update(proof['certificates'])
  paths += [folder/n for n in ['seed.json','maps.json','result.json','certification-protocol.json','certification-ledger.json']]
  entry=next(r for r in search['rows'] if r['id']==ident);m=next(r for r in search['maps'] if r['id']==ident)
  stages=[m['supervision']]+[s['supervision'] for s in entry['stages']]+[s['supervision'] for s in proof['stages']]
  assert all(s['outcome']=='completed' and s['returncode']==0 for s in stages);seconds=sum(s['wall_seconds'] for s in stages);total+=seconds
  j=catalogue.j(seed['curve']);jb=[abs(j.numerator).bit_length(),j.denominator.bit_length()];floor=max(1,-(-(jb[0]-25)//6),-(-(jb[1]-16)//6));maps=cert.read(folder/'maps.json');quartic_bits=[max(abs(int(v)).bit_length() for v in r['discriminant_quartic']) for r in maps['rows']]
  rows.append(dict(id=ident,parameter=row['parameter'],curve=seed['curve'],initial_rank=row['initial_rank'],original_parameter=seed['original_parameter'],original_generic_section_count=17,rank_lower_bound=proof['rank_lower_bound'],discovered_rank_gain=proof['rank_lower_bound']-row['initial_rank'],completed_boxes=49,point_count=proof['point_count'],model_coefficient_bits=seed['model_coefficient_bits'],j_bits=jb,quartic_coefficient_floor=floor,quartic_bits_range=[min(quartic_bits),max(quartic_bits)],finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']},geometry_seconds=m['supervision']['wall_seconds'],point_worker_seconds=next(s['supervision']['wall_seconds'] for s in entry['stages'] if s['name']=='worker'),point_search_cpu_ms=sum(c['search']['search_cpu_ms'] or 0 for c in data['charts']),total_stage_seconds=seconds,catalogue_matches=[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])],inventory_matches=[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])]))
 assert len(rows)==2 and len({catalogue.j(r['curve']) for r in rows})==2
 old=D.with_name('det1092-reduced-chart-point-pilot-v1');failed=cert.read(old/'run-ledger.json');assert failed['status']=='FAILED_OR_CENSORED' and not (old/'ledger.json').exists()
 failure_cost=sum(s['supervision']['wall_seconds'] for s in failed['stages']);total+=failure_cost
 paths += [old/n for n in ['run-ledger.json','protocol.json']]
 inputs.update({str(x.relative_to(ROOT)):cert.hashed(x) for x in paths})
 return dict(schema='elliptic-curves.det1092-reduced-chart-point-pilot.v2',status='PASS',rows=rows,distinct_curves=2,distinct_parent_surfaces=1,completed_boxes=98,retained_point_witnesses=sum(r['point_count'] for r in rows),discovered_rank_gains=sum(r['discovered_rank_gain'] for r in rows),total_stage_seconds=total,failed_intake_seconds=failure_cost,inputs=inputs,sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},boundary='Two unchanged fixed coordinates in an exactly reduced chart of the determinant1092 parent, derived from the generic equation only. All maps precede point searches; the blind-regression factor-free mapping and equal49-box exposures are retained. The unit17-section specialization has exact span15, so gains there are relative to15; outer starts17. The V1 intake failure is preserved and counted. All sections, independent seed subsets, exact point histories, and full-cloud mod2/3/5 certificates replay. Post-search novelty comparison uses pinned626 catalogue equations and201 inventory, never for selection or replacement. Same parent and fibration; no exact whole-curve rank, literature-wide novelty, near-record claim, or automatic parameter sweep.')
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS reduced-chart determinant1092 parent98 boxes:',[(a['parameter'],a['rank_lower_bound'],a['point_count']) for a in r['rows']])
