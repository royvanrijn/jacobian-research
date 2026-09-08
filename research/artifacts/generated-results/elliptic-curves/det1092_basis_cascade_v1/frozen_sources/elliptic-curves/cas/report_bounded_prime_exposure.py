#!/usr/bin/env python3
"""Completed blind regression, exact history gate, point exposure and height audit."""
import argparse
from pathlib import Path
import bounded_prime_point_portfolio as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=batch.ROOT;D=batch.BATCH;ART=batch.ART;OUT=ART/'bounded_prime_point_portfolio_v1.json'
def expected():
 p=batch.protocol();search=cert.read(D/'ledger.json');post=cert.read(D/'post-verification-ledger.json');gate=cert.read(ART/'blind_bounded_prime_28_control_v2.json');height=cert.read(ART/'bounded_prime_history_and_height_v1.json');vis=cert.read(ART/'det1092_refined_visibility_v1.json')
 assert search['status']==post['status']==height['status']==vis['status']==gate['status']=='PASS' and gate['rank_lower_bound']==28 and gate['completed_boxes']==49
 paths=[ART/n for n in ['blind_bounded_prime_28_control_v2.json','bounded_prime_history_and_height_v1.json','det1092_refined_visibility_v1.json']]+[D/n for n in ['protocol.json','ledger.json','post-verification-ledger.json','post-verification-protocol.json','history-height-protocol.json','history-height.supervisor.json','history-height-v2-protocol.json','history-height-v2.supervisor.json','visibility-protocol.json','visibility.supervisor.json']];inputs={};rows=[];total=0
 assert height['source_sha256']==cert.hashed(batch.CAS/'verify_bounded_prime_history_and_height_v2.sage') and all(cert.hashed(ROOT/n)==h for n,h in height['inputs'].items())
 assert all(cert.hashed(ROOT/n)==h for n,h in {**vis['inputs'],**vis['sources']}.items())
 for r in p['rows']:
  ident=r['id'];folder=D/ident;entry=next(e for e in search['rows'] if e['id']==ident);m=next(e for e in search['maps'] if e['id']==ident);h=next(a for a in height['rows'] if a['id']==ident);total+=m['supervision']['wall_seconds'];paths += [folder/'maps.json',folder/'seed.json']
  row=dict(id=ident,initial_rank=r['initial_rank'],point_boxes=entry['point_boxes'],history=h,geometry_seconds=m['supervision']['wall_seconds'])
  if entry['status']=='SKIPPED_NO_PROVED_NEW_HISTORY_BOX':
   assert h['proved_covered_by_prior_boxes']==49 and entry['point_boxes']==0 and not (folder/'result.json').exists();row.update(status='PROVED_PREVIOUSLY_COVERED',rank_lower_bound=r['initial_rank'],discovered_rank_gain=0,retained_points=0)
  else:
   assert entry['status']=='PASS' and h['proved_new_history_boxes']==49
   proof=cert.read(folder/'certification-ledger.json');cp=cert.read(folder/'certification-protocol.json');data=cert.read(folder/'result.json');assert proof['status']=='PASS' and proof['completed_boxes']==49 and proof['initial_rank']==r['initial_rank']
   assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['sources'],**cp['inputs'],**proof['certificates']}.items());inputs.update(proof['certificates'])
   paths += [folder/n for n in ['result.json','certification-ledger.json','certification-protocol.json']]
   total += sum(a['supervision']['wall_seconds'] for a in entry['stages']+proof['stages'])
   row.update(status='COMPLETED_AND_CERTIFIED',rank_lower_bound=proof['rank_lower_bound'],discovered_rank_gain=proof['discovered_rank_gain'],retained_points=proof['point_count'],finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']},returned_finite_points=sum(len(a['search']['finite_curve_points']) for a in data['charts']),point_search_cpu_ms=sum(a['search']['search_cpu_ms'] for a in data['charts']))
  rows.append(row)
 for name in ['history-height','history-height-v2','visibility']:
  supervision=cert.read(D/(name+'.supervisor.json'));total+=supervision['wall_seconds']
  assert supervision['returncode'] is not None
  if name!='history-height':assert supervision['outcome']=='completed' and supervision['returncode']==0
 inputs.update({str(f.relative_to(ROOT)):cert.hashed(f) for f in paths})
 return dict(schema='elliptic-curves.bounded-prime-point-portfolio.v1',status='PASS',rows=rows,completed_point_boxes=sum(r['point_boxes'] for r in rows),skipped_proved_covered_boxes=sum(r['history']['proved_covered_by_prior_boxes'] for r in rows),discovered_rank_gains=sum(r['discovered_rank_gain'] for r in rows),total_stage_seconds=total,blind_control_stage_seconds=gate['total_stage_seconds'],combined_supervised_seconds=total+gate['total_stage_seconds'],inputs=inputs,sources={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))},boundary='No new parameter population. Full49 original27-only restricted-prime control recovered28; two retained27 rosters prove covered by declared completed history and are skipped. Both new-parent49 rosters prove new coordinate exposure and complete without new finite points. All exact maps, history witnesses, point clouds and ranks independently checked. Coefficient bounds concern integral short and completed-square quartic models, not point-search heights or rank upper bounds. Original signed generic-point visibility is retrospective and does not cover all translated representatives. No new curve rank or near-record result.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS bounded-prime refinement:',r['completed_point_boxes'],'new boxes;',r['skipped_proved_covered_boxes'],'covered boxes skipped;',r['discovered_rank_gains'],'directions;',r['combined_supervised_seconds'],'combined seconds')
