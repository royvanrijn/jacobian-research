#!/usr/bin/env python3
"""Fail-closed final ranks, exposure, strata, costs and post-search catalogue check."""
import argparse,gzip,json
from pathlib import Path
import mixed_reduced_parent_exposure as batch
import audit_inventory200_current_catalogue as catalogue
from research_runtime.store import checkpoint
cert=batch.cert;ROOT=batch.ROOT;D=batch.BATCH;ART=batch.ART;OUT=ART/'mixed_reduced_parent_exposure_v2.json';REC=batch.L/'mixed-reduced-parent-relocation-replay-v1'
def expected():
 p=batch.protocol();search=cert.read(D/'ledger.json');post=cert.read(D/'post-ledger.json');driver=cert.read(D/'run-ledger.json');continuation=cert.read(REC/'ledger.json');assert search['status']==continuation['status']=='PASS' and post['status']==driver['status']=='RUNNING'
 assert len(search['rows'])==len(continuation['rows'])==len(search['maps'])==9 and len(post['rows'])==4
 paths=[D/n for n in ['protocol.json','prepare-protocol.json','intake-input.json','intake-result.json','run-protocol.json','run-ledger.json','ledger.json','post-ledger.json','seed-replay/protocol.json','seed-replay/result.json']]
 paths += [REC/'protocol.json',REC/'ledger.json']
 fresh=D/'seed-replay';assert all(cert.hashed(fresh/n)==h for n,h in cert.read(fresh/'protocol.json')['files'].items())
 scorepath=batch.L/'det1092-reduced-score-strata-v1/result.json';sp=batch.L/'det1092-reduced-score-strata-v1/protocol.json';score=cert.read(scorepath);assert score['status']=='PASS' and score['protocol_sha256']==cert.hashed(sp);paths += [scorepath,sp]
 dbpath=ART/'icarm_catalogue_630_b742f700b834.json.gz';idxpath=ART/'new_high_rank_curve_index_v22.json';db=json.loads(gzip.decompress(dbpath.read_bytes()));index=cert.read(idxpath);assert len(db['curves'])==630 and len(index['curves'])==201;paths += [dbpath,idxpath];byj={};local={}
 for r in db['curves']:byj.setdefault(catalogue.j(r['ainvs']),[]).append(r)
 for r in index['curves']:local.setdefault(catalogue.j(r['curve']),[]).append(r)
 interruption=next(cert.read(x) for x in (REC/'interruption').iterdir() if x.name.endswith('certification-driver.supervisor.json'));assert interruption['outcome']=='backend_failure'
 inputs={};rows=[];seen=set();total=interruption['wall_seconds']+sum(s['supervision']['wall_seconds'] for s in driver['stages'])+score['total_score_seconds']
 for row in p['rows']:
  ident=row['id'];folder=D/ident;seed=cert.read(folder/'seed.json');data=cert.read(folder/'result.json');continued=next(r for r in continuation['rows'] if r['id']==ident);pf=ROOT/continued['folder'];proof=cert.read(pf/'certification-ledger.json');cp=cert.read(pf/'certification-protocol.json')
  assert proof['status']=='PASS' and proof['initial_rank']==row['initial_rank'] and proof['completed_boxes']==49 and len(data['charts'])==49
  assert all(c['search']['status']=='bounded_search_complete' for c in data['charts']) and data['rank_lower_bound']<=proof['rank_lower_bound']
  assert all(v==proof['rank_lower_bound'] for v in proof['odd_modulus_ranks'].values())
  assert cert.hashed(folder/'seed.json')==row['seed_sha256'] and cert.hashed(fresh/(ident+'-seed.json'))==row['seed_sha256']
  assert all(cert.hashed(ROOT/n)==h for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items());inputs.update(proof['certificates'])
  paths += [folder/n for n in ['seed.json','maps.json','result.json']]+[pf/n for n in ['certification-protocol.json','certification-ledger.json']]
  assert all(cert.hashed(pf/n)==cert.hashed(folder/n) for n in ['seed.json','maps.json','result.json'])
  entry=next(r for r in search['rows'] if r['id']==ident);m=next(r for r in search['maps'] if r['id']==ident);stages=[m['supervision']]+[s['supervision'] for s in entry['stages']]+[s['supervision'] for s in proof['stages']]
  assert all(s['outcome']=='completed' and s['returncode']==0 for s in stages);seconds=sum(s['wall_seconds'] for s in stages);total+=seconds
  j=catalogue.j(seed['curve']);assert j not in seen;seen.add(j)
  matches=[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])];own=[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])]
  if row['parent']=='X948':assert own==[ident]
  score_row=next((r for r in score['selected'] if ident=='d1092-'+r['id']),None)
  rows.append(dict(id=ident,parent=row['parent'],family=row['family'],parameter=row['parameter'],original_parameter=row.get('original_parameter'),stratum=row['stratum'],curve=seed['curve'],initial_rank=row['initial_rank'],rank_lower_bound=proof['rank_lower_bound'],discovered_rank_gain=proof['rank_lower_bound']-row['initial_rank'],completed_boxes=49,point_count=proof['point_count'],finite_ranks={'2':proof['rank_lower_bound'],**proof['odd_modulus_ranks']},training_score_units=score_row['score_units'] if score_row else None,j_bits=[abs(j.numerator).bit_length(),j.denominator.bit_length()],geometry_seconds=m['supervision']['wall_seconds'],point_worker_seconds=next(s['supervision']['wall_seconds'] for s in entry['stages'] if s['name']=='worker'),point_search_cpu_ms=sum(c['search']['search_cpu_ms'] or 0 for c in data['charts']),total_stage_seconds=seconds,catalogue_matches=matches,inventory_matches=own))
 strata=[]
 for band in ['strong','moderate','lower_fixed','retained26_gap']:
  r=[r for r in rows if r['stratum']==band];g=sum(v['discovered_rank_gain'] for v in r);t=sum(v['point_worker_seconds'] for v in r);strata.append(dict(stratum=band,curves=len(r),initial_ranks=[v['initial_rank'] for v in r],final_bounds=[v['rank_lower_bound'] for v in r],completed_boxes=49*len(r),discovered_rank_gains=g,worker_seconds=t,total_stage_seconds=sum(v['total_stage_seconds'] for v in r),point_search_cpu_ms=sum(v['point_search_cpu_ms'] for v in r),certified_directions_per_worker_hour=3600*g/t,certified_directions_per_total_stage_hour=3600*g/sum(v['total_stage_seconds'] for v in r)))
 inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in paths})
 return dict(schema='elliptic-curves.mixed-reduced-parent-exposure.v2',status='PASS',rows=rows,strata=strata,distinct_parents=2,completed_boxes=441,retained_point_witnesses=sum(r['point_count'] for r in rows),discovered_rank_gains=sum(r['discovered_rank_gain'] for r in rows),retained26_rank_gains=sum(r['discovered_rank_gain'] for r in rows if r['parent']=='X948'),total_stage_seconds=total,relocation_interruption_seconds=interruption['wall_seconds'],certificate_continuation='Five certificate-only replays after repository relocation; original interrupted ledgers preserved. No point worker rerun.',score_worker_seconds=score['total_score_seconds'],inputs=inputs,source_sha256=cert.hashed(Path(__file__)),boundary='Complete fixed441-box portfolio. Six new-parent fibres share a64-bit j-height band, but their certified initial section spans differ17/17/17/17/15/16, so the small score-stratum comparison is descriptive, not causal or a calibrated predictor. Exact discoveries and completed point/verification exposure are recorded separately from scores and point counts. Three retained26 curves close prior source-centre gaps. Selection uses no validation primes or public exceptional points; catalogue matching is post-search against pinned630 equations and201 inventory only. No exact whole-curve rank, literature-wide novelty, record claim or automatic following wave.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();r=expected()
 if a.check:assert cert.read(OUT)==r
 else:assert not OUT.exists();checkpoint(OUT,r)
 print('PASS441 boxes',r['discovered_rank_gains'],'directions;',r['total_stage_seconds'],'supervised seconds',flush=True)
