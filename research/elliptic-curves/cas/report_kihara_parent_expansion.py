#!/usr/bin/env python3
"""Record completed Kihara exposure, portable replay and pinned comparisons."""
import argparse,json,hashlib
from pathlib import Path
import kihara_fresh_point_pilot as batch
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
ROOT=batch.ROOT;ART=batch.ART;LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'kihara_parent_expansion_report_v1.json'
def compute():
 d=batch.BATCH;p=batch.protocol();ledger=cert.read(d/'ledger.json');verified=cert.read(d/'verification-ledger.json')
 assert ledger['status']==verified['status']=='PASS' and len(ledger['rows'])==len(verified['rows'])==4
 dbpath=catalogue.D/'database.json';meta=catalogue.D/'metadata.json';db=cert.read(dbpath);index=ART/'new_high_rank_curve_index_v22.json';inventory=cert.read(index)
 assert cert.hashed(dbpath)==cert.read(meta)['sha256'] and len(db['curves'])==620 and len(inventory['curves'])==201
 paths=[Path(__file__).resolve(),Path(batch.__file__),dbpath,meta,index,d/'protocol.json',d/'ledger.json',d/'verification-ledger.json',ART/'kihara_fresh_point_pilot_sage_replay_v1.json',ART/'kihara_five_parent_distinctness_v1.json']
 byj={};local={}
 for r in db['curves']:byj.setdefault(catalogue.j(r['ainvs']),[]).append(r)
 for r in inventory['curves']:local.setdefault(catalogue.j(r['curve']),[]).append(r)
 def supervision(path):
  r=cert.read(path);assert r['outcome']=='completed' and r['returncode']==0 and r['failure_reason'] is None
  log=Path(r['log']);assert cert.hashed(log)==r['log_sha256'];paths.extend([path,log]);return r['wall_seconds']
 rows=[]
 for i,row in enumerate(p['rows']):
  folder=d/row['id'];seedpath=folder/'seed.json';resultpath=folder/'result.json';seed=cert.read(seedpath);result=cert.read(resultpath);v=verified['rows'][i]
  assert v['id']==row['id'] and result['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(result['charts'])==49
  assert ledger['rows'][i]['result_sha256']==cert.hashed(resultpath)
  cloudpath=ROOT/v['mod2_certificate'];oddpath=ROOT/v['modl_certificate'];cloud=cert.read(cloudpath);odd=cert.read(oddpath)
  assert cert.hashed(cloudpath)==v['mod2_sha256'] and cert.hashed(oddpath)==v['modl_sha256']
  assert cloud['rank_lower_bound']==result['rank_lower_bound']==14 and len(cloud['points'])==14
  assert cloud['curve']==seed['curve'] and cloud['points']==seed['points']
  assert all(a['finite_column_rank']==14 for a in odd['audits'])
  times={stage:supervision(folder/(stage+'.supervisor.json')) for stage in ['maps','worker','replay','mod2-build','mod2-check','modl-build','modl-check']}
  for version in [1,2]:times['seed_intake_v'+str(version)]=supervision(LOCAL/f'kihara-fresh-fibres-v{version}'/f'fibre{i}-supervisor.json')
  oldpath=LOCAL/'kihara-fresh-fibres-v1'/row['id']/'seed.json';old=cert.read(oldpath)
  j=catalogue.j(seed['curve']);completed=sum(c['search']['status']=='bounded_search_complete' for c in result['charts'])
  rows.append({**row,'curve':seed['curve'],'old_model_bits':old['model_coefficient_bits'],'model_bits':seed['model_coefficient_bits'],'old_seed_index_in_corrected':6,
   'rank_lower_bound':14,'discovered_rank_gains':0,'retained_points':14,'completed_boxes':completed,'attempted_boxes':49,'incomplete_boxes':49-completed,
   'stage_seconds':times,'gp_cpu_milliseconds':sum(c['search']['search_cpu_ms'] or 0 for c in result['charts']),
   'catalogue_matches':[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])],
   'inventory_matches':[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])]})
  paths.extend([seedpath,resultpath,folder/'maps.json',oldpath,cloudpath,oddpath])
 shared={name:supervision(d/(name+'-supervisor.json')) for name in ['geometry','independent-rank']}
 # Cloud wrapper contains the four per-row audit jobs above: do not add it twice.
 cloud_wrapper=supervision(d/'cloud-audits-supervisor.json')
 for name in ['global','fixed']:shared['parent_geometry_'+name]=supervision(LOCAL/'kihara-parent-geometry-v1'/(name+'-supervisor.json'))
 shared['parent_distinctness']=supervision(LOCAL/'kihara-parent-distinctness-v1'/'supervisor.json')
 replay=LOCAL/'kihara-parent-standalone-v1';shared['standalone_parent_replay']=supervision(replay/'supervisor.json');transcript=(replay/'replay.log').read_text()
 assert transcript.splitlines()[-1]=='PASS four new Q-distinct parents; index6 seeds; 44/45 comparisons; 10 new good-prime counts'
 for path in [ROOT/'elliptic-curves/cas/verify_kihara_parent_bundle.sage',ART/'kihara_parent_replay_bundle_v1.json']:
  assert path.read_bytes()==(replay/path.name).read_bytes();paths.append(path)
 parents=cert.read(ART/'kihara_five_parent_distinctness_v1.json');assert sum(r['status']=='PROVED_Q_DISTINCT' for r in parents['pairwise_separations'])==44
 assert len(set(catalogue.j(row['curve']) for row in rows))==4
 total=sum(sum(r['stage_seconds'].values()) for r in rows)+sum(shared.values())
 return {'schema':'elliptic-curves.kihara-parent-expansion-report.v1','status':'PASS','rows':rows,'completed_boxes':sum(r['completed_boxes'] for r in rows),'discovered_rank_gains':0,
  'new_Q_distinct_tested_parents':4,'tested_Q_distinct_parent_lower_bound':11,'control_additional_parent_status':'UNKNOWN relative to t=5/2',
  'shared_supervised_seconds':shared,'supervised_seconds_including_failed_v1_seed_intakes_and_independent_audits':total,
  'cloud_wrapper_seconds_not_added_again':cloud_wrapper,'standalone_replay_transcript':transcript,
  'sources':{str(path.relative_to(ROOT)):cert.hashed(path) for path in paths},
  'scope':'Four fixed Kihara path fibres, each from a new Q-distinct K3 parent, tested with identical49-box limits. All models, seeds, histories, full retained clouds and independent finite ranks replay. Index-six seed correction adds no direction. Lower bounds14 only; no exact-rank, complete-saturation, parent-superiority, new geometric NS type or near-record claim. Pinned620 catalogue and201 inventory comparisons are post-search and do not establish literature-wide novelty. No score input or validation-prime selection, and no automatic follow-up campaign.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=compute()
 if a.check:assert r==cert.read(OUT)
 else:
  with OUT.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
 print('PASS',r['completed_boxes'],'boxes;0 gains; four new Q-distinct parents;',r['supervised_seconds_including_failed_v1_seed_intakes_and_independent_audits'],'seconds')
 print('Pinned comparisons:',[(x['parameter'],x['catalogue_matches'],x['inventory_matches']) for x in r['rows']])
