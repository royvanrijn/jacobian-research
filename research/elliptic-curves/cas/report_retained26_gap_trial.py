#!/usr/bin/env python3
"""Bind completed retained26 exposure, full-cloud proofs and portable replay."""
import argparse
from pathlib import Path
import retained26_gap_trial as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
OUT=batch.ART/'retained26_gap_trial_report_v1.json'
def compute():
 p=batch.protocol();d=batch.BATCH;ledger=cert.read(d/'ledger.json');audit=cert.read(d/'verification-ledger.json');post=cert.read(d/'post-verification-ledger.json');score=cert.read(batch.ART/'retained26_gap_scores_v1.json');gap=cert.read(batch.ART/'retained26_source_gap_v1.json')
 assert all(r['status']=='PASS' for r in [ledger,audit,post,score,gap]) and len(ledger['rows'])==len(audit['rows'])==2 and cert.hashed(Path('/usr/bin/gp'))==p['gp_sha256']
 assert [r['id'] for r in p['rows']]==score['selected_ids']
 paths=[Path(__file__).resolve(),Path(batch.__file__),d/'protocol.json',d/'ledger.json',d/'verification-ledger.json',d/'post-verification-ledger.json',batch.ART/'retained26_gap_scores_v1.json',batch.ART/'retained26_source_gap_v1.json'];rows=[]
 def duration(path):
  s=cert.read(path);assert s['outcome']=='completed' and s['returncode']==0 and s['failure_reason'] is None
  assert cert.hashed(Path(s['log']))==s['log_sha256'];paths.extend([path,Path(s['log'])]);return s['wall_seconds']
 for i,row in enumerate(p['rows']):
  folder=d/row['id'];seed=cert.read(folder/'seed.json');result=cert.read(folder/'result.json');maps=cert.read(folder/'maps.json');v=audit['rows'][i]
  assert v['id']==row['id'] and ledger['rows'][i]['result_sha256']==cert.hashed(folder/'result.json')
  assert result['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(result['charts'])==len(maps['rows'])==49 and len(maps['sample'])==2048
  cloudpath=batch.ROOT/v['mod2_certificate'];oddpath=batch.ROOT/v['modl_certificate'];cloud=cert.read(cloudpath);odd=cert.read(oddpath)
  assert cert.hashed(cloudpath)==v['mod2_sha256'] and cert.hashed(oddpath)==v['modl_sha256']
  lower=cloud['rank_lower_bound'];assert lower>=result['rank_lower_bound'] and all(a['finite_column_rank']==lower for a in odd['audits'])
  assert all(c['parity']>>17 for c in maps['centres']) and len({c['parity'] for c in maps['centres']})==49
  times={name:duration(folder/(name+'.supervisor.json')) for name in ['maps','worker','replay','mod2-build','mod2-check','modl-build','modl-check']}
  completed=sum(c['search']['status']=='bounded_search_complete' for c in result['charts']);prior=next(r for r in gap['rows'] if r['id']==row['id']);sc=next(r for r in score['rows'] if r['id']==row['id'])
  rows.append({**row,'training_score_units':sc['score_units'],'initial_rank':26,'rank_lower_bound':lower,'discovered_rank_gains':lower-26,'retained_points':len(cloud['points']),'modular_lower_bounds':{'2':lower,**v['odd_modulus_lower_bounds']},
   'source_generic_centre_boxes':prior['completed_source_boxes'],'attempted_new_own26_boxes':49,'completed_new_own26_boxes':completed,'incomplete_new_boxes':49-completed,'stage_seconds':times,'gp_cpu_ms':sum(c['search']['search_cpu_ms'] or 0 for c in result['charts']),
   'rank_certificate':v['mod2_certificate']})
  paths.extend([folder/'seed.json',folder/'result.json',folder/'maps.json',cloudpath,oddpath])
 fresh=d/'standalone';files=[batch.CAS/'verify_retained26_gap_rank.sage']+[batch.ROOT/r['rank_certificate'] for r in rows]
 for path in files:assert path.read_bytes()==(fresh/path.name).read_bytes();paths.append(path)
 independent=duration(fresh/'supervisor.json');geometry=duration(d/'geometry.supervisor.json');wrapper=duration(d/'cloud-audits.supervisor.json');transcript=(fresh/'replay.log').read_text()
 expected='PASS independent retained26 full clouds '+str([(r['retained_points'],r['rank_lower_bound']) for r in rows]);assert transcript.strip()==expected
 total=score['total_worker_seconds']+sum(sum(r['stage_seconds'].values()) for r in rows)+independent+geometry
 return {'schema':'elliptic-curves.retained26-gap-trial-report.v1','status':'PASS','rows':rows,'eligible_source_gaps':gap['eligible_ids'],'selected_ids':score['selected_ids'],'remaining_source_gap_ids':[i for i in gap['eligible_ids'] if i not in score['selected_ids']],
  'completed_boxes':sum(r['completed_new_own26_boxes'] for r in rows),'discovered_rank_gains':sum(r['discovered_rank_gains'] for r in rows),'score_worker_seconds':score['total_worker_seconds'],'geometry_replay_seconds':geometry,'standalone_rank_seconds':independent,'cloud_audit_wrapper_seconds_not_added_twice':wrapper,'total_supervised_seconds_including_score_and_independent_checks':total,'standalone_transcript':transcript,
  'sources':{str(path.relative_to(batch.ROOT)):cert.hashed(path) for path in paths},
  'scope':'Two fixed already-retained rank26 curves selected by training score with distinct fibration labels from seven exact source-prefix gaps. All98 own26 boxes completed; new centre parities are outside the old generic17 image, certified by the finite-mod2 independent26 seed. Full point clouds, odd-prime quotient bounds, exact maps and independent complete groups replay. This is not a new parameter scan, parent expansion, randomized score efficacy trial, exact-rank bound or proof of point absence. Eleven previously exposed26 curves were excluded; five source gaps remain unscheduled. No automatic next wave or inventory change.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();r=compute()
 if args.check:assert r==cert.read(OUT)
 else:
  assert not OUT.exists();checkpoint(OUT,r)
 print('PASS',r['completed_boxes'],'completed boxes;',r['discovered_rank_gains'],'gains;',r['total_supervised_seconds_including_score_and_independent_checks'],'seconds')
if False:pass
