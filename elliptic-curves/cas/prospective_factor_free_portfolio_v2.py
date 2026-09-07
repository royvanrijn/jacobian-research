#!/usr/bin/env python3
"""Two retained26 curves: fixed own-subgroup exposure after a source-gap audit."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BATCH=ROOT/'artifacts/local/elliptic-curves/prospective-factor-free-portfolio-v2';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['prospective_factor_free_portfolio_v2.py','prepare_factor_free_portfolio_v2.sage']}}
def freeze():
 assert not (BATCH/'protocol.json').exists()
 gate=ART/'blind_factor_free_28_control_v1.json';g=cert.read(gate)
 assert g['status']=='PASS' and g['rank_lower_bound']==28 and g['completed_boxes']==49
 index=cert.read(ART/'new_high_rank_curve_index_v22.json');audit=cert.read(ART/'factor_free_rank27_box_audit_v1.json');scores=cert.read(ART/'retained26_gap_scores_v1.json');gaps=cert.read(ART/'retained26_source_gap_v1.json')
 paths=[gate,ART/'new_high_rank_curve_index_v22.json',ART/'factor_free_rank27_box_audit_v1.json',ART/'retained26_gap_scores_v1.json',ART/'retained26_source_gap_v1.json',ART/'retained26_gap_trial_report_v1.json'];rows=[];inputs={}
 def add(ident,seed,parent,role,extra=None):
  rank=len(seed['points']);proof=seed['rank_certificate']
  assert seed['points'][:len(seed['generic_points'])]==seed['generic_points']
  assert digest(checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime']))==digest(proof)
  dest=BATCH/ident/'seed.json';assert not dest.exists();checkpoint(dest,seed);inputs[str(dest.relative_to(ROOT))]=cert.hashed(dest)
  rows.append(dict(id=ident,family=seed['family'],parameter=seed['parameter'],parent=parent,role=role,initial_rank=rank,generic_dimension=len(seed['generic_points']),**(extra or {})))
 retained27=sorted([r for r in index['curves'] if r['rank_lower_bound']==27 and not r['current_catalogue_matches']],key=lambda r:r['id'])
 assert len(retained27)==7
 for r in retained27:
  prior=next(a for a in audit['rows'] if a['id']==r['id']);folder=ROOT/prior['folder'];seedpath=folder/'seed.json';oldmaps=folder/'maps.json';newmaps=ROOT/'artifacts/local/elliptic-curves/factor-free-rank27-box-audit-v1'/r['id']/'maps.json';oldprotocol=folder.parent/'protocol.json';paths += [seedpath,oldmaps,newmaps,oldprotocol]
  assert prior['counts']=={'PROVED_NEW_COORDINATE_EXPOSURE':49}
  add(r['id'],cert.read(seedpath),'X948','all retained unmatched rank27',dict(previous_maps=str(oldmaps.relative_to(ROOT)),factor_free_maps=str(newmaps.relative_to(ROOT)),sample_domain=cert.read(oldprotocol)['sample_domain']))
 # The prior pair already completed; fill precisely two other source gaps by frozen training order and distinct labels.
 completed=set(scores['selected_ids']);chosen=[];families=set()
 for ident in scores['ordering']:
  if ident in completed or ident not in gaps['eligible_ids']:continue
  r=next(r for r in index['curves'] if r['id']==ident)
  if r['family'] in families:continue
  chosen.append(r);families.add(r['family'])
  if len(chosen)==2:break
 assert [r['id'] for r in chosen]==['new-20260906-42','new-20260906-49']
 for r in chosen:
  source=ART/r['source_certificate'];paths.append(source);q=cert.read(source)['curves'][r['source_curve_index']]
  assert r['rank_lower_bound']==26 and not r['current_catalogue_matches']
  add(r['id'],{k:q[k] for k in ['family','parameter','curve','points','generic_points','rank_certificate']},'X948','remaining rank26 gaps by frozen training score')
 # Higher certified subgroup from previous completed new-parent pilot; no parameter generation.
 for ident,old,cloudname,rank,parent in [('mestre-u11-unit','mestre-parent-calibration-v1/u11-unit','mestre_parent_calibration_u11_unit_mod2_v1.json',12,'Mestre u=11'),('kihara-v3-unit','kihara-positive-point-pilot-v1/v3-unit','kihara_positive_v3_unit_cloud_v1.json',14,'positive Kihara v=3')]:
  seedpath=ROOT/'artifacts/local/elliptic-curves'/old/'seed.json';cloudpath=ART/cloudname;paths += [seedpath,cloudpath];oldseed=cert.read(seedpath);cloud=cert.read(cloudpath)
  if 'independent_points' in cloud:points=cloud['independent_points']
  else:points=[cloud['points'][i] for i in cloud['audits'][0]['independent_indices']]
  assert len(points)==rank and points[:len(oldseed['points'])]==oldseed['points']
  # Finite proof only: no new point search or score, at primes <=997.
  from research_runtime.search_state import raw_state
  from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
  from research_runtime.memory_store import MemoryFactStore
  cache=QuotientOnlyReductionCache(MemoryFactStore());model=tuple(map(cert.F,oldseed['curve']));pts=tuple(tuple(map(cert.F,P)) for P in points);state=raw_state(model,pts,cache=cache,prime_bound=1000)
  proof=checked_rank(model,pts,state.reductions.primes,state.no_two_torsion_prime)
  seed=dict(family=oldseed['family'],parameter=oldseed['parameter'],curve=oldseed['curve'],points=points,generic_points=oldseed['points'],rank_certificate=proof)
  add(ident,seed,parent,'strongest updated own subgroup in each new-parent construction; Mestre tie broken by smaller u')
 assert len(rows)==11
 inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in paths})
 checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.prospective-factor-free-portfolio.v2',sources=sources(),inputs=inputs,rows=rows,sample_size=2048,sample_domain='full11952-specialized-followup-v1',charts=49,height=125000,seconds_per_chart=10,rank_stop=False,target_rank=32,maximum_point_boxes=539,geometry_wall_seconds=180,worker_wall_seconds=600,replay_wall_seconds=300,rss_bytes=2147483648,maximum_workers=1,gp_sha256=cert.hashed(Path('/usr/bin/gp')),
 previous_attempt='prospective-factor-free-portfolio-v1: preparation rejected six-roster sampling-domain mismatch, zero point boxes. V2 explicitly preserves each audited roster domain.',gate='Full49 original27-only factor-free blind control recovered28 and passed independent exact mod2/3/5, point and map replays. Seven retained27 factor-free49 rosters have proved new finite-coordinate exposure. Two unexposed26 curves are selected by existing training scores. Two new-parent fibres use improved own14/12 subgroups outside their earlier centre spans.',
 selection='All seven unmatched rank27; two highest remaining rank26 source gaps from different fibration labels after the previously completed pair; strongest updated own subgroup in the Mestre and positive Kihara constructions, with smaller u breaking the Mestre tie. No public target point, new parameter population, validation score, outcome stopping, adaptive centre reselection or refill.',
 centre_policy='2048 distinct SHA masks nonzero above original generic17 or the earlier new-parent seed (Mestre11/Kihara12).384-bit canonical heights rounded at10^6; LLL/CVP, exact norm/parity transport;49 largest computed norms. Reuse audited factor-free maps for retained27; prepare new maps for26 and updated parents. All539 maps freeze before points. The new-parent generic_points field names the inherited prior specialized subgroup, not a generic-function-field claim.',
 failure_policy='Any failed/censored geometry prevents all point search; preserve stage evidence and stop on worker/replay failure. All539 boxes planned, no rank stop or automatic subsequent wave.',following_campaign=None))
 print('FROZEN11 curves,3 parent surfaces,539 maximum boxes',[(r['id'],r['initial_rank']) for r in rows],flush=True)
def protocol():
 p=cert.read(BATCH/'protocol.json');assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p
def masks(p):
 result=[];i=0
 while len(result)<p['sample_size']:
  m=int(digest([ROW.get('sample_domain',p['sample_domain']),i]),16)%(1<<ROW['initial_rank']);i+=1
  if m>>ROW['generic_dimension'] and m not in result:result.append(m)
 return result
def configure(index):
 global ROW,D,SEED
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);s=run([SAGE,str(CAS/'prepare_factor_free_portfolio_v2.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
  ok=s['outcome']=='completed' and s['returncode']==0;ledger['maps'].append({'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger);print(row['id'],'maps',s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('all maps required before points')
 ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);entry={'id':row['id'],'status':'RUNNING','stages':[]};ledger['rows'].append(entry);checkpoint(out,ledger)
  for stage in ['worker','replay']:
   s=run([sys.executable,str(Path(__file__).resolve()),stage,'--index',str(i)],limits=Limits(p[stage+'_wall_seconds'],p['rss_bytes']),log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'),cwd=ROOT)
   ok=s['outcome']=='completed' and s['returncode']==0;entry['stages'].append({'name':stage,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger);print(row['id'],stage,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
   if not ok:entry['status']=ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed/censored worker or replay')
  result=cert.read(D/'result.json');entry.update(status='PASS',rank_lower_bound=result['rank_lower_bound'],result_sha256=cert.hashed(D/'result.json'));checkpoint(out,ledger)
 ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','launch','worker','replay']);p.add_argument('--index',type=int);args=p.parse_args()
 if args.stage in ['worker','replay']:configure(args.index);getattr(engine,args.stage)()
 else:globals()[args.stage]()
