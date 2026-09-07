#!/usr/bin/env python3
"""Fixed six-fibre score strata plus the three remaining retained26 source gaps."""
import argparse,sys,shutil
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';L=ROOT/'artifacts/local/elliptic-curves';BATCH=L/'mixed-reduced-parent-exposure-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['mixed_reduced_parent_exposure.py','prepare_mixed_reduced_parent_seeds.sage','verify_mixed_reduced_parent_seeds.sage','prepare_mixed_reduced_parent_maps.sage']}}
def prepare():
 assert not (BATCH/'prepare-protocol.json').exists()
 select=L/'det1092-reduced-score-strata-v1/result.json';selected=cert.read(select);assert selected['status']=='PASS' and len(selected['selected'])==6
 source=ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json';parent=cert.read(source);chart=ART/'det1092_reduced_parameter_chart_v1/chart-search.json';index=cert.read(ART/'new_high_rank_curve_index_v22.json');gap=cert.read(ART/'retained26_source_gap_v1.json');scores=cert.read(ART/'retained26_gap_scores_v1.json');old=cert.read(ART/'prospective_factor_free_portfolio_v2.json')
 completed=set(scores['selected_ids'])|{r['id'] for r in old['rows']};remaining=[ident for ident in scores['ordering'] if ident in gap['eligible_ids'] and ident not in completed];assert remaining==['new-20260906-73','new-20260906-75','new-20260906-63']
 rows=[];retained={};paths=[select,source,chart,ART/'new_high_rank_curve_index_v22.json',ART/'retained26_source_gap_v1.json',ART/'retained26_gap_scores_v1.json',ART/'prospective_factor_free_portfolio_v2.json',ART/'blind_factor_free_28_control_v1.json']
 for r in selected['selected']:
  rows.append(dict(id='d1092-'+r['id'],family='det1092-reduced',parent='det1092',parameter=r['parameter'],original_parameter=r['original_parameter'],stratum=r['stratum'],score_model=r['model']))
 for ident in remaining:
  r=next(r for r in index['curves'] if r['id']==ident);assert r['rank_lower_bound']==26
  path=ART/r['source_certificate'];paths.append(path);v=cert.read(path)['curves'][r['source_curve_index']];retained[ident]={k:v[k] for k in ['curve','points','generic_points','rank_certificate']};assert len(v['points'])==26
  rows.append(dict(id=ident,family=r['family'],parent='X948',parameter=r['parameter'],stratum='retained26_gap'))
 inp=dict(parent=parent,parameter_matrix=cert.read(chart)['selected']['state']['parameter_matrix'],rows=rows,retained=retained);checkpoint(BATCH/'intake-input.json',inp)
 checkpoint(BATCH/'prepare-protocol.json',dict(sources={**sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}},input_sha256=cert.hashed(BATCH/'intake-input.json'),maximum_denominator=64,maximum_coefficient=64,seconds=180,rss_bytes=2147483648,workers=1,scope='Fixed nine curves. New fibres take an exact finite-independent subset of their specialized17 generic sections; prove every original section in its rational span using bounded proposals and exact identities. No rank-based replacement; deficient or failed intake stops all point searches. Retained26 inputs unchanged.'))
 print('PREPARED fixed nine curves on two parents',flush=True)
def freeze():
 assert not (BATCH/'protocol.json').exists();q=cert.read(BATCH/'prepare-protocol.json');assert all(cert.hashed(ROOT/n)==h for n,h in q['sources'].items());intake=cert.read(BATCH/'intake-result.json');assert intake['status']=='PASS' and cert.read(BATCH/'seed-replay/result.json')['status']=='PASS'
 rows=intake['rows'];paths=[BATCH/n for n in ['prepare-protocol.json','intake-input.json','intake-result.json','seed-replay/result.json']]+[BATCH/r['id']/'seed.json' for r in rows]
 checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.mixed-reduced-parent-exposure.v1',sources=sources(),inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},rows=rows,sample_size=2048,sample_domain='full11952-specialized-followup-v1',charts=49,height=125000,seconds_per_chart=10,rank_stop=False,target_rank=32,maximum_point_boxes=441,geometry_wall_seconds=180,worker_wall_seconds=1200,replay_wall_seconds=1200,rss_bytes=2147483648,maximum_workers=1,gp_sha256=cert.hashed(Path('/usr/bin/gp')),selection='Six unchanged source-equation-only fibres in the256..319-bit j-numerator band: two strongest training scores, two moderate and a fixed SHA lower-third sample; plus all three remaining retained26 source gaps by their previous frozen order. No validation scores, public target points, catalogue matches, refill or stopping outcomes enter selection.',centre_policy='2048 nonzero parities in the certified specialized section span for each new fibre; parities above generic17 for each retained26 curve.384-bit numerical metric, exact rounded norm/parity transport,49 largest computed norms. All441 factor-free maps freeze before points. Same mapping and point boxes as the successful full original27-only blind28 control. No metric optimality or covering claim.',gate='Reduced parent chart and full blind28 recovery independently verified; exact specialized seed spans prevent generic-rank assumptions. The prior two fixed fibres were not score-selected. This bounded comparison tests score scheduling at comparable arithmetic height, while closing retained high-rank point-exposure gaps.',failure_policy='Preserve all failed or censored stages. No refill, automatic budget extension or following wave.',following_campaign=None))
 print('FROZEN441 boxes',[(r['id'],r['initial_rank']) for r in rows],flush=True)
def protocol():
 p=cert.read(BATCH/'protocol.json');assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p
def masks(p):
 result=[];i=0
 while len(result)<p['sample_size']:
  m=int(digest([p['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
  if m>>ROW['mask_floor'] and m not in result:result.append(m)
 return result
def configure(index):
 global ROW,D,SEED
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks
def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[],'seed_hashes':{r['id']:cert.hashed(BATCH/r['id']/'seed.json') for r in p['rows']}};checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);s=run([SAGE,str(CAS/'prepare_mixed_reduced_parent_maps.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
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
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','freeze','launch','worker','replay']);p.add_argument('--index',type=int);args=p.parse_args()
 if args.stage in ['worker','replay']:configure(args.index);getattr(engine,args.stage)()
 else:globals()[args.stage]()
