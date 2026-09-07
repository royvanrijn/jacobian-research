#!/usr/bin/env python3
"""Two retained26 curves: fixed own-subgroup exposure after a source-gap audit."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BATCH=ROOT/'artifacts/local/elliptic-curves/retained26-gap-trial-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['retained26_gap_trial.py','prepare_retained26_gap.sage']}}
def freeze():
 assert not (BATCH/'protocol.json').exists();scores=cert.read(ART/'retained26_gap_scores_v1.json');gaps=cert.read(ART/'retained26_source_gap_v1.json');index=cert.read(ART/'new_high_rank_curve_index_v22.json')
 assert scores['status']==gaps['status']=='PASS' and len(scores['selected_ids'])==2
 paths=[ART/n for n in ['retained26_gap_scores_v1.json','retained26_source_gap_v1.json','new_high_rank_curve_index_v22.json','factor_free_known28_control_v1.json']];rows=[];inputs={}
 gate=cert.read(paths[-1]);assert gate['status']=='PASS' and gate['rank_lower_bound']==28
 for ident in scores['selected_ids']:
  assert ident in gaps['eligible_ids'];r=next(r for r in index['curves'] if r['id']==ident);source=ART/r['source_certificate'];paths.append(source);q=cert.read(source)['curves'][r['source_curve_index']]
  assert r['rank_lower_bound']==26 and not r['current_catalogue_matches'] and len(q['points'])==26 and q['points'][:17]==q['generic_points']
  proof=q['rank_certificate'];assert digest(checked_rank(tuple(map(cert.F,q['curve'])),[tuple(map(cert.F,p)) for p in q['points']],[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime']))==digest(proof)
  seed={k:q[k] for k in ['family','parameter','curve','points','generic_points','rank_certificate']};path=BATCH/ident/'seed.json';assert not path.exists();checkpoint(path,seed);inputs[str(path.relative_to(ROOT))]=cert.hashed(path)
  rows.append({'id':ident,'family':r['family'],'parameter':r['parameter'],'initial_rank':26,'generic_dimension':17,'source_certificate':r['source_certificate'],'source_curve_index':r['source_curve_index']})
 assert len({r['family'] for r in rows})==2
 inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in paths})
 checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.retained26-gap-trial.v1','sources':sources(),'inputs':inputs,'rows':rows,'sample_size':2048,'sample_domain':'full11952-specialized-followup-v1','charts':49,'height':125000,'seconds_per_chart':10,'rank_stop':False,'target_rank':32,'maximum_point_boxes':98,'geometry_wall_seconds':180,'worker_wall_seconds':600,'replay_wall_seconds':300,'rss_bytes':2147483648,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),
 'gate':'Seven retained rank26 source runs used only17 generic directions. Exclude all eleven retained26 curves with later own26 records in the bounded layout audit. The completed factor-free known28 recovery is a feasibility control, not a yield forecast. Test two highest training-score curves from different fibration labels within the seven.','selection':'Exactly the two preselected retained equations and their existing26-point seeds. No new parameter scan, record input, validation-prime input, outcome stopping, refill or next wave. This is not a parent-diversity or score-efficacy experiment.',
 'centre_policy':'Calibrated2048 distinct SHA256 masks in the26-point group, nonzero above generic17.384-bit heights rounded at10^6; unimodular LLL and numerical CVP;49 largest computed norms with exact parity/norm checks. Factor-free Gauss/hyperellred maps. Both complete map sets must freeze before any point search. New parity classes lie outside the source generic17 subgroup, whose finite-mod2 embedding is injective. No CVP optimality claim.',
 'failure_policy':'Preserve any failed or censored stage and stop; no retry, refill or following campaign. All98 boxes attempted if preparation and workers complete.','following_campaign':None})
 print('FROZEN two retained26 source gaps;98 maximum boxes',rows,flush=True)
def protocol():
 p=cert.read(BATCH/'protocol.json');assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p
def masks(p):
 result=[];i=0
 while len(result)<p['sample_size']:
  m=int(digest([p['sample_domain'],i]),16)%(1<<26);i+=1
  if m>>17 and m not in result:result.append(m)
 return result
def configure(index):
 global ROW,D,SEED
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);s=run([SAGE,str(CAS/'prepare_retained26_gap.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
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
