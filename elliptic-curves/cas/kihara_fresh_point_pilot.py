#!/usr/bin/env python3
"""Fixed four-fibre point visibility test on the index-six enlarged Kihara seeds."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BATCH=ROOT/'artifacts/local/elliptic-curves/kihara-fresh-point-pilot-v1';INTAKE=ROOT/'artifacts/local/elliptic-curves/kihara-fresh-fibres-v2';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['kihara_fresh_point_pilot.py','prepare_kihara_fresh_point_pilot.sage','intake_kihara_fresh_fibres_v2.sage']}}
def freeze():
 assert not (BATCH/'protocol.json').exists();rows=[];inputs={}
 for i,t in enumerate(['3/2','5/2','11/3','1009/101']):
  source=INTAKE/('fibre'+str(i))/'seed.json';seed=cert.read(source);assert seed['status']=='PASS' and len(seed['points'])==14 and seed['parameter']==t
  dest=BATCH/('fibre'+str(i))/'seed.json';assert not dest.exists();checkpoint(dest,seed)
  inputs[str(source.relative_to(ROOT))]=cert.hashed(source);inputs[str(dest.relative_to(ROOT))]=cert.hashed(dest)
  rows.append({'id':'fibre'+str(i),'family':'kihara-rank14','parameter':t,'initial_rank':14})
 checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.kihara-fresh-point-pilot.v1','sources':sources(),'inputs':inputs,'rows':rows,'sample_size':2048,'sample_domain':'kihara-fresh-point-pilot-v1','charts':49,'height':125000,'seconds_per_chart':10,'rank_stop':False,'target_rank':32,'maximum_point_boxes':196,'geometry_wall_seconds':180,'worker_wall_seconds':600,'replay_wall_seconds':300,'rss_bytes':2147483648,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),
  'gate':'Four fixed fresh fibres of the Kihara rank>=14 path have exact14-point proofs after a proved index-six seed enlargement. One fibre has2264 coefficient bits and parameter1009/101. Test practical point visibility in this different construction before any population or score sweep.',
  'selection':'Exactly the four already retained parameters and corrected seeds. No scoring, public target, validation-prime input, replacement or rank stopping. No parent-superiority or generic-rank-only yield inference.',
  'centre_policy':'2048 distinct nonzero SHA256 parity masks in the own14 group;384-bit canonical heights rounded at10^6;unimodular LLL and numerical CVP with exact parity/norm checks;49 largest computed norms. Factor-free Gauss/hyperellred maps. All196 maps must complete before any point search.',
  'failure_policy':'Preserve all failed/censored stages and prior intake costs. No automatic retry, refill or following campaign. Any preparation failure prevents every point attempt in this four-fibre pilot.','following_campaign':None})
 print('FROZEN4 KIHARA FIBRES; MAX196 BOXES',flush=True)
def protocol():
 p=cert.read(BATCH/'protocol.json')
 if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in p['inputs'].items()):raise ArithmeticError('frozen Kihara pilot inputs changed')
 return p
def masks(p):
 result=[];i=0
 while len(result)<p['sample_size']:
  m=int(digest([p['sample_domain'],i]),16)%(1<<14);i+=1
  if m and m not in result:result.append(m)
 return result
def configure(index):
 global ROW,D,SEED
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks
def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);r=run([SAGE,str(CAS/'prepare_kihara_fresh_point_pilot.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
  ok=r['outcome']=='completed' and r['returncode']==0;ledger['maps'].append({'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(out,ledger);print(row['id'],'maps',r['outcome'],r['returncode'],r['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('all maps required before any point search')
 ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);entry={'id':row['id'],'status':'RUNNING','stages':[]};ledger['rows'].append(entry);checkpoint(out,ledger)
  for stage in ['worker','replay']:
   r=run([sys.executable,str(Path(__file__).resolve()),stage,'--index',str(i)],limits=Limits(p[stage+'_wall_seconds'],p['rss_bytes']),log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'),cwd=ROOT)
   ok=r['outcome']=='completed' and r['returncode']==0;entry['stages'].append({'name':stage,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(out,ledger);print(row['id'],stage,r['outcome'],r['returncode'],r['wall_seconds'],flush=True)
   if not ok:entry['status']=ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('point worker or exact history replay failed')
  result=cert.read(D/'result.json');entry.update(status='PASS',rank_lower_bound=result['rank_lower_bound'],result_sha256=cert.hashed(D/'result.json'));checkpoint(out,ledger)
 ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args()
 if a.stage in ('worker','replay'):configure(a.index);getattr(engine,a.stage)()
 else:globals()[a.stage]()
