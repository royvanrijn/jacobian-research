#!/usr/bin/env python3
"""Four fixed curves; restricted-prime coordinate refinement with history gate."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';BATCH=LOCAL/'bounded-prime-point-portfolio-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
PRIMES=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['bounded_prime_point_portfolio.py','prepare_bounded_prime_point_portfolio.sage','bounded_prime_pari_mapping_v2.sage','research_runtime/projective_box_change.py']}}
def freeze():
 assert not (BATCH/'protocol.json').exists();gate=ART/'blind_bounded_prime_28_control_v2.json';g=cert.read(gate);assert g['status']=='PASS' and g['rank_lower_bound']==28 and g['completed_boxes']==49
 entries=[('new-20260906-186','prospective-factor-free-portfolio-v2/new-20260906-186',17,'full11952-specialized-followup-v1','full11952-specialized-followup-v1/new-20260906-186',1000000),('new-20260906-90','prospective-factor-free-portfolio-v2/new-20260906-90',16,'new27-specialized-parity-v1','new27-specialized-parity-six-v1/new-20260906-90',125000),('det1092-unit','det1092-point-pilot-v2/unit',0,'full11952-specialized-followup-v1',None,None),('det1092-outer','det1092-point-pilot-v2/outer',0,'full11952-specialized-followup-v1',None,None)]
 inputs={str(gate.relative_to(ROOT)):cert.hashed(gate)};rows=[]
 for ident,path,floor,domain,prior,height in entries:
  folder=LOCAL/path;seed=cert.read(folder/'seed.json');proof=seed['rank_certificate'];assert digest(checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime']))==digest(proof)
  assert seed['points'][:len(seed['generic_points'])]==seed['generic_points'];dest=BATCH/ident/'seed.json';assert not dest.exists();checkpoint(dest,seed)
  paths=[folder/'seed.json',folder/'maps.json',dest];history=[]
  if prior:
   old=LOCAL/prior/'maps.json';paths.append(old)
   result=LOCAL/('full11952-million-height-v1/new-20260906-186/result.json' if ident.endswith('186') else prior+'/result.json');paths.append(result);r=cert.read(result)
   assert len(r['charts'])==49 and r['maps_sha256']==cert.hashed(old) and all(c['search']['status']=='bounded_search_complete' and c['search']['height_bound']==height for c in r['charts'])
   history.append(dict(maps_path=str(old.relative_to(ROOT)),height=height))
  inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in paths})
  rows.append(dict(id=ident,family=seed['family'],parameter=seed['parameter'],initial_rank=len(seed['points']),generic_dimension=len(seed['generic_points']),mask_floor=floor,sample_domain=domain,base_maps=str((folder/'maps.json').relative_to(ROOT)),prior_same_centre_boxes=history))
 checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.bounded-prime-point-portfolio.v1',sources=sources(),inputs=inputs,rows=rows,sample_size=2048,sample_domain='full11952-specialized-followup-v1',charts=49,height=125000,seconds_per_chart=10,minimization_primes=PRIMES,rank_stop=False,target_rank=32,maximum_point_boxes=196,maximum_workers=1,rss_bytes=2147483648,geometry_wall_seconds=180,worker_wall_seconds=600,replay_wall_seconds=300,gp_sha256=cert.hashed(Path('/usr/bin/gp')),
 selection='Fixed retained27 representatives from R17 (ID186) and MW16 (ID90), and both previously fixed determinant1092 fibres. No score, public exceptional point, validation statistic, new parameter population or outcome refill. All old masks and centres retained. Geometry-only skip: skip a curve if no new-coordinate witness is proved outside every declared same-centre prior box.',
 history_scope='Compare each refined chart with its completed factor-free125000 box and, for retained27, the completed corresponding global-minimal chart at1000000 for ID186 or125000 for ID90. This is a declared same-centre history comparison, not an exhaustive proof against every earlier different-centre chart or rational translation.',
 policy='Explicit25-prime minimisation after the factor-free map, then hyperellred; no full discriminant factorization requested. All196 refined maps freeze before points. A curve with at least one proved new-history coordinate witness receives its complete49-chart roster; counts of new, duplicate and unresolved boxes remain explicit. Otherwise zero boxes and no refill. Exact geometry/history and whole-cloud independent proofs precede promotion.',following_campaign=None))
 print('FROZEN four curves;196 maximum boxes; exact history gate',flush=True)
def protocol():
 p=cert.read(BATCH/'protocol.json');assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p
def configure(index):
 global ROW,D,SEED
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks
def masks(p):
 values=[];i=0
 while len(values)<p['sample_size']:
  m=int(digest([ROW['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
  if m>>ROW['mask_floor'] and m not in values:values.append(m)
 return values
def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger=dict(status='RUNNING_GEOMETRY',maps=[],rows=[]);checkpoint(out,ledger)
 def stage(name,cmd,seconds):
  s=run(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT);ok=s['outcome']=='completed' and s['returncode']==0;print(ROW['id'],name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger.update(status='FAILED_OR_CENSORED',failed=dict(id=ROW['id'],name=name,supervision=s));checkpoint(out,ledger);raise ArithmeticError('preserve failed stage')
  return s
 for i,row in enumerate(p['rows']):
  configure(i);s=stage('maps',[SAGE,str(CAS/'prepare_bounded_prime_point_portfolio.sage'),'--index',str(i)],180);maps=cert.read(D/'maps.json');ledger['maps'].append(dict(id=row['id'],status='PASS',counts=maps['history_counts'],supervision=s));checkpoint(out,ledger)
 ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);maps=cert.read(D/'maps.json');n=maps['history_counts'].get('PROVED_NEW_VS_DECLARED_HISTORY',0)
  if n==0:ledger['rows'].append(dict(id=row['id'],status='SKIPPED_NO_PROVED_NEW_HISTORY_BOX',stages=[],point_boxes=0));checkpoint(out,ledger);print(row['id'],'SKIPPED no new history exposure',flush=True);continue
  entry=dict(id=row['id'],status='RUNNING',stages=[]);ledger['rows'].append(entry);checkpoint(out,ledger)
  for name,seconds in [('worker',600),('replay',300)]:
   s=stage(name,[sys.executable,str(Path(__file__)),name,'--index',str(i)],seconds);entry['stages'].append(dict(name=name,status='PASS',supervision=s));checkpoint(out,ledger)
  result=cert.read(D/'result.json');entry.update(status='PASS',point_boxes=49,rank_lower_bound=result['rank_lower_bound']);checkpoint(out,ledger)
 ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','launch','worker','replay']);parser.add_argument('--index',type=int);args=parser.parse_args()
 if args.stage in ['worker','replay']:configure(args.index);getattr(engine,args.stage)()
 else:globals()[args.stage]()
