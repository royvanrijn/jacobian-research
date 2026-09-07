#!/usr/bin/env python3
"""Two fixed non-anchor fibres of the newly certified determinant1092 MW17 parent."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BATCH=ROOT/'artifacts/local/elliptic-curves/det1092-point-pilot-v2';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
 return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['det1092_point_pilot_v2.py','prepare_det1092_point_pilot_v2.sage','verify_det1092_point_seeds.sage']}}
def freeze():
 assert not (BATCH/'protocol.json').exists()
 parent=ART/'curve302_recovered_mw17_parent_v1.json';gate=ART/'blind_factor_free_28_control_v1.json';p=cert.read(parent);g=cert.read(gate)
 authority=next(e for e in cert.read(ROOT/'MATH_STATUS.json')['entries'] if e['id']=='EC-CURVE302-RECOVERED-MW17-PARENT');assert authority['state']=='proved'
 assert int(p['height_determinant'])==1092
 assert g['status']=='PASS' and g['rank_lower_bound']==28 and g['completed_boxes']==49 and len(p['basis_weierstrass_coordinates'])==17
 # Export only the generic equation and seventeen sections. No public-group words, rank31 points, anchor parameter, or jump labels enter the worker input.
 redacted={k:p[k] for k in ['a_invariants','basis_weierstrass_coordinates']};redacted['source_sha256']=cert.hashed(parent)
 parentfile=BATCH/'parent-sections.json';assert not parentfile.exists();checkpoint(parentfile,redacted)
 rows=[dict(id='unit',family='det1092-MW17',parameter='1',initial_rank=17,generic_dimension=17,mask_floor=0),dict(id='outer',family='det1092-MW17',parameter='1009/101',initial_rank=17,generic_dimension=17,mask_floor=0)]
 checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.det1092-point-pilot.v2',sources=sources(),inputs={str(x.relative_to(ROOT)):cert.hashed(x) for x in [parentfile,gate]},rows=rows,sample_size=2048,sample_domain='full11952-specialized-followup-v1',charts=49,height=125000,seconds_per_chart=10,rank_stop=False,target_rank=32,maximum_point_boxes=98,geometry_wall_seconds=180,worker_wall_seconds=600,replay_wall_seconds=300,rss_bytes=2147483648,maximum_workers=1,gp_sha256=cert.hashed(Path('/usr/bin/gp')),
 previous_attempt='det1092-point-pilot-v1: zero point boxes; both fixed seeds have finite rank17. Intake compared string coordinates with Fraction coordinates and rejected the representation. V2 compares exact rational values; no curve, point or exposure rule changes.',gate='New parallel exact determinant1092 arithmetic MW17 parent with explicit generic sections; full-roster original27-only factor-free known28 regression passed. No extra parameter population, lattice search or recovered record point is required.',
 selection='Two fixed fibres T=1 and1009/101, matching the earlier new-parent calibration coordinates and excluding the known anchor0. No arithmetic score, public exceptional point, anchor jump label, outcome-based replacement or rank stopping. Parent equation and generic sections only enter intake and execution.',
 centre_policy='2048 distinct nonzero SHA parities in the specialized17-point generic subgroup,384-bit numerical heights rounded at10^6, LLL/CVP exact parity/norm transport,49 largest computed norms. Both49-map rosters freeze before points. Since these are new fibres with no additional own directions yet, mask_floor=0; generic_dimension remains17. Same factor-free Gauss/hyperellred mapping and equal boxes as the successful control.',
 failure_policy='Any failed/censored intake or map prevents all point attempts. Preserve failures; no refill, larger parameter scan or automatic following wave.',following_campaign=None))
 print('FROZEN two non-anchor determinant1092 fibres;98 maximum boxes',flush=True)
def intake():
 p=protocol();parent=cert.read(BATCH/'parent-sections.json')
 from research_runtime.search_state import raw_state
 from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
 from research_runtime.memory_store import MemoryFactStore
 def polynomial(c,t):
  v=cert.F(0)
  for a in reversed(c):v=v*t+cert.F(a)
  return v
 def value(r,t):return polynomial(r['numerator'],t)/polynomial(r['denominator'],t)
 for row in p['rows']:
  t=cert.F(row['parameter']);assert t!=0
  a=[value(r,t) for r in parent['a_invariants']];a1,a2,a3,a4,a6=a;b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
  model=(cert.F(0),cert.F(0),cert.F(0),-c4/48,-c6/864);images=[tuple(value(r,t) for r in P) for P in parent['basis_weierstrass_coordinates']]
  assert all(cert.is_on_weierstrass_curve(a,P) for P in images)
  points=[(x+b2/12,y+(a1*x+a3)/2) for x,y in images];cache=QuotientOnlyReductionCache(MemoryFactStore());state=raw_state(model,points,cache=cache,prime_bound=1000)
  assert state.rank==17 and [tuple(map(cert.F,P)) for P in state.basis]==points
  proof=checked_rank(model,points,state.reductions.primes,state.no_two_torsion_prime)
  seed=dict(family=row['family'],parameter=row['parameter'],curve=list(map(str,model)),points=[list(map(str,P)) for P in points],generic_points=[list(map(str,P)) for P in points],rank_certificate=proof,long_curve=list(map(str,a)),long_points=[list(map(str,P)) for P in images],protocol_hash=digest(p),model_coefficient_bits=max(max(abs(v.numerator).bit_length(),v.denominator.bit_length()) for v in model))
  dest=BATCH/row['id']/'seed.json';assert not dest.exists();checkpoint(dest,seed);print(row['id'],'exact17 seed; bits',seed['model_coefficient_bits'],flush=True)
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
 ROW=protocol()['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json';
 if (BATCH/'ledger.json').exists():assert cert.hashed(SEED)==cert.read(BATCH/'ledger.json')['seed_hashes'][ROW['id']]
 engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def launch():
 p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[],'seed_hashes':{r['id']:cert.hashed(BATCH/r['id']/'seed.json') for r in p['rows']}};checkpoint(out,ledger)
 for i,row in enumerate(p['rows']):
  configure(i);s=run([SAGE,str(CAS/'prepare_det1092_point_pilot_v2.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT)
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
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','intake','launch','worker','replay']);p.add_argument('--index',type=int);args=p.parse_args()
 if args.stage in ['worker','replay']:configure(args.index);getattr(engine,args.stage)()
 else:globals()[args.stage]()
