#!/usr/bin/env python3
"""Source-guarded exact support for the every-order off-center cubic family."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import argparse,copy,hashlib,json
if not __debug__:raise SystemExit('Optimized execution is not accepted.')
HERE=Path(__file__).resolve().parent
PN='OFFCENTER_CUBIC_ALL_MULTIPLES_OF_THREE_IMAGE_FAMILY.md'
CN='verify_offcenter_cubic_all_order.py';IN='offcenter_cubic_all_order_inputs.json'
LN='offcenter_cubic_all_order_source_ledger.json';RN='offcenter_cubic_all_order_verification.json'
PROOF_SHA='82c28b517070c5f2d65b98df2f9cc925cf1c5a64bd5b030a2fe31cadc297d39d'
ASSURANCE=dict(external_review=False,formal_verification=False,independent_replay=False,worldwide_priority=False)
PREFIX='research/foundations/extended-geometry/'
SPECS={
 'written_proof':('local',PN,'complete_direct_endpoint_original_lift_and_every_order_fixed_Z_proof'),
 'exact_checker':('local',CN,'bounded_literal_original_and_universal_identity_support_only'),
 'typed_inputs':('local',IN,'exact_j_m_degree_marker_and_antecedent_scope'),
 'balanced_normalization':('repo',PREFIX+'TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md','ONLY_equations3_4_arbitrary_balanced_contraction_to_sphere_NOT_complete_cubic_safety'),
 'image_kernel':('repo',PREFIX+'TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md','ONLY_section2_equations2_1_2_2_direct_charzero_Image_kernel_NOT_old_witness_existence')}
SCOPE=dict(schema=1,reviewed_proof_sha256=PROOF_SHA,parameter_j_min=1,
 characteristic_zero=True,pure_orders='EVERY positive integer m',fixed_marker='Z=xi1*z2',
 fixed_marker_orders='EVERY positive integer m',even_orders_only=False,
 balanced_degree='3*j+3',ordinary_total_degree='6*j+6',
 integer_original_lift=True,radial_point=[1,0,0,1],radial_value=-1,
 original_curve_value='16*tau^3-16*tau',original_curve_derivative_at_zero=-16,
 scalar_nonpower_over_every_extension=True,
 all_degree5_or7_resolved=False,all_balanced_degrees_claimed=False,
 arbitrary_mixed_classification=False,full_conditional_theorem_input_needed=False,
 all_Aut_or_424_needed=False,FPYZ_needed=False,old_existence_family_needed=False,
 all_order_conclusion_from_finite_prefix=False,
 controls=dict(chart_j_max=6,chart_m_max=6,original_cases=[[1,3],[2,3],[3,3]],endpoint_m_max=12,kappa_m_max=40,degree_j_max=32),
 assurance_flags=ASSURANCE)
def need(v,msg):
 if not v:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
 if type(a)is not type(b):return False
 if isinstance(b,dict):return set(a)==set(b) and all(typed(a[k],v) for k,v in b.items())
 if isinstance(b,list):return len(a)==len(b) and all(typed(x,y) for x,y in zip(a,b))
 return a==b
def inputs(v):need(typed(v,SCOPE),'Exact typed every-order family scope')
REVIEWED_RAW={
 'balanced_normalization':'3e3afcf39fbbf3ec1641843cc81b550248d15db531020f4635a2ea0f5759e825',
 'image_kernel':'5dba378cfcde66b1ca7f862533c1ae3277cc2acafd9157c80bec2dd0063384a2'}
def sources(v,repo):
 need(isinstance(v,dict) and set(v)=={'schema','sources','assurance_flags'},'Exact source fields')
 need(type(v['schema'])is int and v['schema']==1 and typed(v['assurance_flags'],ASSURANCE),'Typed schema/assurance')
 rows=v['sources'];need(isinstance(rows,list) and len(rows)==len(SPECS),'Every exact mandatory role')
 out={}
 for row in rows:
  need(isinstance(row,dict) and set(row)=={'role','kind','path','relation','sha256'},'Exact source row')
  role=row['role'];need(isinstance(role,str) and role in SPECS and role not in out,'Unique known role')
  need((row['kind'],row['path'],row['relation'])==SPECS[role],'Exact role/relation/path')
  root=HERE if row['kind']=='local' else repo;p=(root/row['path']).resolve()
  need(p.is_relative_to(root) and p.is_file() and sha(p)==row['sha256'],'Missing/stale source:'+role)
  out[role]=row
 need(out['written_proof']['sha256']==PROOF_SHA,'Exact reviewed proof')
 for role,h in REVIEWED_RAW.items():need(out[role]['sha256']==h,'Exact reviewed raw mathematical body:'+role)
 return out
# Sparse polynomial arithmetic; tuples also permit Laurent x exponents.
def clean(p):return {e:c for e,c in p.items() if c}
def add(*ps):
 r={}
 for p in ps:
  for e,a in p.items():r[e]=r.get(e,0)+a
 return clean(r)
def scale(p,c):return clean({e:c*a for e,a in p.items()})
def mul(p,q):
 r={}
 for e,a in p.items():
  for f,b in q.items():
   k=tuple(x+y for x,y in zip(e,f));r[k]=r.get(k,0)+a*b
 return clean(r)
def power(p,k):
 n=len(next(iter(p)));r={(0,)*n:1}
 for _ in range(k):r=mul(r,p)
 return r
def variable(i,n):
 e=[0]*n;e[i]=1;return {tuple(e):1}
def constant(a,n):return {(0,)*n:a} if a else {}
def mon(i,n,c=1):return {(i,)+(0,)*(n-1):c}
def kappa(m):return sum((Q((-1)**r*comb(m,r),m+2*r+1) for r in range(m+1)),Q(0))
def kappa_product(m):
 den=1
 for r in range(m+1):den*=m+2*r+1
 return Q(2**m*factorial(m),den)
def angular(j):
 need(type(j)is int and j>=1,'j>=1')
 n=3;x,y,t=(variable(i,n) for i in range(n));one=constant(1,n);q=add(one,x);bt=add(t,one)
 d4={(i-1,0,0):comb(4,i) for i in range(1,5)}
 d2={(i-1,0,0):comb(2,i) for i in range(1,3)}
 return add(mul(d4,power(bt,3)),scale(mul(d2,bt),-4),
  scale(mul(power(q,4),mul(mon(j-1,n),power(bt,2))),-3),
  scale(mul(power(q,4),mul(mon(2*j-1,n),bt)),3),
  scale(mul(power(q,4),mon(3*j-1,n)),-1),
  scale(mul(power(q,2),mon(j-1,n)),4),
  scale(mul(y,add(t,constant(3,n))),-2))
def chart(j):
 x,t=variable(0,2),variable(1,2);one=constant(1,2);q=add(one,x)
 z=add(t,one,scale(mon(j,2),-1))
 p=add(mul(power(q,4),power(z,3)),scale(mul(power(q,2),z),-4))
 return {(i-1,k):c for(i,k),c in p.items()}
def sphere(p,marker=0):
 return sum((Q(c,k+1) for(i,k),c in p.items() if i+marker==0 and k%2==0),Q(0))
def original(j):
 G=angular(j);D=3*j+3;out={}
 for(i,l,k),c in G.items():
  r=D-i-l-k;need(r>=0,'Literal degree bound')
  for u in range(r+1):
   for v in range(k+1):
    e=(i+u+v,l+u+v);out[e]=out.get(e,0)+c*2**l*comb(r,u)*comb(k,v)*(-1)**(k-v)
 out=clean(out);need(out.get((D,0))==-1,'ORIGINAL radial point -1')
 return out,D
def curve(F,D):
 # Literal original xi1=z2=1, xi2=(1-tau)/2,z1=(1+tau)/2.
 out={}
 for(i,k),c in F.items():
  for a in range(D-i+1):
   for b in range(k+1):
    e=a+b
    out[e]=out.get(e,Q(0))+Q(c*comb(D-i,a)*(-1)**a*comb(k,b),2**(D-i+k))
 return clean(out)
def universal_controls():
 n=4;x,y,t,a=(variable(i,n) for i in range(n));one=constant(1,n);q=add(one,x);bt=add(t,one)
 p=add(mul(power(q,4),power(add(bt,scale(a,-1)),3)),scale(mul(power(q,2),add(bt,scale(a,-1))),-4))
 xG=add(scale(mul(mul(x,y),add(t,constant(3,n))),-2),
  mul(add(power(q,4),scale(one,-1)),power(bt,3)),
  scale(mul(add(power(q,2),scale(one,-1)),bt),-4),
  scale(mul(power(q,4),mul(a,power(bt,2))),-3),
  scale(mul(power(q,4),mul(power(a,2),bt)),3),
  scale(mul(power(q,4),power(a,3)),-1),scale(mul(power(q,2),a),4))
 defect=scale(mul(add(power(t,2),scale(mul(x,y),2),scale(one,-1)),add(t,constant(3,n))),-1)
 need(add(xG,scale(p,-1))==defect,'Universal regularity/divisibility modulo quadric')
 # Clearing the Q denominator in the exact v-profile transformation.
 n=2;U,v=variable(0,n),variable(1,n);one=constant(1,n)
 left=add(mul(power(U,4),power(scale(v,2),3)),scale(mul(power(U,4),scale(v,2)),-4))
 right=scale(mul(power(U,4),add(power(v,3),scale(v,-1))),8)
 need(left==right,'Universal cleared canonical cubic profile')
 # Universal coefficient-polynomial identity for the m->m+2 kappa ratio.
 m=variable(0,1);one=constant(1,1)
 A=scale(mul(mul(power(add(m,one),2),add(m,constant(2,1))),mul(add(scale(m,3),constant(5,1)),add(scale(m,3),constant(7,1)))),12)
 B=scale(mul(mul(add(m,one),add(m,constant(2,1))),mul(add(scale(m,3),constant(3,1)),mul(add(scale(m,3),constant(5,1)),add(scale(m,3),constant(7,1))))),4)
 need(A==B,'Universal kappa product cancellation')
 R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}
 need(add(power(R,2),scale(power(T,2),-1),scale(mul(Z,W),-2))=={},'Original rank-one quadric')
 return ['regularity_mod_quadric','cleared_endpoint_profile','kappa_product_ratio','original_rank_one_quadric']
def endpoint_kappa_controls():
 endpoints=[]
 for m in range(1,SCOPE['controls']['endpoint_m_max']+1):
  primitive={m+2*r+1:Q((-1)**r*comb(m,r),m+2*r+1) for r in range(m+1)}
  shift={}
  for k,c in primitive.items():
   for b in range(k+1):shift[b]=shift.get(b,Q(0))+c*comb(k,b)
  shift[0]-=kappa(m);shift=clean(shift)
  need(min(primitive)>=m+1 and min(shift)>=m+1,'Both exact endpoint jets vanish through m')
  endpoints.append(m)
 kappas=[]
 for m in range(SCOPE['controls']['kappa_m_max']+1):
  k=kappa(m);need(k==kappa_product(m) and k>0,'Raw interval kappa equals positive product')
  need(3*(3*m+5)*(3*m+7)*kappa(m+2)==4*(m+1)*(m+2)*k,'Exact m+2 recurrence')
  kappas.append(str(k))
 return endpoints,kappas
def chart_degree_controls():
 rows=[];degrees=[]
 for j in range(1,SCOPE['controls']['degree_j_max']+1):
  p=angular(j);D=3*j+3
  need(max(map(sum,p))==D and p.get((D,0,0))==-1,'Exact least-degree top normal form')
  degrees.append(dict(j=j,D=D,pure_x_top=-1))
 for j in range(1,SCOPE['controls']['chart_j_max']+1):
  p=chart(j);pp={(0,0):1}
  for m in range(1,SCOPE['controls']['chart_m_max']+1):
   pp=mul(pp,p);pure=sphere(pp);marked=sphere(pp,1)
   need(pure==0 and marked==(-8)**m*kappa(m),'Independent raw Laurent interval contraction')
   rows.append(dict(j=j,m=m,pure='0',marked=str(marked)))
 return rows,degrees
def original_controls():
 out=[]
 for j,maximum in SCOPE['controls']['original_cases']:
  F,D=original(j);need(curve(F,D)=={1:-16,3:16},'FULL ORIGINAL simple curve 16tau^3-16tau')
  pp={(0,0):1}
  for m in range(1,maximum+1):
   pp=mul(pp,F);n=m*D
   pure=sum(c*factorial(i)*factorial(n-i) for(i,k),c in pp.items() if i==k)
   marked=sum(c*factorial(i+1)*factorial(n-i) for(i,k),c in pp.items() if k==i+1)
   expected=(-8)**m*kappa(m)
   need(pure==0 and Q(marked,factorial(n+2))==expected,'DIRECT ORIGINAL two-pair factorial contraction')
   need(marked!=0,'EVERY order including odd')
   need(Q(marked,factorial((D-1)*m+2))!=expected,'Wrong radial factorial is rejected')
   out.append(dict(j=j,m=m,D=D,original_terms=len(F),power_terms=len(pp),pure='0',
    marked=str(marked),expected_normalized=str(expected),radial_point=-1,curve=[0,-16,0,16]))
 return out
def negatives(ledger,repo):
 out=[]
 def reject(label,fn):
  try:fn()
  except(ValueError,KeyError,TypeError,OSError):out.append(label);return
  raise ValueError('Negative accepted:'+label)
 for role in SPECS:
  bad=copy.deepcopy(ledger);bad['sources']=[r for r in bad['sources'] if r['role']!=role]
  reject('omit:'+role,lambda bad=bad:sources(bad,repo))
  bad=copy.deepcopy(ledger);next(r for r in bad['sources'] if r['role']==role)['sha256']='0'*64
  reject('badpin:'+role,lambda bad=bad:sources(bad,repo))
 for key,value in [('parameter_j_min',0),('fixed_marker_orders','EVERY positive EVEN integer m'),
  ('even_orders_only',True),('fixed_marker','Z^2'),('balanced_degree','3*j'),('radial_value',0),
  ('original_curve_derivative_at_zero',0),('all_Aut_or_424_needed',True),('FPYZ_needed',True),
  ('old_existence_family_needed',True),('all_degree5_or7_resolved',True),('all_balanced_degrees_claimed',True),
  ('all_order_conclusion_from_finite_prefix',True),('schema',True)]:
  bad=copy.deepcopy(SCOPE);bad[key]=value;reject('false_scope:'+key,lambda bad=bad:inputs(bad))
 bad=copy.deepcopy(ledger);bad['schema']=True;reject('Boolean_ledger_schema',lambda:sources(bad,repo))
 reject('j0_wrong_parameter',lambda:angular(0))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
 ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args()
 repo=ns.repo_root.resolve();v=load(HERE/IN);inputs(v);L=load(HERE/LN);pins=sources(L,repo)
 universal=universal_controls();endpoints,ks=endpoint_kappa_controls()
 chartrows,degrees=chart_degree_controls();originals=original_controls();rejects=negatives(L,repo)
 report=dict(schema=1,status='PASS_OFFCENTER_CUBIC_EVERY_ORDER_IMAGE_FAMILY',scope=v,
 proof_sha256=sha(HERE/PN),checker_sha256=sha(HERE/CN),source_ledger_sha256=sha(HERE/LN),
 sources=[pins[k] for k in sorted(pins)],mandatory_source_roles=len(pins),
 universal_identity_families=universal,endpoint_orders=endpoints,kappa_controls=ks,
 raw_chart_controls=chartrows,degree_controls=degrees,original_four_variable_controls=originals,
 source_scope_controls=rejects,assurance_flags=ASSURANCE,
 meaning='The written endpoint congruence proves every j,m. Finite controls reconstruct ORIGINAL lifts/contractions and authenticate exact narrow sources; no old witness or broader classification is imported.')
 ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',roles=len(pins),universal=len(universal),chart=len(chartrows),
 original=len(originals),endpoints=len(endpoints),kappas=len(ks),source_controls=len(rejects))))
if __name__=='__main__':main()

