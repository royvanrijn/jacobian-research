#!/usr/bin/env python3
"""Bounded exact identities and source authentication for mixed-cubic reductions."""
from fractions import Fraction as F
from pathlib import Path
import copy
HERE=Path(__file__).parent
import hashlib, json, argparse
def clean(p): return {e:F(a) for e,a in p.items() if a}
def zero(n): return {}
def con(a,n): return {(0,)*n:F(a)} if a else {}
def var(i,n):
 e=[0]*n;e[i]=1;return {tuple(e):F(1)}
def add(*ps):
 r={}
 for p in ps:
  for e,a in p.items():r[e]=r.get(e,F(0))+a
 return clean(r)
def scale(p,a):return clean({e:a*v for e,v in p.items()})
def mul(p,q):
 r={}
 for e,a in p.items():
  for f,b in q.items():
   k=tuple(i+j for i,j in zip(e,f));r[k]=r.get(k,F(0))+a*b
 return clean(r)
def pw(p,k,n):
 r=con(1,n)
 for _ in range(k):r=mul(r,p)
 return r
def same(p,q):return clean(p)==clean(q)
def avg_t(p):
 r={}
 for e,a in p.items():
  if not e[2]%2:
   k=e[:2]+(0,);r[k]=r.get(k,F(0))+a/F(e[2]+1)
 return clean(r)
def reduce_c2(p,q):
 r={}
 for e,a in p.items():
  term=mul(con(a,3),pw(scale(q,-F(1,28)),e[1]//2,3))
  mon={(e[0],e[1]%2,e[2]):F(1)}
  r=add(r,mul(term,mon))
 return r
def rem_b(p,q):
 p=clean(p);q=clean(q);d=max(e[0] for e in q)
 assert all(e[1:]==(0,0) for e in p)
 assert all(e[1:]==(0,0) for e in q)
 lc=q[(d,0,0)]
 while p and max(e[0] for e in p)>=d:
  j=max(e[0] for e in p);a=p[(j,0,0)]/lc
  p=add(p,scale(mul({(j-d,0,0):F(1)},q),-a))
 return p
def eval_b(p,a):return sum(v*a**e[0] for e,v in p.items())
def need(ok,label):
 if not ok:raise AssertionError(label)
 return label

INPUT_EXPECTED={'schema': 1, 'scope': {'characteristic_zero': True, 'nonnegative_angular_all_heights': True, 'zero_center_cubic': True, 'cubic_prefix': 4, 'prefix_three_suffices': False, 'conditional_invariant_zero': True, 'conditional_E0_nonzero': True, 'conditional_offcenter_C0_over_B0': [-3, 3], 'conditional_degree_ceiling': 7, 'conditional_safe_degrees': [5, 7], 'purity_implies_invariant_zero': False, 'unrestricted_mixed_classification': False, 'global_degree5_7_resolution': False}, 'analytic_boundary': {'interval_rigidity': 'FPYZ2011 Corollary4.1, printed27; complex coefficients and distinct interval endpoints', 'field_transfer': 'Only the finitely generated coefficient field embeds in C', 'cubic_prefix': 'Direct rational identities; no external finite-prefix theorem imported', 'conditional_branch': 'Written formal endpoint recurrence and rational-function valuations; finite controls do not prove it'}, 'assurance': {'external_review': False, 'formal_verification': False, 'independent_replay': False, 'worldwide_priority': False}}
PROOF_HASH='c4fee78a3b99b286e22818f63c90f8da450293bde974b354fb22bd2866bf1dc6'
ROLES={
 'proof':('MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md',PROOF_HASH,'mathematical_input'),
 'checker':('verify_mixed_cubic_hopf.py',None,'checker'),
 'inputs':('mixed_cubic_hopf_inputs.json',None,'typed_scope'),
 'classical_interval_rigidity':('mixed_cubic_sources/francoise_pakovich_yomdin_zhao_2011.pdf','6dcafd3fa8920502e6cb3e45e691b0e328bc364deccdb85aa81b2377a1d4d67f','mathematical_input'),
 'classical_stronger_context':('mixed_cubic_sources/muger_tuset_2019.pdf','39292887e1ba9d6a44a0909cb64136e3b8dab21c45457eb08b7b55b9de9f1b56','context')}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strict(a,b):
 if type(a) is not type(b):return False
 if isinstance(a,dict):return set(a)==set(b) and all(strict(a[k],b[k]) for k in a)
 if isinstance(a,list):return len(a)==len(b) and all(strict(x,y) for x,y in zip(a,b))
 return a==b
def authenticate(inp,L):
 if not strict(inp,INPUT_EXPECTED):raise ValueError('typed mathematical scope changed')
 if not isinstance(L,dict) or set(L)!={'schema','files','assurance'}:raise ValueError('ledger structure')
 if type(L['schema']) is not int or L['schema']!=1:raise ValueError('schema')
 if not strict(L['assurance'],INPUT_EXPECTED['assurance']):raise ValueError('assurance')
 if not isinstance(L['files'],list) or len(L['files'])!=len(ROLES):raise ValueError('mandatory records')
 got={}
 for row in L['files']:
  if set(row)!={'role','path','sha256','relation','locator'}:raise ValueError('source row')
  role=row['role']
  if role in got or role not in ROLES:raise ValueError('role duplicate/unknown')
  got[role]=row
 for role,(name,pinned,relation) in ROLES.items():
  row=got[role]
  if row['path']!=name or row['relation']!=relation:raise ValueError('path/relation')
  p=HERE/name
  if row['sha256']!=digest(p):raise ValueError('source mismatch')
  if pinned is not None and row['sha256']!=pinned:raise ValueError('reviewed source mismatch')
  if not isinstance(row['locator'],str) or not row['locator']:raise ValueError('locator')
def source_negatives(inp,L):
 controls=[]
 def reject(name,a,b):
  try:authenticate(a,b)
  except (ValueError,KeyError,TypeError,FileNotFoundError):controls.append(name);return
  raise AssertionError('accepted corrupted source/scope: '+name)
 for j,row in enumerate(L['files']):
  bad=copy.deepcopy(L);bad['files'].pop(j);reject('omit:'+row['role'],inp,bad)
  bad=copy.deepcopy(L);bad['files'][j]['sha256']='0'*64;reject('bad_hash:'+row['role'],inp,bad)
 for key,value in inp['scope'].items():
  bad=copy.deepcopy(inp)
  if type(value) is bool:bad['scope'][key]=not value
  elif type(value) is int:bad['scope'][key]=value+1
  else:bad['scope'][key]=[]
  reject('false_scope:'+key,bad,L)
 bad=copy.deepcopy(L);bad['schema']=True;reject('boolean_schema',inp,bad)
 bad=copy.deepcopy(inp);bad['scope']['cubic_prefix']=True;reject('boolean_integer_domain',bad,L)
 bad=copy.deepcopy(L);bad['files'][0]['relation']='context';reject('logical_input_downgraded',inp,bad)
 return controls

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--output',type=Path,default=HERE/'mixed_cubic_hopf_verification.json')
 ap.add_argument('--repo-root',type=Path)
 args=ap.parse_args()
 if not __debug__:raise RuntimeError('optimized execution is not accepted')
 inp=json.loads((HERE/'mixed_cubic_hopf_inputs.json').read_text())
 ledger=json.loads((HERE/'mixed_cubic_hopf_source_ledger.json').read_text())
 authenticate(inp,ledger)
 source_controls=source_negatives(inp,ledger)
 b,c,t=(var(i,3) for i in range(3));one=con(1,3)
 b2=pw(b,2,3);b3=pw(b,3,3);b4=pw(b,4,3);c2=pw(c,2,3);c4=pw(c,4,3)
 g=add(pw(t,3,3),mul(b,t),mul(c,add(pw(t,2,3),con(-F(1,3),3))))
 mu=[avg_t(pw(g,m,3)) for m in range(1,5)]
 q=add(scale(b2,105),scale(b,126),con(45,3))
 P=add(scale(b2,56),scale(b,105),con(45,3))
 R=add(scale(b4,28028),scale(b3,84084),scale(b2,93275),scale(b,46410),con(8955,3))
 N4=add(scale(b4,27027),scale(b3,77220),scale(mul(b2,c2),37752),scale(b2,90090),
         scale(mul(b,c2),61776),scale(b,49140),scale(c4,2288),scale(c2,26520),con(10395,3))
 checks=[]
 checks.append(need(not mu[0],'normalized_first_moment'))
 checks.append(need(same(mu[1],scale(add(q,scale(c2,28)),F(1,315))),'second_moment'))
 checks.append(need(same(mu[2],scale(mul(c,add(scale(b2,63),scale(b,108),scale(c2,4),con(45,3))),F(4,945))),'third_moment'))
 checks.append(need(same(mu[3],scale(N4,F(1,135135))),'fourth_moment'))
 checks.append(need(same(reduce_c2(mu[2],q),scale(mul(c,P),F(8,2205))),'third_reduction'))
 checks.append(need(same(reduce_c2(mu[3],q),scale(R,-F(16,735735))),'fourth_reduction'))
 rq=scale(add(scale(b,117),con(-110,3)),-F(84,125))
 rp=scale(add(scale(b,299),con(159,3)),-F(315,128))
 checks.append(need(same(rem_b(R,q),rq),'q_remainder'))
 checks.append(need(same(rem_b(R,P),rp),'P_remainder'))
 checks.append(need(eval_b(q,F(110,117))==F(1169375,4563),'q_contradiction'))
 checks.append(need(eval_b(P,-F(159,299))==F(446976,89401),'P_contradiction'))
 g2=add(pw(t,2,3),mul(b,t),con(-F(1,3),3))
 m2=avg_t(pw(g2,2,3));m3=avg_t(pw(g2,3,3))
 checks.append(need(same(add(m3,scale(m2,-F(4,5))),con(-F(256,4725),3)),'quadratic_contradiction'))
 # Independent polynomial identities for the depressed cubic.
 n=5;e,s,h,T,v=(var(i,n) for i in range(n))
 C=mul(e,add(scale(pw(s,3,n),-1),mul(pw(h,2,n),s)))
 B=mul(e,add(scale(pw(s,2,n),3),scale(pw(h,2,n),-1)))
 D=scale(mul(e,s),-3)
 invariant=add(scale(mul(C,pw(e,2,n)),27),scale(mul(mul(B,D),e),-9),scale(pw(D,3,n),2))
 checks.append(need(not invariant,'depressed_cubic_invariant'))
 poly=add(C,mul(B,T),mul(D,pw(T,2,n)),mul(e,pw(T,3,n)))
 dep=mul(e,add(pw(add(T,scale(s,-1)),3,n),scale(mul(pw(h,2,n),add(T,scale(s,-1))),-1)))
 checks.append(need(same(poly,dep),'depressed_cubic_expansion'))
 Tv=add(s,mul(h,v))
 sub=add(C,mul(B,Tv),mul(D,pw(Tv,2,n)),mul(e,pw(Tv,3,n)))
 profile=scale(mul(mul(e,pw(h,3,n)),add(v,scale(pw(v,3,n),-1))),-1)
 checks.append(need(same(sub,profile),'endpoint_profile'))
 # Nonzero cubic with first THREE moments zero over Q[b]/(q):
 odd=add(pw(t,3,3),mul(b,t))
 oddmom=[avg_t(pw(odd,m,3)) for m in range(1,5)]
 need(all(not rem_b(p,q) for p in oddmom[:3]),'three_prefix_control')
 need(bool(rem_b(oddmom[3],q)),'fourth_does_not_vanish_control')
 negatives=['first_three_moments_do_not_force_zero']
 need(not same(rem_b(R,q),add(rq,con(1,3))),'corrupt_remainder_control')
 negatives.append('corrupt_remainder')
 corrupt=add(scale(mul(C,pw(e,2,n)),27),scale(mul(mul(B,D),e),-8),scale(pw(D,3,n),2))
 need(bool(corrupt),'corrupt_invariant_control')
 negatives.append('corrupt_depressed_invariant')
 # Actual safe quintic: arbitrary-alpha shear of y(t+alpha).
 x,alpha,T=b,c,t
 z=add(T,pw(x,2,3))
 raw_shift=scale(mul(add(con(1,3),scale(pw(z,2,3),-1)),add(z,alpha)),F(1,2))
 C5=scale(add(alpha,pw(x,2,3),scale(mul(alpha,pw(x,4,3)),-1),scale(pw(x,6,3),-1)),F(1,2))
 B5=scale(add(con(1,3),scale(mul(alpha,pw(x,2,3)),-2),scale(pw(x,4,3),-3)),F(1,2))
 D5=scale(add(alpha,scale(pw(x,2,3),3)),-F(1,2));E5=con(-F(1,2),3)
 P5=add(C5,mul(B5,T),mul(D5,pw(T,2,3)),mul(E5,pw(T,3,3)))
 checks.append(need(same(P5,raw_shift),'safe_quintic_literal_numerator'))
 I5=add(scale(mul(C5,pw(E5,2,3)),27),scale(mul(mul(B5,D5),E5),-9),scale(pw(D5,3,3),2))
 expectedI=scale(mul(alpha,add(con(9,3),scale(pw(alpha,2,3),-1))),F(1,4))
 checks.append(need(same(I5,expectedI),'safe_quintic_invariant'))
 # Set alpha=1 and divide by x; negative angular exponents are exact.
 f5={}
 for ex,a0 in P5.items():
  key=(ex[0]-1,0,ex[2]);f5[key]=f5.get(key,F(0))+a0
 f5=clean(f5)
 top=max(sum(ex) for ex in f5)
 checks.append(need(top==5 and {ex:a0 for ex,a0 in f5.items() if sum(ex)==top}=={(5,0,0):-F(1,2)},'safe_quintic_exact_top_degree'))
 J3={}
 for ex,a0 in f5.items():
  if ex[2]>=3:
   key=(ex[0]+3,0,ex[2]-3)
   J3[key]=J3.get(key,F(0))+a0*ex[2]*(ex[2]-1)*(ex[2]-2)
 checks.append(need(same(J3,{(2,0,0):F(-3)}),'safe_quintic_exact_height'))
 rawmom=[avg_t(pw(f5,m,3)) for m in range(1,5)]
 need(all(not any(ex[0]==0 for ex in value) for value in rawmom),'safe_quintic_raw_pure_controls')
 need(sum(I5.values())==2 and all(ex[0]==ex[2]==0 for ex in I5),'purity_does_not_imply_invariant_zero')
 negatives.append('purity_does_not_imply_invariant_zero_even_quintic')

 proof=HERE/'MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md'
 report={'status':'PASS_SOURCE_GUARDED_MIXED_CUBIC_REDUCTIONS','proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest(),
         'source_ledger_sha256':digest(HERE/'mixed_cubic_hopf_source_ledger.json'),
         'scope':inp['scope'],'analytic_boundary':inp['analytic_boundary'],
         'mandatory_source_roles':len(ROLES),'source_scope_controls':source_controls,
         'universal_identities':checks,'negative_controls':negatives,
         'assurance':{'external_review':False,'formal_verification':False,'independent_replay':False}}
 args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':'PASS','identities':len(checks),'negatives':len(negatives),'output':str(args.output)}))
if __name__=='__main__':main()

