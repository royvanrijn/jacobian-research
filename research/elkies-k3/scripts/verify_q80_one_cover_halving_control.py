#!/usr/bin/env python3
"""Bounded exact one-cover control: finite quotients and independent characters."""
import sys,json,resource,argparse
from hashlib import sha256
from math import isqrt
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_q80_branch_trace_specialization import quotient,rank
from audit_q80_branch_trace_specialization import prime
root=Path(__file__).resolve().parents[2]; s=json.load(open(root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json')); b=json.load(open(root/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json'))['bisections'][0]
q=list(map(int,b['branch']['numerator_coefficients'])); rows=[]; records=[]
def ev(cs,t,p):
 v=0
 for c in reversed(cs):
  z=F(c);v=(v*t+z.numerator*pow(z.denominator,-1,p))%p
 return v
def pt(sec,t,p):
 return tuple(ev(sec[k]['numerator_coefficients_low_to_high'],t,p)*pow(ev(sec[k]['denominator_coefficients_low_to_high'],t,p),-1,p)%p for k in ['X','Y'])
for p in range(5,998):
 if not prime(p) or q[2]%p==0 or (q[1]**2-4*q[0]*q[2])%p==0:continue
 for t in range(p):
  if ev(q,t,p):continue
  try:
   a=ev(s['weierstrass_model']['A_coefficients_low_to_high'],t,p);c=ev(s['weierstrass_model']['B_coefficients_low_to_high'],t,p)
   pts=[pt(sec,t,p) for sec in s['sections']['records']]
  except ValueError:continue
  if (4*a**3+27*c*c)%p==0:continue
  assert all((y*y-x*x*x-a*x-c)%p==0 for x,y in pts)
  labels,dim,n,nd=quotient(a,c,p)
  rr=[sum(((labels[P]>>i)&1)<<j for j,P in enumerate(pts)) for i in range(dim)]
  cc=[]
  for e in range(p):
   if (e**3+a*e+c)%p:continue
   vv=[(x-e)%p if x!=e else (3*e*e+a)%p for x,y in pts]
   assert all(vv)
   cc.append(sum((pow(v,(p-1)//2,p)==p-1)<<j for j,v in enumerate(vv)))
  assert rank(cc)==rank(rr)==rank(cc+rr)
  records.append(dict(p=p,t=t,A=a,B=c,rows=rr,character_rows=cc,order=n,doubles=nd));rows+=rr
  if rank(rows)>=16:break
 if rank(rows)>=16:break
ker=[v for v in range(1<<17) if all((v&r).bit_count()%2==0 for r in rows)]

assert rank(rows)==16 and ker==[0,65538]
assert sum((n%2)<<i for i,n in enumerate(b['section_basis_w']))==65538
# Exact polynomial identities for the retained section; arithmetic in Q[t,w]/(w²-q).
def plus(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def scale(a,k):return [k*x for x in a]
def equal(a,b):return not any(plus(a,scale(b,-1)))
model=s['weierstrass_model'];A=list(map(F,model['A_coefficients_low_to_high']));B=list(map(F,model['B_coefficients_low_to_high']))
lift=b['lifted_section'];x,z,y,v=[list(map(F,lift[k+'_coefficients'])) for k in ['x0','x1','y0','y1']]
assert equal(plus(mul(y,y),mul(q,mul(v,v))),plus(plus(mul(mul(x,x),x),scale(mul(q,mul(x,mul(z,z))),3)),plus(mul(A,x),B)))
assert equal(scale(mul(y,v),2),plus(plus(scale(mul(mul(x,x),z),3),mul(q,mul(mul(z,z),z))),mul(A,z)))
assert any(z) # Q differs from sigma(Q); torsion-freeness is inherited.
disc=q[1]**2-4*q[0]*q[2]
assert disc!=0 and (disc<0 or isqrt(disc)**2!=disc)
# Degree-one reduction witnesses no branch-field rational 2-torsion and smoothness.
witness=next(r for r in records if r['p']==23 and r['t']==0)
assert all((e**3+witness['A']*e+witness['B'])%23 for e in range(23))
pairpath=root/'artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json'
pair=json.loads(pairpath.read_text())['pairs'][0];point=pair['v4_base_point']
u=F(point['u']);w=F(point['left_square_root'])
assert w*w==sum(F(c)*u**i for i,c in enumerate(q))
paths=[root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json',root/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json',pairpath,Path(__file__),Path(__file__).with_name('verify_q80_branch_trace_specialization.py'),Path(__file__).with_name('audit_q80_branch_trace_specialization.py')]
result=dict(schema='q80-one-cover-halving-control-v1',label=b['label'],rank=rank(rows),kernel=ker,branch_quadratic=q,records=records,rational_base_point={'t':str(u),'w':str(w)},exact_lift_identity=True,noninvariant_x=True,branch_torsion_witness=witness,inputs={str(p.relative_to(root)):sha256(p.read_bytes()).hexdigest() for p in paths},proof_boundary='Exact rank18 uses the retained Q80 saturated rank17 and unramified Kummer theorem; no second gain constructed.',prime_pool=[5,997],cpu_limit_seconds=20,memory_limit_bytes=1024**3)
parser=argparse.ArgumentParser();parser.add_argument('--record',type=Path,required=True);args=parser.parse_args();args.record.parent.mkdir(parents=True,exist_ok=True)
with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({'rank':16,'kernel':ker,'exact_lift_identity':True,'record':str(args.record)}))
