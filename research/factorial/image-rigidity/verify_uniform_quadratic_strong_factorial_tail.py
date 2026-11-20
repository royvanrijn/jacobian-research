#!/usr/bin/env python3
"""Exact Bell/Laurent and boundary support for the uniform quadratic tail.

The written analytic proof supplies the all-order contour argument. Literal
coefficient controls are not used as a substitute for that proof.
"""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: exact checks must execute.')
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse,json,math
HERE=Path(__file__).resolve().parent
ROWS=[]
def add(a,b):
 out=dict(a)
 for k,v in b.items():out[k]=out.get(k,Q(0))+v
 return{k:v for k,v in out.items()if v}
def scale(a,c):return{k:v*c for k,v in a.items()if v*c}
def mul(a,b):
 out={}
 for i,u in a.items():
  for j,v in b.items():
   k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,Q(0))+u*v
 return{k:v for k,v in out.items()if v}
def power(a,n,dim):
 result={(0,)*dim:Q(1)}
 for _ in range(n):result=mul(result,a)
 return result
def ck(name,a,b,op='<'):
 a,b=Q(a),Q(b);assert(a<b if op=='<'else a==b),(name,a,op,b)
 ROWS.append(dict(name=name,left=str(a),relation=op,right=str(b)))
def eu(x):
 x=Q(x);assert 0<=x<65
 term=total=Q(1)
 for j in range(1,65):term*=x/j;total+=term
 return total+term*x/65/(1-x/66)
def bell_Q(d,m,omit_secondary=False):
 maxj=d*m;kap=[{} for _ in range(maxj+1)]
 kap[1]={(1,):Q(1,d*m)}
 for k in range(2,maxj+1):
  fac=math.factorial(k-1)
  kap[k]={(k-2,):Q(fac,d**(k-2)*m**(k-1))}
  if not omit_secondary:kap[k][(k,)]=Q(fac,d**k*m**k)
 moments=[{(0,):Q(1)}]
 for n in range(1,maxj+1):
  value={}
  for k in range(1,n+1):value=add(value,scale(mul(kap[k],moments[n-k]),math.comb(n-1,k-1)))
  moments.append(value)
 h={(j-1,j):Q(math.comb(d,j),d**j)for j in range(1,d+1)}
 hp={(0,0):Q(1)};out={}
 for r in range(m+1):
  for (e,j),c in hp.items():
   for (f,),u in moments[j].items():
    key=(e+f,m-r);out[key]=out.get(key,Q(0))+math.comb(m,r)*c*u
  hp=mul(hp,h)
 return{k:v for k,v in out.items()if v}
def direct_Q(d,m,start=1):
 # Independent rising-factorial construction is Laurent before cancellation.
 a={(-2,0):Q(d*d)};am=scale(a,m);C=(m*d*d)**d
 shift={(-2*d,0):Q(-C),(-2*d+1,1):Q(C)}
 out={}
 for k in range(m+1):
  rise={(0,0):Q(1)}
  for j in range(start,start+d*k):rise=mul(rise,add(am,{(0,0):Q(j)}))
  term=scale(mul(rise,power(shift,m-k,2)),Q(math.comb(m,k),C**m))
  term={(e+(2*d-1)*m,z):v for(e,z),v in term.items()}
  out=add(out,term)
 return out
def gaussian_Q(m):
 out={}
 for j in range(m//2+1):
  out[(0,m-2*j)]=Q(math.comb(m,2*j)*math.factorial(2*j),2**j*math.factorial(j)*m**j)
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'uniform_quadratic_strong_factorial_tail_verification.json');args=ap.parse_args()
 delta=Q(1,10**12);eps=Q(1,10**7);m0=10**16;R=1+Q(1,10**6)
 crit=R**3/(6*(1-delta*R/2)**2)+(R**3/2)*eu(delta*R)/(1-delta*R)
 ck('real-a finite contour right-half-plane',0,1-Q(16,5)*delta/2)
 ck('uniform bulk stationary delta derivative',crit,1)
 ck('vertical bulk action barrier',2,m0*(Q(59,100)*eps-2*delta))
 ck('log6 below2',6,sum(Q(2**j,math.factorial(j))for j in range(6)))
 transfer=Q(49,100)*delta/(Q(1,10**5)-delta)
 ck('homotopy endpoint zeta transfer',transfer,Q(5,10**8))
 ck('Hermite endpoint imaginary barrier',Q(97,100)*eps-Q(5,10**8),Q(47,10**9),'=')
 ck('actual zero band below boundary barrier',Q(5200,m0),Q(47,10**9))
 ck('homotopy delta inside uniform phase',delta,Q(1,10**10))
 ck('boundary orders exceed outside threshold',10**12,m0)
 ck('quadratic actual delta relation',10**24*delta**2,1,'=')
 ck('minimum hard lower order at d2',m0,10**24)
 ck('coupled radial exponential factor',2/(1-2*delta**2),3)
 ck('coupled radial bulk real containment',eps+3*delta+9*delta**2,2*eps)
 ck('coupled radial endpoint real containment',eps+3*delta+9*delta**2,Q(1,10**5))
 ck('coupled radial bulk imaginary containment',Q(19,10)*(1+3*delta**2),Q(191,100))
 ck('coupled radial endpoint imaginary containment',Q(21,10)*(1+3*delta**2),Q(21001,10000))
 ck('quadratic threshold endpoint separation',15600,10**12)
 ck('actual bulk packet ratio below1 via exp(x)>=1+x',40,1+Q(2,5)*2*10**12)
 # Universal cumulant identity in affine exponent vectors. Each coordinate
 # pair gives coefficient of formal k and its constant part for delta,d,m.
 first={(1,-2,-1,2,-1,1):Q(1)}
 second={(1,0,-1,0,-1,0):Q(1)}
 multiply_by={(-1,0,1,0,1,0):Q(1)}
 assert mul(add(first,second),multiply_by)=={(0,-2,0,2,0,1):Q(1),(0,0,0,0,0,0):Q(1)}
 # kappa1 cancellation, after multiplication by m*sqrt(a).
 am={(1,1):Q(1)};one={(0,0):Q(1)}
 assert add(add(am,one),scale(am,-1))==one
 # Universal coefficient recurrence of monic He_n:
 # (n+1-2j)+2j=n+1 after common factorial denominator.
 n={(1,0):Q(1)};j={(0,1):Q(1)}
 assert add(add(add(n,one),scale(j,-2)),scale(j,2))==add(n,one)
 cases=[];stream=[]
 shapes=[(d,m)for d in range(2,7)for m in range(1,6)]+[(8,3),(9,2)]
 for d,m in shapes:
  b=bell_Q(d,m);q=direct_Q(d,m)
  assert b==q and all(e>=0 for e,z in q)
  assert{(e,z):v for(e,z),v in q.items()if z==m}=={(0,m):Q(1)}
  atzero={(e,z):v for(e,z),v in q.items()if e==0};assert atzero==gaussian_Q(m)
  row=dict(d=d,m=m,terms=len(q),largest_delta_degree=max(e for e,z in q),coefficients=[dict(delta=e,z=z,value=str(v))for(e,z),v in sorted(q.items())]);stream.append(row)
  cases.append({k:v for k,v in row.items()if k!='coefficients'})
 # Four controls reject correct-looking but false alternatives.
 assert direct_Q(2,2,start=0)!=bell_Q(2,2) # a*m vs a*m+1 shape.
 assert bell_Q(2,2,omit_secondary=True)!=direct_Q(2,2) # Lost +1 in Gamma shape.
 assert gaussian_Q(2)=={(0,2):Q(1),(0,0):Q(1,2)}
 assert gaussian_Q(2)!={(0,2):Q(1),(0,0):Q(1)} # Wrong unit variance.
 assert Q(49,100)*Q(1,10**8)/(Q(1,10**5)-Q(1,10**8))>Q(97,100)*eps
 proof=HERE/'UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md'
 result=dict(schema=1,status='PASS_EXACT_BELL_AND_BOUNDARY_SUPPORT',scope='Exact support for the written monic Gamma/Bell homotopy, real-a boundary proof and uniform quadratic tail. The analytic family/contour/argument-principle proof is written, not inferred from these finite controls.',proof_sha256=sha256(proof.read_bytes()).hexdigest(),checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),rational_guards=ROWS,universal_identity_count=3,exact_laurent_bell_cases=cases,coefficient_stream_sha256=sha256(json.dumps(stream,sort_keys=True,separators=(',',':')).encode()).hexdigest(),coefficient_stream=stream,negative_controls=4)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='PASS',rational_guards=len(ROWS),universal_identities=3,exact_laurent_bell_cases=len(cases),negative_controls=4)))
if __name__=='__main__':main()
