#!/usr/bin/env python3
"""Exact finite unit-ideal certificate for a specified Hopf quintic chart."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from functools import lru_cache
from math import comb
import argparse,json,hashlib
CERT=Path(__file__).resolve().with_name('hopf_dirichlet_quintic_certificate.json')
EXPECTED_SHA='fe4bc9c0a5479c33f4f96d4b5ae61e77de156d376acac32d435dc399d8f756cb'
def add(P,Q):
 z=dict(P)
 for e,c in Q.items():z[e]=z.get(e,0)+c
 return {e:c for e,c in z.items() if c}
def mul(P,Q):
 z=defaultdict(int)
 for a,c in P.items():
  for b,d in Q.items():z[tuple(x+y for x,y in zip(a,b))]+=c*d
 return {e:c for e,c in z.items() if c}
def power(P,n):
 z={(0,)*len(next(iter(P))):1}
 for _ in range(n):z=mul(z,P)
 return z
def scale(P,c):return {e:c*v for e,v in P.items() if c*v}
@lru_cache(None)
def rf(a,n):
 z=F(1)
 for i in range(n):z*=a+i
 return z
@lru_cache(None)
def mu(i,j):return rf(F(1,2),i)*rf(F(1,2),j)/rf(F(3,2),i+j)
def load(rows):return {tuple(e['powers']):int(e['coefficient']) for e in rows}
def exact_certificate(G,M,N):
 z={}
 for a,b in zip(G,M):z=add(z,mul(a,b))
 assert z=={(0,0,0):N},'nonzero integer Bezout identity'
def centered_dirichlet_moments():
 # Keys: u,v,A,B,C. D=-A/5-1/15-(B+C)/3 is already substituted.
 Q={(2,0,1,0,0):F(1),(1,1,0,0,0):F(1),(1,0,0,1,0):F(1),(0,1,0,0,1):F(1),(0,0,1,0,0):F(-1,5),(0,0,0,1,0):F(-1,3),(0,0,0,0,1):F(-1,3),(0,0,0,0,0):F(-1,15)}
 P={(0,0,0,0,0):F(1)};out={}
 for m in range(1,6):
  P=mul(P,Q);z=defaultdict(F)
  for e,c in P.items():z[e[2:]]+=c*mu(e[0],e[1])
  out[m]={e:c for e,c in z.items() if c}
 assert out[1]=={}
 return out
def angular_integral(P):
 out=defaultdict(F)
 for (x,t,a,b,c),v in P.items():
  if x==0 and t%2==0:out[(a,b,c)]+=v/F(t+1)
 return {e:c for e,c in out.items() if c}
def raw_angular_moments():
 # Keys: x,t,A,B,C. Direct native functions B0,C0 are expanded literally.
 B0={(-1,0,0,0,0):1,(-1,2,0,0,0):-1,(0,2,0,0,0):-2,(1,2,0,0,0):-1}
 C0={(-1,0,0,0,0):1,(-1,2,0,0,0):-1,(0,2,0,0,0):-1}
 u=scale(C0,-1);v=add(C0,scale(B0,-1))
 aa={(0,0,1,0,0):1};bb={(0,0,0,1,0):1};cc={(0,0,0,0,1):1}
 d=add(add(scale(aa,F(-1,5)),scale(bb,F(-1,3))),add(scale(cc,F(-1,3)),{(0,0,0,0,0):F(-1,15)}))
 Q=add(add(mul(aa,power(u,2)),mul(u,v)),add(add(mul(bb,u),mul(cc,v)),d))
 P={(0,0,0,0,0):1};out={}
 for m in range(1,6):P=mul(P,Q);out[m]=angular_integral(P)
 return out

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,help='Optional new JSON receipt path; existing files are preserved');a=ap.parse_args()
 if a.output is not None and a.output.exists():raise FileExistsError('receipt output already exists: '+str(a.output))
 data=CERT.read_bytes();assert hashlib.sha256(data).hexdigest()==EXPECTED_SHA
 c=json.loads(data);G=list(map(load,c['integer_generators']));M=list(map(load,c['integer_multipliers']));N=int(c['nonzero_rhs']);assert N>0
 moments=centered_dirichlet_moments();raw=raw_angular_moments();assert moments==raw,'literal Hopf/Dirichlet moment coefficients disagree'
 for m,g,k in zip(range(2,6),G,c['generator_scalings_to_centered_moments']):
  assert g==scale(moments[m],F(k)),('moment numerator scaling',m)
 exact_certificate(G,M,N)
 bad=[dict(z) for z in M];bad[0][(0,0,0)]=bad[0].get((0,0,0),0)+1
 try:exact_certificate(G,bad,N)
 except AssertionError:pass
 else:raise AssertionError('altered multiplier accepted')
 # Further raw B0,C0 mixed moments, independently evaluated by the direct coefficient sum.
 mixed=0
 for m in range(11):
  for n in range(11-m):
   actual=sum((F((-1)**(i+j)*comb(m,i)*comb(n,j)*comb(2*i+j,m+n),2*(i+j)+1) for i in range(m+1) for j in range(n+1) if 2*i+j>=m+n),F(0))
   expected=(-1)**(m+n)*rf(F(1,2),n)*__import__('math').factorial(m+n)/(F(__import__('math').factorial(n))*rf(F(3,2),m+n))
   assert actual==expected,(m,n,actual,expected)
   mixed+=1
 enc=lambda P:[{'powers':list(e),'coefficient':str(v)} for e,v in sorted(P.items())]
 receipt={'schema':1,'status':'PASS','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'certificate_sha256':EXPECTED_SHA,'moments':{str(m):enc(z) for m,z in moments.items()},'generator_scalings':c['generator_scalings_to_centered_moments'],'multiplier_degrees':c['multiplier_degrees'],'multiplier_terms':c['multiplier_terms'],'nonzero_rhs':str(N),'raw_mixed_moment_controls':mixed,'corrupted_multiplier_rejected':True,'scope':'Exact necessary equations at orders1--5 for this specific normalized five-term chart, checked against literal native Hopf expansion. The nonzero integer identity excludes all simultaneous solutions over every characteristic-zero field. No global quintic classification or finite-prefix extrapolation. The all-orders Dirichlet law is separately written.','formal_verification':False,'independent_replay':False,'external_review':False}
 if a.output is not None:a.output.write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':'PASS','mixed_controls':mixed,'certificate_sha256':EXPECTED_SHA}))
if __name__=='__main__':main()
