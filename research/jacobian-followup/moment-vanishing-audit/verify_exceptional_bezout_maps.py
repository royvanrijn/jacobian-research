#!/usr/bin/env python3
"""Construct two specified regular Bezout maps and check weighted Haar controls.

Universal slice/topology/transfer proofs remain written, not machine verified.
"""
from pathlib import Path
from fractions import Fraction
from math import comb
import hashlib,json,sys
import sympy as s
import verify_haar_shear_witnesses as H
ROOT=Path(__file__).resolve().parent
x,y,z,w=s.symbols('x y z w')
def add_term(out,e,v):
 out[e]=out.get(e,Fraction(0))+v
 if not out[e]:del out[e]
def sparse(poly):return {e:Fraction(int(c.p),int(c.q)) for e,c in s.Poly(poly,x,y,z,w).terms() if c}
def D(poly):
 out={}
 for (a,b,c,d),v in poly.items():
  if a:add_term(out,(a-1,b,c+1,d),a*v)
  if b:add_term(out,(a,b-1,c,d+1),b*v)
 return out
def scale(poly,q):return {e:v*q for e,v in poly.items() if v*q}
def add(a,b):
 out=dict(a)
 for e,v in b.items():add_term(out,e,v)
 return out
def encode(poly):return [[list(e),str(v)] for e,v in sorted(poly.items())]
def case(name,f,g,n,m):
 N=n+m-1
 fc={b:Fraction(c) for (a,b),c in s.Poly(f,z,w).terms()};gc={b:Fraction(c) for (a,b),c in s.Poly(g,z,w).terms()}
 matrix=s.zeros(N+1)
 for j in range(m):
  for i,v in fc.items():matrix[i+j,j]=v
 for j in range(n):
  for i,v in gc.items():matrix[i+j,m+j]=v
 inverse=matrix.inv();assert matrix*inverse==s.eye(N+1)
 U={};V={}
 for i in range(N+1):
  factor=(-1)**i*comb(N,i)
  for j in range(m):
   q=-inverse[j,N-i]*factor
   if q:V[(N-i,i,m-1-j,j)]=Fraction(q)
  for j in range(n):
   q=inverse[m+j,N-i]*factor
   if q:U[(N-i,i,n-1-j,j)]=Fraction(q)
 F=sparse(f);G=sparse(g)
 assert add(H.mul(U,G),scale(H.mul(V,F),-1))==sparse((x*w-y*z)**N)
 du=D(U);dv=D(V)
 point=(1,0) if f.subs({z:1,w:0})!=0 else (1,1)
 value=Fraction(f.subs({z:point[0],w:point[1]}));htop={}
 for (a,b,c,d),v in du.items():add_term(htop,(a,b,0,0),v*point[0]**c*point[1]**d/value)
 assert du==H.mul(F,htop) and dv==H.mul(G,htop)
 c=s.expand((s.diff(f,z)*s.diff(g,w)-s.diff(f,w)*s.diff(g,z))/(n*m))
 cbottom=sparse(c);ctop=sparse(c.subs({z:x,w:y},simultaneous=True))
 first=next(iter(ctop));kappa=htop[first]/ctop[first]
 assert htop==scale(ctop,kappa)
 norm=H.haar(H.mul(cbottom,ctop));assert kappa*norm==1
 W=H.mul(cbottom,htop);assert H.haar(W)==1
 af=sparse((x*s.diff(f,z)+y*s.diff(f,w))/n)
 ag=sparse((x*s.diff(g,z)+y*s.diff(g,w))/m)
 nu=add(H.mul(V,af),scale(H.mul(U,ag),-1))
 determinant=sparse(x*w-y*z)
 assert D(nu)==add(H.mul(W,determinant),scale(sparse((x*w-y*z)**N),-1))
 tau=add(H.mul(U,G),H.mul(V,F));tau2=H.mul(tau,tau)
 htau2=H.haar(tau2);weighted=H.haar(H.mul(W,tau2));assert weighted==Fraction(1,3)
 assert H.haar(H.mul(W,H.mul(U,G)))==Fraction(1,2)
 # Original torus weights certify extraction of the first Hopf moment.
 assert {a+b-c-d for (a,b,c,d) in U}=={m}
 assert {a+b-c-d for (a,b,c,d) in V}=={n}
 assert {a+b-c-d for (a,b,c,d) in tau}=={0}
 assert 1-3*htau2!=0
 row={'n':n,'m':m,'r':n+m-2,'sharp_minimum_D_depth':N,'U_terms':len(U),'V_terms':len(V),'nu_terms':len(nu),
  'D_nu_identity_checked_mod_determinant':True,
  'kappa':str(kappa),'c_top_bottom_norm':str(norm),'h_W':str(H.haar(W)),
  'h_tau2':str(htau2),'h_weighted_tau2':str(weighted),'h_weighted_UG':'1/2',
  'h_pulled_Hopf_phase':str(1-3*htau2),'h_weighted_Hopf_phase':'0',
  'h_first_weighted_Hopf_marked_value':str(1-weighted)}
 artifact={'U':encode(U),'V':encode(V),'nu':encode(nu),'H_top':encode(htop),'source_f':encode(F),'source_g':encode(G),'source_c':encode(cbottom)}
 return row,artifact
def main():
 cases={};maps={}
 for name,f,g,n,m in [
  ('2O',z**8+14*z**4*w**4+w**8,z**2*w**2*(z**4-w**4)**2,8,12),
  ('2I',z*w*(z**10+11*z**5*w**5-w**10),z**20-228*z**15*w**5+494*z**10*w**10+228*z**5*w**15+w**20,12,20)]:
  cases[name],maps[name]=case(name,f,g,n,m)
 frozen=ROOT/'exceptional_bezout_maps.json'
 if '--freeze' in sys.argv:frozen.write_text(json.dumps(maps,indent=2)+'\n')
 else:assert maps==json.loads(frozen.read_text())
 result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'maps_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),'sympy_version':s.__version__,
  'imported_haar_arithmetic_sha256':hashlib.sha256((ROOT/'verify_haar_shear_witnesses.py').read_bytes()).hexdigest(),
  'cases':cases,'universal_slice_topology_and_weighted_transfer_proofs_machine_verified':False,
  'canonical_b0_purity':'EXCLUDED by exact nonzero first moment on both models',
  'other_b_purity':'OPEN:no family sweep performed','all_exceptional_MZ_decisions':'OPEN',
  'independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'exceptional_bezout_verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
