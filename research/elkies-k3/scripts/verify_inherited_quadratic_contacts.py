#!/usr/bin/env python3
"""Replay a fixed signed-basis contact attempt with exact polynomial witnesses."""
import json,hashlib,resource
from pathlib import Path
from math import isqrt
import sympy as S
resource.setrlimit(resource.RLIMIT_CPU,(30,30));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
root=Path(__file__).resolve().parents[2]
source=root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
packet=root/'artifacts/generated-results/elkies-k3-inherited-quadratic-contact-v1'
data=json.loads(source.read_text());inp=json.loads((packet/'input.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest()==inp['source_sha256']
t=S.symbols('t')
def pol(cs):return S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(cs)),t)
def coord(i,key):
 c=data['sections']['records'][i][key]
 n,d=pol(c['numerator_coefficients_low_to_high']),pol(c['denominator_coefficients_low_to_high'])
 assert S.gcd(n,d).degree()==0
 return n,d
def power(a,n,f):
 v=S.Poly(1,t,modulus=f.get_modulus())
 while n:
  if n&1:v=(v*a)%f
  a=(a*a)%f;n//=2
 return v
def exclude(f,p):
 assert S.isprime(p)
 cs=[int(c.p%p)*pow(int(c.q),-1,p)%p for c in f.all_coeffs()]
 assert cs[0]!=0
 ff=S.Poly.from_list(cs,t,modulus=p);z=S.Poly(t,t,modulus=p)
 assert S.gcd(ff,power(z,p*p,ff)-z).degree()==0
A=pol(data['weierstrass_model']['A_coefficients_low_to_high']);B=pol(data['weierstrass_model']['B_coefficients_low_to_high']);H=data['sections']['height_gram']
assert {(c['trace'],c['basis_section']) for c in inp['cases']}=={(i,j) for i in [1,2] for j in range(17)}
assert len(inp['cases'])==34
rows=[];candidates=[];low=high=0
for c in inp['cases']:
 i,j=c['trace'],c['basis_section'];X,D=coord(i,'X');Y,DY=coord(i,'Y');x,dx=coord(j,'X');y,dy=coord(j,'Y')
 NN=(3*x*x+A*dx*dx)**2*dy*dy-8*x*dx**3*y*y;DD=4*y*y*dx**4
 common=S.gcd(NN,DD);NN=NN.exquo(common);DD=DD.exquo(common)
 cp=S.gcd(D,DD).monic()
 if cp.degree()>0:assert i==j and cp==D.monic()
 N=X*DD-D*NN
 while True:
  g=S.gcd(N,D*DD)
  if g.degree()==0:break
  N=N.exquo(g)
 assert N.degree()==c['degree']
 heights=[int(H[i][i]+4*H[j][j]-4*sgn*H[i][j]) for sgn in [1,-1]]
 low+=sum(h<20 for h in heights);high+=sum(h>=20 for h in heights)
 if 'prime' in c:exclude(N,c['prime'])
 else:
  prod=S.Poly(1,t)
  for f in c['factors']:
   q=pol(f['coefficients']);assert q.LC()==1 and q.degree()==f['degree']
   prod*=q**f['multiplicity']
   if q.degree()>2:exclude(q,f['prime']);continue
   assert q.degree()==2 and i in [1,2] and j==6
   disc=S.discriminant(q);assert disc!=0
   if disc>0:assert isqrt(int(disc.p))**2!=disc.p or isqrt(int(disc.q))**2!=disc.q
   def red(n,d):return (n*S.invert(d,q))%q
   xx,yy,xt,yt=red(x,dx),red(y,dy),red(X,D),red(Y,DY)
   lm=red(3*xx*xx+A,2*yy);xd=(lm*lm-2*xx)%q;yd=(lm*(xx-xd)-yy)%q
   assert xt==xd and yt==yd and not yd.is_zero
   assert heights[0]==8 and heights[1]==40
   assert S.gcd(q,4*A**3+27*B**2).degree()==0
   candidates.append({'trace':i,'inherited_section':6,'sign':1,'quadratic':[str(z) for z in reversed(q.all_coeffs())],'distance_height':8,'outcome':'Exact half, but contact intersection forces reducibility'})
  assert prod==N.monic()
 rows.append({'trace':i,'basis_section':j,'signed_distance_heights':heights,'nonpole_degree':N.degree(),'common_pole_degree':cp.degree()})
# Original trace poles are quadratic and have half O. Their contact image has square0.
for i in [1,2]:
 X,D=coord(i,'X');q=S.sqf_part(D).monic();assert q.degree()==2 and D.monic()==q*q
 disc=S.discriminant(q);assert disc<0 or isqrt(int(disc.p))**2!=disc.p or isqrt(int(disc.q))**2!=disc.q
 assert H[i][i]==8
assert len(candidates)==2 and (low,high)==(32,36)
# C.C_R = h(T-2R)/4-1 for arithmetic genus3; two distinct contacts need >=4.
assert S.Rational(8,4)-1==1<4 and S.Rational(20,4)-1==4
out={'schema':'inherited-quadratic-contact-v1','status':'PASS','source_sha256':inp['source_sha256'],
 'input_sha256':hashlib.sha256((packet/'input.json').read_bytes()).hexdigest(),
 'traces':[1,2],'signed_cases':68,'below_height_cutoff':low,'at_or_above_cutoff':high,
 'cases':rows,'actual_quadratic_halves':candidates,
 'boundary':'Only quadratic closed contact points on these two trace-halving curves, with halves supplied by O or one signed member of the retained17-section basis. Diagonal contacts, higher MW combinations, other trace parities and other image genera are not excluded. Intersection interpretation is written; finite algebra is replayed. No new rank-gain cover is constructed.'}
path=packet/'result.json'
if '--write' in __import__('sys').argv:path.write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads(path.read_text())
print('PASS: 68 signed cases; two genuine quadratic halves, both geometrically degenerate for this construction')
