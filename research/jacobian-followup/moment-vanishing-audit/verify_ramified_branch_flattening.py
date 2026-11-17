#!/usr/bin/env python3
"""Exact algebraic controls; global topology/completion arguments are written."""
from pathlib import Path
from math import factorial
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parent
X,Y,F,G,eta,k,l,a,b=s.symbols('X Y F G eta k l a b')
def haar(poly):
 out=0
 for (i,j,p,q),c in s.Poly(s.expand(poly),X,Y,F,G).terms():
  if p==j and q==i:out+=c*(-1)**j*s.Rational(factorial(i)*factorial(j),factorial(i+j+1))
 return s.simplify(out)
def chart(poly):return s.cancel(poly.subs(Y,(X*G-1)/F))
def delta(poly):return F*s.diff(poly,X)+G*s.diff(poly,Y)
def main():
 cases={}
 for name,n,m,K,lam in [('2O',8,12,F**3*G/36-3*G**3,-s.Rational(1,3)),('2I',12,20,12*F**5+G**3/144,s.Integer(144))]:
  r=X*G;E=lam*K*X**3
  J=Y*(r**2+r+1)-X**3*F**2*G/108 if name=='2O' else Y*(r**2+r+1)+1728*X**3*F**4
  assert chart(F*J-E+1)==0
  assert 3*m==2*(n+m-2)
  for power in [1,2]:
   S=sum(E**j for j in range(power));xp=X*E**power;yp=Y+X*G*J*S
   assert chart(xp*G-yp*F-1)==0
   assert chart(delta(xp)-(3*power+1)*E**power*F)==0
   assert chart(delta(yp)-(3*power+1)*E**power*G)==0
  controls=[]
  for aa,bb in [(0,0),(0,1),(1,0),(1,1),(2,0),(2,2)]:
   value=haar(E**aa*r**bb);assert value==s.Rational(1,3*aa+bb+1)
   controls.append([aa,bb,str(value)])
  # F=0 permits a single generic degree-three bottom monomial in K.
  assert s.expand(lam*K).subs(F,0)==G**3
  cases[name]={'n':n,'m':m,'lambda':str(lam),'regular_J_identity':True,
   'bounded_determinant_and_derivative_controls_k':[1,2],
   'bounded_scalar_moment_controls':controls,
   'generic_degree':'3*k+1','H3_scalar':1,
   'height_branch_order':'at least2*k','normalized_source_weight_branch_order':'at least2*k+1'}
 d=3*k+1
 assert s.expand((3*k+1)*(3*l+1)-(3*(k+l+3*k*l)+1))==0
 assert s.expand(3*(k+d*a+k*b)+b+1-d*(3*a+b+1))==0
 mean=s.factor(1-3*(1-4/(3*k+2)+4/(6*k+3)))
 assert s.simplify(mean+2*k*(6*k+1)/((3*k+2)*(2*k+1)))==0
 # Exact polynomial parameter identity, with Fj=e-1 and UG-VF=1.
 U,V,e,j,S,B=s.symbols('U V e j S B')
 # e^k-1=(e-1)S is the finite geometric identity used universally.
 z=s.symbols('z')
 assert s.expand((z-1)*(1+z+z**2)-(z**3-1))==0
 # In the formal limit,Fj=e-1 gives b+U_b*j/(1-e)=-U/F.
 formal=s.cancel(B+(U+F*B)*j/(1-e)+U/F)
 assert s.cancel(formal.subs(j,(e-1)/F))==0
 # The new weight is (3k+1)*e^k*W;it has unit normalized Haar
 # by the written global weighted-transfer proof,not by a finite jet.
 receipt={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':s.__version__,
  'cases':cases,'semigroup_index':'k+l+3*k*l','weighted_scalar_identity_denominator_checked':True,
  'unweighted_ambient_Hopf_first_moment':str(mean),
  'ambient_Hopf_purity_for_every_k_positive':'EXCLUDED by the displayed nonzero first moment',
  'admissible_etaX_image':'lambda^k*eta^(2k+1)*X^(3k+1),violates derivative order3k+1 for k>0',
  'formal_parameter_limit':'-U/F,independent of the starting b','formal_limit_has_genuine_global_pole':True,
  'source_unweighted_Hopf_purity':'OPEN:not inferred from ambient or weighted moments',
  'universal_Haar_pairing_topology_and_completion_proofs_machine_verified':False,
  'all_exceptional_MZ_decisions':'OPEN','independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'ramified_branch_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
