#!/usr/bin/env python3
"""Exact controls for the written all-degree clearing/automorphism proofs.

No universal topology, automorphism exhaustion or MZ decision is certified.
"""
from pathlib import Path
from math import factorial,gcd
import hashlib,json
import sympy as s

ROOT=Path(__file__).resolve().parent
X,Y,F,G,eta,t,b=s.symbols('X Y F G eta t b')
V=(X,Y,F,G)

def haar(poly):
 out=0
 for (a,c,d,e),coef in s.Poly(s.expand(poly),*V).terms():
  if d==c and e==a:out+=coef*(-1)**c*s.Rational(factorial(a)*factorial(c),factorial(a+c+1))
 return s.simplify(out)

def delta(poly):return s.expand(F*s.diff(poly,X)+G*s.diff(poly,Y))

def weight_set(poly,weights):
 return sorted({sum(a*w for a,w in zip(exp,weights)) for exp,coef in s.Poly(poly,X,Y,F,G,eta).terms() if coef})

def cyclic_genus(degree,orders):
 rhs=-2*degree+sum(degree-gcd(degree,abs(e)) for e in orders)
 assert rhs%2==0
 return 1+rhs//2

def main():
 images=ROOT/'quadratic_trace_polynomial_images.json'
 frozen=json.loads(images.read_text());cases={}
 for name,n,m in [('2O',8,12),('2I',12,20)]:
  r=n+m-2;weights=(m,n,-n,-m,-r)
  K=F**3*G/36-3*G**3 if name=='2O' else 12*F**5+G**3/144
  assert weight_set(X*G-Y*F,weights)==[0]
  assert weight_set(eta**2-K,weights)==[-2*r]
  assert weight_set(F,weights)[0]-weight_set(X,weights)[0]==-(n+m)
  assert weight_set(G,weights)[0]-weight_set(Y,weights)[0]==-(n+m)
  assert m-r==2-n and n-r==2-m
  # A mixed-parity expression, plus the extremal sharpness monomial.
  a=X**3+eta*Y**2+X*F*G+4
  derivs=[a]
  while derivs[-1]!=0:derivs.append(delta(derivs[-1]))
  assert len(derivs)-2==3
  assert max(weight_set(a,weights))==3*m
  assert max(weight_set(eta**3*a,weights))==-3*(n-2)
  for j in range(1,4):
   assert s.expand(delta(eta**3*derivs[j-1])-eta**j*(eta**(3-j)*derivs[j]))==0
  assert s.expand(delta(delta(delta(X**3)))-factorial(3)*F**3)==0
  # All terms of the actual full highest images have ORIGINAL weight+n.
  image=sum(s.Rational(c)*s.prod(v**e for v,e in zip(V,exps)) for exps,c in frozen[name]['highest_image_even'])
  assert frozen[name]['highest_image_odd']==[]
  assert weight_set(image,weights)==[n]
  if name=='2O':
   assert s.expand(K.subs({F:3*t**2,G:t**3/2},simultaneous=True))==0
   assert s.cancel((K-b**2).subs(F**3,108*G**2+36*b**2/G))==0
   assert s.discriminant(108*G**3+36*b**2,G)!=0
   orders=[1,1,1,-1,-2];genus=cyclic_genus(3,orders)
   assert genus==3
  else:
   assert s.expand(K.subs({F:t**3,G:-12*t**5},simultaneous=True))==0
   assert s.expand((K-b**2).subs(G**3,144*b**2-1728*F**5))==0
   assert s.discriminant(144*b**2-1728*F**5,F)!=0
   orders=[1,1,1,1,1,-5];genus=cyclic_genus(3,orders)
   assert genus==4
  assert s.gcd(K,F)==1
  cases[name]={'n':n,'m':m,'r':r,'weights':list(weights),
   'generic_eta_base_curve_genus':genus,'cyclic_ramification_orders':orders,
   'sharp_universal_exponent':'depth d','sharp_negative_weight_gap_per_depth':n-2,
   'mixed_expression_depth':3,'mixed_expression_cleared_max_weight':-3*(n-2),
   'actual_highest_image_positive_weight':n,'actual_highest_image_term_count':len(frozen[name]['highest_image_even'])}
 # Nonzero Haar pairings used in the period-pairing argument.
 polys=[X**2+Y*F+s.I*G,X+2*Y+F*G,1+X*G]
 norms=[]
 for a in polys:
  adj=s.conjugate(a).subs({s.conjugate(X):G,s.conjugate(Y):-F,s.conjugate(F):-Y,s.conjugate(G):X},simultaneous=True)
  value=haar(a*adj);assert value.is_positive
  norms.append(str(value))
 # The odd-coefficient period T and even Haar trace L obey L(a)=T(eta*a).
 # This checks the coefficient relation with nonconstant K and both parities.
 for K in [F**3*G/36-3*G**3,12*F**5+G**3/144]:
  a0=X+Y*F;a1=F+G*X
  product=s.rem(eta*(a0+eta*a1),eta**2-K,eta)
  assert s.Poly(product,eta).coeff_monomial(eta)==a0
  assert s.Poly(product,eta).coeff_monomial(1)==s.expand(K*a1)
 result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'sympy_version':s.__version__,'quadratic_images_sha256':hashlib.sha256(images.read_bytes()).hexdigest(),
  'cases':cases,'nonzero_Haar_pairing_norm_controls':norms,
  'all_degree_filtration_topology_and_automorphism_proofs_machine_verified':False,
  'exceptional_MZ_decisions':'OPEN on both models','independent_replay':False,'external_review':False,'formal_verification':False}
 (ROOT/'branch_clearing_verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))

if __name__=='__main__':main()
