#!/usr/bin/env python3
"""Exact coefficient/weight controls; universal density proof is written."""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parent
u,T,q,X,Y,F,G=s.symbols('u T q X Y F G')
def CT(poly):return s.expand(poly).coeff(u,0)
def flat(poly):return s.integrate(poly,(T,0,1))
def main():
 P=(1+u)/u*(1-(1+u)**2*T**2)
 controls=[]
 for m in range(1,7):
  C=sum((-1)**j*comb(m,j)*comb(m+2*j,m)*T**(2*j) for j in range(m+1))
  B=sum((-1)**j*comb(m,j)*comb(m+2*j,m-1)*T**(2*j) for j in range(m+1))
  assert s.expand(CT(s.expand(P**m))-C)==0
  assert s.expand(CT(s.expand(u*P**m))-B)==0
  assert s.degree(C,T)==2*m and C.coeff(T,2*m)==(-1)**m*comb(3*m,m)
  cm=s.Rational(4**m*factorial(m)**2,factorial(2*m+1))
  assert flat(C)==0 and flat(B)==cm
  A=s.integrate((1-q**2)**m,(q,0,q))
  a=s.Poly(s.expand(A.subs(q,1+u)),u)
  assert all(a.coeff_monomial(u**j)==0 for j in range(1,m+1))
  R=s.Poly(s.expand((1+u)**(m-1)*A.subs(q,1+u)),u)
  assert R.coeff_monomial(u**m)==0 and R.coeff_monomial(u**(m-1))==cm
  controls.append({'m':m,'C_coefficients':[str(C.coeff(T,2*j)) for j in range(m+1)],'flat_C':'0','flat_B':str(cm)})
 v=2*F*G;uu=-2*X*Y;tau=X*G+Y*F
 pieces=[v,uu*v-2*tau**2,-3*uu*tau**2,-uu**2*tau**2]
 weights=[-1,0,1,2];depths=[]
 for part,k in zip(pieces,weights):
  der=s.expand(part);depth=0
  assert {a+b-c-d for (a,b,c,d),coef in s.Poly(der,X,Y,F,G).terms()}=={2*k}
  while True:
   der=s.expand(F*s.diff(der,X)+G*s.diff(der,Y))
   if not der:break
   depth+=1
  depths.append(depth)
 assert depths==[0,2,4,6]
 fractions=sorted({Fraction(j-k,j) for k,d in zip(weights,depths) for j in range(1,d+1)})
 assert fractions==[Fraction(-1),Fraction(0),Fraction(1,3),Fraction(1,2),Fraction(3,5),Fraction(2,3),Fraction(3,4),Fraction(1)]
 resonances={}
 for name,L,expected in [('2O',20,[-20,0,10,12,20]),('2I',32,[-32,0,16,24,32])]:
  vals=[int(L*f) for f in fractions if (L*f).denominator==1 and int(L*f)%2==0]
  assert vals==expected;resonances[name]=vals
 # Local critical-point lower-bound exponents are strictly below1.
 assert Fraction(2,3)<1 and Fraction(1,2)+Fraction(1,3)==Fraction(5,6)<1
 e=s.symbols('e',positive=True)
 # Area{|xy|<=e} in [-1,1]^2 for0<e<1.
 area=4*(e+s.integrate(e/q,(q,e,1)))
 assert s.simplify(area-4*e*(1-s.log(e)))==0
 assert s.limit(area/e,e,0,dir='+')==s.oo
 # The canonical height is not silently treated as compact-real.
 maps_path=ROOT/'exceptional_bezout_maps.json';maps=json.loads(maps_path.read_text())
 star_differences={}
 def dec(rows):return {tuple(a):Fraction(c) for a,c in rows}
 def add(a,b):
  out=dict(a)
  for k,v in b.items():
   out[k]=out.get(k,Fraction(0))+v
   if not out[k]:del out[k]
  return out
 import verify_haar_shear_witnesses as H
 for name,row in maps.items():
  t=add(H.mul(dec(row['U']),dec(row['source_g'])),H.mul(dec(row['V']),dec(row['source_f'])))
  star={(d,c,b,a):coef*(-1)**(b+c) for (a,b,c,d),coef in t.items()}
  diff=add(t,{e:-c for e,c in star.items()});assert diff
  star_differences[name]=len(diff)
 assert star_differences=={'2O':80,'2I':174}
 receipt={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':s.__version__,
  'bounded_coefficient_controls':controls,'universal_phase_piece_delta_depths':depths,
  'first_moment_resonance_ratios':[str(f) for f in fractions],'homogeneous_candidate_weights':resonances,
  'critical_mass_lower_bound_exponents':['2/3','5/6'],'saddle_band_mass':'4*epsilon*(1-log(epsilon)) on the unit square',
  'canonical_compact_star_difference_terms':star_differences,'frozen_maps_sha256':hashlib.sha256(maps_path.read_bytes()).hexdigest(),
  'imported_haar_arithmetic_sha256':hashlib.sha256((ROOT/'verify_haar_shear_witnesses.py').read_bytes()).hexdigest(),
  'all_order_algebra_and_universal_density_topology_proofs_machine_verified':False,
  'real_compact_height_branch':'EXCLUDED for noncyclic rotation symmetry','complex_height_branch':'OPEN',
  'five_candidate_weights_are_sufficient_for_purity':False,'mixed_weight_b_classification':'OPEN','exceptional_MZ_decisions':'OPEN on both models',
  'independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'hopf_height_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
