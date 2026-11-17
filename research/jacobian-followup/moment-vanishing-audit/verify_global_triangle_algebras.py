#!/usr/bin/env python3
"""Exact global-chart controls and two complete lowest-window linear spaces.

The universal algebra/valuation/filtration arguments remain written proofs.
"""
from pathlib import Path
from fractions import Fraction
from math import factorial
import hashlib,json,sys
import sympy as s
import verify_haar_shear_witnesses as H
ROOT=Path(__file__).resolve().parent
r,j,X,Y,F,G,t=s.symbols('r j X Y F G t')
def laurent(poly):
 out={}
 for term in s.expand(poly).as_ordered_terms():
  k=int(term.as_powers_dict().get(j,0));out[k]=out.get(k,0)+term/j**k
 return {k:s.expand(v) for k,v in out.items() if v!=0}
def member(poly,p,q):
 rows=laurent(poly)
 if not rows:return True
 for k,a in rows.items():
  if k>0 and s.rem(a,r**(p*k),r)!=0:return False
  if k<0 and s.rem(a,(r-1)**(-q*k),r)!=0:return False
 degree=max(s.degree(v,r) for v in rows.values());low=min(rows)
 for l in range(1,int(degree)+1):
  a=s.expand(sum(v.coeff(r,l)*j**k for k,v in rows.items())*j**(-low))
  if s.rem(a,(j-1)**((l+1)//2),j)!=0:return False
 return True
def haar(poly):return s.integrate(laurent(poly).get(0,0),(r,0,1))
def source_haar(poly):
 return H.haar({e:Fraction(int(c.p),int(c.q)) for e,c in s.Poly(s.expand(poly),X,Y,F,G).terms() if c})
def linear_space(p,q,a,b):
 degree=2*(a+b);basis=[]
 for k in range(-a,b+1):
  factor=r**(p*k) if k>0 else (r-1)**(-q*k) if k<0 else s.Integer(1)
  for l in range(degree-int(s.degree(factor,r))+1):basis.append(s.expand(factor*r**l*j**k))
 matrix=s.Matrix([[s.diff(sum(v.coeff(r,l)*j**k for k,v in laurent(poly).items()),j,z).subs(j,1) for poly in basis] for l in range(1,degree+1) for z in range((l+1)//2)])
 null=matrix.nullspace();assert all(matrix*v==s.zeros(matrix.rows,1) for v in null)
 return len(basis),len(null),matrix.rank()
def case(name,n,m,A,B,kappa,p,q,expected_norm,frozen):
 high=sum(s.Rational(c)*X**a*Y**b*F**c1*G**d for (a,b,c1,d),c in frozen[name]['highest_image_even'])
 high=s.expand(high.subs({X:r/G,Y:(r-1)/F},simultaneous=True));depth=n//2
 K=F**3*G/36-3*G**3 if name=='2O' else 12*F**5+G**3/144
 mid=s.cancel((F*G)**depth*s.diff(high,r,depth)/(factorial(depth)*K**(depth//2)))
 numerator,den=s.fraction(mid);coef,fac=s.factor_list(den);shifts={F:0,G:0}
 for v,e in fac:assert v in shifts;shifts[v]+=e
 poly=0
 for (a,b,c1),v in s.Poly(numerator,F,G,r).terms():
  aa=a-shifts[F];bb=b-shifts[G]
  assert aa%A==0 and bb==-B*(aa//A)
  kk=aa//A;poly+=v/coef/kappa**kk*r**c1*j**kk
 poly=s.expand(poly);assert member(poly,p,q) and member(poly**2,p,q)
 rows=laurent(poly);degree=max(s.degree(a,r) for a in rows.values());window=[min(rows),max(rows)]
 assert degree==n//2==2*(window[1]-window[0])
 assert haar(poly)==0 and haar(poly**2)==expected_norm
 # Independent direct monomial control in the ORIGINAL source SL2 ring.
 f_top=sum(s.Rational(c)*X**(n-i)*Y**i for i,c in frozen[name]['f_coefficients'])
 source_mid=s.Poly(s.expand(f_top.subs({X:X+t*F,Y:Y+t*G},simultaneous=True)),t).coeff_monomial(t**depth)
 assert source_haar(source_mid)==0 and source_haar(source_mid**2)==expected_norm
 a,b=(1,1) if name=='2O' else (1,2)
 columns,nullity,rank=linear_space(p,q,a,b);assert nullity==2
 assert window==[-a,b] and degree>0 # 1 and the actual mid span the full kernel.
 if name=='2O':assert (columns,rank)==(8,6)
 else:
  assert (columns,rank)==(14,12)
  assert linear_space(p,q,1,1)[1]==1
 old_long=1-2*r+r*j-(1-r)/j
 assert not member(old_long,p,q) and not member(2*r-1,p,q)
 return {'n':n,'m':m,'p':p,'q':q,'middle_radial_degree':int(degree),'middle_Laurent_window':window,
  'middle_Haar_mean':'0','middle_Haar_second_moment':str(expected_norm),'direct_source_mean_and_norm_match':True,
  'complete_lowest_window_linear_columns':columns,'complete_lowest_window_rank':rank,'complete_lowest_window_nullity':nullity,
  'old_unrestricted_Long_and_radial_height_admitted':False},[[k,[[list(e)[0],str(c)] for e,c in s.Poly(a,r).terms()]] for k,a in sorted(rows.items())]
def main():
 source=ROOT/'quadratic_trace_polynomial_images.json';data=json.loads(source.read_text());cases={};polys={}
 for name,n,m,A,B,kappa,p,q,norm in [('2O',8,12,3,2,s.Rational(1,108),4,3,s.Rational(112,3)),('2I',12,20,5,3,s.Integer(-1728),3,5,s.Rational(275,13))]:
  cases[name],polys[name]=case(name,n,m,A,B,kappa,p,q,norm,data)
 assert member(r**2*(j-1),2,3) and not member(r**2*(j-1),4,3)
 frozen=ROOT/'global_triangle_lowest_polarizations.json'
 if '--freeze' in sys.argv:frozen.write_text(json.dumps(polys,indent=2)+'\n')
 else:assert polys==json.loads(frozen.read_text())
 result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':s.__version__,
  'quadratic_images_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'lowest_polarizations_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),
  'imported_haar_arithmetic_sha256':hashlib.sha256((ROOT/'verify_haar_shear_witnesses.py').read_bytes()).hexdigest(),
  'cases':cases,'extra_octahedral_G0_condition_detected':True,
  'universal_global_regular_algebra_marker_and_filtration_proofs_machine_verified':False,
  'lowest_window_pure_null_locus':'zero only on both models',
  'larger_windows_and_complex_logarithmic_height_identity':'OPEN','all_original_phases_reduced_to_T_invariant':False,
  'all_exceptional_MZ_decisions':'OPEN','independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'global_triangle_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
