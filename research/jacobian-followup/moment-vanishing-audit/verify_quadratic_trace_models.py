#!/usr/bin/env python3
"""Exact Klein charts, quadratic trace images and nonzero Haar controls.

The all-degree topology and torsor arguments remain written proofs.
Default replays frozen polynomial images; --freeze reconstructs that artifact.
"""
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
import hashlib,json,sys
import sympy as s
import verify_haar_shear_witnesses as H
from verify_degree12_four_weights import groups

ROOT=Path(__file__).resolve().parent
z,w,t,F,G,Q,X,Y=s.symbols('z w t F G Q X Y')
VARIABLES=(X,Y,F,G)

def encode(poly):
 return [[list(e),str(c)] for e,c in s.Poly(poly,*VARIABLES).terms() if c]

def haar(poly):
 return H.haar({e:Fraction(int(c.p),int(c.q)) for e,c in s.Poly(poly,*VARIABLES).terms() if c})

def regular(poly):
 """Recover a polynomial in the target SL2 ring from its F!=0 chart."""
 poly=s.cancel(poly);num,den=s.fraction(poly)
 coef,fac=s.factor_list(den);assert all(v==F for v,e in fac)
 shift=sum(e for v,e in fac);out=0
 for (k,),c in s.Poly(num,F).terms():
  e=k-shift;c/=coef
  if e>=0:out+=c*F**e
  else:
   q,r=s.div(c,(X*G-1)**(-e),X,G);assert r==0
   out+=q*Y**(-e)
 return s.expand(out)

def reduce_target(poly):return s.cancel(poly.subs(Y,(X*G-1)/F))

def case(name,f,g,n,m):
 field,gens,values=groups(name)
 def coefficients(poly):return {b:Fraction(int(c.p),int(c.q)) for (a,b),c in s.Poly(poly,z,w).terms()}
 H.invariant_form(field,gens,coefficients(f),n)
 H.invariant_form(field,gens,coefficients(g),m)
 assert s.gcd(f.subs(w,1),g.subs(w,1))==1
 assert f.subs({z:1,w:0})!=0 or g.subs({z:1,w:0})!=0
 h=s.expand((s.diff(f,z)*s.diff(g,w)-s.diff(f,w)*s.diff(g,z))/(n*m));d=n+m-2
 K=F**3*G/36-3*G**3 if name=='2O' else G**3/144+12*F**5
 assert s.expand(h**2-K.subs({F:f,G:g},simultaneous=True))==0
 assert n*m==2*len(values) and n%4==m%4==0 and d%4==2
 # Euler and determinant give the precise additive-torsor transition.
 af=(X*s.diff(f,z)+Y*s.diff(f,w))/n
 ag=(X*s.diff(g,z)+Y*s.diff(g,w))/m
 assert s.expand(af*g-ag*f-h*(X*w-Y*z))==0
 assert s.expand(z*s.diff(af,X)+w*s.diff(af,Y)-f)==0
 assert s.expand(z*s.diff(ag,X)+w*s.diff(ag,Y)-g)==0
 assert s.cancel((1+t**d)/((1-t**n)*(1-t**m))-(1-t**(2*d))/((1-t**n)*(1-t**m)*(1-t**d)))==0
 def invariant(poly):
  poly=s.Poly(poly,z,w);N=poly.total_degree();basis=[];names=[]
  for c in [0,1]:
   for a in range(N//n+1):
    b=(N-a*n-c*d)//m
    if b<0 or a*n+b*m+c*d!=N:continue
    basis.append(s.Poly(f**a*g**b*h**c,z,w));names.append(F**a*G**b*Q**c)
  exponents=sorted(set(poly.monoms())|{e for b in basis for e in b.monoms()})
  matrix=s.Matrix([[b.coeff_monomial(e) for b in basis] for e in exponents])
  rhs=s.Matrix([poly.coeff_monomial(e) for e in exponents])
  coefficients,parameters=matrix.gauss_jordan_solve(rhs);assert parameters.rows==0
  assert matrix*coefficients==rhs
  return sum(a*b for a,b in zip(coefficients,names))
 # Express the entire highest invariant f(top) on the regular quadratic model.
 numerator=s.Poly(s.expand(f.subs({z:n*f*t*z+s.diff(f,w),w:n*f*t*w-s.diff(f,z)},simultaneous=True)),t)
 image=0
 for (k,),coefficient in numerator.terms():image+=invariant(coefficient)*Q**k*X**k/(F**k*(n*F)**n)
 num,den=s.together(image).as_numer_denom()
 image=s.Poly(s.cancel(s.rem(num,Q**2-K,Q)/den),Q)
 a0=regular(image.coeff_monomial(1));a1=regular(image.coeff_monomial(Q))
 assert a1==0  # f has degree divisible by4; it is even under the deck involution.
 # Verify EVERY derivative-divisibility condition for this full polynomial.
 derivative=a0
 for j in range(1,n+1):
  derivative=s.expand(F*s.diff(derivative,X)+G*s.diff(derivative,Y))
  quotient=regular(reduce_target(derivative/K**((j+1)//2)))
  assert reduce_target(derivative-K**((j+1)//2)*quotient)==0
 assert s.expand(F*s.diff(derivative,X)+G*s.diff(derivative,Y))==0
 source=H.haar(H.mul({(n-i,i,0,0):c for i,c in coefficients(f).items()},
                     {(0,0,n-i,i):c for i,c in coefficients(f).items()}))
 target=haar(F*a0);assert target==source
 assert target==(Fraction(8,15) if name=='2O' else Fraction(25,1092))
 return {'n':n,'m':m,'h_degree':d,'group_order':len(values),'branch_polynomial':str(K),
         'f_coefficients':[[i,str(c)] for i,c in coefficients(f).items()],
         'g_coefficients':[[i,str(c)] for i,c in coefficients(g).items()],
         'highest_image_even':encode(a0),'highest_image_odd':encode(a1),
         'all_highest_image_derivative_conditions_verified':n,
         'nonzero_norm_direct_and_quadratic_trace':str(target)}

def main():
 frozen=ROOT/'quadratic_trace_polynomial_images.json'
 cases={
  '2O':case('2O',z**8+14*z**4*w**4+w**8,z**2*w**2*(z**4-w**4)**2,8,12),
  '2I':case('2I',z*w*(z**10+11*z**5*w**5-w**10),z**20-228*z**15*w**5+494*z**10*w**10+228*z**5*w**15+w**20,12,20)}
 if '--freeze' in sys.argv:frozen.write_text(json.dumps(cases,indent=2)+'\n')
 else:assert cases==json.loads(frozen.read_text())
 # A normalization control with Gamma trivial,f=z,g=w^2,h=w:
 # source matrix [[QX,Y],[F,Q]],Q^2=G, has determinant1.
 assert s.expand(Q*X*Q-Y*F-1-(X*G-Y*F-1)-X*(Q**2-G))==0
 # Both the residue formula and the quadratic trace formula are checked on
 # small monomials; the note proves the displayed formulas for all exponents.
 residue_controls=trace_controls=0
 for a in range(4):
  for b in range(4):
   for c in range(4):
    for d in range(4):
     expected=H.haar({(a,b,c,d):Fraction(1)})
     beta=sum(Fraction((-1)**(b-j)*comb(b,j),a+j+1) for j in range(b+1)) if c==b and d==a else Fraction(0)
     assert beta==expected;residue_controls+=1
     traced=Fraction(0) if (a+d)%2 else haar(X**a*Y**b*F**c*G**((a+d)//2))
     assert traced==expected;trace_controls+=1
 # The old circuits fail the first new regularity condition in every unary degree.
 v=2*F*G;uu=-2*(F*t)*(G*t-1/F);tau=2*F*G*t-1
 hopf=(1+uu)*(v-(2+uu)*tau**2)
 U,V,T=uu**2,v**2,tau**2
 normalizer=(1+U)*(V-T*(U+2)*(2-T*(U**2+2*U+2)))
 hp=s.Poly(s.expand(hopf),t);np=s.Poly(s.expand(normalizer),t)
 assert hp.degree()==6 and s.expand(hp.LC()+16*F**4*G**4)==0
 assert np.degree()==20 and s.expand(np.LC()-4096*F**12*G**12)==0
 for branch in [F**3-108*G**2,G**3+1728*F**5]:assert s.gcd(branch,F*G)==1
 receipt={
  'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'polynomial_images_sha256':hashlib.sha256(frozen.read_bytes()).hexdigest(),
  'sympy_version':s.__version__,'cases':{name:{k:v for k,v in row.items() if not k.startswith('highest_image_') and not k.endswith('_coefficients')} for name,row in cases.items()},
  'residue_monomial_controls':residue_controls,'trivial_group_quadratic_trace_controls':trace_controls,
  'old_circuit_slice_orders_and_leading_coefficients':{'Hopf':[6,'-16*F^4*G^4'],'normalizer':[20,'4096*F^12*G^12']},
  'imported_implementations_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['verify_haar_shear_witnesses.py','verify_degree12_four_weights.py','verify_octahedral_minimal_type.py']},
  'completed':True,'universal_torsor_topology_and_trace_proof_machine_verified':False,
  'all_exceptional_MZ_decisions':'OPEN','independent_replay':False,'formal_verification':False,'external_review':False}
 (ROOT/'quadratic_trace_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
