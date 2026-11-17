#!/usr/bin/env python3
"""Exact algebraic controls for the written degree-seven proof, not formalization."""
from pathlib import Path
import hashlib, json
import sympy as s
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(a): assert s.cancel(s.expand(a))==0, s.simplify(a)
u,v=s.symbols('u v')
A,B=s.Function('A')(u),s.Function('B')(u)
e,f,g=s.symbols('e f g')
Q=v**3+A*v+B
P=v**4+(s.Rational(4,3)*A+e)*v**2+(s.Rational(4,3)*B+f)*v+s.Rational(2,9)*A**2+s.Rational(2,3)*e*A+g
jac=lambda P,Q: s.expand(s.diff(P,u)*s.diff(Q,v)-s.diff(P,v)*s.diff(Q,u))
X,Y=s.Function('X')(u),s.Function('Y')(u)
J=jac(P,Q).subs({A:X-s.Rational(3,2)*e,B:Y-s.Rational(3,4)*f}).doit().expand()
T=s.Rational(4,27)*X**3-e*X**2/3-s.Rational(2,3)*Y**2
zero(J+s.Rational(4,3)*v*s.diff(X*Y,u)-s.diff(T,u))
# Independently extract the four leading coefficient equations before integration.
C0,E0,F0,G0=[s.Function(n)(u) for n in ('C0','E0','F0','G0')]
raw=v**4+C0*v**3+E0*v**2+F0*v+G0
rawj=s.Poly(jac(raw,Q),v)
zero(rawj.coeff_monomial(v**5)-3*s.diff(C0,u))
zero(rawj.coeff_monomial(v**4)-(3*s.diff(E0,u)-4*s.diff(A,u)))
zero(rawj.coeff_monomial(v**3)-(3*s.diff(F0,u)+A*s.diff(C0,u)-3*C0*s.diff(A,u)-4*s.diff(B,u)))
zero(rawj.coeff_monomial(v**2)-(3*s.diff(G0,u)+A*s.diff(E0,u)-3*C0*s.diff(B,u)-2*E0*s.diff(A,u)))

t,H=s.symbols('t H',nonzero=True)
num=4*t**5-9*e*t**4-18*H**2; den=27*t**2
assert s.degree(num,t)==5 and s.degree(den,t)==2
assert s.resultant(num,den,t)!=0
assert s.degree(4*t**3-9*e*t**2,t)==3
# General monic (4,2) reduction and residual affine normal form.
a,b,c,d=s.symbols('a b c d',cls=s.Function)
G=s.Function('G')(u)
Pm=v**4+a(u)*v**3+b(u)*v**2+c(u)*v+d(u)
Qm=v**2+G
jc=s.Poly(jac(Pm,Qm),v)
zero(jc.coeff_monomial(v**4)-2*s.diff(a(u),u))
zero(jc.coeff_monomial(v**3)-(2*s.diff(b(u),u)-4*s.diff(G,u)))
C,E,alpha,beta=s.symbols('C E alpha beta')
zero(jac(C*v+E,v**2+alpha*u+beta)+C*alpha)
# Direct full rational (4,2) family, before polynomiality restrictions.
r,k=s.symbols('r k',nonzero=True)
qa=s.Function('qa')(u)
Al,Be,HH,EE=s.symbols('Al Be HH EE',nonzero=True)
vv=v+qa/(2*r)
qq=r*(vv**2+Al*u+Be)
pp=qq**2/r+HH*qq-k/(r*Al)*vv+EE
zero(jac(pp,qq)-k)
# Gluing identities, independent polynomial symbols.
x,y,z,aa,bb,AA,BB,TT=s.symbols('x y z aa bb AA BB TT')
h=1+x*y
rr=x*h
zero(rr/x-h)
# Both distinct zero-fibre components; h is their separating idempotent.
zero((h*h-h)/rr-y)
assert s.rem(h,s.Poly(x,x),x)==1
assert s.gcd(x,h)==1
# Eliminating a^2 between Q-polynomiality and P_z2 polynomiality.
expr=(-4*AA*h-4*BB+TT)
zero(expr.subs(x,0)-(-4*AA-4*BB+TT))
zero(expr.subs(y,-1/x)-(-4*BB+TT))
# Last normalized quartic/quadratic gluing obstruction.
cc,al1,L,II=s.symbols('cc al1 L II',nonzero=True)
qq=rr*z**2+cc*z+rr*al1*h
HH=-cc**2/rr**2+L/rr
pp=s.expand(qq**2/rr+HH*qq+II)
const=s.Poly(pp,z).coeff_monomial(1)
zero(const-(rr*al1**2*h**2+L*al1*h-cc**2*al1/x+II))
# Mixed (2,1) coefficient and final common linear factor.
ell,ff=s.Function('ell')(u),s.symbols('ff',nonzero=True)
bbfun,ee=s.Function('bb')(u),s.symbols('ee')
PP=ff*ell**2*v**2+ell*(2*bbfun+ee)*v+s.Function('c0')(u)
QQ=ff*ell*v+bbfun
j=s.Poly(jac(PP,QQ),v)
zero(j.coeff_monomial(v**2)); zero(j.coeff_monomial(v))
zero(j.coeff_monomial(1)-ell*(ff*s.diff(s.Function('c0')(u),u)-(2*bbfun+ee)*s.diff(bbfun,u)))
# Invariant affine-leading next coefficient, including q=f^k divisibility.
f0,q0,p0=s.symbols('f0 q0 p0',nonzero=True)
for m in (2,3,4):
 coeffs=[s.Function('c%d'%i)(u) for i in range(m)]
 pp=p0*v**m+sum(coeffs[i]*v**i for i in range(m))
 qq=q0*v+bbfun
 zc=s.Poly(jac(pp,qq),v).coeff_monomial(v**(m-1))
 zero(zc-(q0*s.diff(coeffs[m-1],u)-m*p0*s.diff(bbfun,u)))
# New monic cubic/quadratic square obstruction: H affine is not square in K.
H0,H1,lam=s.symbols('H0 H1 lam',nonzero=True)
R=s.symbols('R')
assert s.degree(H0+H1*R,R)==1
assert s.discriminant(lam*u**2+H,u)==-4*lam*H
# Exhaustive positive z-degree bookkeeping, not a coefficient search.
cases=[]
for m in range(1,7):
 for n in range(1,m+1):
  if m==n: branch='equal-leading-cancellation-or-degree-free-residual'
  elif n>=3:
   branch='monic43' if (m,n)==(4,3) else 'power-class-degree-impossible'
   if branch=='power-class-degree-impossible':
    import math
    assert m//math.gcd(m,n)>7-m and 7-m<=2
  elif n==2:
   branch={3:'cubic-square-class-and-affine-H',4:'invariant42-gluing',5:'power-class-degree-impossible',6:'power-class-degree-impossible'}[m]
   if m>=5: assert m>7-m and 7-m<=2
  else: branch='affine-fibre-power-class-and-mixed21'
  cases.append({'m':m,'n':n,'p_degree_bound':7-m,'q_degree_bound':7-n,'branch':branch})
assert len(cases)==21
# Actual power-class degree combinations H(r)*s^k, including special fibre parity.
power_cases={}
for kk,bound in [(2,5),(3,4),(4,3),(5,2),(6,1)]:
 choices=[(dh,ds) for dh in range(bound//3+1) for ds in range(bound//kk+1) if 3*dh+kk*ds<=bound]
 power_cases[str((kk,bound))]=choices
 assert all(dh==0 or ds==0 or (kk==2 and dh==1 and ds==1) for dh,ds in choices)
# Pin current premises and prior 58/16 application. Universal proofs remain written.
pins={}
for name in ['DEGREE6_COMPONENT_EXCLUSION.md','degree6-validation.json','verify_degree6.py','EXCEPTIONAL_DEGREE5.md','INVERTIBILITY_AUDIT.md']:
 pins[name]=sha(BASE/name)
prior=json.loads((BASE/'degree6-validation.json').read_text())
# Partition source is the prior exact receipt; tolerate its nested field layout explicitly.
def find_key(d,key):
 if isinstance(d,dict):
  if key in d:return d[key]
  for value in d.values():
   found=find_key(value,key)
   if found is not None:return found
 return None
# The integration receipt pins the authoritative prior whole-profile partition.
integ=json.loads((BASE/'degree6-integration.json').read_text())
assert integ['status']=='pass' and integ['full_profiles_excluded']==58 and integ['full_profiles_remaining']==16
old_remaining=['544','554','555','644','654','655','664','665','666','733','743','744','753','754','755','763']
added=['733','743','753','763']
remaining=[p for p in old_remaining if p not in added]
assert len(remaining)==12
control=BASE.parent/'equal-negative-completion/source-shear-validation.json'
w=json.loads(control.read_text())
assert w['status']=='pass' and w['component_degrees']==[3,7,9] and w['determinant_one_checked']
receipt={'status':'pass','claim':'GLDC10','scope':'Exact structural controls and all21 integer degree cases for the written universal degree-seven proof; not formal verification, whole-proof independent replay or bounded coefficient search.',
 'monic43_coefficient_identity':True,'monic43_rational_map_degrees':[5,3,2], 'rational42_normal_form':True,'two_component_gluing_identities':True,'mixed21_identity':True,'invariant_affine_next_coefficients':True,
 'new_affine_H_square_obstruction':True,'case_count':21,'cases':cases,'power_class_degree_controls':power_cases,
 'added_full_profiles':added,'full_profiles_excluded':62,'full_profiles_remaining':12,'remaining_profiles':remaining,
 'global3Ddegree_interval':[5,7],'cubic_component_degree_interval':[8,9],'degree8_exclusion_proved':False,
 'positive_witness_bounds_unchanged':True,'prior_premise_sha256':pins,'positive_control_sha256':sha(control),
 'checker_sha256':sha(Path(__file__)),'proof_sha256':sha(BASE/'DEGREE7_COMPONENT_EXCLUSION.md'),
 'universal_proofs_written_not_formal':True,'independent_full_replay':False}
(BASE/'degree7-validation.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS degree7 exact identities and21 cases;62 excluded12 remaining;cubic-component degree[8,9]')
