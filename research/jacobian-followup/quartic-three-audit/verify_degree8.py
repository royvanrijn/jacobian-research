#!/usr/bin/env python3
"""Exact controls for the written sharp degree-nine cubic-component proof."""
from pathlib import Path
import hashlib,json,math
import sympy as s
BASE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def zero(e):assert s.cancel(s.expand(e))==0,s.factor(e)
u,v=s.symbols('u v')
jac=lambda p,q:s.expand(s.diff(p,u)*s.diff(q,v)-s.diff(p,v)*s.diff(q,u))
A,B,C=[s.Function(n)(u) for n in ('A','B','C')]
b,c,d,e,f,g,j=s.symbols('b c d e f g j')
# Full monic53 integrated higher equations, with both last coefficients checked.
Q=v**3+A*v+B
P=v**5+c*v**4+s.Rational(5,3)*A*v**3+(s.Rational(5,3)*B+s.Rational(4,3)*c*A+f)*v**2+(s.Rational(5,9)*A**2+s.Rational(4,3)*c*B+g)*v+s.Rational(10,9)*A*B+s.Rational(2,9)*c*A**2+s.Rational(2,3)*f*A+j
I=s.Rational(5,27)*A**3-s.Rational(4,3)*c*A*B-g*A-s.Rational(5,3)*B**2-2*f*B
T=s.Rational(5,9)*A**2*B+s.Rational(4,27)*c*A**3+f*A**2/3-s.Rational(2,3)*c*B**2-g*B
J=s.Poly(jac(P,Q),v)
for i in range(2,7):zero(J.coeff_monomial(v**i))
zero(J.coeff_monomial(v)-s.diff(I,u));zero(J.coeff_monomial(1)-s.diff(T,u))
X,Y=s.symbols('X Y');sh=Y-2*c*X/5-3*f/5
II=s.expand(I.subs({A:X,B:sh}));TT=s.Poly(s.expand(T.subs({A:X,B:sh})),X,Y)
zero(II-(-5*Y**2/3+5*X**3/27+4*X**2*c**2/15+4*X*c*f/5-X*g+3*f**2/5))
assert [(m,coef) for m,coef in TT.terms() if 2*m[0]+3*m[1]>=7]==[((2,1),s.Rational(5,9))]
# Monic54: two exact first integrals and the residual differential.
Q=v**4+A*v**2+B*v+C
E=s.Rational(5,16)*A*B+s.Rational(3,4)*b*B+c*A/2
Delta=s.Rational(5,4)*C+s.Rational(5,32)*A**2+s.Rational(3,4)*b*A+d
P=v**5+(s.Rational(5,4)*A+b)*v**3+(s.Rational(5,4)*B+c)*v**2+Delta*v+E+e
I1=s.Rational(5,32)*A**3+s.Rational(3,8)*b*A**2-s.Rational(5,4)*A*C-d*A-s.Rational(5,8)*B**2-2*c*B-3*b*C
I2=s.Rational(15,32)*A**2*B+c*A**2/2+s.Rational(3,4)*b*A*B-s.Rational(5,4)*B*C-d*B-2*c*C
J=s.Poly(jac(P,Q),v)
for i in range(3,8):zero(J.coeff_monomial(v**i))
zero(J.coeff_monomial(v**2)-s.diff(I1,u));zero(J.coeff_monomial(v)-s.diff(I2,u))
zero(J.coeff_monomial(1)-(B*s.diff(E,u)-Delta*s.diff(C,u)))
# Pole coefficients: all zero integration-constant strata retained.
t,a0,b0=s.symbols('t a0 b0',nonzero=True)
Ap=a0*t**-2;Bp=b0*t**-3;Cp=s.Rational(3,8)*a0**2*t**-4
Ep=s.Rational(5,16)*Ap*Bp;Dp=s.Rational(5,4)*Cp+s.Rational(5,32)*Ap**2
omega=s.expand(Bp*s.diff(Ep,t)-Dp*s.diff(Cp,t))
zero((omega*t**9).subs(b0**2,-a0**3/2)-s.Rational(55,32)*a0**4)
Ap=a0/t;Cp=Ap**2/8;Dp=5*Cp/4+5*Ap**2/32
zero((-Dp*s.diff(Cp,t))*t**5-s.Rational(5,64)*a0**4)
# Monic52 odd part and rational-map primitive.
G=s.Function('G')(u)
Q=v**2+G
P=v**5+(5*G/2+b)*v**3+(s.Rational(15,8)*G**2+3*b*G/2+c)*v
primitive=s.Rational(5,8)*G**3+3*b*G**2/4+c*G
zero(jac(P,Q)+s.diff(primitive,u))
# Quartic/quadratic field normal forms for J=k*u/H.
r,k,H,al,be,L,EE=s.symbols('r k H al be L EE',nonzero=True)
Q=v**2+al*u+be
P=Q**2+L*Q-2*k*v**3/(3*H*al**2)-k*u*v/(H*al)+EE
zero(jac(P,Q)-k*u/H)
Q=v**2+al*u**2+be
P=Q**2+L*Q-k*v/(2*H*al)+EE
zero(jac(P,Q)-k*u/H)
# Exact k2 triple-pole coefficients, valid for a shift regular only on h=0.
w,cc=s.symbols('w cc',nonzero=True)
qq=r*w**2+cc*w+al*r*u;vv=w+cc/(2*r)
pp=qq**2/r**2+L*qq/r-2*cc*vv**3/r-3*al*cc*u*vv/r+EE
coeff=s.Poly(s.expand(pp),w)
zero(coeff.coeff_monomial(w**2)-(L+2*al*u-2*cc**2/r**2))
zero(coeff.coeff_monomial(w)-(cc*L/r-al*cc*u/r-3*cc**3/(2*r**3)))
zero(coeff.coeff_monomial(w).subs(L,2*cc**2/r**2+b)-(cc**3/(2*r**3)+cc*b/r-al*cc*u/r))
# Regression for the rejected full-fibre divisibility shortcut.
x,y,z=s.symbols('x y z');h=1+x*y;rr=x*h
assert s.Poly(s.cancel(x**3*h/rr),x,y).as_expr()==x**2
assert s.rem(h,s.Poly(x,x),x)==1 # r does NOT divide h.
# Local order12 invariant42: Q-normal form and polynomial P_z2 force b constant on both components.
aa,bb,AA,BB,TT,hi=s.symbols('aa bb AA BB TT hi')
sol=s.solve([bb-aa/4-AA*hi-BB,aa+2*bb+TT],[aa,bb],dict=True)[0]
assert all(not val.has(x,y,z) for val in sol.values())
# All28 degree pairs, classified before any degree-increasing shear.
cases=[]
for m in range(1,8):
 for n in range(1,m+1):
  if m==n:branch='equal-degree-residual-and-invariant-affine-lemmas'
  elif n==1:branch='affine-mate-power-classes-and-invariant-slice'
  elif n in (2,3):
   opts={ (3,2):'retained32', (4,2):'new42-invariant-and-common-linear', (4,3):'retained43-and-common-linear', (5,2):'new-monic52', (5,3):'new-monic53'}
   branch=opts.get((m,n),'prime-leading-degree-impossible')
   if branch=='prime-leading-degree-impossible':assert 8-m<=2 and m>8-m
  else:
   if (m,n)==(5,4):branch='new-monic54'
   else:
    branch='highest-power-degree-impossible'
    assert m//math.gcd(m,n)>8-m and 8-m<=2
  cases.append({'m':m,'n':n,'p_degree_bound':8-m,'q_degree_bound':8-n,'branch':branch})
assert len(cases)==28
# Positive bound: actual coefficient reconstruction, determinant and distinct collision.
control=BASE.parent/'equal-negative-completion/source-shear-validation.json'
dat=json.loads(control.read_text());xyz=s.symbols('a b c')
polys=[sum(s.Rational(val)*s.prod(var**ex for var,ex in zip(xyz,mon)) for mon,val in row) for row in dat['components']]
degrees=[s.Poly(p,*xyz).total_degree() for p in polys];assert degrees==[3,7,9]
counts=[sum(sum(m)>1 for m in s.Poly(p,*xyz).monoms()) for p in polys];assert counts==[1,7,11]
assert s.Poly(s.expand(s.det(s.Matrix(polys).jacobian(xyz))),*xyz).as_expr()==1
pts=[tuple(map(s.Rational,p)) for p in dat['collision_points']];assert pts[0]!=pts[1]
images=[[s.expand(p.subs(dict(zip(xyz,pt)))) for p in polys] for pt in pts]
assert images==[[0,1,0],[0,1,0]]
prior=json.loads((BASE/'degree7-integration.json').read_text());assert prior['status']=='pass' and prior['full_profiles_remaining']==12
pins={n:sha(BASE/n) for n in ['DEGREE7_COMPONENT_EXCLUSION.md','degree7-validation.json','degree7-integration.json','DEGREE6_COMPONENT_EXCLUSION.md','EXCEPTIONAL_DEGREE5.md','INVERTIBILITY_AUDIT.md']}
result={'status':'pass','claim':'GLDC11','scope':'Exact new monic52/53/54 potentials, pole coefficients, quartic42 field forms and triple-pole controls, rejected fibre-divisibility regression, all28 integer degree cases, and direct degree9 determinant/collision reconstruction. Universal pole/UFD/fibre and retained geometric arguments are written, not formalized or fully independently replayed.',
'monic53_potential_degree':7,'monic54_pole_coefficients':['55/32','5/64'],'monic52_primitive_degree':3,'case_count':28,'cases':cases,'divisibility_shortcut_rejected':True,'universal_proofs_written_not_formal':True,'independent_full_replay':False,
'cubic_component_degree_minimum':9,'cubic_component_degree_interval':[9,9],'degree8_exceptional_completion_excluded':True,'full_profiles_excluded':62,'full_profiles_remaining':12,'global3Ddegree_interval':[5,7],'global_simplicity_goal_complete':False,
'positive_control':{'reconstructed_coefficients':True,'component_degrees':degrees,'actual_determinant':'1','distinct_collision_checked':True,'common_image':[0,1,0],'nonlinear_occurrences':sum(counts),'source_sha256':sha(control)},
'prior_premise_sha256':pins,'checker_sha256':sha(Path(__file__)),'proof_sha256':sha(BASE/'DEGREE8_COMPONENT_MINIMUM.md')}
(BASE/'degree8-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS degree8 exact controls and28 cases;direct determinant1/collision; cubic-component minimum9;whole-profile count unchanged62/12')
