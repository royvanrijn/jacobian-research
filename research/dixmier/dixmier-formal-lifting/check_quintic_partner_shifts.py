#!/usr/bin/env python3
"""Exact controls for QPS1; all rational shifts are covered by the written proof."""
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE=Path(__file__).resolve().parent
H,h,c,cv,z,lam,mu,w,C=s.symbols("H h c c_v kappa lambda mu w C")
def clean(op):
    return {n:s.factor(a) for n,a in op.items() if s.factor(a)!=0}
def add(*ops):
    out={}
    for op in ops:
        for n,a in op.items():out[n]=out.get(n,0)+a
    return clean(out)
def scale(op,a):return clean({n:a*b for n,b in op.items()})
def mul(left,right):
    out={}
    for n,a in left.items():
        for m,b in right.items():
            out[n+m]=out.get(n+m,0)+a.subs(H,H+m*h)*b
    return clean(out)
def inverse(op):
    assert len(op)==1
    n,a=next(iter(op.items()))
    return {-n:s.factor(1/a.subs(H,H-n*h))}
def power(op,n):
    if n<0:return power(inverse(op),-n)
    out={0:s.Integer(1)}
    for _ in range(n):out=mul(out,op)
    return out
def comm(left,right):
    return add(mul(left,right),scale(mul(right,left),-1))

# Normal kappa coefficients: the pole split has the displayed leading terms.
R=-(z-3)*(z-2)/z**2
U=(z-3)**3/(z-2)**2
g_cases=[
 ("regular_constant",lam*h),
 ("regular_polynomial",lam*h*z**3),
 ("simple_pole",lam*h/z),
 ("double_pole",lam*h/z**2),
 ("triple_pole",lam*h/z**3),
 ("other_pole",lam*h/(z-3)**2),
 ("mixed_poles",lam*h/z**2+mu*h/(z-2)),
]
def order_zero(expr):
    num,den=s.fraction(s.cancel(expr))
    a=min(k[0] for k,_ in s.Poly(num,z).terms())
    b=min(k[0] for k,_ in s.Poly(den,z).terms())
    return a-b
pole_controls=[]
for label,g in g_cases:
    go=order_zero(g)
    Jg=s.factor(c*R*(H-5*z*g))
    actual_order=order_zero(Jg)
    if go>=-1:
        residue=s.limit(z*g,z,0)
        expected_order=-2
        expected=-6*c*(H-5*residue)
    else:
        ell=-go
        a=s.limit(z**ell*g,z,0)
        expected_order=-(ell+1)
        expected=30*c*a
    leading=s.factor(s.limit(z**(-actual_order)*Jg,z,0))
    assert actual_order==expected_order
    assert s.factor(leading-expected)==0
    pole_controls.append({"case":label,"g":str(g),
                          "kappa_order":actual_order,"leading_coefficient":str(leading)})

# Actual skew-Laurent leading coefficients for the first four pole orders.
Fs=H**2
kap={5:cv/Fs}
Rs=add({0:-s.Integer(1)},scale(power(kap,-1),5),scale(power(kap,-2),-6))
J0=scale(mul(Rs,{0:H}),c)
for ell in range(4):
    Jg=add(J0,scale(mul(Rs,power(kap,1-ell)),-5*c*lam*h))
    if ell<=1:
        theta=5*lam*h if ell==1 else 0
        expected=-6*c/cv**2*Fs.subs(H,H-5*h)*Fs.subs(H,H-10*h)*(H-theta)
        assert min(Jg)==-10 and s.factor(Jg[-10]-expected)==0
    else:
        q=ell+1
        expected=30*c*lam*h/cv**q
        for j in range(1,q+1):expected*=Fs.subs(H,H-5*j*h)
        assert min(Jg)==-5*q and s.factor(Jg[-5*q]-expected)==0

# A genuine signed rational norm yields actual ambient Weyl coordinates.
B=(H-2*h)/(H+3*h)
norm=s.factor(s.prod(B.subs(H,H+i*h) for i in range(5)))
wm={1:1/B}
pm=mul(inverse(wm),{0:H})
assert power(wm,5)=={5:s.factor(1/norm)}
assert comm(pm,wm)=={0:h}
assert mul(wm,{0:B})=={1:s.Integer(1)}
assert mul(wm,pm)=={0:H}
source_Fs=C*norm
source_kappa={5:cv/source_Fs}
assert source_kappa==scale(power(wm,5),cv/C)

# Integer signed divisors distinguish norm and nonnorm cases.
T=s.symbols("T")
phi5=s.Poly(T**4+T**3+T**2+T+1,T)
def mod_phi5(expr):
    num,den=s.fraction(s.cancel(expr))
    return s.rem(s.Poly(num,T)*s.invert(s.Poly(den,T),phi5),phi5).as_expr()
def integer_orbit_divisor(expr):
    num,den=s.fraction(s.cancel(expr))
    out=s.Integer(0)
    for p,sign in [(num,1),(den,-1)]:
        for factor,e in s.factor_list(p,H)[1]:
            fp=s.Poly(factor,H)
            if fp.degree()==1:
                root=s.cancel(-fp.all_coeffs()[1]/fp.all_coeffs()[0])
                if root.is_Integer:out+=sign*e*T**int(root)
    return s.expand(out)
norm_divisor=integer_orbit_divisor(norm.subs(h,1))
assert mod_phi5(norm_divisor)==0
nonnorm_divisor=integer_orbit_divisor(H**2)
assert mod_phi5(nonnorm_divisor)==2
for q in [3,4,5,7]:
    d=sum(T**(5*j) for j in range(1,q+1))*nonnorm_divisor
    assert mod_phi5(d)==2*q

# Rational ambient momentum translations in the norm coordinates.
def omul(left,right):
    out={}
    for i,a in left.items():
        for j,b in right.items():
            for n in range(i+1):
                out[i+j-n]=out.get(i+j-n,0)+a*comb(i,n)*h**n*s.diff(b,w,n)
    return clean(out)
def ocomm(A,B):return add(omul(A,B),scale(omul(B,A),-1))
kw=cv*w**5
Rw=s.factor(R.subs(z,kw))
Uw=s.factor(U.subs(z,kw))
ambient_controls=[]
for label,g in g_cases:
    gw=s.factor(g.subs(z,kw))
    correction=-5*kw*gw/w
    shifted_p={1:s.Integer(1),0:correction}
    assert ocomm(shifted_p,{0:w})=={0:h}
    transported=omul({0:c*Rw*w},shifted_p)
    expected=clean({1:c*Rw*w,0:-5*c*Rw*kw*gw})
    assert transported==expected
    assert ocomm(transported,{0:Uw})=={0:s.factor(-5*c*h*Uw)}
    ambient_controls.append({"case":label,"canonical_translation":True,
                             "J0_to_Jg":True,"U_eigenvalue":"-5*c*h"})
for n in [-2,-1,0,1,2]:
    assert s.factor(5*c*h*z*R*s.diff(U**n,z)+5*c*h*n*U**n)==0

names=[
 "QUINTIC_PARTNER_SHIFT_RIGIDITY.md",
 "check_quintic_partner_shifts.py",
 "QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md",
 "CLASSICAL_QUINTIC_COVER_SQUARE.md",
]
receipt={
 "claim":"QPS1","status":"PASS",
 "environment":{"python":platform.python_version(),"sympy":s.__version__},
 "pole_split":pole_controls,
 "actual_first_pole_orders":"ell=0,1,2,3 verified in finite skewLaurent algebra",
 "norm_coordinate":{
  "B":str(B),"Fs":"C*product B(H+i*h), i=0..4",
  "w":"v*B(H)^-1","p":"w^-1*H",
  "w_fifth_power":"C*v^5/Fs","canonical_commutator":"[p,w]=h",
  "recovers_original_v_H":True,
  "signed_norm_divisor":str(norm_divisor),"norm_character_mod_Phi5":"0",
  "nonnorm_H_squared_character":"2",
 },
 "ambient_translations":ambient_controls,
 "written_theorem_scope":{
  "Fs":"any nonzero rational function of H",
  "g":"any rational function of kappa",
  "scalar_eigenvalue_set":"5*c*h*Z",
  "Euler_factorization":"no actual P,Q with [P,Q]=c*h and Jg=Q*P+scalar",
  "target_fifth_roots":"no actual root P^5=U*G(Jg), G nonzero rational",
  "h_smallness_required":False,
  "formal_completion_covered":False,
 },
 "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
           for name in names},
 "coverage_boundary":(
  "Exact leading coefficients for regular/simple/higher partner poles, "
  "actual skewLaurent first pole cases, signed rational norm and Weyl "
  "coordinate identities, rational ambient momentum translations and U "
  "eigenrelations. All rational shifts are exhausted by the written "
  "divisor-character/norm dichotomy, fixed-field leading-coefficient "
  "argument and QQL1. No mixed source/partner correction or Long lift."
 ),
 "independent_replay":False,"formal_verification":False,"external_review":False,
}
out=Path(tempfile.mkdtemp(prefix="qps1-check-"))/"quintic_partner_shift_controls.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"claim":"QPS1","status":"PASS","pole_cases":len(pole_controls),
                  "ambient_translations":len(ambient_controls),"receipt":str(out)},indent=2))
