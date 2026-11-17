#!/usr/bin/env python3
"""Exact controls for QCS1; no quantum quintic extension is constructed."""
from pathlib import Path
from math import comb
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE = Path(__file__).resolve().parent
h,r,y,w,pp,qq,u,zeta,eta,aa,ss,tt = s.symbols(
    "h r y w p q u zeta eta a s t")

def clean(op):
    return {i:s.factor(c) for i,c in op.items() if s.factor(c)!=0}
def add(*ops):
    out={}
    for op in ops:
        for i,c in op.items(): out[i]=out.get(i,0)+c
    return clean(out)
def scale(op,c): return clean({i:c*a for i,a in op.items()})
def mul(left,right):
    out={}
    for i,a in left.items():
        for j,b in right.items():
            for k in range(i+1):
                out[i+j-k]=out.get(i+j-k,0)+a*comb(i,k)*h**k*s.diff(b,y,k)
    return clean(out)
def comm(A,B): return add(mul(A,B),scale(mul(B,A),-1))

f=-y**3/(y+1)
fp=s.factor(s.diff(f,y))
F={0:f}
T={1:1/fp}
Z=scale(T,-1)
Eta={1:s.Integer(-1)}
assert comm(T,F)=={0:h}
assert comm(F,Z)=={0:h}
assert comm({0:y},Eta)=={0:h}
assert mul({0:fp},Z)==Eta
assert mul({0:fp},T)=={1:s.Integer(1)}
assert s.factor(y**3+f*y+f)==0
dc=s.factor(s.discriminant(y**3+u*y+u,y))
assert dc==-u**2*(4*u+27)
assert s.Poly(4*u+27,u).degree()==1
assert s.gcd(s.Poly(4*u+27,u),s.Poly(4*u+27,u).diff()).degree()==0

# The actual ordered source and target quantum radicands agree.
radicand_target=scale(mul({0:1/f},Z),20/r**2)
radicand_source=scale(mul({0:1/(f*fp)},Eta),20/r**2)
assert radicand_target==radicand_source

# Depressing the retained Long cubic.
p_of_st=ss/r-s.Rational(4,3)/r**2
q_of_st=4*tt/r+2*ss/(3*r**2)-s.Rational(16,27)/r**3
zz=s.symbols("zz")
long_cubic=r*aa**3-2*aa**2+ss*aa+4*tt
depressed=s.factor(long_cubic.subs(aa,zz+2/(3*r))/r)
assert s.factor(depressed-(zz**3+p_of_st*zz+q_of_st))==0

def poisson_st(A,B):
    return s.factor(s.diff(A,ss)*s.diff(B,tt)-s.diff(A,tt)*s.diff(B,ss))
assert poisson_st(p_of_st,q_of_st)==4/r**2

# Rational base change and canonical quintic partner.
u_of_pq=pp**3/qq**2
w_of_pq=qq/pp
p_of_uw=u*w*w
q_of_uw=u*w**3
assert s.factor(u_of_pq.subs({pp:p_of_uw,qq:q_of_uw})-u)==0
assert s.factor(w_of_pq.subs({pp:p_of_uw,qq:q_of_uw})-w)==0
jac=s.det(s.Matrix([p_of_uw,q_of_uw]).jacobian([u,w]))
assert s.factor(jac-u*w**4)==0

def poisson_pq(A,B):
    return s.factor(4/r**2*(s.diff(A,pp)*s.diff(B,qq)-
                            s.diff(A,qq)*s.diff(B,pp)))
bracket_uw=poisson_pq(u_of_pq,w_of_pq)
assert s.factor(bracket_uw.subs({pp:p_of_uw,qq:q_of_uw})-
                4/(r*r*u*w**4))==0
zeta_uw=r*r*u*w**5/20
zeta_pq=s.factor(zeta_uw.subs({u:u_of_pq,w:w_of_pq}))
assert poisson_pq(u_of_pq,zeta_pq)==1
assert s.factor(w**5-20*zeta_uw/(r*r*u))==0

# Both source and target columns are cotangent canonical pairs.
zeta_yw=zeta_uw.subs(u,f)
eta_yw=fp*zeta_yw
pref=4/(r*r*f*fp*w**4)
def poisson_yw(A,B):
    return s.factor(pref*(s.diff(A,y)*s.diff(B,w)-s.diff(A,w)*s.diff(B,y)))
assert poisson_yw(f,zeta_yw)==1
assert poisson_yw(y,eta_yw)==1
assert s.factor(w**5-20*eta_yw/(r*r*f*fp))==0
assert s.factor((w*y)**3+p_of_uw*w*y+q_of_uw-
                w**3*(y**3+u*y+u))==0

# Full source rational inverse and correct relative Poisson form.
xi=w*y
a_yw=xi+2/(3*r)
v_yw=r*w*w*y*y*(2*y+3)/(2*(y+1))
kappa=s.factor(2*v_yw/(r*xi*xi))
assert s.factor((3-kappa)/(kappa-2)-y)==0
assert s.factor(xi/y-w)==0
assert s.factor(poisson_yw(a_yw,v_yw)-1/v_yw)==0
s_yw=r*p_of_uw.subs(u,f)+s.Rational(4,3)/r
assert s.factor((s_yw-4*a_yw+3*r*a_yw*a_yw)/2-v_yw)==0

names=[
    "CLASSICAL_QUINTIC_COVER_SQUARE.md",
    "check_cubic_quintic_square.py",
    "QUANTUM_DIVISION_AND_COMPLETION.md",
    "QUANTUM_RADICAL_TOWERS_AND_LONG_OBSTRUCTION.md",
]
receipt={
    "claim":"QCS1","status":"PASS",
    "environment":{"python":platform.python_version(),"sympy":s.__version__},
    "actual_finite_rational_pair":{
        "f":str(f),"f_prime":str(fp),
        "T_coefficient_of_b":str(s.factor(1/fp)),
        "commutator_T_f":"h","commutator_f_Z":"h",
        "cubic_discriminant":str(dc),
        "nonconstant_discriminant_factor":"4*u+27",
        "normal_closure":"S3, written degree/discriminant argument",
        "spatial_poles":["y=0","y=-1","2*y+3=0","infinity"],
        "Long_global_admission":False,
    },
    "classical_cover_square":{
        "depressed_p":str(p_of_st),"depressed_q":str(q_of_st),
        "bracket_p_q":"4/r^2",
        "u":"p^3/q^2","w":"q/p",
        "inverse_p":"u*w^2","inverse_q":"u*w^3",
        "Jacobian_pq_uw":"u*w^4",
        "zeta":"r^2*u*w^5/20",
        "bracket_u_zeta":"1","eta":"f_prime(y)*zeta",
        "bracket_y_eta":"1",
        "horizontal_degrees":[5,5],"vertical_degrees":[3,3],
        "degrees_proof":"Written cubic function-field and quintic Eisenstein arguments",
        "fixed_u_all_rational_conjugates":"zeta+g(u), written derivative argument",
    },
    "Long_source_inverse":{
        "a":"w*y+2/(3*r)",
        "v":"r*w^2*y^2*(2*y+3)/(2*(y+1))",
        "kappa":"2*v/(r*(a-2/(3*r))^2)",
        "y":"(3-kappa)/(kappa-2)","w":"(a-2/(3*r))/y",
        "bracket_a_v":"1/v",
    },
    "quantum_common_radicand":{
        "target":"20/r^2 * U^-1 * Z",
        "source":"20/r^2 * (f*f_prime)^-1 * eta",
        "normal_operator":{str(i):str(c) for i,c in radicand_target.items()},
        "ordered_equality":"PASS",
        "fifth_root_constructed":False,
        "Weyl_division_realization":False,
    },
    "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
              for name in names},
    "coverage_boundary":(
        "Exact finite rational Weyl pair and source/target ordered radicand; "
        "classical cubic depression, discriminant, rational inverse, Poisson "
        "brackets and quintic identities. All-field-degree and all-conjugate "
        "statements use written proofs. No actual quantum fifth root, Weyl "
        "division model, global polynomial admission or bounded Long lift."
    ),
    "independent_replay":False,"formal_verification":False,"external_review":False,
}
out=Path(tempfile.mkdtemp(prefix="qcs1-check-"))/"cubic_quintic_square_controls.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"claim":"QCS1","status":"PASS",
                  "actual_cubic_pair":True,"quantum_fifth_root":False,
                  "receipt":str(out)},indent=2))

