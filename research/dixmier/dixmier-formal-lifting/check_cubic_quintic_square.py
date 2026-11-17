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


# Explicit actual homogeneous quintic Weyl division models.
H, alpha, gamma, beta, kap = s.symbols("H alpha gamma beta kappa")
c = 4/r**2
cv = 2/r
def hom_mul(left,right,shift):
    out={}
    for i,a0 in left.items():
        for j,b0 in right.items():
            out[i+j]=out.get(i+j,0)+a0.subs(H,H+j*shift)*b0
    return clean(out)
def hom_inv(monomial,shift):
    assert len(monomial)==1
    n,coefficient=next(iter(monomial.items()))
    return {-n:s.factor(1/coefficient.subs(H,H-n*shift))}
def hom_power(op,n,shift):
    out={0:s.Integer(1)}
    for _ in range(n): out=hom_mul(out,op,shift)
    return out
def hom_comm(A0,B0,shift):
    return add(hom_mul(A0,B0,shift),scale(hom_mul(B0,A0,shift),-1))
target_p={1:s.Integer(1)}
target_q={-1:H+c*h}
assert hom_comm(target_p,target_q,-c*h)=={0:c*h}
target_U={5:1/(H*(H-beta*c*h))}
target_Z=scale(hom_mul(hom_inv(target_U,-c*h),{0:H},-c*h),1/(5*c))
assert hom_comm(target_U,target_Z,-c*h)=={0:h}
assert hom_mul(target_U,{0:H*(H-beta*c*h)},-c*h)=={5:s.Integer(1)}
target_U1=clean({n:a0.subs(beta,1) for n,a0 in target_U.items()})
assert hom_mul(hom_power(target_p,3,-c*h),
               hom_power(hom_inv(target_q,-c*h),2,-c*h),-c*h)==target_U1
source_kappa={5:cv/((H+alpha*h)*(H+gamma*h))}
source_K=scale(hom_mul(hom_inv(source_kappa,h),{0:H},h),-s.Rational(1,5))
assert hom_comm(source_kappa,source_K,h)=={0:h}
assert scale(hom_mul(source_kappa,{0:(H+alpha*h)*(H+gamma*h)},h),r/2)=={5:s.Integer(1)}
source_b={-1:H}
source_k21=clean({n:a0.subs({alpha:2,gamma:1}) for n,a0 in source_kappa.items()})
assert scale(hom_mul({3:s.Integer(1)},hom_power(hom_inv(source_b,h),2,h),h),cv)==source_k21

# The cubic transport and exact fifth-power equation in k(kappa)[H;5h*kappa*d/dkappa].
def delta_kap(a0,n):
    for _ in range(n): a0=5*h*kap*s.diff(a0,kap)
    return a0
def kap_mul(left,right):
    out={}
    for i,a0 in left.items():
        for j,b0 in right.items():
            for k in range(i+1):
                out[i+j-k]=out.get(i+j-k,0)+a0*comb(i,k)*delta_kap(b0,k)
    return clean(out)
def kap_comm(A0,B0):
    return add(kap_mul(A0,B0),scale(kap_mul(B0,A0),-1))
ym=(3-kap)/(kap-2)
Km={1:-1/(5*kap)}
Etam=scale(Km,-(kap-2)**2)
Um={0:s.factor(f.subs(y,ym))}
Zm=kap_mul({0:s.factor(1/fp.subs(y,ym))},Etam)
assert kap_comm({0:ym},Etam)=={0:h}
assert kap_comm(Um,Zm)=={0:h}
Rm=-(kap-3)*(kap-2)/kap**2
Jm=scale(kap_mul(Um,Zm),5*c)
assert Jm=={1:s.factor(c*Rm)}
Abeta=kap_mul(Um,kap_mul(Jm,add(Jm,{0:-beta*c*h})))
Aexpected=clean({
    2:c*c*Um[0]*Rm*Rm,
    1:c*c*h*Um[0]*(5*kap*Rm*s.diff(Rm,kap)-beta*Rm)
})
assert Abeta==Aexpected
Ahalf=clean({n:a0.subs(beta,s.Rational(5,2)) for n,a0 in Abeta.items()})
Aclosed=clean({
    2:c*c*(kap-3)**5/kap**4,
    1:c*c*s.Rational(5,2)*h*(kap-3)**4*(kap+12)/kap**4,
})
assert Ahalf==Aclosed
residue_beta=s.factor(s.limit((kap-2)*Abeta[1],kap,2))
assert s.factor(residue_beta-c*c*h*(2*beta-5)/8)==0
assert s.factor(residue_beta.subs(beta,s.Rational(5,2)))==0
VV,BB=s.symbols("v b")
classical_A=Ahalf[2].subs(kap,cv*VV**3/BB**2)*(VV*BB)**2
classical_p=cv*VV-3*BB**2/VV**2
assert s.factor(classical_A-classical_p**5)==0



# Conditional generic fourth connection: exact affine variations and quadratic curvature.
pvar,qvar=s.symbols("P Q")
Sexpr=r*pvar+s.Rational(4,3)/r
Texpr=r*qvar/4-Sexpr/(6*r)+s.Rational(4,27)/r**2
def Ds_derivative(a0):
    return s.factor(s.diff(a0,r)-4*pvar/(5*r)*s.diff(a0,pvar)-
                    6*qvar/(5*r)*s.diff(a0,qvar))
assert s.factor(Ds_derivative(Sexpr)-(Sexpr/(5*r)-8/(5*r**2)))==0
assert s.factor(Ds_derivative(Texpr)-(-Texpr/(5*r)+Sexpr/(10*r**2)))==0

def st_clean(op):
    return {key:s.factor(a0) for key,a0 in op.items() if s.factor(a0)!=0}
def st_add(*ops):
    out={}
    for op in ops:
        for key,a0 in op.items():out[key]=out.get(key,0)+a0
    return st_clean(out)
def st_mul(left,right):
    out={}
    for (i,j),a0 in left.items():
        for (k,l),b0 in right.items():
            for m in range(min(j,k)+1):
                factor=s.prod(k-n for n in range(m))
                key=(i+k-m,j+l-m)
                out[key]=out.get(key,0)+a0*b0*comb(j,m)*factor*(-h)**m
    return st_clean(out)
def st_comm(A0,B0):
    return st_add(st_mul(A0,B0),{key:-a0 for key,a0 in st_mul(B0,A0).items()})
Sop={(1,0):s.Integer(1)}
Top={(0,1):s.Integer(1)}
Fop={(1,1):1/(5*r),(0,1):-8/(5*r**2),(2,0):-1/(20*r**2)}
assert st_comm(Sop,Top)=={(0,0):h}
assert st_comm(Fop,Sop)=={(1,0):-h/(5*r),(0,0):8*h/(5*r**2)}
assert st_comm(Fop,Top)=={(0,1):h/(5*r),(1,0):-h/(10*r**2)}

srca,srcv=s.symbols("source_a source_v")
Sclass=4*srca-3*r*srca**2+2*srcv
Tclass=(r*srca**3-srca**2-srcv*srca)/2
Ds_minus_d=srcv**2*srca/(5*r)+srcv**2/(5*r**2)
Fclass=Sclass*Tclass/(5*r)-8*Tclass/(5*r**2)-Sclass**2/(20*r**2)
Dlong_minus_d=srca**4/4-3*r*srca**5/10+srcv*srca**3/2
assert s.factor(Ds_minus_d+Fclass-Dlong_minus_d)==0


names=[
    "CLASSICAL_QUINTIC_COVER_SQUARE.md",
    "check_cubic_quintic_square.py",
    "QUANTUM_DIVISION_AND_COMPLETION.md",
    "QUANTUM_RADICAL_TOWERS_AND_LONG_OBSTRUCTION.md",
    "QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md",
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
    "actual_quantum_quintic_models":{
        "target_U":"p^5 * [J*(J-beta*c*h)]^-1",
        "target_Z":"U^-1*J/(5*c)",
        "target_commutator":"[U,Z]=h",
        "source_kappa":"(2/r)*v^5*[(H+alpha*h)*(H+gamma*h)]^-1",
        "source_K":"-kappa^-1*H/5",
        "source_commutator":"[kappa,K]=h",
        "standard_ordering":"alpha=2,gamma=1 gives kappa=(2/r)*v^3*b^-2",
        "cubic_transport_J":"c*R(kappa)*H",
        "compatibility":"P^5=U*J*(J-beta*c*h)",
        "classical_root":"2*v/r-3*b^2/v^2",
        "beta_pole_residue":str(residue_beta),
        "pole_canceling_beta":"5/2",
        "A_at_beta_5_over_2":{str(i):str(a0) for i,a0 in Ahalf.items()},
        "individual_Weyl_division_models":True,
        "individual_quintic_indices":"5 on both sides, written finite-basis/deck proof",
        "conditional_relative_index":"3 on both sides, written dimension composition",
        "compatible_root_P_constructed":False,
        "displayed_root_gate":"All rational target roots and actual Euler pairs excluded by QQL1 exact spectrum",
        "finite_global_admission":False,
    },
    "conditional_fourth_connection":{
        "positive_coefficient_field":"L(r), L=K0((h)), not K0(r)((h))",
        "ordering_parameter_scope":"hregular and r-independent for full connection",
        "Ds":"d+v^2*a/(5*r)+v^2/(5*r^2)",
        "delta_S":"S/(5*r)-8/(5*r^2)",
        "delta_T":"-T/(5*r)+S/(10*r^2)",
        "F":"S*T/(5*r)-8*T/(5*r^2)-S^2/(20*r^2)",
        "exact_quadratic_commutators":"PASS",
        "classical_D_minus_d":"a^4/4-3*r*a^5/10+v*a^3/2",
        "classical_fourth_image_exact_match":True,
        "conditional_formula_only":True,
        "displayed_root_nonexistence":"QQL1 written proof",
    },
    "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
              for name in names},
    "coverage_boundary":(
        "Exact finite rational Weyl pair and source/target ordered radicand; "
        "classical cubic depression, discriminant, rational inverse, Poisson "
        "brackets, quintic identities, actual homogeneous Weyl division models, "
        "cubic transported radicand and beta=5/2 pole cancellation. "
        "All-field-degree and all-conjugate "
        "statements use written proofs. No actual quantum fifth root, Weyl "
        "compatible full Long division model, global polynomial admission or bounded Long lift."
    ),
    "independent_replay":False,"formal_verification":False,"external_review":False,
}
out=Path(tempfile.mkdtemp(prefix="qcs1-check-"))/"cubic_quintic_square_controls.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"claim":"QCS1","status":"PASS",
                  "actual_cubic_pair":True,"quantum_fifth_root":False,
                  "receipt":str(out)},indent=2))
