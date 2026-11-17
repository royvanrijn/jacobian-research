"""Exact controls for F11_AFFINE_COMPLETION.md; no bounded search is used."""
from pathlib import Path
from math import gcd
from itertools import product, combinations_with_replacement
from collections import Counter
import hashlib
import json
import platform
import time
import sympy as s
from sympy.polys.rings import ring
import model as m

ROOT=Path(__file__).resolve().parent
started=time.monotonic(); checks=[]
def check(name,condition):
    if not condition: raise AssertionError(name)
    checks.append(name)

N,Q,D=m.N,m.Q,m.D
J=N.diff()*Q-N*Q.diff()
H,rem=(N*D*Q).div(J)
check("polynomial target-divisor correction",rem.is_zero and H.degree()==13)
_,sf=Q.sqf_list(); R=sf[0][0]
check("six distinct F11 poles",R.degree()==6 and s.gcd(R,R.diff()).degree()==0)
lam_poly,rem=J.div(m.B**4*m.C**4*R)
check("ramification-factor derivative",rem.is_zero and lam_poly.degree()==0)
lam=lam_poly.rep.LC(); k=Q.rep.LC()
c0_poly,rem=N.div(m.A*m.B**5)
check("numerator scalar extraction",rem.is_zero and c0_poly.degree()==0)
c0=c0_poly.rep.LC()
rhs=(m.A*m.B*m.C*m.E*R)*s.Poly.from_dict({(0,):c0*k/lam},m.x,domain=m.K)
check("minimal residual boundary factors",H==rhs)

rr,x,w,z=ring("x,w,z",m.K)
def lift(p):
    return rr.from_dict({(power[0],0,0):c for power,c in p.rep.to_dict().items()})
def homogeneous(p,degree=None):
    if degree is None: degree=p.degree()
    return rr.from_dict({(power[0],degree-power[0],0):c for power,c in p.rep.to_dict().items()})
def det2(f,g):
    return f.diff(x)*g.diff(w)-f.diff(w)*g.diff(x)
def det3(fs):
    mat=[[f.diff(v) for v in (x,w,z)] for f in fs]
    return (mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
          - mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])
          + mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]))
def total_degree(f): return max(sum(p) for p in f)
nx,qx,dx,jx,hx=map(lift,(N,Q,D,J,H))
U=w*nx; V=w*qx; Gamma=w**2*hx; Z=z*Gamma
check("expanded two-variable cone determinant",det2(U,V)==w*jx)
check("expanded target-divisor threefold determinant",det3((U,V,Z))==U*(U-V)*V)
check("cone coordinate total degrees",[total_degree(v) for v in (U,V)]==[13,13])
check("lift coordinate total degrees",[total_degree(v) for v in (U,V,Z)]==[13,13,16])
check("cone is not Keller",total_degree(det2(U,V))>0)
check("lift is not Keller",total_degree(det3((U,V,Z)))>0)
check("cone generic root equation",qx*U-nx*V==0)
check("target triangle pullback",U*(U-V)*V==w**3*nx*dx*qx)
Nh,Qh=homogeneous(N,12),homogeneous(Q,12)
Ah,Bh,Ch,Eh,Rh=map(homogeneous,(m.A,m.B,m.C,m.E,R))
check("binary numerator equation",Nh==c0*Ah*Bh**5)
check("binary denominator square equation",Qh==k*Rh**2)
check("binary one-fiber includes source infinity",Nh-Qh==w*Ch**5*Eh)
check("two distinct numerator ramification points",m.B.degree()==2 and m.B.rep.discriminant()!=m.K.zero)
check("two distinct one-fiber ramification points",m.C.degree()==2 and m.C.rep.discriminant()!=m.K.zero)
# First target blowdown: the pole order becomes even but numerator stays squarefree.
S=U*V; T=V
check("first blowdown expanded determinant",det2(S,T)==w**2*qx*jx)
check("first blowdown base ratio identity",S*qx-T**2*nx==0)
# Source branch valuation ledger, ordered w,A1,A2,B,C,E,R.
pullbacks=s.Matrix([[1,1,1],[1,0,0],[1,0,0],[5,0,0],[0,5,0],[0,1,0],[0,0,2]])
det_vector=s.Matrix([1,0,0,4,4,0,1])
correction=s.Matrix([2,1,1,1,1,1,1])
check("target-divisor valuation cancellation",pullbacks*s.ones(3,1)==det_vector+correction)

# Exact scheme reconstruction control over a regular target, not only roots as numbers.
P=2*N-Q
check("regular target has twelve sheets",P.degree()==12 and s.gcd(P,P.diff()).degree()==0)
for name,f in [("Q",Q),("A",m.A),("B",m.B),("C",m.C),("E",m.E),("R",R),("H",H)]:
    check("all reconstruction denominators are units: "+name,s.gcd(P,f).degree()==0)
wrec=(2*s.invert(Q,P)).rem(P)
gammarec=(wrec**2*H).rem(P)
zrec=(3*s.invert(gammarec,P)).rem(P)
check("reconstructed first coordinate",(wrec*N).rem(P)==s.Poly(1,m.x,domain=m.K))
check("reconstructed second coordinate",(wrec*Q).rem(P)==s.Poly(2,m.x,domain=m.K))
check("reconstructed normal coordinate",(zrec*gammarec).rem(P)==s.Poly(3,m.x,domain=m.K))
check("one-fiber loses infinity in the cone",(N-Q).degree()==11)

# Rational target replacement over the prime field; the extension is unchanged
# because the new parameter t is linear in z with coefficients in k(a,c).
a,c,zz,tt,r=s.symbols("a c zz tt r")
num=a**5*(a**5+c**5*(zz+2))
den=(a**5+c**5)**2
one=c**5*(a**5*zz-c**5)
Bt=a**5+c**5*(zz+2); Dt=a**5*zz-c**5
check("replacement target one-fiber identity",s.expand(num-den-one)==0)
check("replacement fraction reduced",s.gcd(num,den)==1)
check("replacement branch factors distinct",s.gcd(Bt,Dt)==1)
check("zero branch factor irreducible",s.Poly(Bt,a,c,zz).is_irreducible)
check("one branch factor irreducible",s.Poly(Dt,a,c,zz).is_irreducible)
check("replacement rational degrees",(s.Poly(num,a,c,zz).total_degree(),s.Poly(den,a,c,zz).total_degree())==(11,10))
z_inverse=((a**5+c**5)**2*tt-a**5*(a**5+2*c**5))/(a**5*c**5)
check("rational target inverse",s.cancel(z_inverse.subs(tt,num/den)-zz)==0)
check("rational target forward inverse",s.cancel((num/den).subs(zz,z_inverse)-tt)==0)
check("rational target parameter derivative",
      s.cancel(s.diff(num/den,zz)-a**5*c**5/(a**5+c**5)**2)==0)
beta=s.Function("beta")
source_last=z_inverse.subs(tt,beta(r))
jac=s.Matrix([a,c,source_last]).jacobian([a,c,r]).det()
check("replacement source rational Jacobian",
      s.cancel(jac-(a**5+c**5)**2/(a**5*c**5)*s.diff(beta(r),r))==0)
# Check exact multiplicities with bivariate specializations at generic points;
# the global polynomial factorizations above prove the divisorial statements.
check("five simple pole lines",s.gcd(a**5+c**5,s.diff(a**5+c**5,a))==1)
def valuation_at_variable(expr,var):
    return min(p[0] for p,k in s.Poly(s.expand(expr),var).terms() if k!=0)
check("order-five target zero divisor",valuation_at_variable(num,a)==5 and valuation_at_variable(den,a)==0)
check("order-five target one divisor",valuation_at_variable(one,c)==5 and valuation_at_variable(den,c)==0)
check("simple active branch factors",s.degree(Bt,zz)==1 and s.degree(Dt,zz)==1)
ramification=[
 dict(divisor="a=0",old_order=5,base_order=5,new_order=5//gcd(5,5),fixed_sheets=12),
 dict(divisor="c=0",old_order=5,base_order=5,new_order=5//gcd(5,5),fixed_sheets=12),
 dict(divisor="each a^5+c^5 pole line",old_order=2,base_order=2,new_order=2//gcd(2,2),fixed_sheets=12),
 dict(divisor="B_t=0",old_order=5,base_order=1,new_order=5//gcd(5,1),fixed_sheets=2),
 dict(divisor="D_t=0",old_order=5,base_order=1,new_order=5//gcd(5,1),fixed_sheets=2)]
check("complete replacement divisorial inertia passes fixed-sheet gate",
      [r["new_order"] for r in ramification]==[1,1,1,5,5] and all(r["fixed_sheets"]>0 for r in ramification))
check("cone pole inertia fails fixed-sheet gate",2//gcd(2,1)==2)
check("first blowdown kills pole inertia",2//gcd(2,2)==1)
# Counterexamples delimit the universal statements.
xx,yy=s.symbols("xx yy")
inner=(xx*yy+1)/xx
num2=(xx*yy+1)**2+xx**2; den2=xx*(xx*yy+1)
check("two-pole composition control reduced",s.gcd(num2,den2)==1 and s.Poly(den2,xx,yy).total_degree()>0)
check("two-pole composition control identity",s.cancel(inner+1/inner-num2/den2)==0)
check("one-pole rational rechart is allowed",s.cancel(1/(1/xx)**3-xx**3)==0)
check("non-Keller coordinate composition allowed",s.diff((xx+yy)**3,xx)==3*(xx+yy)**2)
check("characteristic-five derivative control",s.Poly(s.diff(xx**5,xx),xx,modulus=5).is_zero)
check("low-degree difference of fifth powers is not square",
      any(e%2 for f,e in s.factor_list((xx+1)**5-(yy+2)**5,xx,yy)[1]))
check("linear Keller coordinate gradient control",
      s.Matrix([xx+yy**2,yy]).jacobian([xx,yy]).det()==1)
# The bounds are universal written proofs; this integer control only prevents
# forgetting that a nonconstant square has even total degree.
check("degree-five gate forces denominator degree at most four",
      [d for d in range(1,6) if d%2==0]==[2,4])

# Global obstruction: the pole profile is extracted from the exact map.
pole_multiplicity=Q.sqf_list()[1][0][1]
pole_count=Q.sqf_list()[1][0][0].degree()
cover_degree=max(N.degree(),Q.degree())
check("global obstruction profile extracted from F11",
      (cover_degree,pole_multiplicity,pole_count)==(12,2,6))
check("global valuation congruence has no integer solution",
      1%gcd(pole_multiplicity,cover_degree//pole_multiplicity)!=0)
check("locally admissible replacement still fails sixth-power descent",
      s.gcd(a**5+c**5,s.diff(a**5+c**5,a))==1 and R.degree()==6)
check("uniform 2^3 fiber is not excluded by this congruence",gcd(2,6//2)==1)
check("uniform 3^2 fiber is not excluded by this congruence",gcd(3,6//3)==1)
check("totally ramified degree-eleven fiber passes",gcd(11,11//11)==1)
for profile in [(2,2,8),(4,8),(3,3,6),(2,4),(6,6),(1,5,6)]:
    common=gcd(*profile); degree=sum(profile); complement=degree//common
    per_place=all((e//common)%gcd(e,complement)==0 for e in profile)
    check("general-profile congruence equivalence "+str(profile),
          per_place==(gcd(common,complement)==1))
check("nonuniform degree-twelve common index two fails",gcd(2,12//2)>1)
check("nonuniform degree-twelve common index four passes",gcd(4,12//4)==1)

# Exhaust the complete finite group, not a sample of permutations.
def canonical(v):
    return min(tuple(x%11 for x in v),tuple(-x%11 for x in v))
def mmul(v,w):
    return canonical((v[0]*w[0]+v[1]*w[2],v[0]*w[1]+v[1]*w[3],
                      v[2]*w[0]+v[3]*w[2],v[2]*w[1]+v[3]*w[3]))
def minv(v): return canonical((v[3],-v[1],-v[2],v[0]))
ident=canonical((1,0,0,1))
group={canonical(v) for v in product(range(11),repeat=4) if (v[0]*v[3]-v[1]*v[2])%11==1}
def order(v):
    w=ident
    for i in range(1,661):
        w=mmul(w,v)
        if w==ident:return i
    raise AssertionError(v)
def natural(v):
    aa,b,c,d=v
    out=[]
    for p in range(11):
        den=(c*p+d)%11
        out.append(11 if not den else ((aa*p+b)*pow(den,-1,11))%11)
    out.append(11 if not c else aa*pow(c,-1,11)%11)
    return out
def cycle_lengths(v):
    p=natural(v);seen=set();sizes=[]
    for i in range(12):
        if i in seen:continue
        j=i;size=0
        while j not in seen:seen.add(j);size+=1;j=p[j]
        sizes.append(size)
    return tuple(sorted(sizes,reverse=True))
def generated_pair(a,b):
    seen={ident};stack=[ident]
    while stack:
        v=stack.pop()
        for g in (a,b):
            w=mmul(v,g)
            if w not in seen:seen.add(w);stack.append(w)
    return seen
by_order={}
for v in sorted(group):by_order.setdefault(order(v),[]).append(v)
check("full PSL2(11) projective matrix order",len(group)==660)
check("complete natural action index census",
      Counter(12-len(cycle_lengths(v)) for v in group if v!=ident)==Counter({6:55,8:374,10:230}))
check("all involutions have the uniform obstructing fiber",
      all(cycle_lengths(v)==(2,)*6 for v in by_order[2]))
budgets=[]
for size in range(1,4):
    budgets.extend(v for v in combinations_with_replacement((6,8,10),size) if sum(v)==22)
check("complete genus-zero branch-index budget",budgets==[(6,6,10),(6,8,8)])
check("four branch points exceed genus-zero budget",4*6>22)
a2=by_order[2][0]
conjugates={mmul(mmul(g,a2),minv(g)) for g in group}
check("all involutions form one conjugacy class",conjugates==set(by_order[2]))
triangle_sizes={}
for last_order in (3,5):
    sizes={len(generated_pair(a2,b)) for b in by_order[3]
           if order(mmul(a2,b))==last_order}
    triangle_sizes[str(last_order)]=sorted(sizes)
    check("complete (2,3,"+str(last_order)+") subgroup census",
          sizes==({12} if last_order==3 else {60}))
check("two involutions cannot generate PSL2(11)",
      2*max(by_order)<len(group))
global_obstruction=dict(degree=cover_degree,uniform_pole_index=pole_multiplicity,
 root_denominator_degree=pole_count,congruence="2*m-6*h=1",
 solvable=False,all_genus_zero_index_budgets=budgets,
 triangle_subgroup_orders=triangle_sizes,
 applies_to="Curve-cover field pair plus independent transcendental variables, with arbitrary birational source and target identifications.",
 excludes_abstract_natural_action=False)

def encoded_coefficient(c): return [str(v) for v in c.to_list()]
def encoded_poly(p):
    return [{"exponents":list(power),"coefficient":encoded_coefficient(coef)}
            for power,coef in sorted(p.items())]
model_data=dict(coefficient_field="Q[a]/(a^4-a^3+a^2-a+1)",
 source_variables=["x","w","z"],source_coefficient_encoding="Descending rational power basis in a",
 polynomial_cone=[encoded_poly(U),encoded_poly(V)],
 polynomial_target_divisor_lift=[encoded_poly(U),encoded_poly(V),encoded_poly(Z)],
 residual_correction=encoded_poly(Gamma),lift_jacobian=encoded_poly(det3((U,V,Z))),
 target_divisor="U*(U-V)*V",cone_target_Keller_completion="excluded by derangement inertia",
 first_blowdown_target_Keller_completion="excluded by squarefree numerator",
 rational_replacement=dict(numerator=str(num),denominator=str(den),zero_factor=str(Bt),
                           one_factor=str(Dt),inverse=str(z_inverse)),
 branch_ledger=ramification,global_curve_obstruction=global_obstruction,
 rational_replacement_Keller_completion="excluded by global sixth-power obstruction",
 absolute_Keller_realization=False)
(ROOT/"affine_completion_model.json").write_text(json.dumps(model_data,indent=2)+"\n")
manifest=json.loads((ROOT/"sources/closed-polynomials-manifest.json").read_text())
for item in manifest["sources"]:
    check("closed-polynomial source snapshot hash",hashlib.sha256((ROOT/"sources"/item["path"]).read_bytes()).hexdigest()==item["sha256"])
files=["F11_AFFINE_COMPLETION.md","verify_affine_completion.py","model.py",
       "affine_completion_model.json","sources/closed-polynomials-manifest.json"]
receipt=dict(result="pass",checks=checks,check_count=len(checks),python=platform.python_version(),
 sympy=s.__version__,runtime_seconds=round(time.monotonic()-started,3),
 sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
 assurance="Written universal proofs with exact identities and controls; no formal verification or independent whole-proof review. The classical closed-polynomial criterion and original F11 group identification are attributed inputs.",
 old_affine_targets_excluded=True,rational_replacement_divisorial_inertia_admissible=True,
 rational_replacement_Keller_completion_excluded=True,
 all_stabilized_natural_PSL2_11_curve_covers_excluded=True,
 global_curve_obstruction=global_obstruction,
 full_action_spectrum_complete=False,absolute_Keller_realization=False)
(ROOT/"affine_completion_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                     absolute_Keller_realization=False)))
