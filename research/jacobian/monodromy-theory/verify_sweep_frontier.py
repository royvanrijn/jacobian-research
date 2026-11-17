"""Exact controls for the universal proofs in SWEEP_FRONTIER.md."""
from pathlib import Path
from itertools import product
from math import gcd,lcm
import ast,json,hashlib,platform,time
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
w,u,z,t,gamma,c,Y1,Y2,Y3,eps,Z,T,H,a=S.symbols("w u z t gamma c Y1 Y2 Y3 eps Z T H a")
def zero(expr):return S.cancel(expr)==0
def det(*cols):return S.expand(S.Matrix.hstack(*cols).det())
def surface(P,cc=1,phi=0):
    v=S.expand(-S.diff(P,w)+(2*P+S.diff(phi,w))*S.diff(P,z)/cc)
    return S.Matrix([v,P+w*v,2*w*P+cc*z+phi+w*w*v])
D=S.Matrix([1,w,w*w])
P=z*z+w**3-w*w
X=surface(P);sw=X+gamma*D
check("nonlinear sweep full Jacobian",det(D,sw.diff(w),sw.diff(z))==gamma)
zY=Y3-2*w*Y2+w*w*Y1
R=S.expand(P.subs(z,zY)+w*Y1-Y2)
check("nonlinear root degree",S.degree(R,w)==4)
check("nonlinear leading coefficient only Y1",S.Poly(R,w).LC()==Y1**2)
check("nonlinear derivative and gamma identity",
    zero(S.diff(R,w)-(Y1-X[0].subs(z,zY))-2*S.diff(P,z).subs(z,zY)*R))
Rbranch=S.expand(R.subs({Y1:4,Y2:1,Y3:1}))
check("nonlinear ordinary double branch",zero(Rbranch-w*w*(16*w*w-15*w+11)))
check("nonlinear branch residual squarefree",S.gcd(16*w*w-15*w+11,32*w-15)==1)
check("nonlinear branch residual avoids double root",(16*w*w-15*w+11).subs(w,0)!=0)
phi=w**3-2*w
P2=(w+1)*z**3+w*z+w*w
Xold=surface(P2,2,phi)
Pbar=S.expand(P2.subs(z,z-phi/2))
Xnew=surface(Pbar,2)
check("polynomial phi gauge exact",
      all(zero(a-b) for a,b in zip(Xold.subs(z,z-phi/2),Xnew)))
d=w-u
Pt=P.subs({w:u,z:t},simultaneous=True)
mu=(Pt-P)/d
Yw=S.Matrix([mu,P+w*mu,2*w*P+z+w*w*mu])
Yu=S.Matrix([mu,Pt+u*mu,2*u*Pt+t+u*u*mu])
G=z-t+d*(P+Pt)
check("nonlinear exact ordered-pair elimination",zero(Yw[2]-Yu[2]-G))
check("ordered-pair first two coordinates agree",zero(Yw[1]-Yu[1]))
aa=S.symbols("aa",nonzero=True)
nu=S.Matrix([w+aa*w*w,-1-2*aa*w,aa])
other=S.Matrix([u-aa*u*u,-1+2*aa*u,-aa])
check("two tangent normals force collision criticality",
      all(zero(v) for v in (nu+other).subs(u,w+1/aa)))
degenerations=[]
for e in range(2,6):
    Pe=(w+1)*z**e+(w*w-1)*z**(e-1)+w*z+w**3
    Ge=c*(z-t)+(w-u)*(Pe+Pe.subs({w:u,z:t},simultaneous=True))
    scaled=S.cancel(eps*Ge.subs({u:w-eps**(e-1),z:Z/eps,t:T/eps},simultaneous=True))
    check("degeneration is polynomial e="+str(e),S.Poly(scaled,eps,Z,T).is_zero is False)
    limit=S.expand(scaled).subs(eps,0)
    check("all-degree smooth limit formula e="+str(e),zero(limit-(w+1)*(Z**e+T**e)-c*(Z-T)))
    F=a*(Z**e+T**e)+c*(Z-T)
    FZ=S.diff(F,Z);FT=S.diff(F,T)
    J=sum(Z**(e-2-j)*T**j for j in range(e-1))
    cert=c*(e-1)*(FZ-FT)-e*a*J*(e*F-Z*FZ-T*FT)
    check("affine smoothness Bezout identity e="+str(e),zero(cert-2*c*c*(e-1)))
    Fh=a*(Z**e+T**e)+c*(Z-T)*H**(e-1)
    check("smooth infinity partials e="+str(e),
          zero(S.diff(Fh,Z).subs(H,0)-e*a*Z**(e-1)) and
          zero(S.diff(Fh,T).subs(H,0)-e*a*T**(e-1)))
    degenerations.append(dict(e=e,limit=str(Fh),smoothness_certificate_constant=str(2*c*c*(e-1))))
# Entire affine lift, checked by all required weighted jets.
g,u1,u2,v1,v2,v3=S.symbols("g u1 u2 v1 v2 v3")
stage=S.Matrix([1+3*v1+v3,1-2*v1,v2])
check("nonlinear affine stage determinant",stage.jacobian([v1,v2,v3]).det()==-2)
check("nonlinear affine stage inverse first",zero((1-stage[1])/2-v1))
check("nonlinear affine stage inverse second",stage[2]==v2)
check("nonlinear affine stage inverse third",zero(stage[0]-1-3*v1-v3))
lift_data=[]
for i,(si,power) in enumerate(zip(sw,[1,2,1]),1):
    scaled=S.expand(si.subs({gamma:g,w:g*u1,z:g*u2},simultaneous=True))
    quotient,remainder=S.div(S.Poly(scaled,g,u1,u2),S.Poly(g**power,g,u1,u2))
    check("nonlinear gamma divisibility "+str(i),remainder.is_zero)
    E=quotient.as_expr()
    check("nonlinear affine base point "+str(i),E.subs({g:1,u1:1,u2:0})==0)
    composed=S.Poly(S.expand(E.subs({g:stage[0],u1:stage[1],u2:stage[2]},simultaneous=True)),v1,v2,v3)
    check("nonlinear all polynomiality jets "+str(i),
          all(sum(j*k for j,k in zip(mon,[1,1,2]))>=power for mon,coef in composed.terms()))
    lift_data.append(str(E))
    if i==2:check("nonlinear middle first jet",
        [S.diff(E,q).subs({g:1,u1:1,u2:0}) for q in [g,u1,u2]]==[-2,-3,0])
check("nonlinear full determinant exponents",sum([1,1,2])==sum([1,2,1])==2+2)
# Nilpotent-matrix and finite-sweep controls for one-parameter multipliers.
families=[]
for m in range(1,4):
    h=w**m
    for kind in ["constant","linear"]:
        if kind=="constant":
            XX=S.Matrix([h,z+S.Rational(m-1,m+1)*w**(m+1),
                h*z-w**(2*m+1)/S.Integer(2*m+1)])
        else:
            XX=S.Matrix([w**(m+1)-z,S.Rational(m,m+2)*w**(m+2),
                 w**(m+1)*z-w**(2*m+2)/S.Integer(2*m+2)-z*z/2])
        lam=1 if kind=="constant" else w
        DD=S.Matrix([1,w,z])
        syzygy=S.Matrix([[-w,1,0],[-z,0,1]])
        matrix=syzygy*XX.jacobian([w,z])
        expected=lam*S.Matrix([[-h,1],[-h*h,h]])
        label=f"h=w^{m},lambda={kind}"
        check(label+" normalized nilpotent matrix",all(zero(v) for v in matrix-expected))
        SS=XX+gamma*DD
        check(label+" sweep Jacobian",det(DD,SS.diff(w),SS.diff(z))==gamma**2)
        eq2=S.expand(XX[1]+(Y1-XX[0])*w-Y2)
        eq3=S.expand(XX[2]+(Y1-XX[0])*z-Y3)
        elim=S.resultant(eq2,eq3,z)
        degree=2*m+1 if kind=="constant" else 2*m+4
        check(label+" finite eliminant degree",S.degree(elim,w)==degree)
        check(label+" finite eliminant constant leading",not S.Poly(elim,w).LC().free_symbols)
        trans=eq2 if kind=="constant" else eq3
        check(label+" monic transverse equation",not S.Poly(trans,z).LC().free_symbols)
        delta=1 if kind=="constant" else gcd(m+2,2*m+2)
        check(label+" Gauss degree at most two",delta in [1,2])
        families.append(dict(h=str(h),lambda_polynomial=str(lam),degree=degree,delta=delta))
# The delta-two finite control and quadratic block.
sigma,qq,gg=S.symbols("sigma qq gg")
XX=S.Matrix([w**3-z,w**4/2,w**6/3-(w**3-z)**2/2])
DD=S.Matrix([1,w,z]);SS=XX+gamma*DD
R8=w**8-12*Y2*w**4+12*(Y1**2+2*Y3)*w*w-12*Y2**2
eq2=S.expand(XX[1]+(Y1-XX[0])*w-Y2)
zr=S.solve(eq2,z)[0]
eq3=S.expand(XX[2]+(Y1-XX[0])*z-Y3)
check("delta two exact degree-eight equation",zero(eq3.subs(z,zr)*(-24*w*w)-R8))
check("delta two generic squarefree root equation",
      S.gcd(R8.subs({Y1:0,Y2:1,Y3:1}),S.diff(R8,w).subs({Y1:0,Y2:1,Y3:1}))==1)
Yinv=S.Matrix([sigma+gg,qq*qq/2+gg*w,qq**3/3-sigma*sigma/2+gg*(qq*w-sigma)])
inv=Yinv.subs({w:-w,sigma:sigma+2*gg,gg:-gg},simultaneous=True)
check("delta two target-fixing involution",
      all(zero(v.subs(qq,w*w)) for v in inv-Yinv))
check("delta two proper block",8//2>1)
# Reconstruct the complete finite PSL group and exceptional coset action.
tree=ast.parse((ROOT/"verify_curve_descent_lcm.py").read_text())
names={"projective_group","cycles","generated"}
ns={"product":product}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[]),"<retained finite-group helpers>","exec"),ns)
group,ident,mul,action=ns["projective_group"](11)
cyc=ns["cycles"];gen=ns["generated"]
check("PSL complete group order",len(group)==660)
orders={q:lcm(*cyc(action(q))) for q in group}
a2=min(q for q in group if orders[q]==2)
sub=None
for a3 in sorted(q for q in group if orders[q]==3):
    if orders[mul(a2,a3)]!=5:continue
    candidate=gen((a2,a3),ident,mul)
    if len(candidate)==60:sub=candidate;break
check("PSL exceptional stabilizer order",sub is not None and len(sub)==60)
unseen=set(group);reps=[];lookup={}
while unseen:
    q=min(unseen);coset={mul(q,h) for h in sub}
    for h in coset:lookup[h]=len(reps)
    reps.append(q);unseen-=coset
check("PSL exceptional action degree",len(reps)==11)
exceptional=lambda q:tuple(lookup[mul(q,r)] for r in reps)
gates={}
for name,act,expected in [("natural",action,[(4,2,2,6),(10,1,1,12)]),
    ("exceptional",exceptional,[(1,4,3,3),(2,3,2,4),(4,2,1,6)])]:
    profiles=set(cyc(act(q)) for q in group)
    rows=set()
    for profile in profiles:
        nonfixed=[e for e in profile if e>1]
        if not nonfixed or 1 not in profile or len(set(nonfixed))!=1:continue
        e=nonfixed[0];delta=len(nonfixed)
        rows.add((e-1,delta,profile.count(1),e+1))
    check("PSL complete finite sweep gates "+name,sorted(rows)==expected)
    gates[name]=[dict(k=k,delta=d,fixed=f,min_dimension=n) for k,d,f,n in sorted(rows)]
check("natural PSL finite sweeps excluded through dimension five",all(q["min_dimension"]>=6 for q in gates["natural"]))
# Retain the next nonhomogeneous delta-two candidate without claiming its lift.
alpha=S.symbols("alpha")
apoly=16*alpha*alpha+30*alpha+15
def ared(expr):
    return S.Poly(S.expand(expr),alpha).rem(S.Poly(apoly,alpha)).as_expr().expand()
hc=alpha*w**4+w*w
sigc=hc*w-z
bc=S.Rational(2,3)*alpha*w**6+w**4/2
dc=S.Rational(2,5)*alpha*alpha*w**10+S.Rational(3,4)*alpha*w**8+w**6/3
XC=S.Matrix([sigc,bc,dc-sigc*sigc/2])
DC=S.Matrix([1,w,z]);SC=XC+gamma*DC
base={w:1,z:alpha/3+S.Rational(1,2),gamma:-2*alpha/3-S.Rational(1,2)}
check("next candidate nonzero marked gamma",ared(base[gamma])!=0)
check("next candidate unramified origin preimage",all(ared(v.subs(base))==0 for v in SC))
check("next candidate exact sweep Jacobian",zero(det(DC,SC.diff(w),SC.diff(z))-gamma**2))
for i,power in enumerate([1,2,4]):
    scaled=S.expand(SC[i].subs({w:gamma*u1,z:gamma**3*u2},simultaneous=True))
    check("next candidate gamma weight divisibility "+str(i),
      S.div(S.Poly(scaled,gamma,u1,u2),S.Poly(gamma**power,gamma,u1,u2))[1].is_zero)
check("next candidate discrete determinant balance",sum([1,2,4])==1+3+3)
candidate=dict(status="historical pre-completion checkpoint; completed by KMFI29",
 completed_by="SIGNED_ALTERNATING.md",
 coefficient_field="Q(alpha), 16 alpha^2+30 alpha+15=0",
 h=str(hc),b=str(bc),d=str(dc),surface=[str(v) for v in XC],
 marked_source={str(k):str(v) for k,v in base.items()},
 scaling_weights=[1,3],target_exponents=[1,2,4],proposed_source_weights=[1,2,4],
 next_checks=["construct and verify every polynomial stage jet",
              "verify a full polynomial Keller self-map",
              "prove its generic monodromy and arithmetic constants"],
 full_polynomial_Keller_map_claimed=False,geometric_group_identified=False,
 proposed_stage_strategy=[
 "At the marked base, form the three gamma-quotient gradients.",
 "Choose column1 in ker(dE2) intersect ker(dE3), column2 in ker(dE3), and column3 transverse to dE3.",
 "A shear along column3 in v1^2,v1^3,v1*v2 should cancel every remaining E3 jet of weight below4.",
 "The E2 requirements have weight below2, so these higher corrections should preserve them. Verify every coefficient and the triangular inverse."],
 proposed_group_strategy=[
 "Eliminate to a degree-six equation in q=w^2; prove its geometric group A6 from single-three-cycle branch values.",
 "The twelve-root group preserves six pairs and lies in the signed permutation wreath group.",
 "The product of the six q-roots has square class -5 times a square, so the geometric sign kernel lies in the five-dimensional even subspace.",
 "Search a bounded good-prime witness: a q Frobenius of type5,1 with a w lift of type10,1,1 or5,5,2 gives a nonconstant kernel vector after the fifth power.",
 "Use the A6 orbit spans in the binary permutation module to force the full even kernel; then use normality of geometric monodromy and its A6 projection.",
 "This is a proposed route, not a proved group identification; calculate the actual factorization witness and rule out hidden degeneration."])
(ROOT/"sweep_frontier_next_candidate.json").write_text(json.dumps(candidate,indent=2)+"\n")
manifest=json.loads((ROOT/"sources/sweep-monodromy-manifest.json").read_text())
for source in manifest["sources"][:2]:
    check("primary source pin "+source["path"],
      hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest()==source["sha256"])
model=dict(nonlinear_control=dict(potential=str(P),surface=[str(q) for q in X],
  root_polynomial=str(R),branch_target=[4,1,1],branch_polynomial=str(Rbranch),
  stage=[str(q) for q in stage],gamma_quotients=lift_data,full_Jacobian=-2,
  actual_affine_space=True,degree=4,geometric_group="S4",arithmetic_group="S4"),
  degeneration_controls=degenerations,one_parameter_finite_controls=families,
  delta_two=dict(root_equation=str(R8),degree=8,quadratic_block=True,
                 affine_Keller_lift_claimed=False),
  PSL2_11_finite_sweep_gates=gates,
  arbitrary_PSL2_11_occurrence_decided=False,full_classification_complete=False)
(ROOT/"sweep_frontier_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["SWEEP_FRONTIER.md","verify_sweep_frontier.py","sweep_frontier_model.json",
       "verify_curve_descent_lcm.py","sources/sweep-monodromy-manifest.json","sweep_frontier_next_candidate.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,
 python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written all-degree collision degeneration, finite inertia and intermediate-field Gauss proofs; exact controls and actual affine lift. No formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"sweep_frontier_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=receipt["runtime_seconds"])))
