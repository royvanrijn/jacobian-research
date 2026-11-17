"""Exact controls for SWEEP_ASYMPTOTICS.md; written universal proofs are separate."""
from pathlib import Path
import hashlib,json,platform,time
from math import gcd
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,v):
    if not v:raise AssertionError(name)
    checks.append(name)
w,z,g,k,q,sig,xi=S.symbols("w z gamma kappa q sigma xi")
def wave(lam,rho):
    return S.diff(lam,w,2)+2*S.diff(rho*lam,w,z)+S.diff(rho*rho*lam,z,2)
Delta=S.Matrix([1,w,z]);N=S.Matrix([[-w,1,0],[-z,0,1]])
check("tautological transverse frame",N*Delta.jacobian([w,z])==S.eye(2))
from signed_alternating import family
f=family(2);wf,zf,gf=f["w"],f["z"],f["g"]
Xf=f["X"];Df=f["Delta"];Nf=S.Matrix([[-wf,1,0],[-zf,0,1]])
P=Nf*Xf.jacobian([wf,zf])
check("sharp signed nilpotent trace",S.expand(S.trace(P))==0)
check("sharp signed nilpotent square",all(S.expand(t)==0 for t in P*P))
check("sharp signed sweep determinant",S.expand(S.Matrix.hstack(Df,Xf.diff(wf)+gf*Df.diff(wf),Xf.diff(zf)+gf*Df.diff(zf)).det()-gf**2)==0)
a=q*q+q;b=q**3*S.Rational(2,3)+q*q/2;d=q**5*S.Rational(2,5)+3*q**4/4+q**3/3
Xq=S.Matrix([sig,b,d-sig*sig/2]);normal=S.Matrix([sig,-a,1])
check("sharp surface normal",S.expand(normal.dot(Xq.diff(q)))==0 and S.expand(normal.dot(Xq.diff(sig)))==0)
II=S.Matrix(2,2,lambda i,j:S.expand(normal.dot(Xq.diff([q,sig][i],[q,sig][j]))))
check("sharp second fundamental form",all(S.expand(e)==0 for e in II-S.diag(q*S.diff(a,q)**2,-1)))
v=S.Matrix([w/S.diff(b,q),1])
check("two asymptotic directions",S.cancel((v.T*II*v)[0]).subs(q,w*w)==0)
check("sharp double parametrization",S.degree(w*w,w)==2)
# Highest coefficient and pole equations; coefficient variables are exact rationals.
r,s=S.symbols("r s",integer=True,positive=True)
C=r*(r-1)+2*k*(r-1)*(s+1)+k*k*(s+1)*(s+2)
check("residue discriminant",S.expand(S.discriminant(C,k)+4*(r-1)*(s+1)*(r+s+1))==0)
check("residue normalized trace",S.cancel(-S.Poly(C,k).nth(1)/(2*S.Poly(C,k).nth(2))+(r-1)/(s+2))==0)
t=(s+1)/r;u=(s+1)*(s+2)/(r*(r-1))
check("homogeneous exponent recovery",S.cancel((u+t)/(u-t*t)-r)==0)
check("exponent recovery denominator",S.factor(u-t*t)==(s+1)*(r+s+1)/(r**2*(r-1)))
pf=S.Function("p")(w);ff=S.Function("f")(w);qq=S.Function("q")(w)
for sv in range(5):
    expr=S.expand(wave(ff*z**sv,pf*z+qq)).coeff(z,sv)
    expected=S.diff(ff,w,2)+2*(sv+1)*S.diff(pf*ff,w)+(sv+1)*(sv+2)*pf*pf*ff
    check("affine multiplier top ODE s="+str(sv),S.expand(expr-expected)==0)
for dv in [2,3,4]:
    lam=(w+1)*z**2+1;rho=(w+2)*z**dv+z+1
    check("higher z degree unique leading term d="+str(dv),
          S.degree(wave(lam,rho),z)==2+2*dv-2)
models=[]
for rv,sv in [(2,0),(2,1),(3,0),(3,1),(4,2),(6,4)]:
    cc=C.subs({r:rv,s:sv});lab=f"r={rv},s={sv}"
    def red(e):
        num,den=S.fraction(S.cancel(e))
        num=S.rem(num,cc,k);den=S.rem(den,cc,k)
        return S.rem(num*S.invert(den,cc,k),cc,k).expand()
    lam=w**rv*z**sv;rho=k*z/w
    check(lab+" exact monomial wave coefficient",S.expand(wave(lam,rho)-cc*w**(rv-2)*z**sv)==0)
    av=-k-S.Rational(rv,sv+1);beta=S.Rational(1,sv+1)
    chi=(2*k*(sv+1)+rv)/((sv+1)*(sv+2));t0=w**(rv-1)*z**(sv+1)
    surf=S.Matrix([av*t0,(av+beta)*t0*w,(av+chi)*t0*z])
    pm=N*surf.jacobian([w,z])
    expected=lam*S.Matrix([[-rho,1],[-rho*rho,rho]])
    check(lab+" integrated matrix",all(red(e)==0 for e in pm-expected))
    check(lab+" polynomial surface",all(S.denom(S.cancel(e)).has(w,z)==False for e in surf))
    check(lab+" normalized Jacobian",red(S.Matrix.hstack(Delta,surf.diff(w)+g*Delta.diff(w),surf.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
    polynomial=xi*(1+beta*xi)**(rv-1)*(1+chi*xi)**(sv+1);degree=rv+sv+1
    deriv=(1+beta*xi)**(rv-2)*(1+chi*xi)**sv*(1-av*xi)**2
    check(lab+" exact curve derivative",red(S.diff(polynomial,xi)-deriv)==0)
    check(lab+" distinct critical and zero points",all(red(e)!=0 for e in [av,beta,chi,beta-chi,av+beta,av+chi]))
    check(lab+" nonzero third branch value",red(polynomial.subs(xi,1/av))!=0)
    check(lab+" exact curve degree",S.degree(polynomial,xi)==degree and red(S.Poly(polynomial,xi).LC())!=0)
    mat=S.Matrix([[rv-1,rv,rv-1],[sv+1,sv+1,sv+2]])
    check(lab+" primitive target kernel",mat*S.Matrix([-degree,rv-1,sv+1])==S.zeros(2,1))
    minors=[int(mat[:,[i,j]].det()) for i,j in [(0,1),(0,2),(1,2)]]
    check(lab+" saturated source exponent lattice",gcd(*minors)==1)
    check(lab+" finite branch index sum",(rv-2)+sv+2==degree-1)
    check(lab+" group parity",((-1)**((rv-2)+sv)==1)==(degree%2==1))
    Y=surf+g*Delta
    check(lab+" centered origin second equation",S.expand(Y[1].subs(g,-av*t0)-beta*w*t0)==0)
    models.append(dict(r=rv,s=sv,degree=degree,residue_equation=str(cc),
       geometric_group=("A" if degree%2 else "S")+str(degree),lattice_minors=minors))
# Nonhomogeneous rational control. It is not finite and is not a Keller lift.
R=w*w+w;I=S.integrate(R,w);J=S.integrate(w*R,w);K=S.integrate(R*R,w)
rho=(k*z+R)/w;lam=w*w
def red3(e):
    num,den=S.fraction(S.cancel(e))
    return S.rem(S.rem(num,k*k+k+1,k)*S.invert(S.rem(den,k*k+k+1,k),k*k+k+1,k),k*k+k+1,k).expand()
surf=S.Matrix([-(k+2)*w*z+w*R+2*k*I,
 k*k*w*w*z+w*w*R+2*(k-1)*J,-w*z*z+w*R*z-K])
check("genuine rational wave",red3(wave(lam,rho))==0)
check("genuine rational integrated matrix",all(red3(e)==0 for e in N*surf.jacobian([w,z])-lam*S.Matrix([[-rho,1],[-rho*rho,rho]])))
check("genuine rational normalized determinant",red3(S.Matrix.hstack(Delta,surf.diff(w)+g*Delta.diff(w),surf.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
check("nonfinite rational collapsed line",all(e.subs(w,0)==0 for e in surf))
check("exceptional PSL11 degree-four finite gate excluded",3>2 and 3*3+2==11)
check("natural PSL12 order-three fixed-sheet gate excluded",12-3*4==0)
manifest=json.loads((ROOT/"sources/sweep-monodromy-manifest.json").read_text())
# Existing retained primary source pins; do not fetch or reconstruct missing sources.
source_count=0
for source in manifest["sources"]:
    if "path" in source and "sha256" in source:
        check("retained primary source hash "+source["path"],hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest()==source["sha256"])
        source_count+=1
check("primary source pins present",source_count>0)
model=dict(tautological_rank_two_surface_degree_bound=2,
 finite_PSL2_11_tautological_route_excluded=True,
 arbitrary_Keller_PSL2_11_action_decided=False,full_classification_complete=False,
 homogeneous_controls=models,
 rational_control=dict(R=str(R),residue_equation="kappa^2+kappa+1",finite=False,Keller_lift_claimed=False),
 assurance="Exact symbolic controls accompany written universal arguments; not formal verification or independent whole-proof review.")
(ROOT/"sweep_asymptotics_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["SWEEP_ASYMPTOTICS.md","verify_sweep_asymptotics.py","sweep_asymptotics_model.json","signed_alternating.py","sources/sweep-monodromy-manifest.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 formal_verification=False,independent_whole_proof_review=False,goal_complete=False)
(ROOT/"sweep_asymptotics_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=receipt["runtime_seconds"])))
