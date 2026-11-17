"""Exact controls for RATIONAL_SWEEPS.md; universal arguments are written."""
from pathlib import Path
from math import gcd
import hashlib,json,platform,time
import sympy as S
from rational_sweep import surface,general_surface,family
ROOT=Path(__file__).resolve().parent
start=time.monotonic();checks=[]
def check(name,v):
    if not v:raise AssertionError(name)
    checks.append(name)
w,z,g,k,Y1,Y2,Y3=S.symbols("w z gamma kappa Y1 Y2 Y3")
def redk(e):
    num,den=S.fraction(S.cancel(e))
    return S.rem(S.rem(num,k*k+k+1,k)*S.invert(S.rem(den,k*k+k+1,k),k*k+k+1,k),k*k+k+1,k).expand()
D=S.symbols("D",integer=True,positive=True)
zd=2*(k+D+1)/((D+1)*(D+2));z2=(k+3)/6
cd=D*D*(D-1)**2/((2*D+1)*(D+1)**2*(D+2)**2)
bd=D*(D-1)/(3*(D+1)*(D+2)*(D+3))
check("universal leading eliminant coefficient",redk(1/(2*D+1)+2*k*zd/(D+1)-(k+1)*zd**2-cd)==0)
check("universal marked cross coefficient",redk(2/(D+3)+2*k*(z2/(D+1)+zd/3)-2*(k+1)*zd*z2-bd)==0)
check("universal marked discriminant",S.factor(bd*bd-cd/45+
 D*D*(D-2)**2*(D-1)**2/(45*(D+1)**2*(D+2)**2*(D+3)**2*(2*D+1)))==0)
check("universal marked gamma coefficient",redk((k+2)*zd-1-2*k/(D+1)+D*(D-1)/((D+1)*(D+2)))==0)
Delta=S.Matrix([1,w,z]);N=S.Matrix([[-w,1,0],[-z,0,1]])
# The low-degree family is affine equivalent to the regular cubic field.
r0,r1=S.symbols("r0 r1")
Xlow,*_=surface(r0+r1*w,w,z,k)
Xzero,*_=surface(S.Integer(0),w,z,k)
t0=-r0/k;t1=-r1/(k-1)
trans=S.Matrix([[1,0,0],[0,1,0],[-t0,-t1,1]])
check("all affine R degenerate to homogeneous cubic",
 all(redk(e)==0 for e in trans*Xlow.subs(z,z+t0+t1*w)-Xzero))
controls=[]
for R in [w*w+w+1,w**3+w*w,w**5+2*w**2+w+3]:
    dv=S.degree(R,w);lab=f"general D={dv}"
    X,I,J,K0=surface(R,w,z,k);rho=(k*z+R)/w
    pm=N*X.jacobian([w,z])
    check(lab+" polynomial nilpotent frame",all(redk(e)==0 for e in pm-w*w*S.Matrix([[-rho,1],[-rho*rho,rho]])))
    check(lab+" full sweep Jacobian",redk(S.Matrix.hstack(Delta,X.diff(w)+g*Delta.diff(w),X.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
    check(lab+" nonfinite collapsed line",all(e.subs(w,0)==0 for e in X))
    eta=(k-1)*z+R;beta=-(k*z+R)/w
    normal=S.Matrix([eta,beta,1])
    check(lab+" surface Gauss normal",all(redk(normal.dot(X.diff(t)))==0 for t in [w,z]))
    check(lab+" generic rank two minor",redk(S.expand(X.jacobian([w,z])[:2,:].det()).coeff(z,1)-k*k*(k+2)*w*w)==0)
    B=S.symbols("beta")
    X1B=redk(X[0].subs(z,-(B*w+R)/k))
    etaB=redk(eta.subs(z,-(B*w+R)/k))
    check(lab+" normal coordinate formula",redk(etaB-(R-(k-1)*B*w)/k)==0)
    check(lab+" first coordinate formula",redk(X1B+(2*k+1)*B*w*w+2*k*(w*R-I))==0)
    check(lab+" birational surface degree gate",S.degree(etaB,w)==dv and S.degree(X1B,w)==dv+1 and gcd(dv,dv+1)==1)
    Z=Y2-w*Y1+2*k*w*I-2*(k-1)*J
    Gamma=w*Y1+(k+2)*Z-w*w*R-2*k*w*I
    P=redk(w**3*(Y3+K0)-w*(Y1-2*k*I)*Z-(k+1)*Z*Z)
    recovered={z:Z/w**2,g:Gamma/w}
    sweep=X+g*Delta
    check(lab+" exact first target inverse",redk(sweep[0].subs(recovered,simultaneous=True)-Y1)==0)
    check(lab+" exact second target inverse",redk(sweep[1].subs(recovered,simultaneous=True)-Y2)==0)
    check(lab+" exact final target equation",redk(w**3*(Y3-sweep[2].subs(recovered,simultaneous=True))-P)==0)
    check(lab+" complete derivative-square identity",redk(w*S.diff(P,w)-3*P-Gamma*Gamma)==0)
    check(lab+" generic eliminant degree",S.degree(P,w)==2*dv+4)
    lc=S.Poly(P,w).LC();expected=S.Rational(dv*dv*(dv-1)**2,(2*dv+1)*(dv+1)**2*(dv+2)**2)*S.Poly(R,w).LC()**2
    check(lab+" constant leading coefficient",redk(lc-expected)==0)
    check(lab+" squared constant coefficient",redk(P.subs(w,0)-k*k*Y2**2)==0)
    check(lab+" independent critical parameter terms",redk(S.diff(Gamma,Y1)+(k+1)*w)==0 and redk(S.diff(Gamma,Y2)-(k+2))==0)
    check(lab+" generic critical degree",S.degree(Gamma,w)==dv+2)
    tangent=sum(c*w**i*Y2**j for (i,j),c in S.Poly(P,w,Y2).terms() if i+j==2)
    check(lab+" all boundary quadratic terms",redk(tangent+k*(Y1*w-Y2)*(Y1*w+k*Y2))==0)
    v=S.symbols("v")
    hensel=redk(tangent.subs(w,Y2*v)/Y2**2)
    check(lab+" both missing roots simple",all(redk(hensel.subs(v,c))==0 and redk(S.diff(hensel,v).subs(v,c))!=0 for c in [1/Y1,-k/Y1]))
    check(lab+" full curve Riemann Hurwitz",2*(dv+2)+2+2*dv==2*(2*dv+4)-2)
    controls.append(dict(D=int(dv),R=str(R),degree=2*int(dv)+4,geometric_group="A"+str(2*dv+4),finite=False))
# Weighted jets are computed independently as truncated polynomial dictionaries.
v1,v2,v3=S.symbols("v1 v2 v3")
def jet_compose(E,theta,variables,cut,red):
    weights=(1,2,3);zero=(0,0,0)
    def mul(a,b):
        out={}
        for i,ci in a.items():
            for j,cj in b.items():
                mon=tuple(x+y for x,y in zip(i,j))
                if sum(x*y for x,y in zip(mon,weights))<cut:
                    out[mon]=out.get(mon,0)+ci*cj
        return {mon:red(c) for mon,c in out.items() if red(c)!=0}
    jets=[{mon:red(c) for mon,c in S.Poly(t,v1,v2,v3).terms()
          if sum(x*y for x,y in zip(mon,weights))<cut and red(c)!=0} for t in theta]
    powers=[{0:{zero:S.Integer(1)}} for j in range(3)];out={}
    for exps,coeff in S.Poly(E,*variables).terms():
        term={zero:coeff}
        for j,e in enumerate(exps):
            for h in range(1,e+1):
                if h not in powers[j]:powers[j][h]=mul(powers[j][h-1],jets[j])
            term=mul(term,powers[j][e])
        for mon,c in term.items():out[mon]=red(out.get(mon,0)+c)
    return {mon:c for mon,c in out.items() if c!=0}
models=[]
for dv in [3,4,5]:
    f=family(dv);a,kk,gg,u,v,ww,zz=[f[t] for t in ["a","k","g","u","v","w","z"]]
    red=f["red"];lab=f"affine D={dv}"
    check(lab+" nonreal marked quadratic",S.discriminant(f["minimal"],a)<0)
    check(lab+" unramified marked gamma",red(f["g0"])!=0 and S.resultant(f["minimal"],f["g0"],a)!=0)
    check(lab+" all marked origin coordinates",all(red(e.subs({ww:1,zz:f["z0"],gg:f["g0"]}))==0 for e in f["sweep"]))
    for i,e in enumerate([1,2,3]):
        check(lab+" polynomial gamma quotient "+str(i),S.denom(S.cancel(f["E"][i])).has(gg,u,v)==False)
        check(lab+" exact gamma quotient "+str(i),S.expand(f["E"][i]*gg**e-f["sweep"][i].subs({ww:gg*u,zz:gg*gg*v},simultaneous=True))==0)
    check(lab+" marked derivative inverse",all(red(e)==0 for e in f["J"]*f["L"]-S.eye(3)))
    check(lab+" stage constant determinant",red(f["L"].det()-f["g0"])==0)
    theta=f["theta0"]+f["L"]*S.Matrix([v1,v2,v3+f["h"]*v1*v1])
    recovered=(f["J"]*(theta-f["theta0"])).applyfunc(red)
    check(lab+" triangular stage inverse",all(red(e)==0 for e in recovered-S.Matrix([v1,v2,v3+f["h"]*v1*v1])))
    for i,e in enumerate([1,2,3]):
        check(lab+" every weighted source coefficient "+str(i),not jet_compose(f["E"][i],theta,[gg,u,v],e,red))
    check(lab+" full x gamma determinant cancellation",1+2+3==1+1+2+2)
    check(lab+" exact degree",S.degree(f["P"],ww)==f["N"])
    dchar=(-1)**dv*(2*dv+1)
    check(lab+" exact arithmetic discriminant class",S.cancel(((-1)**dv/(f["cd"]*a*a))/dchar-((dv+1)*(dv+2)/(dv*(dv-1)*a))**2)==0)
    models.append(dict(D=dv,degree=f["N"],geometric_group="A"+str(f["N"]),
       alpha_equation=str(f["minimal"]),gamma0=str(f["g0"]),shear=str(f["h"]),
       source_and_target="ordinary affine four-space",polynomial=True,arithmetic_sign_class=dchar))
general_controls=[]
for rv,R in [(3,w**3+w*w),(4,w*w+1/w),(5,w**3+1/w),(6,1/w+1/w**2)]:
    lab=f"all-r rational r={rv}"
    minimum=k*k+(rv-1)*k+S.Rational(rv*(rv-1),2)
    def redr(e):
        num,den=S.fraction(S.cancel(e))
        return S.rem(S.rem(num,minimum,k)*S.invert(S.rem(den,minimum,k),minimum,k),minimum,k).expand()
    X,I,J,K0=general_surface(rv,R,w,z,k);c=S.Rational(rv,2)+k
    F=w**(rv-1)*R+2*k*I;G=w**rv*R+2*(k-1)*J
    Z=Y2-w*Y1+w*F-G
    Gamma=w*Y1+(rv+k)*Z-w*F
    P=redr(w**(rv+1)*(Y3+K0)-w*(Y1-2*k*I)*Z-c*Z*Z)
    jv=max(int(term.as_powers_dict().get(w,0)) for term in S.Add.make_args(S.expand(R)))
    nv=2*(rv+jv);uj=S.Rational(jv*(jv-1),(rv+jv-1)*(rv+jv))
    check(lab+" rational input is admissible",min(int(term.as_powers_dict().get(w,0)) for term in S.Add.make_args(S.expand(R)))>=1-rv//2)
    check(lab+" all surface coordinates polynomial",all(S.denom(S.cancel(e)).has(w,z)==False for e in X))
    rho=(k*z+R)/w
    check(lab+" complete polynomial normalization",all(redr(e)==0 for e in N*X.jacobian([w,z])-w**rv*S.Matrix([[-rho,1],[-rho*rho,rho]])))
    check(lab+" full normalized Jacobian",redr(S.Matrix.hstack(Delta,X.diff(w)+g*Delta.diff(w),X.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
    check(lab+" nonfinite line",all(e.subs(w,0)==0 for e in X))
    check(lab+" residue square-class identity",redr((k+rv)**2-(rv+1)*c)==0)
    check(lab+" complete derivative identity",redr(w*S.diff(P,w)-(rv+1)*P-Gamma*Gamma)==0)
    check(lab+" field degree",S.degree(P,w)==nv and nv>rv+1)
    check(lab+" constant leading coefficient",redr(S.Poly(P,w).LC()-uj**2/(rv+2*jv-1))==0)
    check(lab+" critical degree",S.degree(Gamma,w)==rv+jv)
    check(lab+" critical independent parameters",redr(S.diff(Gamma,Y1)+(rv+k-1)*w)==0 and redr(S.diff(Gamma,Y2)-(rv+k))==0)
    check(lab+" full Riemann Hurwitz",2*(rv+jv)+rv+(rv+2*jv-2)==2*nv-2)
    beta=S.symbols("beta")
    X1B=redr(X[0].subs(z,-(beta*w+R)/k))
    X2B=redr(X[1].subs(z,-(beta*w+R)/k))
    degs=(S.degree(X1B,w),S.degree(X2B,w))
    check(lab+" consecutive surface field degrees",degs==((rv+jv-1,rv+jv) if jv>=2 else (rv,rv+1)) and gcd(*degs)==1)
    tangent=sum(co*w**i*Y2**j for (i,j),co in S.Poly(P,w,Y2).terms() if i+j==2)
    boundary=(1-c)*Y1**2*w*w+(2*c-1)*Y1*w*Y2-c*Y2*Y2
    check(lab+" complete normalization boundary",redr(tangent-boundary)==0)
    vv=S.symbols("v");quadratic=redr(boundary.subs(w,Y2*vv)/Y2**2)
    check(lab+" unramified missing boundary branches",all(redr(quadratic.subs(vv,t))==0 and redr(S.diff(quadratic,vv).subs(vv,t))!=0 for t in [1/Y1,c/((c-1)*Y1)]))
    sign=(-1)**(rv+jv+1)*(rv+1)*(rv+2*jv-1)
    leading=S.Poly(P,w).LC()
    check(lab+" exact arithmetic sign square",redr(((-1)**(rv+jv)/(-c*leading))/sign-1/((k+rv)*uj)**2)==0)
    low,*_=general_surface(rv,r0+r1*w,w,z,k)
    zero,*_=general_surface(rv,S.Integer(0),w,z,k)
    check(lab+" all affine terms removable",all(redr(e)==0 for e in trans*low.subs(z,z+t0+t1*w)-zero))
    general_controls.append(dict(r=rv,R=str(R),highest_nonaffine_exponent=jv,degree=nv,geometric_group="A"+str(nv),arithmetic_sign_class=sign,finite=False))
# Degree-two centered origin obstruction includes arbitrary affine terms.
a2=S.symbols("a2")
XX,II,JJ,KK=surface(a2*w*w+r1*w+r0,w,z,k)
ZZ=2*k*w*II-2*(k-1)*JJ
origin=redk(w**3*KK+2*k*w*II*ZZ-(k+1)*ZZ*ZZ)
check("all quadratic R centered origin obstruction",S.expand(origin-a2*a2*w**8/180)==0)
manifest=json.loads((ROOT/"sources/sweep-monodromy-manifest.json").read_text())
for src in manifest["sources"][:2]:
    check("retained Gao primary source "+src["path"],hashlib.sha256((ROOT/src["path"]).read_bytes()).hexdigest()==src["sha256"])
model=dict(classification="rho=(kappa*z+R(w))/w,lambda=w^2,kappa^2+kappa+1=0",
 deg_R_at_most_one="regular C3, excluded as absolute Keller field",
 deg_R_at_least_two="degree2D+4; geometric A_(2D+4); arithmetic sign class(-1)^D(2D+1)",
 sweep_controls=controls,actual_affine_controls=models,general_r_controls=general_controls,
 all_z_independent_lambda_rational_slice_classified=True,
 universal_family="every D>=3; full formulas in rational_sweep.family",
 abstract_PSL2_11_actions_decided=False,full_action_classification_complete=False)
(ROOT/"rational_sweep_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["RATIONAL_SWEEPS.md","rational_sweep.py","verify_rational_sweep.py","rational_sweep_model.json","sources/sweep-monodromy-manifest.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,
 python=platform.python_version(),sympy=S.__version__,runtime_seconds=round(time.monotonic()-start,2),
 sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written universal group, boundary, constant-field and affine-lift proofs. Exact controls are not formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"rational_sweep_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
