"""Exact controls for the all-parameter proof in RATIONAL_CHART_CLASSIFICATION.md."""
from pathlib import Path
from math import gcd
import hashlib,json,platform,time
import sympy as S
from rational_wave import rational_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(n,v):
    if not v:raise AssertionError(n)
    checks.append(n)
w,z,k,t,g,Y1,Y2,Y3=S.symbols("w z kappa t gamma Y1 Y2 Y3")
rr,ss,jj,nn=S.symbols("r s j n",integer=True)
def C(a,b):return a*(a-1)+2*k*(a-1)*(b+1)+k*k*(b+1)*(b+2)
AA=1+(k-jj)*t
BB=2*(nn-1+2*k)+(-2*jj*nn+jj*jj+jj+2*k*nn-4*k*jj-2*k+4*k*k)*t
EE=nn*AA+(jj-1)+2*k+k*(2*k-jj-1)*t
t0=-1/(k-jj)
check("universal leading-wave regular coefficient",S.cancel(BB.subs(t,t0)+jj*(jj-1)/(k-jj))==0)
check("universal surface Euler coefficient",S.cancel(EE.subs(t,t0)+jj*(jj-1)/(k-jj))==0)
leading=(k-jj)**2*ss*(ss-1)+S.diff(BB,t)*ss+C(nn,0)
check("universal top residue equation",S.expand(leading.subs(nn,rr+ss*jj)-C(rr,ss))==0)
M=rr+(ss+1)*jj;H=rr+(ss+2)*jj-1
raw=(ss+2)*M+(ss+1)*H-jj*(ss+1)*(ss+2)
E=(ss+1)*(rr-1)
check("universal resultant degree",S.expand(raw-E-(ss+2)*M)==0)
check("universal resultant zero valuation",S.expand(rr*(ss+2)+(rr-1)*(ss+1)-rr*(ss+2)-E)==0)
b=(rr+2*k*(ss+1))/((ss+1)*(ss+2));beta0=1/(ss+1)
at,TT=S.symbols("a T")
J=S.Matrix([at*(Y1+beta0*TT)-1,Y1+b*TT]).jacobian([at,TT])
check("universal boundary Hensel determinant",S.cancel(J.det().subs(TT,-Y1/b)-Y1*(b-beta0))==0)
# Complete parity argument is written; these controls cover its four parity cases.
for rv,sv in [(2,0),(2,1),(3,0),(3,1),(4,2),(5,3),(10,5),(11,6)]:
    gv=gcd(sv+1,rv-1)
    check(f"boundary parity r={rv},s={sv}",((sv+1-gv)%2==1)==(rv%2==0 and sv%2==1))
Delta=S.Matrix([1,w,z]);Nframe=S.Matrix([[-w,1,0],[-z,0,1]])
cases=[(2,0,w*w+w**3),(2,1,w*w),(3,1,w*w+w**3),(2,2,w*w),
       (3,2,w**3+w*w),(6,1,1/w),(8,1,1/w**2),(10,2,1/w**2+1/w),
       (7,1,1/w+w*w),(2,3,w*w)]
models=[]
for rv,sv,R in cases:
    f=rational_chart(rv,sv,R,w,z,k);red=f["red"];lam=f["lambda_polynomial"];rho=f["rho"];X=f["X"]
    label=f"r={rv},s={sv},R={R}"
    check(label+" complete wave",red(S.diff(lam,w,2)+2*S.diff(rho*lam,w,z)+S.diff(rho*rho*lam,z,2))==0)
    check(label+" top coefficient",red(S.Poly(lam,z).LC()-w**rv)==0)
    check(label+" polynomial frame",all(not S.denom(S.cancel(e)).has(w,z) for e in [lam,rho*lam,rho*rho*lam]))
    check(label+" exact integrated frame",all(red(e)==0 for e in Nframe*X.jacobian([w,z])-lam*S.Matrix([[-rho,1],[-rho*rho,rho]])))
    check(label+" full Jacobian",red(S.Matrix.hstack(Delta,X.diff(w)+g*Delta.diff(w),X.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
    check(label+" centered collapsed line",all(e.subs(w,0)==0 for e in X))
    check(label+" exact first fiber equation",red((X+g*Delta)[1].subs(g,Y1-X[0])-f["V"]-w*Y1)==0)
    check(label+" exact second fiber equation",red((X+g*Delta)[2].subs(g,Y1-X[0])-f["W"]-z*Y1)==0)
    scalar=red((rv+2*k*(sv+1))/((sv+1)*(sv+2)))
    check(label+" first leading z coefficient",red(S.Poly(f["V"],z).LC()-w**rv/(sv+1))==0)
    check(label+" second leading z coefficient",red(S.Poly(f["W"],z).LC()-scalar*w**(rv-1))==0)
    jv=f["j"];mv=rv+(sv+1)*jv;hv=rv+(sv+2)*jv-1
    vp=sum(co*w**i*z**j for (i,j),co in S.Poly(f["V"],w,z).terms() if i+jv*j==mv)
    wp=sum(co*w**i*z**j for (i,j),co in S.Poly(f["W"],w,z).terms() if i+jv*j==hv)
    va=red(vp.subs({w:1,z:t}));wb=red(wp.subs({w:1,z:t}))
    leading_resultant=red(S.resultant(va,wb,t))
    check(label+" coprime weighted leading equations",leading_resultant!=0)
    # The resultant is calculated at one ordinary target; universal degree and
    # coefficient arguments are in the proof, rather than inferred from this point.
    res=red(S.resultant(f["V"]+2*w-3,f["W"]+2*z-5,z))
    extraneous=(sv+1)*(rv-1);poly=S.Poly(res,w)
    valuation=min(mon[0] for mon in poly.monoms())
    check(label+" exact extraneous w valuation",valuation==extraneous)
    Q=S.cancel(res/w**extraneous)
    check(label+" exact generic-degree control",S.degree(Q,w)==f["N"])
    check(label+" leading coefficient equals weighted resultant",red(S.Poly(Q,w).LC()-leading_resultant)==0)
    check(label+" exact nonzero constant coefficient",red(Q.subs(w,0)-(-1)**(sv+2)*scalar**(sv+1)*3**(sv+2))==0)
    check(label+" composite degree",f["N"]==(sv+2)*mv and sv+2>1 and mv>1)
    # Full Newton initial forms at the omitted normalization branches.
    uu,a0,z0=S.symbols("u a0 z0")
    gg=gcd(sv+1,rv-1);ee=(sv+1)//gg;hh=(rv-1)//gg
    sub={w:uu**ee*a0,z:uu**(-hh)*z0}
    eq1=S.expand((f["V"]+w*Y1).subs(sub,simultaneous=True)-uu**ee)
    eq2=S.expand((f["W"]+z*Y1-Y3).subs(sub,simultaneous=True))
    def coefficient(expr,power):
        return red(sum(term/uu**power for term in S.Add.make_args(expr)
                    if term.as_powers_dict().get(uu,0)==power))
    init1=coefficient(eq1,ee);init2=coefficient(eq2,-hh)
    T=a0**(rv-1)*z0**(sv+1)
    check(label+" first complete Newton initial equation",red(init1-a0*(Y1+T/(sv+1))+1)==0)
    check(label+" second complete Newton initial equation",red(init2-z0*(Y1+scalar*T))==0)
    check(label+" boundary sheet count and fixed sheet",ee*gg==sv+1 and f["N"]-(sv+1)>=1)
    model=rational_chart(rv,sv,w**jv,w,z,k)
    LL=red(model["lambda_polynomial"].subs({w:1,z:t}))
    check(label+" leading lambda avoids singular point",red(LL.subs(t,-1/(k-jv)))!=0)
    if sv:
        check(label+" leading lambda has no multiple zero",red(S.discriminant(LL,t))!=0)
    check(label+" full group parity control",(f["geometric_group"].startswith("S"))==((sv+1-gg)%2==1))
    models.append(dict(r=rv,s=sv,R=str(R),degree=f["N"],geometric_group=f["geometric_group"],
       boundary_index=ee,boundary_branches=gg,boundary_fixed_sheets=f["N"]-(sv+1),
       resultant_degree=int(S.degree(res,w)),extraneous_w_power=extraneous,
       lambda_polynomial=str(lam)))
# Generic target coefficients are also retained in one positive-z example.
f=rational_chart(2,1,w*w,w,z,k);red=f["red"]
res=red(S.resultant(f["V"]+w*Y1-Y2,f["W"]+z*Y1-Y3,z))
Q=S.cancel(res/w**2);scalar=(2+4*k)/6
check("generic-target resultant degree",S.degree(Q,w)==18)
check("generic-target exact plane constant",red(Q.subs(w,0)+scalar**2*Y2**3)==0)
# The sharper pole bound rejects cases allowed by the old top-only denominator bound.
for rv,sv,exponent in [(2,1,-1),(4,1,-1),(8,2,-2)]:
    try:rational_chart(rv,sv,w**exponent,w,z,k)
    except ValueError as e:check(f"exact pole-bound rejection {rv},{sv},{exponent}","pole" in str(e))
    else:raise AssertionError("inadmissible pole accepted")
for rv,sv in [(2,0),(3,0),(2,1),(2,2)]:
    f=rational_chart(rv,sv,S.Integer(0),w,z,k)
    check(f"homogeneous limit {rv},{sv}",f["lambda_polynomial"]==w**rv*z**sv and f["N"]==rv+sv+1)
# No classification claim is based on the finite examples above.
model=dict(classification="all z-dependent rho in k(w)[z] with polynomial nilpotent normalization matrix",
 exact_integrability_bound="r+(s+2)*m>=2 for a negative least Laurent exponent m",
 homogeneous_degree="r+s+1",nonhomogeneous_degree="(s+2)*(r+(s+1)*j), composite",
 nonhomogeneous_group="S_N iff r even and s odd; otherwise A_N",
 affine_boundary="only Y2=0 can contain omitted normalization branches",
 principal_branch="birational image X; full inertia one isolated3-cycle",
 controls=models,generic_target_control=dict(r=2,s=1,R="w^2",degree=18),
 specified_PSL2_11_field_pairs_in_this_chart="excluded for both actions",
 general_positive_z_arithmetic_constants_complete=False,
 abstract_PSL2_11_actions_decided=False,full_action_classification_complete=False)
(ROOT/"rational_chart_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["RATIONAL_CHART_CLASSIFICATION.md","rational_wave.py","verify_rational_chart.py","rational_chart_model.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written universal recursion, weighted resultant, finite-normalization and full-group proofs. Exact controls are not formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"rational_chart_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
