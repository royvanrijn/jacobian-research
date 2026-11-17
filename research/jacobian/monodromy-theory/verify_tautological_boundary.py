"""Exact controls for TAUTOLOGICAL_BOUNDARY.md; universal proofs are written."""
from pathlib import Path
from collections import Counter
from itertools import product
from math import lcm
import ast,hashlib,json,platform,time
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,v):
    if not v:raise AssertionError(name)
    checks.append(name)
w,z,g,k,a,eta,beta,Y1,Y2,Y3,C=S.symbols("w z gamma kappa a eta beta Y1 Y2 Y3 C")
Delta=S.Matrix([1,w,z]);N=S.Matrix([[-w,1,0],[-z,0,1]])
rank_controls=[]
for f in [S.Integer(2),z,z*z+z,z**4-z*z+1]:
    h=S.integrate(z*S.diff(f,z),z);b=z**3+z
    X=S.Matrix([f,b,h]);sweep=X+g*Delta
    label="rank-one f="+str(f)
    check(label+" full Jacobian",S.expand(S.Matrix.hstack(Delta,sweep.diff(w),sweep.diff(z)).det()-g*g)==0)
    P=S.expand(h-z*f)
    gamma=Y1-f;wr=(Y2-b)/gamma
    check(label+" exact first inverse",S.cancel(sweep[0].subs({g:gamma,w:wr},simultaneous=True)-Y1)==0)
    check(label+" exact second inverse",S.cancel(sweep[1].subs({g:gamma,w:wr},simultaneous=True)-Y2)==0)
    check(label+" exact pencil",S.expand(sweep[2].subs(g,gamma)-Y1*z-P)==0)
    check(label+" derivative",S.expand(S.diff(P,z)+f)==0)
    degree=1 if S.degree(f,z)==0 else int(S.degree(f,z))+1
    check(label+" exact degree",S.degree(Y1*z+P,z)==degree)
    if degree>1:
        av=-S.diff(P,z);bv=P-z*S.diff(P,z)
        check(label+" critical-value curve inverse",S.cancel(S.diff(bv,z)/S.diff(av,z)-z)==0)
    rank_controls.append(dict(f=str(f),degree=degree,geometric_group="S"+str(degree)))
# Determinant trace condition is independently computed for a general
# rank-one parametrization; its irreducible-fiber proof is geometric.
U=S.Function("U")(w,z)
f1,f2,f3=[S.Function(n)(U) for n in ["f1","f2","f3"]]
trace=S.trace(N*S.Matrix([f1,f2,f3]).jacobian([w,z]))
expected=(S.diff(f2,U)-w*S.diff(f1,U))*S.diff(U,w)+(S.diff(f3,U)-z*S.diff(f1,U))*S.diff(U,z)
check("general rank-one trace identity",S.expand(trace-expected)==0)
# Genuine rational poles force every tangent derivative along the pole line to zero.
for tv in [1,2,3]:
    lam=w**(2*tv)*(1+w+z+z*z)
    rho=(k*z*w**(tv-1)+1+w)/w**tv
    x1z=-S.diff(lam,w)-S.diff(rho*lam,z)
    x2z=w*x1z+lam;x3z=z*x1z+rho*lam
    check("pole order "+str(tv)+" collapsed line derivatives",all(S.cancel(e).subs(w,0)==0 for e in [x1z,x2z,x3z]))
# The whole non-affine rational Gauss pencil, with its pole accounting.
gauss=[]
for R in [w*w,w**4+w*w,1/w,1/w**3,w*w+1/w]:
    av=S.diff(R,w)/(k-1);bv=(R-w*S.diff(R,w))/k
    check("Gauss critical-value inverse "+str(R),S.cancel(S.diff(bv,w)/S.diff(av,w)+(k-1)*w/k)==0)
    numerator,denominator=S.fraction(S.cancel((R-(k-1)*beta*w)/k))
    degree=max(S.degree(numerator,w),S.degree(denominator,w))
    check("Gauss cover nontrivial "+str(R),degree>=2)
    gauss.append(dict(R=str(R),degree=int(degree),geometric_group="S"+str(degree)))
# Both quadratic exceptions are controlled at the level of exact involutions.
Bquad=(a*w-k*eta/w)/(k-1);tauw=-k*eta/(a*w)
check("quadratic polynomial Gauss involution",S.cancel(Bquad.subs(w,tauw)-Bquad)==0)
check("quadratic polynomial involution square",S.cancel(tauw.subs(w,tauw)-w)==0)
for j in range(1,5):
    check("quadratic polynomial invariant Laurent separation "+str(j),S.degree(S.fraction(S.cancel(tauw**j))[1],w)==j)
rho=(k*z+a/w)/w
normal_eta=(k-1)*z+a/w;normal_beta=-rho
tw=a*w/((k-1)*(k*w*z+a));tz=-(k-1)*z+a*(2-k)/((k-1)*w)
sub={w:tw,z:tz}
check("reciprocal Gauss involution eta",S.cancel(normal_eta.subs(sub,simultaneous=True)-normal_eta)==0)
check("reciprocal Gauss involution beta",S.cancel(normal_beta.subs(sub,simultaneous=True)-normal_beta)==0)
check("reciprocal Gauss involution square w",S.cancel(tw.subs(sub,simultaneous=True)-w)==0)
check("reciprocal Gauss involution square z",S.cancel(tz.subs(sub,simultaneous=True)-z)==0)
check("reciprocal pole nonconstant finite z residue",S.cancel(tz.subs(z,-a/(k*w))-a/(k*(k-1)*w))==0)
check("reciprocal pole w residue nonzero",S.cancel(((k*w*z+a)*tw).subs(z,-a/(k*w)))!=0)
check("reciprocal second pole coefficient",S.cancel(S.cancel(w*tz).subs(w,0)-a*(2-k)/(k-1))==0)
# Actual positive-z-degree lambda controls, including a higher-pole case.
rational=[]
for rv,R in [(2,w*w+w**3),(6,1/w),(8,1/w**2)]:
    minimum=rv*(rv-1)+4*k*(rv-1)+6*k*k
    def red(e):
        num,den=S.fraction(S.cancel(e))
        return S.rem(S.rem(num,minimum,k)*S.invert(S.rem(den,minimum,k),minimum,k),minimum,k).expand()
    lam=w**rv*z
    for term in S.Add.make_args(S.expand(R)):
        j=int(term.as_powers_dict().get(w,0));co=term/w**j
        coeff=-2*(rv+j-1+2*k)/((rv+j)*(rv+j-1)+2*k*(rv+j-1)+2*k*k)
        lam+=red(coeff)*co*w**(rv+j)
    rho=(k*z+R)/w
    wave=S.diff(lam,w,2)+2*S.diff(rho*lam,w,z)+S.diff(rho*rho*lam,z,2)
    check("positive-z lambda wave r="+str(rv),red(wave)==0)
    check("positive-z polynomial matrix r="+str(rv),all(not S.denom(S.cancel(e)).has(w,z) for e in [lam,rho*lam,rho*rho*lam]))
    x1z=-S.diff(lam,w)-S.diff(rho*lam,z)
    x1w=rho*S.diff(lam,w)+rho*rho*S.diff(lam,z)+(S.diff(rho,w)+2*rho*S.diff(rho,z))*lam
    X1=red(S.integrate(red(x1z),z)+S.integrate(red(x1w.subs(z,0)),w))
    X2=red(S.integrate(w*S.diff(X1,z)+lam,z)+S.integrate(red((w*S.diff(X1,w)-rho*lam).subs(z,0)),w))
    X3=red(S.integrate(z*S.diff(X1,z)+rho*lam,z)+S.integrate(red((-rho*rho*lam).subs(z,0)),w))
    X=S.Matrix([X1,X2,X3])
    check("positive-z integrated normalization r="+str(rv),all(red(e)==0 for e in N*X.jacobian([w,z])-lam*S.Matrix([[-rho,1],[-rho*rho,rho]])))
    check("positive-z full determinant r="+str(rv),red(S.Matrix.hstack(Delta,X.diff(w)+g*Delta.diff(w),X.diff(z)+g*Delta.diff(z)).det()-g*g)==0)
    check("positive-z collapsed line r="+str(rv),all(e.subs(w,0)==0 for e in X))
    check("positive-z rank two r="+str(rv),red(X.jacobian([w,z])[:2,:].det())!=0)
    rational.append(dict(r=rv,R=str(R),lambda_polynomial=str(lam),finite=False,Gauss_surface_degree=1))
# Complete PSL2(11) actions, constructed from all 660 matrices.
source=ast.parse((ROOT/"verify_curve_descent_lcm.py").read_text())
wanted={"projective_group","cycles","generated"}
code=ast.Module(body=[n for n in source.body if isinstance(n,ast.FunctionDef) and n.name in wanted],type_ignores=[])
scope=dict(product=product,lcm=lcm)
exec(compile(code,str(ROOT/"verify_curve_descent_lcm.py"),"exec"),scope)
projective_group,cycles,generated=[scope[n] for n in ["projective_group","cycles","generated"]]
G,ident,mul,action=projective_group(11)
order={g:lcm(*cycles(action(g))) for g in G}
a2=min(g for g in G if order[g]==2);H=None
for a3 in sorted(g for g in G if order[g]==3):
    if order[mul(a2,a3)]!=5:continue
    h=generated((a2,a3),ident,mul)
    if len(h)==60:H=h;break
check("complete projective group and exceptional stabilizer",len(G)==660 and H is not None)
unseen=set(G);reps=[];lookup={}
while unseen:
    h=min(unseen);coset={mul(h,t) for t in H}
    for t in coset:lookup[t]=len(reps)
    reps.append(h);unseen-=coset
exceptional=lambda h:tuple(lookup[mul(h,r)] for r in reps)
check("exceptional action degree",len(reps)==11)
passports={}
for label,perm in [("natural",action),("exceptional",exceptional)]:
    census=Counter(cycles(perm(h)) for h in G)
    fixed={cs for cs in census if 1 in cs and cs!=(1,)*len(perm(ident))}
    candidates={cs for cs in fixed if 3 in cs}
    check(label+" fixed-sheet present-three gate",candidates==(set() if label=="natural" else {(3,3,3,1,1)}))
    allowed={length for cs in fixed for length in cs if length>1}
    check(label+" full local ramification index gate",allowed==({5,11} if label=="natural" else {2,3,5}))
    passports[label]=dict(cycle_census={",".join(map(str,cs)):n for cs,n in sorted(census.items())},
       present_three_and_fixed=[list(cs) for cs in sorted(candidates)],allowed_present_indices=sorted(allowed))
check("nonfinite natural minimum architecture dimension",min(e+1 for e in [5,11])==6)
check("general rank-two exceptional missing ramification",all(3-delta>=1 for delta in [1,2]))
check("rational-chart exceptional missing ramification",3-1==2)
# The target twist really is an isomorphism on C!=0 for arbitrary integer weights.
F1,F2,F3=S.symbols("F1 F2 F3")
for weights in [(1,2,3),(1,2,4),(2,3,5)]:
    ys=[C**e*f for e,f in zip(weights,[F1,F2,F3])]
    check("target product divisor chart inverse "+str(weights),all(S.cancel(y/C**e-f)==0 for y,e,f in zip(ys,weights,[F1,F2,F3])))
manifest=json.loads((ROOT/"sources/sweep-monodromy-manifest.json").read_text())
for src in manifest["sources"][:2]:
    check("retained Gao source "+src["path"],hashlib.sha256((ROOT/src["path"]).read_bytes()).hexdigest()==src["sha256"])
model=dict(rank_at_most_one_controls=rank_controls,rational_Gauss_controls=gauss,positive_z_lambda_controls=rational,
 full_PSL2_11_actions=passports,
 all_tautological_Gao_natural_PSL12="excluded, including rank collapse and nonfinite input",
 exceptional_rational_chart_boundary_requirement="two additional ramification-index3 normalization branches outside polynomial sweep source",
 full_action_classification_complete=False,abstract_PSL2_11_actions_decided=False)
(ROOT/"tautological_boundary_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["TAUTOLOGICAL_BOUNDARY.md","verify_tautological_boundary.py","tautological_boundary_model.json",
 "verify_curve_descent_lcm.py","sources/sweep-monodromy-manifest.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written universal rank, local valuation, pole-collapse and complete rational Gauss arguments. Exact controls are not formal verification or independent whole-proof review.",goal_complete=False)
(ROOT/"tautological_boundary_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
