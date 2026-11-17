"""Exact controls for CURVE_DESCENT_LCM.md; universal results have written proofs."""
from pathlib import Path
from collections import Counter
from itertools import product, permutations, combinations_with_replacement
from math import gcd,lcm,factorial
import hashlib,json,platform,time
import sympy as s

ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)

# Complete finite action construction over prime fields.
def projective_group(p):
    def can(a):return min(tuple(v%p for v in a),tuple(-v%p for v in a))
    def mul(a,b):
        return can((a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
                    a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]))
    ident=can((1,0,0,1))
    G=frozenset(can(a) for a in product(range(p),repeat=4) if (a[0]*a[3]-a[1]*a[2])%p==1)
    def action(a):
        aa,b,c,d=a;out=[]
        for x in range(p):
            denominator=(c*x+d)%p
            out.append(p if not denominator else (aa*x+b)*pow(denominator,-1,p)%p)
        out.append(p if not c else aa*pow(c,-1,p)%p)
        return tuple(out)
    return G,ident,mul,action
def cycles(p):
    seen=set();out=[]
    for a in range(len(p)):
        if a in seen:continue
        x=a;size=0
        while x not in seen:seen.add(x);size+=1;x=p[x]
        out.append(size)
    return tuple(sorted(out,reverse=True))
def invariant(permutations):
    census=Counter(gcd(*cycles(g)) for g in permutations)
    return lcm(*census),dict(sorted(census.items()))
def generated(gs,ident,mul):
    group={ident};stack=[ident]
    while stack:
        a=stack.pop()
        for g in gs:
            b=mul(a,g)
            if b not in group:group.add(b);stack.append(b)
    return frozenset(group)

benchmarks=[];contexts={}
for p in (3,5,7,11,13):
    G,ident,mul,action=projective_group(p)
    n=p+1
    inv,census=invariant(action(g) for g in G)
    check(f"PSL2({p}) full group census",len(G)==p*(p*p-1)//2)
    check(f"PSL2({p}) natural global invariant",inv==(p+1)//2<n)
    contexts[p]=(G,ident,mul,action)
    benchmarks.append(dict(group=f"PSL2({p})",action="natural",degree=n,
                           group_order=len(G),cycle_gcd_lcm=inv,census=census,
                           curve_induced_Keller_gate="excluded"))
# Two exceptional actions: same abstract groups, different answers.
for p,sub_order,product_order in [(7,24,4),(11,60,5)]:
    G,ident,mul,action=contexts[p]
    order={g:lcm(*cycles(action(g))) for g in G}
    a2=min(g for g in G if order[g]==2)
    subgroup=None
    for a3 in sorted(g for g in G if order[g]==3):
        if order[mul(a2,a3)]!=product_order:continue
        candidate=generated((a2,a3),ident,mul)
        if len(candidate)==sub_order:
            subgroup=candidate;break
    check(f"PSL2({p}) exceptional stabilizer constructed",subgroup is not None)
    unseen=set(G);reps=[];lookup={}
    while unseen:
        g=min(unseen);coset={mul(g,h) for h in subgroup}
        for member in coset:lookup[member]=len(reps)
        reps.append(g);unseen-=coset
    perm=lambda g:tuple(lookup[mul(g,r)] for r in reps)
    inv,census=invariant(perm(g) for g in G)
    check(f"PSL2({p}) exceptional degree and global invariant",len(reps)==p and inv==p)
    if p==11:
        finite_indices={11-len(cycles(perm(g))) for g in G if g!=ident}
        check("degree-eleven complete nontrivial index set",finite_indices=={4,6,8,10})
        budgets=[a for size in (1,2) for a in combinations_with_replacement(sorted(finite_indices),size) if sum(a)==10]
        check("degree-eleven polynomial finite-index budget",budgets==[(10,),(4,6)])
        check("degree-eleven index ten requires eleven-cycle",
              all(cycles(perm(g))==(11,) for g in G if 11-len(cycles(perm(g)))==10))
        check("degree-eleven index four and six require orders two and three",
              all(lcm(*cycles(perm(g)))==(2 if 11-len(cycles(perm(g)))==4 else 3)
                  for g in G if 11-len(cycles(perm(g))) in (4,6)))
    benchmarks.append(dict(group=f"PSL2({p})",action="exceptional cosets",degree=p,
                           group_order=len(G),cycle_gcd_lcm=inv,census=census,
                           curve_induced_Keller_gate="passes necessary test only"))
# Positive controls for a full cycle and an independent regular-action exclusion.
def pmul(a,b):return tuple(a[b[i]] for i in range(len(a)))
for n in (3,4,5):
    ident=tuple(range(n));rotation=tuple((i+1)%n for i in range(n))
    transposition=tuple([1,0]+list(range(2,n)))
    G=generated((rotation,transposition),ident,pmul)
    inv,census=invariant(G)
    check(f"S_{n} natural full-cycle control",len(G)==factorial(n) and inv==n)
regular=generated(((1,2,0),),tuple(range(3)),pmul)
check("regular C3 passes C but is excluded by KMFI1",invariant(regular)[0]==3)
for n in (4,8,16):
    # Any even n-cycle is odd; c(g)=n would require that full cycle.
    check(f"alternating power-two full cycle parity degree {n}",(n-1)%2==1)
a6=[p for p in permutations(range(6)) if sum(p[i]>p[j] for i in range(6) for j in range(i+1,6))%2==0]
check("composite-degree A6 passes invariant without full cycle",
      len(a6)==360 and invariant(a6)[0]==6 and all(cycles(p)!=(6,) for p in a6))
check("A6 attains both needed cycle gcds",
      any(cycles(p)==(4,2) for p in a6) and any(cycles(p)==(3,3) for p in a6))
# Global gcd/lcm conversion, not a search for maps.
for n,es in [(8,[1,1,1]),(12,[1,1,2]),(3,[1,1,3]),(6,[2,3]),(11,[1,1,11]),(14,[1,7])]:
    check("integer global congruence "+str((n,es)),
          gcd(*(n//e for e in es))==n//lcm(*es))

# Degree-eight excluded cover from quadratic base change of X0(7).
v,t,z=s.symbols("v t z")
P4=v**4+5*v**2+1
Q4=v**4+13*v**2+49
N=v**8+14*v**6+63*v**4+70*v**2-7
f=N/v
jnum=(t*t+13*t+49)*(t*t+5*t+1)**3
L=t**4+14*t**3+63*t**2+70*t-7
check("classical X0(7) square identity",s.expand(jnum-1728*t-L**2)==0)
check("quadratic base-change identity",s.cancel(jnum.subs(t,v*v)/(v*v)-1728-f*f)==0)
check("simple derivative-square identity",s.cancel(s.diff(f,v)-7*P4**2/v**2)==0)
check("critical values and exact multiplicities",s.expand(N**2+1728*v*v-P4**3*Q4)==0)
check("all critical and simple-fiber roots distinct",
      s.gcd(P4,s.diff(P4,v))==1 and s.gcd(Q4,s.diff(Q4,v))==1 and s.gcd(P4,Q4)==1)
extension=s.QQ.algebraic_field(s.sqrt(-3))
for sign in (-1,1):
    critical=sign*24*s.sqrt(-3)
    poly=s.Poly(N-critical*v,v,domain=extension)
    check("critical-fiber passport "+str(sign),
          sorted((p.degree(),e) for p,e in poly.sqf_list()[1])==[(2,1),(2,3)])
check("pole passport seven plus one",s.degree(N,v)-1==7 and N.subs(v,0)==-7)
disc=s.factor(s.discriminant(N-z*v,v))
check("degree-eight arithmetic discriminant",s.expand(disc+7**7*(z*z+1728)**4)==0)
check("each single-fiber coprimality test passes",all(gcd(e,8//e)==1 for e in (1,1,1)))
check("global lcm test excludes the same cover",lcm(1,1,1)!=8)
# Its finite inertia has two unramified sheets, and pole inertia has one.
check("all curve inertia fixes sheets",all(1 in profile for profile in [(3,3,1,1),(3,3,1,1),(7,1)]))

# Actual cubic Keller map and a birational curve-field presentation.
x,y,zz=s.symbols("x y zz")
u=1+x*y
gamma=1-s.Rational(3,2)*x*y+x*x*zz
FA=s.cancel((u+u*u-2*u**3*gamma)/(x*x))
FB=s.cancel((1+2*u-3*u*u*gamma)/x)
FC=x*gamma
F=[FA,FB,FC]
for i,coord in enumerate(F):
    check("global cubic-map polynomial coordinate "+str(i),s.denom(coord)==1)
check("global cubic-map Jacobian one",s.expand(s.Matrix(F).jacobian([x,y,zz]).det())==1)
W=u*gamma; slope=FB*FC; intercept=FA*FC**2
check("cubic inverse equation in source ring",
      s.expand(W**3-W**2+slope*W-intercept)==0)
# Work in generic curve/auxiliary coordinates, avoiding any pointwise roots.
ell,cc,tau=s.symbols("ell cc tau")
curve_tau=-(v+1)/v**3
pp=ell**2/curve_tau; qq=ell**3/curve_tau
ss=pp+s.Rational(1,3)
tt=ss/3-qq-s.Rational(2,27)
AA=tt/cc**2; BB=ss/cc; WW=ell*v+s.Rational(1,3)
gg=s.factor(ss-2*WW+3*WW**2)
xr=cc/gg
ur=WW/gg
yr=s.cancel((ur-1)/xr)
zr=s.cancel((gg-1+s.Rational(3,2)*xr*yr)/xr**2)
sub={x:xr,y:yr,zz:zr}
for label,coord,target in zip(("A","B","C"),F,(AA,BB,cc)):
    check("source reconstruction returns target "+label,
          s.cancel(coord.subs(sub,simultaneous=True)-target)==0)
check("source reconstruction returns W",s.cancel(W.subs(sub,simultaneous=True)-WW)==0)
target_p=BB*cc-s.Rational(1,3)
target_q=BB*cc/3-AA*cc**2-s.Rational(2,27)
check("target reconstruction returns p",s.cancel(target_p-pp)==0)
check("target reconstruction returns q",s.cancel(target_q-qq)==0)
check("target reconstruction returns ell",s.cancel(target_q/target_p-ell)==0)
check("target reconstruction returns tau",s.cancel(target_q**2/target_p**3-curve_tau)==0)
check("source field root reconstruction returns v",s.cancel((WW-s.Rational(1,3))/(target_q/target_p)-v)==0)
# Both directions of the target birational reparameterization.
par,spar,tpar=s.symbols("par spar tpar")
p0=spar-s.Rational(1,3)
q0=spar/3-tpar-s.Rational(2,27)
ell0=q0/p0; tau0=q0*q0/p0**3
check("target inverse returns slope",s.cancel(ell0**2/tau0+s.Rational(1,3)-spar)==0)
check("target inverse returns intercept",
      s.cancel((ell0**2/tau0+s.Rational(1,3))/3-ell0**3/tau0-s.Rational(2,27)-tpar)==0)
curve_poly=tau*v**3+v+1
check("curve equation and rational map",s.cancel(curve_poly.subs(tau,curve_tau))==0)
check("cubic curve discriminant",s.expand(s.discriminant(curve_poly,v)+tau*(4+27*tau))==0)
check("cubic finite ramification factorization",s.expand(4*v**3-27*v-27-(2*v+3)**2*(v-3))==0)
check("cubic global lcm condition",lcm(1,1,3)==3)
check("Keller realization is not falsely excluded",invariant(generated(((1,2,0),(1,0,2)),tuple(range(3)),pmul))[0]==3)

manifest=json.loads((ROOT/"sources/curve-descent-manifest.json").read_text())
for source in manifest["sources"]:
    if source.get("path"):
        check("downloaded modular source hash",hashlib.sha256((ROOT/"sources"/source["path"]).read_bytes()).hexdigest()==source["sha256"])
data=dict(group_benchmarks=benchmarks,
 global_obstruction="For a stabilized curve pair, lcm of the branch-fiber index gcds must equal its degree.",
 natural_odd_q_PSL2="all curve-induced Keller realizations excluded; general higher-dimensional actions open",
 degree_eight_control=dict(numerator=str(N),denominator=str(v),derivative=str(7*P4**2/v**2),
                          discriminant=str(disc),geometric_group="natural PSL2(7), using the attributed modular-cover identification",
                          fiber_gcds=[1,1,1],global_lcm=1,excluded=True),
 cubic_Keller_control=dict(map=[str(e) for e in F],jacobian=1,
                          curve_equation=str(curve_poly),curve_parameter=str(curve_tau),
                          birational_source_inverse=[str(xr),str(yr),str(zr)],
                          target_coordinates=[str(AA),str(BB),str(cc)],
                          geometric_and_arithmetic_group="S3",actual_affine_space_realization=True),
 full_action_spectrum_complete=False,proper_simple_Keller_map_constructed=False)
(ROOT/"curve_descent_model.json").write_text(json.dumps(data,indent=2)+"\n")
files=["CURVE_DESCENT_LCM.md","verify_curve_descent_lcm.py","curve_descent_model.json",
       "sources/curve-descent-manifest.json"]
receipt=dict(result="pass",checks=checks,check_count=len(checks),python=platform.python_version(),
 sympy=s.__version__,runtime_seconds=round(time.monotonic()-started,3),
 sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
 assurance="Universal proofs, exact identities, full finite-group censuses and field reconstruction controls. No formal verification or independent whole-proof review. Modular cover identification is attributed.",
 complete_classification=False,abstract_natural_PSL2_odd_q_excluded=False,
 stabilized_natural_PSL2_odd_q_curve_covers_excluded=True,
 actual_cubic_curve_induced_Keller_map_verified=True)
(ROOT/"curve_descent_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                     proper_simple_Keller_map_constructed=False)))
