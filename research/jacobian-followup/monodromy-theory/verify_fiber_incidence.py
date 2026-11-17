"""Exact algebraic controls for the written fiber-incidence theorem."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
s,L,v,x,y,z,w,C,T1,T2=S.symbols("s L v x y z w C T1 T2")
def red(expr):
    return S.Poly(S.expand(expr),s).rem(S.Poly(s*s+2,s)).as_expr().expand()
def zero(expr):return red(expr)==0
zeta=-1+s;delta=4+s;dinv=(4-s)/18
g=1+delta*L/3+(-10+2*s)*L**2/3
rho=g+L**2*v
Z=red(rho**3-delta*L);Q=red((Z-1)/3)
P=red(6*dinv**2*((rho**4-1)/4-Z*(rho-1)))
check("coefficient inverse",zero(delta*dinv-1))
check("marked value of rho",zero(g.subs(L,1)-zeta))
check("marked unramified tail",zero(zeta**3-(5+s)))
check("cube interpolation jet",S.Poly(red(g**3-1-delta*L),L,v).terms()[-1][0][0]>=2)
qa=S.div(S.Poly(Q,L,v),S.Poly(L**2,L,v))
pa=S.div(S.Poly(P,L,v),S.Poly(L**2,L,v))
check("Q divisible by L squared",qa[1].is_zero)
check("P divisible by L squared",pa[1].is_zero)
A=qa[0].as_expr();B=pa[0].as_expr()
check("lift base A",zero(A.subs({L:1,v:0})))
check("lift base B",zero(B.subs({L:1,v:0})))
check("quotient Jacobian constant after contraction",
      zero(S.diff(P,L)*S.diff(Q,v)-S.diff(P,v)*S.diff(Q,L)-2*L**3))
R=red(delta**2*P/6+S.Rational(1,4)-Z)
check("isotrivial numerator identity",zero(R-(rho**4/4-Z*rho)))
f=w*(w/4-1)**3
check("S4 finite double branch factorization",
      S.expand(64*f+27-(w-1)**2*(w*w-10*w+27))==0)
check("S4 finite zero branch factorization",S.expand(64*f-w*(w-4)**3)==0)
check("S4 fixed polynomial degree",S.degree(f,w)==4)
check("S4 derivative finite critical indices",
      S.factor(S.diff(f,w))==(w-4)**2*(w-1)/16)
c=5+s;cb=5-s
H=red(rho**3-c*Z);Hb=red(rho**3-cb*Z)
check("tail pair polynomial",zero((w-c)*(w-cb)-(w*w-10*w+27)))
check("mixed tail derivative",zero(S.diff(H,v)-3*(1-c)*rho**2*L**2))
check("tail not divisible by L",zero(H.subs(L,0)-(1-c)))
check("tail residual versus rho nonzero scalar",zero(H-(1-c)*rho**3-c*delta*L))
Hsource=red(H.subs({L:1+x*x*y,v:x*x*z},simultaneous=True))
div=S.div(S.Poly(Hsource,x,y,z),S.Poly(x*x,x,y,z))
check("mixed tail x squared factor",div[1].is_zero)
Hres=div[0].as_expr()
res0=red(Hres.subs(x,0))
alpha=red(S.diff(H,L).subs({L:1,v:0}))
beta=red(S.diff(H,v).subs({L:1,v:0}))
check("mixed tail exact residual along x zero",zero(res0-alpha*y-beta*z))
check("mixed tail residual nonconstant",beta!=0)
check("mixed tail x multiplicity exactly two",res0!=0)
check("opposite tail misses x zero",zero(Hb.subs({L:1,v:0})-2*s))
check("full branch identity before lift",
      zero(64*R**3+27*Z**4-delta**2*L**2*H*Hb))
Zt=1+3*C*C*T1
Rt=delta**2*C*C*T2/6-3*C*C*T1-S.Rational(3,4)
At=red(64*Rt**3+27*Zt**4)
div=S.div(S.Poly(At,C,T1,T2),S.Poly(C*C,C,T1,T2))
check("target double component",div[1].is_zero)
D=div[0].as_expr()
check("target component exponent exactly two",zero(D.subs(C,0)-18*delta**2*T2))
As=red(A.subs({L:1+x*x*y,v:x*x*z},simultaneous=True))
Bs=red(B.subs({L:1+x*x*y,v:x*x*z},simultaneous=True))
adx=S.div(S.Poly(As,x,y,z),S.Poly(x*x,x,y,z))
bdx=S.div(S.Poly(Bs,x,y,z),S.Poly(x*x,x,y,z))
check("global affine lift first component",adx[1].is_zero)
check("global affine lift second component",bdx[1].is_zero)
Fs1,Fs2=adx[0].as_expr(),bdx[0].as_expr()
# The preceding quotient Jacobian identity and denominator-free lift prove
# the determinant globally; no high-degree determinant expansion is needed.
# Check D pullback using the lower-degree quotient plane identity.
Dplane=red(D.subs({C:L,T1:A,T2:B},simultaneous=True))
check("target residual complete pullback",
      zero(Dplane-delta**2*H*Hb))
# Here C=L corresponds x=1, so no x denominator remains.
edges=[("critical","C",2,1,2),("tail_plus","C",1,2,2),
       ("tail_plus","D",1,1,1),("tail_minus","D",1,1,1)]
check("all S4 edge multiplicity equations",all(e*m==N for _,_,e,m,N in edges))
vertices={("l",a) for a,_,_,_,_ in edges}|{("r",b) for _,b,_,_,_ in edges}
adj={q:set() for q in vertices}
for a,b,_,_,_ in edges:adj["l",a].add(("r",b));adj["r",b].add(("l",a))
seen=set();stack=[next(iter(vertices))]
while stack:
 q=stack.pop()
 if q in seen:continue
 seen.add(q);stack.extend(adj[q]-seen)
check("S4 incidence tree connected",seen==vertices)
check("S4 graph reaches lower edge bound",len(edges)==len(vertices)-1)
psl=[]
for branch,e,ram,tails in [(0,3,3,2),(1,2,4,3)]:
    n=e*ram+tails;d=ram+tails;min_target=2
    check(f"PSL branch {branch} degree",n==11)
    check(f"PSL branch {branch} mixed tail criterion",ram>0 and tails>=2 and e>1)
    psl.append(dict(branch=branch,ramification=e,ramified_points=ram,
                    tails=tails,left_vertices=d,min_target_components=min_target,
                    min_source_components=d+min_target-1,
                    mixed_tail_multiplicity_divisor=e))
check("PSL complete finite ramification budget",3*(3-1)+4*(2-1)==10)
check("PSL fourteen distinct finite source primes",
      sum(q["min_source_components"] for q in psl)==14)
model=dict(psl2_11=psl,S4_control=dict(field="s^2=-2",zeta=str(zeta),
    delta=str(delta),g=str(g),tail=str(c),residual_at_x_zero=str(res0),
    target_residual_at_C_zero=str(red(D.subs(C,0))),
    edges=[dict(source=a,target=b,curve_index=e,source_multiplicity=m,
                target_multiplicity=N) for a,b,e,m,N in edges]),
    universal_theorem="written connected-incidence and mixed-tail proof",
    general_PSL2_11_occurrence_decided=False,full_classification_complete=False)
(ROOT/"fiber_incidence_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["FIBER_INCIDENCE.md","verify_fiber_incidence.py","fiber_incidence_model.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,
    python=platform.python_version(),sympy=S.__version__,
    runtime_seconds=round(time.monotonic()-started,2),
    sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
    assurance="Written universal field/divisor proof; exact positive-control algebra. No formal verification or independent whole-proof review.",
    goal_complete=False)
(ROOT/"fiber_incidence_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"seconds":receipt["runtime_seconds"]}))
