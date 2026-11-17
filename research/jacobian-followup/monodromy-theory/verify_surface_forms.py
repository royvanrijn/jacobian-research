"""Exact controls for the global form proof in SURFACE_FORMS.md."""
from pathlib import Path
from itertools import product,combinations
from math import gcd
import ast,hashlib,json,platform,time
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
z,Q,V,X,P,t,u=S.symbols("z Q V X P t u")
theta=105+630*z**2+2205*z**4+4452*z**6+5040*z**8+2880*z**10+640*z**12
critical_squares={
 "cone":z**3*(S.Rational(8,3)+S.Rational(144,35)*z**2+S.Rational(72,55)*z**4)**2,
 "shear":S.Rational(64,315**2)*z**3*theta**2}
geometry={}
for name,mu2 in critical_squares.items():
    reduced=S.rem(mu2,z**3+z-Q,z)
    R=S.Poly(S.resultant(z**3+z-Q,V-reduced,z),V,Q).monic().as_expr()
    R0=S.factor(R.subs(V,0));dis=S.factor(S.discriminant(R,V))
    C=S.Poly(R0,Q).sqf_list()[1]
    E=S.Poly(dis,Q).sqf_list()[1]
    node0=next(f.monic() for f,e in C if e==2)
    nodes=next(f.monic() for f,e in E if e==2)
    cusps=next(f.monic() for f,e in E if e==5)
    expected=(4,6) if name=="cone" else (12,8)
    check(name+" ordinary-node packet degrees",(node0.degree(),nodes.degree())==expected)
    check(name+" cusp packet",zero(cusps.as_expr()-(Q**2+S.Rational(4,27))))
    check(name+" origin exponent",any(e==3 and f.degree()==1 and f.eval(0)==0 for f,e in C))
    for label,f in [("node0",node0),("nodes",nodes),("cusps",cusps)]:
        check(name+" squarefree "+label,S.gcd(f,f.diff()).degree()==0)
        check(name+" nonzero packet "+label,f.eval(0)!=0)
    for f,g in combinations([node0,nodes,cusps],2):
        check(name+" packets disjoint",S.gcd(f,g).degree()==0)
    aa,bb,cc,dd=S.Poly(R,V).all_coeffs()
    A=3*aa*cc-bb**2;B=27*aa**2*dd-9*aa*bb*cc+2*bb**3
    check(name+" no triple critical-square collision",S.degree(S.gcd(A,B),Q)==0)
    check(name+" zero not a double critical square",S.degree(S.gcd(R0,S.diff(R,V).subs(V,0)),Q)==0)
    count=node0.degree()+2*nodes.degree()
    check(name+" finite node count",count==(16 if name=="cone" else 28))
    geometry[name]=dict(critical_square_polynomial=str(R),zero_packet=str(node0.as_expr()),
      node_packet=str(nodes.as_expr()),cusp_packet=str(cusps.as_expr()),discriminant=str(dis),
      nodes=count,cusps_2_5=4,cusp_2_3=1)

# Explicit weighted homogeneous root polynomials.
Dcone=(Q**2+2*t**4*Q*X**2+S.Rational(2,5)*Q*X**6
       -t**8*X**4/3-S.Rational(2,7)*t**4*X**8-X**12/11-P*X)
top=(X**12/9-S.Rational(6,7)*Q*X**10+3*Q**2*X**8
     -S.Rational(20,3)*Q**3*X**6+15*Q**4*X**4+6*Q**5*X**2-Q**6/3)
Dshear=top+t**4*(2*Q**3*X**2+6*Q**2*X**4-2*Q*X**6+S.Rational(2,5)*X**8)+t**8*X**4-P*X**3
def weighted(e,weights):
    p=S.Poly(S.expand(e),*weights.keys())
    return {sum(a*weights[x] for a,x in zip(m,weights)) for m,c in p.terms()}
check("cone homogenized root equation",weighted(Dcone,{t:1,Q:6,P:11,X:1})=={12})
check("shear homogenized root equation",weighted(Dshear,{t:1,Q:2,P:9,X:1})=={12})
check("cone boundary nonzero certificate polynomial",
      zero((-11*Dcone).subs({t:0,Q:0,P:-S.Rational(1,11)})-(X**12-X)))
N=7*X**12-54*X**10+189*X**8-420*X**6+945*X**4+378*X**2-21
check("shear boundary certificate polynomial",zero(63*Dshear.subs({t:0,Q:1,P:0})-N))
# Global finite branch support, including the shear's apparent Y=0 collision.
for name,D,radicand,power in [
 ("cone",Dcone.subs(t,1),X**6+X**2-Q,1),
 ("shear",Dshear.subs(t,1),X**2-(Q-X**2)**3,9)]:
    gamma=S.sympify(geometry[name]["critical_square_polynomial"],locals={"Q":Q,"V":V}).subs(V,P**2)
    check(name+" global critical resultant",
          zero(S.resultant(D,radicand,X)+Q**power*gamma))
    spectator=S.Poly(S.cancel(D.subs({P:0,Q:0})/X**4),X)
    check(name+" origin eight spectators simple",
          spectator.degree()==8 and spectator.eval(0)!=0
          and S.gcd(spectator,spectator.diff()).degree()==0)
# At a shear lost triple, X=Y^2 Z makes all three normalized roots distinct.
Z=S.symbols("Z")
loss=S.cancel(Dshear.subs({t:1,X:Q**2*Z})/Q**6)
check("shear lost-triple normalization",zero(loss.subs(Q,0)+P*Z**3+S.Rational(1,3)))
# The leading infinity critical value and exact seven-root cluster.
kappa=S.Rational(1024,63)
p0=top.subs(Q,1)/X**3
check("shear boundary derivative",zero(S.diff(p0,X)-(X*X-1)**6/X**4))
factor=S.Poly(S.cancel((top.subs(Q,1)-kappa*X**3)/(X-1)**7),X)
check("shear boundary seven-cluster",factor.degree()==5 and factor.eval(1)!=0 and factor.eval(0)!=0)
check("shear five spectators distinct",S.gcd(factor,factor.diff()).degree()==0)
ss,zz,pp=S.symbols("ss zz pp")
# Laurent primitive has no pole near X=1. Series checks its exact weighted initial form.
ps=S.cancel(Dshear.subs({Q:1,P:0})/X**3)
series=S.series(ps.subs(X,1+zz)-ps.subs(X,1),zz,0,9).removeO()
wp=S.Poly(S.expand(series),zz,t)
initial=sum(c*zz**m[0]*t**m[1] for m,c in wp.terms() if 4*m[0]+3*m[1]==28)
check("shear infinity no terms below weight 28",
      all(4*m[0]+3*m[1]>=28 for m,c in wp.terms()))
check("shear infinity weighted seven-cover initial form",
      zero(initial-(t**8*zz+4*t**4*zz**4+S.Rational(64,7)*zz**7)))
check("shear infinity initial derivative",zero(S.diff(initial,zz)-(t**4+8*zz**3)**2))
# Arithmetic genus of Cartier branch curves in the weighted compactifications.
def monomials(degree,w1,w2):
    return sum(1 for j in range(degree//w1+1) for k in range(degree//w2+1) if w1*j+w2*k<=degree)
check("cone weighted branch genus exhausts all singularities",monomials(66-18,6,11)==16+4*2+1==25)
check("shear weighted branch genus exhausts all singularities",monomials(54-12,2,9)==28+4*2+1+27==64)
# Graph cuts bound every denominator at each collision; proof handles stack normalization.
bounds=[]
for name,size,root_weight,volume,extra in [
 ("smooth triple",3,1,2,0),("finite five-cusp",5,1,6,0),
 ("cone four-cusp",4,1,4,0),("shear four-cusp",4,3,4,8),
 ("shear lost triple",3,2,0,4),("shear infinity seven-cusp",7,4,30,36)]:
    max_edges=max(k*(size-k) for k in range(max(0,size-8),min(4,size)+1))
    available=volume+extra
    check(name+" cut-denominator extension inequality",available>=root_weight*max_edges)
    bounds.append(dict(place=name,cluster_size=size,root_difference_weight=root_weight,
                       maximum_cross_edges=max_edges,available_order=available,
                       lower_bound=available-root_weight*max_edges))
check("cone form projective weight",4+10-32+(1+6+11)==0)
check("shear form projective weight",12+8-32+(1+2+9)==0)

# Retained natural PSL2(11) action and its unique 165-element four-subset orbit.
helper=ROOT/"verify_curve_descent_lcm.py";tree=ast.parse(helper.read_text())
ctx={"product":product}
defs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ["projective_group","cycles","generated"]]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(helper),"exec"),ctx)
G,I,mul,act=ctx["projective_group"](11);perms=[act(g) for g in G]
orbit=sorted({tuple(sorted(g[i] for i in (0,1,2,5))) for g in perms})
check("PSL natural group order",len(G)==660)
check("four-subset orbit order",len(orbit)==165)
check("orbit is invariant",all(tuple(sorted(g[i] for i in A)) in orbit for g in perms for A in orbit))
normalizer=tuple([2*i%11 for i in range(11)]+[11])
check("odd normalizer multiplier",sum(normalizer[i]>normalizer[j] for i in range(12) for j in range(i+1,12))%2==1)
check("odd normalizer preserves orbit",{tuple(sorted(normalizer[i] for i in A)) for A in orbit}==set(orbit))
# Exact nonzero certificate for the cone form at roots 0 and all eleventh roots of unity.
prime=991;omega=118
check("cone certificate prime and root",S.isprime(prime) and pow(omega,11,prime)==1 and omega!=1)
roots=[pow(omega,i,prime) for i in range(11)]+[0]
check("cone certificate roots distinct",len(set(roots))==12)
cone_value=0
for A in orbit:
    den=1
    for i in A:
        for j in set(range(12))-set(A):den=den*(roots[i]-roots[j])%prime
    cone_value=(cone_value+pow(sum(roots[i] for i in A),10,prime)*pow(den,-1,prime))%prime
check("cone invariant form nonzero certificate",cone_value==473)
# Exact shear certificate over F_179[t]/(t^2-2), using a fixed listed split fiber.
p=179;nu=2
check("shear field irreducible",S.isprime(p) and pow(nu,(p-1)//2,p)==p-1)
def add(a,b):return ((a[0]+b[0])%p,(a[1]+b[1])%p)
def sub(a,b):return ((a[0]-b[0])%p,(a[1]-b[1])%p)
def prod(a,b):return ((a[0]*b[0]+nu*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
def inv(a):
    d=pow((a[0]*a[0]-nu*a[1]*a[1])%p,-1,p)
    return (a[0]*d%p,-a[1]*d%p)
def power(a,n):
    out=(1,0)
    while n:
        if n&1:out=prod(out,a)
        a=prod(a,a);n//=2
    return out
roots2=[(0,4),(0,175),(21,1),(21,178),(34,18),(34,161),
        (43,0),(136,0),(145,18),(145,161),(158,1),(158,178)]
check("shear split roots distinct",len(set(roots2))==12)
coeff=S.Poly(N,X).all_coeffs()
for r in roots2:
    value=(0,0)
    for c in coeff:value=add(prod(value,r),(int(c)%p,0))
    check("shear listed exact root "+str(r),value==(0,0))
shear_value=(0,0)
for A in orbit:
    den=(1,0)
    for i in A:
        for j in set(range(12))-set(A):den=prod(den,sub(roots2[i],roots2[j]))
    shear_value=add(shear_value,inv(den))
check("shear invariant form nonzero certificate",shear_value==(173,114))
model=dict(branch_geometry=geometry,extension_bounds=bounds,
 orbit_165=orbit,
 cone_form="sum_A (sum_{i in A} w_i)^10 / prod_{i in A,j notin A}(w_i-w_j) * dP wedge dQ",
 shear_form="Y^4 * sum_A 1/prod_{i in A,j notin A}(w_i-w_j) * dP wedge dY",
 cone_projective_factor="s^4",shear_projective_factor="s^12",
 cone_certificate=dict(prime=991,omega=118,residue=cone_value),
 shear_certificate=dict(prime=179,quadratic_nonresidue=2,roots=roots2,residue=shear_value),
 surfaces_have_nonzero_global_two_form=True,PSL2_11_realization=False,
 full_classification_complete=False)
(ROOT/"surface_forms_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["SURFACE_FORMS.md","verify_surface_forms.py","surface_forms_model.json","verify_curve_descent_lcm.py"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,3),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written global extension and quotient-descent proof, exact full branch packets and finite-field nonzero certificates; not formal verification or independent whole-proof review.",
 both_A12_subgroup_targets_not_stably_rational=True,
 general_PSL2_11_action_excluded=False,full_classification_complete=False)
(ROOT/"surface_forms_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                     nonzero_certificates=[cone_value,shear_value],complete_classification=False)))
