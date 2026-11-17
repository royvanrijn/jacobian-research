"""Bounded exact controls for the written cost and resolvent theorems."""
from pathlib import Path
import hashlib
import json
import platform
import time
from itertools import product
import sympy as s

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
started=time.monotonic();checks=[];modules=[];fibres=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
def rank(rows,p):
    if not rows:return 0
    if not p:return s.Matrix(rows).rank()
    a=[[int(v)%p for v in row] for row in rows]
    pivot=0
    for j in range(len(a[0])):
        k=next((i for i in range(pivot,len(a)) if a[i][j]),None)
        if k is None:continue
        a[pivot],a[k]=a[k],a[pivot]
        inv=pow(a[pivot][j],-1,p)
        a[pivot]=[(v*inv)%p for v in a[pivot]]
        for i in range(len(a)):
            if i!=pivot and a[i][j]:
                c=a[i][j];a[i]=[(v-c*w)%p for v,w in zip(a[i],a[pivot])]
        pivot+=1
        if pivot==len(a):break
    return pivot

# Actual permutation quotients, including p|N and the N=2,p=2 exception.
for n,p in product(range(2,7),[0,2,3,5]):
    Q=s.eye(n-1).row_join(-s.ones(n-1,1))
    J=s.eye(n)[:,:n-1]
    quotient_rows=[];full_rows=[]
    for i in range(n-1):
        perm=s.eye(n);perm.row_swap(i,i+1)
        induced=Q*perm*J
        check(f"N{n} char{p} quotient action well-defined",Q*perm==induced*Q)
        quotient_rows.extend((induced-s.eye(n-1)).tolist())
        full_rows.extend((perm-s.eye(n)).tolist())
    full_fixed=n-rank(full_rows,p)
    quotient_fixed=n-1-rank(quotient_rows,p)
    points=[tuple(int(Q[i,j])%p if p else Q[i,j] for i in range(n-1)) for j in range(n)]
    distinct=len(set(points))==n
    exceptional=n==2 and p==2
    check(f"N{n} char{p} permutation fixed line",full_fixed==1)
    check(f"N{n} char{p} quotient orbit distinctness",distinct!=exceptional)
    check(f"N{n} char{p} quotient fixed space",quotient_fixed==(1 if exceptional else 0))
    check(f"N{n} char{p} full-module minimum architecture",n+1-(full_fixed+1)==n-1)
    if distinct:
        check(f"N{n} char{p} quotient minimum architecture",n-1+2==n+1 and quotient_fixed==0)
    differences=(Q[:,1:]-Q[:,0]*s.ones(1,n-1)).tolist()
    adim=rank(differences,p)
    check(f"N{n} char{p} modular affine-span boundary",
          adim==(n-2 if p and n%p==0 else n-1))
    modules.append({"N":n,"characteristic":p,"full_fixed_dimension":full_fixed,
                    "quotient_fixed_dimension":quotient_fixed,"quotient_points_distinct":distinct,
                    "affine_span_dimension":adim})

# The all-signature plane classification: exact identity and scope controls.
x,y,z,lam=s.symbols("x y z lam")
a,b=s.symbols("a b",positive=True,integer=True)
f=1+2*z+3*z*z;g=2+z
mono=x**2*y**3 # opposite primitive weights (3,-2)
plane=s.Matrix([x*f.subs(z,mono),y*g.subs(z,mono)])
expected=(f*g+2*z*s.diff(f,z)*g+3*z*f*s.diff(g,z)).subs(z,mono)
check("opposite-weight full Jacobian formula",
      s.expand(plane.jacobian([x,y]).det()-expected)==0)
top=s.Poly(f*g+2*z*s.diff(f,z)*g+3*z*f*s.diff(g,z),z)
check("opposite-weight leading positive factor",
      top.LC()==s.Poly(f,z).LC()*s.Poly(g,z).LC()*(1+2*2+3*1))
positive=s.Matrix([x,y+x*x])
check("same-sign positive-weight nonlinear automorphism control",
      positive.jacobian([x,y]).det()==1
      and positive.subs({x:x,y:y-x*x},simultaneous=True)==s.Matrix([x,y]))
for p in [3,5]:
    modular=s.Matrix([x+x**p*y**(p-1),y])
    check(f"char{p} opposite-weight nonlinear countercontrol",
          s.Poly(modular.jacobian([x,y]).det()-1,x,y,modulus=p).is_zero
          and s.Poly(modular[0],x).degree()==p)
    cov=modular[0].subs({x:lam*x,y:y/lam},simultaneous=True)-lam*modular[0]
    check(f"char{p} opposite-weight countercontrol covariance",s.cancel(cov)==0)

# Rank-one canonical Frobenius covariance, including composite degree four.
X,U,V,A=s.symbols("X U V A")
for p,q in [(2,2),(3,3),(2,4),(3,9),(5,5)]:
    phi=X-X**q
    check(f"char{p} q{q} canonical line derivative one",
          s.Poly(s.diff(phi,X)-1,X,modulus=p).is_zero)
    check(f"char{p} q{q} canonical zero fibre whole rank",
          s.Poly(phi,X,modulus=p).degree()==q
          and s.gcd(s.Poly(phi,X,modulus=p),s.Poly(s.diff(phi,X),X,modulus=p)).degree()==0)
for p in [2,3,5]:
    for aa,bb in product(range(1,p),range(p)):
        expr=(aa*X+bb)-(aa*X+bb)**p-aa*(X-X**p)
        check(f"F{p} affine {aa},{bb} line covariance",s.Poly(expr,X,modulus=p).is_zero)
alpha=s.symbols("alpha")
for aa,bb in product([1,alpha,alpha+1],[0,1,alpha,alpha+1]):
    expr=(aa*X+bb)-(aa*X+bb)**4-aa*(X-X**4)
    rem=s.reduced(s.expand(expr),[alpha**2+alpha+1],alpha,X,modulus=2)[1]
    check(f"F4 affine {aa},{bb} canonical line covariance",rem==0)
for p,e in [(2,2),(2,3),(3,2),(5,2)]:
    q=p**e
    T=s.Poly(sum(X**(p**j) for j in range(e)),X,modulus=p)
    outer=s.Poly(X-X**p,X,modulus=p)
    composed=outer.compose(T)
    check(f"char{p} exponent{e} exact nontrivial Keller decomposition",
          composed==s.Poly(X-X**q,X,modulus=p)
          and T.diff()==s.Poly(1,X,modulus=p)
          and outer.diff()==s.Poly(1,X,modulus=p)
          and T.degree()>1 and outer.degree()>1)

# A finite two-variable Keller map has the same constant sign resolvent
# on all complete fibres; source reconstruction is checked in each algebra.
for target in [0,1,2]:
    P=s.Poly(X**3-X-target,X,modulus=3)
    xx=s.Poly(X,X,modulus=3)
    uu=s.Poly(X**3-X,X,modulus=3).rem(P)
    vv=2
    yy=s.Poly(vv-target**2,X,modulus=3)
    out2=(yy+uu**2).rem(P)
    check(f"finite plane fibre {target} complete quotient algebra",
          uu==s.Poly(target,X,modulus=3) and out2==s.Poly(vv,X,modulus=3)
          and P.degree()==3 and s.gcd(P,P.diff()).degree()==0)
    disc=s.discriminant(P.as_expr(),X)
    check(f"finite plane fibre {target} constant discriminant",int(disc)%3==1)
    fibres.append({"characteristic":3,"target_u":target,"target_v":vv,"rank":3,
                    "discriminant_square_class":1})
finite=s.Matrix([x**3-x,y+(x**3-x)**2])
check("finite plane model full Keller determinant",
      s.Poly(finite.jacobian([x,y]).det()+1,x,y,modulus=3).is_zero)
disc=s.discriminant(X+A*X**3-V,X)
check("variable coefficient is constant in the univariate target resolvent",
      s.Poly(disc+A,A,V,modulus=3).is_zero)
check("the promoted coefficient has odd resolvent valuation",
      s.Poly(-A,A,modulus=3).degree()==1)
for p,q in [(2,4),(3,3),(3,9)]:
    inverse=s.Poly(X+(U-1)*X**q-V,X,domain=s.GF(p).frac_field(U,V))
    collapsed=s.Poly(inverse.as_expr().subs(U,1),X,domain=s.GF(p).frac_field(V))
    check(f"q{q} exact nonproper degree-drop divisor",
          inverse.degree()==q and collapsed.degree()==1)

# The modular omitted-point control has an exact polynomial unit-ideal receipt.
for p in [2,3,5]:
    P=lambda value:s.Poly(value,x,y,modulus=p)
    aa=P(1+x)
    ff=P(x)+aa**p*P(y)**p
    gg=P(y-1)-aa**(p-1)*P(y)**p
    check(f"char{p} omitted-point plane full Jacobian",
          ff.diff(x)*gg.diff(y)-ff.diff(y)*gg.diff(x)==P(1))
    cA=(P(y)-gg-P(1))-ff**(p-1)*(P(y)-gg)
    cB=-(aa**p)*(gg**(p-1))*(P(y)-gg)
    check(f"char{p} omitted-point exact Bezout identity",ff*cA+gg*cB==P(1))

manifest=json.loads((HERE/"sources.json").read_text())
for source in manifest["local_sources"]:
    check("source pin "+source["path"],
          hashlib.sha256((RESEARCH/source["path"]).read_bytes()).hexdigest()==source["sha256"])
receipt={
 "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
 "permutation_module_controls":modules,"complete_fibre_controls":fibres,
 "assurance":"Bounded exact controls support written universal arguments. No finite matrix or polynomial sample is used to infer an all-rank representation theorem, a Kummer classification, generic monodromy, torsor versality or a minimum dimension. Not formal verification or independent whole-proof review.",
 "python":platform.python_version(),"sympy":s.__version__,
 "runtime_seconds":round(time.monotonic()-started,3),
 "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [HERE/"PROOF.md",HERE/"LITERATURE_AND_GAPS.md",Path(__file__),HERE/"sources.json"]}
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"module_controls":len(modules),
                  "whole_fibres":len(fibres),"seconds":receipt["runtime_seconds"]}))
