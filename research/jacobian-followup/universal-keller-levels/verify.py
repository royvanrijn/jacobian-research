"""Bounded controls for the written coordinate-level and Frobenius theorems."""
from pathlib import Path
from itertools import product
from collections import Counter
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
started=time.monotonic()
checks=[];counts=[];fibres=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)

# Direct full-source counts, without using the proof's strata to count.
for p,J in [(3,[]),(3,[6]),(3,[4]),(3,[4,5]),
            (5,[]),(5,[5]),(5,[4]),(5,[4,5]),
            (7,[7]),(7,[4])]:
    r=len(J)
    actual=Counter();strata=Counter()
    for aj in product(range(p),repeat=r):
        for x,y,z in product(range(p),repeat=3):
            t=(1+x*y)%p
            q=(t*t*z+y*y*(1+3*t))%p
            B=(y+3*x*q+sum(j*a*t*t*pow(x,j-2,p)*pow(q,j,p)
                           for j,a in zip(J,aj)))%p
            actual[B]+=1
            tag="x0" if x==0 else "t0" if t==0 else "q0" if q==0 else "torus"
            strata[B,tag]+=1
    active=any(j%p for j in J)
    exceptional=p==3 and not active
    for beta in range(p):
        if exceptional:pred=p**(r+2)
        elif active:
            pred=(p**(r+2)+p**r-p**(r-1) if beta else
                  p**(r+2)-p**(r+1)+2*p**r-p**(r-1))
        else:
            pred=p**r*(p*p+1 if beta else p*p-p+1)
        check(f"F{p} J{J} level{beta} exact full-source count",actual[beta]==pred)
        if active:
            expected={
                "x0":p**(r+1),"t0":p**(r+1) if beta else 0,
                "torus":(p-1)**3*p**(r-1),
                "q0":(p-2 if beta else p-1)*p**r}
            check(f"F{p} J{J} level{beta} four stratum sizes",
                  all(strata[beta,tag]==value for tag,value in expected.items()))
    counts.append({"field_size":p,"retained_indices":J,"level_counts":[actual[i] for i in range(p)],
                   "active_seed":active,"affine_exception":exceptional})

L=s.symbols("L")
for r in [1,2,4]:
    nonzero=2*L**(r+1)+(L-1)**3*L**(r-1)+(L-2)*L**r
    zero=L**(r+1)+(L-1)**3*L**(r-1)+(L-1)*L**r
    check(f"r{r} nonzero motive simplifies",
          s.expand(nonzero-(L**(r+2)+L**r-L**(r-1)))==0)
    check(f"r{r} zero motive simplifies",
          s.expand(zero-(L**(r+2)-L**(r+1)+2*L**r-L**(r-1)))==0)

x,y,z,u,v,alpha,S=s.symbols("x y z u v alpha S")
tt=1+x*y;qq=tt**2*z+y**2*(1+3*tt)
pi0,c0=s.symbols("pi0 c0")
pi_curve={x:1,y:u-1,z:(pi0/u-(u-1)**2*(1+3*u))/u**2}
check("Pi nonzero level has exact varying unit t",
      s.cancel((tt*qq).subs(pi_curve)-pi0)==0 and tt.subs(pi_curve)==u)
C=x*(5-3*tt)-x**3*z
c_curve={x:u,y:0,z:(2*u-c0)/u**3}
check("C nonzero level exact varying unit x",s.cancel(C.subs(c_curve)-c0)==0)
check("C zero factors are comaximal",s.cancel(C/x).subs(x,0)==2)
check("Pi zero factors are comaximal certificate",
      s.expand(1-((1-x*y)*tt+x*x*y*y))==0)

seed=s.symbols("seed")
Bcontrol=y+3*x*qq+4*seed*tt**2*x**2*qq**4
check("B zero unit certificate xB=2 modulo t",
      s.cancel((x*Bcontrol-2).subs(y,-1/x))==0)

# Frobenius polynomial covariance on the entire planes, not only field points.
for p,q in [(2,2),(3,3),(2,4),(5,5),(3,9)]:
    K=[u,x+(u-1)*x**q]
    determinant=s.Matrix(K).jacobian([u,x]).det()
    check(f"char{p} q{q} full Jacobian one",
          s.Poly(determinant-1,u,x,modulus=p).is_zero)
    root=s.Poly(x-x**q,x,modulus=p)
    check(f"char{p} q{q} marked whole fibre degree and etaleness",
          root.degree()==q and s.gcd(root,root.diff()).degree()==0)
for p in [2,3,5]:
    for a,b in product(range(1,p),range(p)):
        lhs=(a*x+b)+(u-1)*(a*x+b)**p
        rhs=a*(x+(u-1)*x**p)+b*u
        check(f"F{p} affine element {a},{b} global covariance",
              s.Poly(lhs-rhs,u,x,modulus=p).is_zero)
    poly=s.prod(x-a for a in range(p))-(x**p-x)
    check(f"F{p} marked fibre has exactly its p labels",s.Poly(poly,x,modulus=p).is_zero)
    orbit={(b,(a+b)%p) for a in range(1,p) for b in range(p)}
    check(f"F{p} affine group doubly transitive",
          orbit=={(a,b) for a in range(p) for b in range(p) if a!=b})

# F4, represented by F2[alpha]/(alpha^2+alpha+1).
elements=[0,1,alpha,alpha+1]
def reduce4(expr,*gens):
    return s.reduced(s.expand(expr),[alpha**2+alpha+1],alpha,*gens,modulus=2)[1]
for a,b in product(elements[1:],elements):
    lhs=(a*x+b)+(u-1)*(a*x+b)**4
    rhs=a*(x+(u-1)*x**4)+b*u
    check(f"F4 affine {a},{b} polynomial covariance",reduce4(lhs-rhs,u,x)==0)
check("F4 marked fibre splits as its four labels",
      reduce4(s.prod(x-a for a in elements)-(x**4-x),x)==0)
def mult4(a,b):
    out=0
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return out
orbit={(b,mult4(a,1)^b) for a in range(1,4) for b in range(4)}
check("F4 affine group double transitivity, composite degree four",
      orbit=={(a,b) for a in range(4) for b in range(4) if a!=b})
check("affine source has no translation-fixed origin",s.expand((x+1)-x)==1)
check("target linear affine-group action fixes marked origin",
      all((a*v+b*u).subs({u:0,v:0})==0 for a,b in product(elements[1:],elements)))

# Complete cubic algebras of all three splitting types over F3.
for expression,label in [
    (S**3-S,"split 1+1+1"),
    (S**3+S,"type 1+2"),
    (S**3-S-1,"field degree 3")]:
    P=s.Poly(expression,S,modulus=3)
    d1=int(P.nth(1))%3;d0=int(P.nth(0))%3
    zval=pow(d1,-1,3);uval=(1+zval)%3;vval=(-d0*zval)%3
    E=s.Poly(S+(uval-1)*S**3-vval,S,modulus=3)
    check(label+" normalized entire fibre algebra",E==P.mul_ground(zval))
    check(label+" all roots etale",s.gcd(P,P.diff()).degree()==0)
    root_source=s.Poly(S+(uval-1)*S**3,S,modulus=3).rem(P)
    check(label+" source reconstructs target",root_source==s.Poly(vval,S,modulus=3))
    if label.startswith("field"):check(label+" irreducible",P.is_irreducible)
    fibres.append({"label":label,"field":"F3","rank":3,"target_u":uval,"target_v":vval})
check("split cubic has three distinct rational labels",
      all(s.Poly(S**3-S,S,modulus=3).eval(a)==0 for a in range(3)))
check("1+2 residual quadratic irreducible",s.Poly(S*S+1,S,modulus=3).is_irreducible)

# Fixed discriminant, checked with polynomial discriminants before reduction.
aa,bb,cc,V=s.symbols("aa bb cc V")
disc=s.discriminant(aa*S**3+bb*S+cc-V,S)
check("univariate char3 discriminant independent of target",
      s.Poly(disc+aa*bb**3,aa,bb,cc,V,modulus=3).is_zero)
for p,n,poly in [(3,6,S**6+S**3+S-V),(5,5,S**5+S-V)]:
    disc=s.discriminant(poly,S)
    expected=(-1)**(n*(n-1)//2)
    check(f"univariate char{p} degree{n} constant discriminant",
          s.Poly(disc-expected,V,modulus=p).is_zero)
check("characteristic-three positive one-dimensional fibre degree is allowed",
      s.Poly(s.diff(S**3-S,S),S,modulus=3)==s.Poly(-1,S,modulus=3))
# The permutation representation still has its fixed line in wild characteristic.
line=s.Matrix([u,u,u])
generators=[s.Matrix([[0,1,0],[1,0,0],[0,0,1]]),
            s.Matrix([[1,0,0],[0,0,1],[0,1,0]])]
check("S3 wild permutation fixed-line control",all(g*line==line for g in generators))
check("natural three-label set has no full S3 fixed label",
      all(any(g*s.eye(3)[:,j]!=s.eye(3)[:,j] for g in generators) for j in range(3)))

manifest=json.loads((HERE/"sources.json").read_text())
for item in manifest["local_sources"]:
    check("source pin "+item["path"],
          hashlib.sha256((RESEARCH/item["path"]).read_bytes()).hexdigest()==item["sha256"])
receipt={
    "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
    "full_source_count_controls":counts,"complete_fibre_controls":fibres,
    "assurance":"Bounded exact controls supplement written all-degree and all-field proofs. Point-count specialization, Hilbert 90, associated-bundle descent and generic monodromy are not inferred from these finite examples and are not formally verified.",
    "python":platform.python_version(),"sympy":s.__version__,
    "runtime_seconds":round(time.monotonic()-started,3),
    "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
              [HERE/"PROOF.md",HERE/"LITERATURE_AND_GAPS.md",Path(__file__),HERE/"sources.json"]}
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"full_source_counts":len(counts),
                  "whole_fibres":len(fibres),"seconds":receipt["runtime_seconds"]}))
