#!/usr/bin/env python3
"""QGR1 exact controls; all-degree group exhaustion uses the written proof."""
from pathlib import Path
from math import comb
import hashlib
import json
import sympy as sp
import check_relative_period as rel

HERE=Path(__file__).resolve().parent
x,y,b,p,w,h,t=sp.symbols("x y b p w h t")
r=x*(1+x*y)
u=1+x*y

def zero(F):
    assert sp.cancel(F)==0,sp.cancel(F)

def theta(F):
    return sp.expand(x**2*sp.diff(F,x)-(1+2*x*y)*sp.diff(F,y))

def poisson(F,G):
    return sp.expand(sp.diff(F,p)*sp.diff(G,x)-sp.diff(F,x)*sp.diff(G,p)
                     +sp.diff(F,w)*sp.diff(G,y)-sp.diff(F,y)*sp.diff(G,w))

def projection(G):
    result=sp.expand(G)
    for (a,c),v in sp.Poly(G,x,y).terms():
        if a==2*c:
            result-=v*r**c
    return sp.expand(result)

def degree(G):
    return int(sp.Poly(G,x,y).total_degree()) if G!=0 else -1

def decompose(F):
    # F=theta(G)+Q, with a normalized primitive and the global obstruction.
    work={m:c for m,c in sp.Poly(sp.expand(F),x,y).terms() if c!=0}
    G=sp.S.Zero
    Q=sp.S.Zero
    while work:
        i,j=max(work,key=lambda m:(sum(m),m))
        c=work.pop((i,j))
        if c==0:
            continue
        denominator=i-1-2*j
        if i==0 or denominator==0:
            Q+=c*x**i*y**j
            continue
        coefficient=c/denominator
        G+=coefficient*x**(i-1)*y**j
        if j:
            monomial=(i-1,j-1)
            work[monomial]=sp.cancel(work.get(monomial,0)+j*coefficient)
            if work[monomial]==0:
                del work[monomial]
    return sp.expand(G),sp.expand(Q)

def ore_mul(A,B):
    result=sp.S.Zero
    for (i,),a in sp.Poly(A,b).terms():
        for (j,),c in sp.Poly(B,b).terms():
            deriv=c
            for k in range(i+1):
                result+=comb(i,k)*h**k*a*deriv*b**(i-k+j)
                deriv=theta(deriv)
    return sp.expand(result)

def main():
    B=x**2*p-(1+2*x*y)*w
    D=(1-2*x*y)*p/2+2*y**2*w
    zero(2*(1+2*x*y)*D+4*y**2*B-p)
    zero(2*x**2*D-(1-2*x*y)*B-w)
    zero(poisson(D,r)-sp.Rational(1,2))
    zero(poisson(D,B)+3*y*B)
    zero(ore_mul(b,x)-ore_mul(x,b)-h*x**2)
    zero(ore_mul(b,y)-ore_mul(y,b)+h*(1+2*x*y))
    zero(ore_mul(b,u)-ore_mul(u,b)+h*r)
    A=(x+y)*b**3+x*b+y
    C=(x*y+1)*b**2+x
    product=ore_mul(A,C)
    assert sp.Poly(product,b).degree()==5
    zero(sp.Poly(product,b).coeff_monomial(b**5)-(x+y)*(x*y+1))

    tests=[
        x, y, x*y, x**7*y**2+3*r**6+r**5*x+y**4+3,
        r**8+x**3*y**7-2*r**3+y,
        r**10+sp.Rational(2,3)*x**5*y**2+sp.Rational(5,7)*y**6
    ]
    controls=[]
    for G in tests:
        normalized=projection(sp.expand(G))
        F=theta(G)
        zero(theta(normalized)-F)
        reconstructed,remainder=decompose(F)
        zero(reconstructed-normalized)
        assert remainder==0
        zero(projection(normalized)-normalized)
        assert degree(normalized)<=degree(G)
        assert F!=0 and degree(normalized)==degree(F)-1
        assert rel.normal(F*rel.alpha[0],F*rel.alpha[1])=={}
        shift_p=sp.diff(normalized,x)
        shift_w=sp.diff(normalized,y)
        zero(sp.diff(shift_w,x)-sp.diff(shift_p,y))
        zero(x**2*shift_p-(1+2*x*y)*shift_w-F)
        controls.append({"input_primitive_degree":degree(G),
                         "normalized_primitive_degree":degree(normalized),
                         "correction_degree":degree(F)})
    rejection_controls=[]
    for F in (sp.S.One,x,y,x*y,x**2*y,r*y,r*y**3):
        G,Q=decompose(F)
        zero(F-theta(G)-Q)
        assert projection(G)==G
        assert Q!=0
        assert degree(G)<=degree(F)-1 or G==0
        assert degree(Q)<=degree(F)
        assert rel.normal(F*rel.alpha[0],F*rel.alpha[1])==rel.normal(Q*rel.alpha[0],Q*rel.alpha[1])
        rejection_controls.append({"F":str(F),"obstruction":str(Q)})
    for k in range(8):
        assert projection(sp.expand(r**k))==0
    # A bounded scalar-series potential with infinitely many h coefficients.
    scalar_G=h**2*x/(1-h)
    zero(theta(scalar_G)-h**2*x**2/(1-h))
    reconstructed,remainder=decompose(h**2*x**2/(1-h))
    zero(reconstructed-scalar_G)
    assert remainder==0

    X=x/(1-t*x)
    Y=(1-t*x)*(y-t*u)
    zero(1+X*Y-u*(1-t*x))
    zero(X*(1+X*Y)-r)
    zero(sp.diff(X,t)-X**2)
    zero(sp.diff(Y,t)+(1+2*X*Y))
    zero(theta(X)-X**2)
    zero(theta(Y)+(1+2*X*Y))
    expansion=sp.series(X.subs(t,h**2),h,0,9).removeO()
    for j in range(5):
        assert sp.expand(expansion).coeff(h,2*j)==x**(j+1)

    # Generic irreducibility control: x^2*y+x+h^2 is primitive and linear in y.
    assert sp.gcd(x**2,x+h**2,x)==1
    assert sp.Poly(r+h**2,y).degree()==1
    zero(poisson(D,r+h**2)-sp.Rational(1,2))

    receipt={
        "claim":"QGR1","status":"PASS",
        "ambient_Ore_coordinate_identities":True,
        "Ore_product_order_control":5,
        "primitive_degree_controls":controls,
        "projection_kernel_samples":8,
        "constructive_rejection_controls":rejection_controls,
        "bounded_infinite_scalar_series_control":True,
        "R_preserving_unbounded_formal_flow_control_mod_h":9,
        "scalar_shift_primitive_linear_polynomial_and_bounded_mate":True,
        "coverage_boundary":"Exact finite symbolic controls only. Exhaustion of all center-preserving finite endomorphisms, the full fixed-R stabilizer, sharp primitive bounds and the admission equivalence for every bounded centralizer change use written proofs. Irreducibility gives a written all-degree obstruction to bounded target normalization of r+h^2; a bounded source translation removes the scalar shift in a tuple. No bounded reachability between arbitrary lifts, general two-sided R-fixing or bounded Long lift is established.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("BOUNDED_GAUGE_RIGIDITY.md","check_bounded_gauge.py")},
        "retained_input":{
            "path":"dixmier/dixmier-formal-lifting/check_relative_period.py",
            "sha256":hashlib.sha256(Path(rel.__file__).read_bytes()).hexdigest(),
            "scope":"Imported exact relative-cohomology normal form; this is not an independent replay."}
    }
    (HERE/"bounded_gauge_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
