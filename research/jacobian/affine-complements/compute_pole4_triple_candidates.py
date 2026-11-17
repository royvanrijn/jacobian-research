"""Exact normalized triple-root candidates and their contact quartics.

Normalization coverage is a written geometric question until the atlas
proof is recorded; this producer alone does not claim full coverage.
"""
from pathlib import Path
from time import monotonic
from itertools import combinations
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
delta,b,k,z,w,x=S.symbols("delta b k z w x")
relation=6*delta**2-12*delta+5
field=S.QQ.frac_field(b,k)
f=S.Poly(relation,delta,domain=field)
def red(expr):
    num,den=S.fraction(S.cancel(expr))
    a=S.Poly(num,delta,domain=field)
    d=S.Poly(den,delta,domain=field)
    return S.rem(a*S.invert(d,f),f).as_expr()
def red_poly(expr,*gens):
    p=S.Poly(S.expand(expr),*gens)
    return S.expand(sum(red(coef)*S.prod(v**power for v,power in zip(gens,powers))
                        for powers,coef in p.terms()))
def W(F,G):
    return S.expand(S.diff(F,z)*S.diff(G,w)-S.diff(F,w)*S.diff(G,z))
g=[w/z,w*w/z,(z+w*w)*w/z,(z+w*w)**2/z]
t0=red(3*b*(delta-1)/(3*delta-2))

def case(marked):
    correction=0 if marked else -4*delta+S.Rational(8,3)
    x0=red(-(2*b*b-4*b*t0+(2*delta-1)*k+correction)/(4*(delta-1)))
    ell=red(-2*b*k-4*b*x0+4*delta-S.Rational(4,3)) if marked else red(-2*b*k+4*b-4*b*x0)
    if marked:
        T=g[3]+t0;X=g[2]+x0;Y=b*g[2]+delta*g[1]
        Z=3*b*g[1]+(3*delta-1)*g[0]
    else:
        T=g[0]+t0;X=g[1]+x0;Y=delta*g[2]+b*g[1]
        Z=(3*delta-1)*g[3]+3*b*g[2]
    coeffs=[
        (T,X,12*X*X-12*T*Y+6*k*X-3*ell*T),
        (T,Y,-12*T*X-12*Y-6*k*T-3*ell),
        (T,Z,4*T*T+4*X),
        (X,Y,-12*X-6*k),(X,Z,4*T),(Y,Z,S.Integer(4))]
    # Formula uses Q's +(4/3)Z, already included in exterior coefficients.
    raw=S.expand(sum(v*W(F,G) for F,G,v in coeffs)*z**4)
    assert S.denom(S.cancel(raw))==1
    numerator=red_poly(raw,z,w)
    assert numerator.subs(z,0)==0
    volume=S.expand(numerator/z)
    plane=S.expand(volume.subs(z,x-w*w))
    assert S.Poly(plane,x,w).total_degree()<=3
    quartic=S.cancel(volume.subs(z,0)/w**2)
    assert S.denom(quartic)==1
    q=S.Poly(quartic,w)
    assert q.degree()<=4
    c=[q.coeff_monomial(w**i)/S.binomial(4,i) for i in range(5)]
    M=S.Matrix(3,3,lambda i,j:c[i+j])
    minors=[red(S.det(M.extract(ii,jj))) for ii in combinations(range(3),2)
            for jj in combinations(range(3),2)]
    variables=(delta,b,k)
    polys=[S.together(a).as_numer_denom()[0] for a in minors if a!=0]+[relation]
    gb=S.groebner(polys,*variables,order="lex")
    return dict(name="triple_at_marked" if marked else "triple_away_marked",
                t0=str(t0),x0=str(x0),ell=str(ell),
                plane_volume_numerator=str(plane),contact_quartic=str(quartic),
                contact_minors=list(map(str,minors)),
                groebner_basis=list(map(str,gb.polys)),
                groebner_expressions=[str(p.as_expr()) for p in gb.polys],
                contact_ideal_is_unit=gb.is_zero_dimensional is False and list(gb)==[1])

start=monotonic()
cases=[case(True),case(False)]
for row in cases:
    print(row["name"],"quartic",row["contact_quartic"])
    print("contact ideal basis",row["groebner_expressions"])
payload=dict(as_of="2026-11-14",outcome="EXACT_CANDIDATES_COMPUTED",cases=cases,
             quadratic_relation=str(relation),
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             elapsed_seconds=round(monotonic()-start,6))
(HERE/"pole4_triple_candidates.json").write_text(json.dumps(payload,indent=2)+"\n")
print("seconds",payload["elapsed_seconds"])

