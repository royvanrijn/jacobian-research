"""Small exact Bezout certificates excluding fourth-power contact in both triple charts."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
d,b,k,w=S.symbols("delta b k w")
rel=6*d*d-12*d+5
field=S.QQ.frac_field(b,k)
f=S.Poly(rel,d,domain=field)
def red(expr):
    n,q=S.fraction(S.cancel(expr))
    return S.rem(S.Poly(n,d,domain=field)*S.invert(S.Poly(q,d,domain=field),f),f).as_expr()
def coeff_red(expr,var):
    return S.expand(sum(red(c)*var**power[0] for power,c in S.Poly(expr,var).terms()))
start=monotonic()
data=json.loads((HERE/"pole4_triple_candidates.json").read_text())
out=[]
for row in data["cases"]:
    G=S.Poly(S.sympify(row["contact_quartic"]),w)
    values=[G.coeff_monomial(w**i)/S.binomial(4,i) for i in range(5)]
    if row["name"]=="triple_away_marked":
        values=list(reversed(values))
    c0,c1,c2,c3,c4=values
    F=[red(c0*c2-c1*c1),red(c0*c0*c3-c1**3),red(c0**3*c4-c1**4)]
    a=red(S.diff(F[0],k))
    assert not a.has(k)
    ks=red(-F[0].subs(k,0)/a)
    R=[red(v.subs(k,ks)) for v in F[1:]]
    print(row["name"],"kstar",ks)
    print("residuals",list(map(S.factor,R)))
    # Residuals have the structure A*b^3+C and D*b^4+E*b.
    rp=S.Poly(R[0],b);rq=S.Poly(R[1],b)
    if rp.degree()==0 and R[0]!=0:
        B3=red(1/R[0]);B4=S.Integer(0)
    elif rp.degree()==1 and rq.degree()==2:
        A=rp.coeff_monomial(b);D=rq.coeff_monomial(b*b);E=rq.coeff_monomial(1)
        assert S.expand(R[0]-A*b)==0 and S.expand(R[1]-D*b*b-E)==0
        B3=red(-D*b/(A*E));B4=red(1/E)
    else:
        assert rp.degree()==3 and rq.degree()==4
        A,C=rp.coeff_monomial(b**3),rp.coeff_monomial(1)
        D,E=rq.coeff_monomial(b**4),rq.coeff_monomial(b)
        assert S.expand(R[0]-A*b**3-C)==0 and S.expand(R[1]-D*b**4-E*b)==0
        J=red(E-D*C/A)
        assert J!=0 and C!=0
        B3=red(1/C+D*b**3/(C*J))
        B4=red(-A*b*b/(C*J))
    assert red(B3*R[0]+B4*R[1]-1)==0
    Q3=coeff_red(S.div(S.Poly(F[1]-R[0],k),S.Poly(k-ks,k))[0].as_expr()/a,k)
    Q4=coeff_red(S.div(S.Poly(F[2]-R[1],k),S.Poly(k-ks,k))[0].as_expr()/a,k)
    assert red(F[1]-R[0]-Q3*F[0])==0
    assert red(F[2]-R[1]-Q4*F[0])==0
    B2=red(-B3*Q3-B4*Q4)
    remainder=S.expand(1-B2*F[0]-B3*F[1]-B4*F[2])
    q,r=S.div(S.Poly(remainder,d,b,k,domain=S.QQ),S.Poly(rel,d,b,k,domain=S.QQ))
    assert r.is_zero
    multipliers=[q.as_expr(),B2,B3,B4]
    equations=[rel]+F
    assert S.expand(sum(u*v for u,v in zip(multipliers,equations))-1)==0
    out.append(dict(name=row["name"],reverse_contact_coefficients=row["name"]=="triple_away_marked",
                    equations=list(map(str,equations)),multipliers=list(map(str,multipliers)),
                    solved_k=str(ks),residuals=list(map(str,R)),
                    constant_coefficient=str(c0),
                    identity="sum multipliers[i]*equations[i] == 1"))
payload=dict(as_of="2026-11-14",outcome="PASS",certificates=out,
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             elapsed_seconds=round(monotonic()-start,6))
(HERE/"pole4_contact_certificates.json").write_text(json.dumps(payload,indent=2)+"\n")
print("both exact identities verified; seconds",payload["elapsed_seconds"])

