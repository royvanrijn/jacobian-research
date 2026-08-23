#!/usr/bin/env sage -python
"""
Audit the missing marking in the modular q32 -> D12 compiler.

The exact q32 divisor has degree one on precisely the two resolved I9* spinor
components E10a,E10b, so those components are sections of the NEW q32
fibration.  The binary quartic produced by the compiler has old-base
coordinate U.  Both spinor components lie over the old I9* point U=alpha.

Therefore a first, marking-sensitive test is:
    quartic(U=alpha) is a square in GF(p)(V).
If nonzero, the two rational points
    (U,W)=(alpha,+sqrt), (alpha,-sqrt)
are exactly the natural candidate pair E10a/E10b.  If zero, the common
branch point is rational and the terminal tangent direction must distinguish
the two sections.

This script uses the already-saved q32 signatures, so it does not rerun the
expensive resolved RR compiler.
"""

import json
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    p for p in q8_candidates
    if p.exists()
    and json.loads(p.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
),None)
if Q8 is None:
    raise SystemExit("missing exact q8/D13 child")

q8=json.loads(Q8.read_text())
i9=next(x for x in q8["child"]["finite_fibres"] if x["kodaira"]=="I9*")

RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
assert fQ.degree()==1
alphaQ=-fQ[0]/fQ[1]

records=[]
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        d=json.loads(path.read_text())
        p=ZZ(d["prime"])
    except Exception:
        continue
    if d.get("status")!="PASS_Q32_MODP_SIGNATURE":
        continue
    records.append((p,d,path))
records.sort(key=lambda z:int(z[0]))

if not records:
    raise SystemExit("no q32 signature artifacts")

def redq(q,p,F):
    q=QQ(q)
    den=ZZ(q.denominator())%p
    if not den:
        raise ZeroDivisionError
    return F(ZZ(q.numerator())%p)/F(den)

def poly_sqrt(P):
    R=P.parent()
    F=R.base_ring()
    if not P:
        return R.zero()
    fac=P.factor()
    unit=F(fac.unit())
    if not unit.is_square():
        return None
    out=R(unit.sqrt())
    for f,e in fac:
        if int(e)%2:
            return None
        out*=f**(int(e)//2)
    assert out**2==P
    return out

def rat_sqrt(value,R,K):
    value=K(value)
    if not value:
        return K.zero()
    n=R(value.numerator())
    d=R(value.denominator())
    nr=poly_sqrt(n)
    dr=poly_sqrt(d)
    if nr is None or dr is None or not dr:
        return None
    out=K(nr)/K(dr)
    assert out**2==value
    return out

good=[]
for p,sig,path in records:
    F=GF(p)
    RV=PolynomialRing(F,"V")
    V=RV.gen()
    K=RV.fraction_field()

    alpha=redq(alphaQ,p,F)

    coeffs=[]
    for rec in sig["quartic_coefficients"]:
        n=RV([F(v) for v in rec["num"]])
        d=RV([F(v) for v in rec["den"]])
        coeffs.append(K(n)/K(d))

    qa=K.zero()
    pw=F.one()
    for c in coeffs:
        qa += c*K(pw)
        pw *= alpha

    root=rat_sqrt(qa,RV,K)
    square=(root is not None)
    zero=(not qa)

    if root is None:
        rn=rd=-1
    else:
        rn=RV(root.numerator()).degree()
        rd=RV(root.denominator()).degree()

    qn=RV(qa.numerator()).degree() if qa else -1
    qd=RV(qa.denominator()).degree() if qa else 0

    print(
        "Q32SPINOR_QUARTIC|"
        f"prime={p}|alpha={int(alpha)}|"
        f"value_deg={qn}/{qd}|zero={int(zero)}|square={int(square)}|"
        f"sqrt_deg={rn}/{rd}|"
        f"status={'PASS_RATIONAL_PAIR' if square and not zero else 'PASS_BRANCH_POINT' if zero else 'FAIL_NOT_SQUARE'}",
        flush=True,
    )
    if square:
        good.append((int(p),zero,rn,rd))

status=(
    "PASS_SPINOR_SECTION_ANCHOR"
    if len(good)==len(records)
    else "PARTIAL_SPINOR_SECTION_ANCHOR"
)
print(
    "Q32SPINOR_RESULT|"
    f"compatible={len(good)}/{len(records)}|"
    f"nonzero_pairs={sum(not z for _,z,_,_ in good)}|"
    f"branch_points={sum(z for _,z,_,_ in good)}|"
    f"status={status}",
    flush=True,
)

out=LOCAL/"q32-spinor-quartic-section-anchor.json"
out.write_text(json.dumps({
    "schema":"elkies-k3.h3-q32-spinor-quartic-section-anchor.v1",
    "status":status,
    "old_I9star_root_QQ":str(alphaQ),
    "interpretation":(
        "The exact q32 divisor has degree one on E10a and E10b. "
        "A square quartic specialization at the old I9* base point provides "
        "the rational-point pair needed to mark the binary quartic/Jacobian "
        "before choosing a child zero section."
    ),
    "primes":[
        {
            "prime":p,
            "quartic_value_zero":bool(z),
            "sqrt_numerator_degree":rn,
            "sqrt_denominator_degree":rd,
        }
        for p,z,rn,rd in good
    ],
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
