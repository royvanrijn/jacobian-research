#!/usr/bin/env sage -python
"""
Mark the modular q24/orbit85 D12 child using the two old-I9* spinor components.

Because q24 has degree one on C10a,C10b, they are sections of the new D12
fibration.  Evaluate the q24 binary quartic at the old I9* point, map the two
quartic points through the repository's covariant 2-cover to the Jacobian,
then solve 2*S = phi(C10b)-phi(C10a).  The rational half-points are the
primitive spinor MW direction candidates.
"""
import argparse, json
from pathlib import Path
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3/scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
CORE=SCRIPTS/"elliptic_neighbor_compiler.sage"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)
F=GF(p)

SIG=LOCAL/f"q24-orbit85-d12-signature-mod-{p}.json"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((x for x in q8_candidates if x.exists()),None)
for path in (SIG,CORE):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")
if Q8 is None:
    raise SystemExit("missing exact q8/D13 child")

sig=json.loads(SIG.read_text())
q8=json.loads(Q8.read_text())
assert sig["status"] in (
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
)

# Old I9* point in the q8 base coordinate U.
i9=next(x for x in q8["child"]["finite_fibres"] if x["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
assert fQ.degree()==1
alphaQ=-fQ[0]/fQ[1]
den=ZZ(alphaQ.denominator())%p
assert den
alpha=F(ZZ(alphaQ.numerator())%p)/F(den)

RV=PolynomialRing(F,"V")
V=RV.gen()
KV=RV.fraction_field()
RU=PolynomialRing(KV,"U")
U=RU.gen()

def rf(rec):
    return KV(RV([F(v) for v in rec["num"]]))/KV(
        RV([F(v) for v in rec["den"]])
    )

quartic=RU([rf(rec) for rec in sig["quartic_coefficients"]])
assert quartic.degree()==4
jacA=rf(sig["jacobian_A"])
jacB=rf(sig["jacobian_B"])

def poly_sqrt(P):
    P=RV(P)
    if not P:
        return RV.zero()
    fac=P.factor()
    unit=F(fac.unit())
    if not unit.is_square():
        return None
    out=RV(unit.sqrt())
    for fac0,e in fac:
        if int(e)%2:
            return None
        out*=fac0**(int(e)//2)
    assert out**2==P
    return out

def rat_sqrt(value):
    value=KV(value)
    if not value:
        return KV.zero()
    nr=poly_sqrt(value.numerator())
    dr=poly_sqrt(value.denominator())
    if nr is None or dr is None or not dr:
        return None
    out=KV(nr)/KV(dr)
    assert out**2==value
    return out

qa=KV(quartic(alpha))
root=rat_sqrt(qa)
if root is None:
    raise SystemExit("spinor quartic specialization is not a square")
if not qa:
    raise SystemExit(
        "spinor quartic specialization is a branch point; tangent-direction "
        "marking is required instead of the nonzero covariant point map"
    )

exec(compile(CORE.read_text(),str(CORE),"exec"))

plus=transport_binary_quartic_point_to_jacobian(
    quartic,alpha,F(1),root,minimalizing_unit=1
)
minus=transport_binary_quartic_point_to_jacobian(
    quartic,alpha,F(1),-root,minimalizing_unit=1
)
assert KV(plus["standard_a"])==jacA and KV(plus["standard_b"])==jacB
assert KV(minus["standard_a"])==jacA and KV(minus["standard_b"])==jacB

E=EllipticCurve(KV,[0,0,0,jacA,jacB])
Pplus=E(KV(plus["standard_x"]),KV(plus["standard_y"]))
Pminus=E(KV(minus["standard_x"]),KV(minus["standard_y"]))
Dcov=Pminus-Pplus
assert not Dcov.is_zero()
xd,yd=Dcov.xy()

# x-coordinate equation for 2Q=Dcov on y^2=x^3+A*x+B.
XR=PolynomialRing(KV,"Xh")
Xh=XR.gen()
half_poly=(
    Xh**4
    - 4*KV(xd)*Xh**3
    - 2*jacA*Xh**2
    - (8*jacB+4*jacA*KV(xd))*Xh
    + (jacA**2-4*jacB*KV(xd))
)

halves=[]
for fac,e in half_poly.factor():
    if fac.degree()!=1:
        continue
    xq=-fac[0]/fac[1]
    rhs=xq**3+jacA*xq+jacB
    yq=rat_sqrt(rhs)
    if yq is None:
        continue
    for yy in (yq,-yq):
        Q=E(xq,yy)
        if 2*Q==Dcov:
            halves.append(Q)

# Deduplicate.
uniq={}
for Q in halves:
    xq,yq=Q.xy()
    uniq[(str(xq),str(yq))]=Q
halves=list(uniq.values())

def enc(v):
    v=KV(v)
    return {
        "num":[int(x) for x in RV(v.numerator()).list()],
        "den":[int(x) for x in RV(v.denominator()).list()],
    }

def enc_point(P):
    x,y=P.xy()
    return {"x":enc(x),"y":enc(y)}

status="PASS_Q24_D12_SPINOR_PRIMITIVE_DIRECTION" if halves else "NO_RATIONAL_HALF"
payload={
    "schema":"elkies-k3.h3-q24-d12-spinor-anchor-modp.v1",
    "status":status,
    "prime":int(p),
    "old_I9star_root_QQ":str(alphaQ),
    "old_I9star_root_mod_p":int(alpha),
    "quartic_value":enc(qa),
    "covariant_plus":enc_point(Pplus),
    "covariant_minus":enc_point(Pminus),
    "covariant_difference":enc_point(Dcov),
    "primitive_half_points":[enc_point(Q) for Q in halves],
    "half_count":len(halves),
    "interpretation":(
        "C10a,C10b are q24-degree-one sections. The covariant map is a "
        "2-cover, so their image difference is twice the primitive MW "
        "difference. Rational half-points are the explicit primitive spinor "
        "direction candidates for a geometric-zero D12 marking."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-d12-spinor-anchor-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(
    "Q24D12SPINOR|"
    f"prime={p}|alpha={int(alpha)}|quartic_square=1|"
    f"half_count={len(halves)}|status={status}",
    flush=True,
)
print(f"OUTPUT|{OUT}",flush=True)
