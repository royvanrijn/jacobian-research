#!/usr/bin/env sage -python
"""
Use the rational spinor point on each modular q32 binary quartic to turn the
genus-one curve into a POINTED elliptic curve, then identify it explicitly
with the short Jacobian already emitted by the q32 compiler.

For a shifted quartic
    w^2 = a*u^4 + b*u^3 + c*u^2 + d*u + r^2
with chosen zero (u,w)=(0,r), the classical quartic->Weierstrass map gives
    a1=d/r
    a2=c-d^2/(4r^2)
    a3=2*r*b
    a4=-4*r^2*a
    a6=a2*a4.

The opposite spinor (0,-r) maps to
    Q=(-a2, a1*a2-a3).

After completing square/cube and applying the invariant normalization used by
binary_quartic_invariants(), the compiler's short Jacobian is obtained by
    x_jac = 9*x_short
    y_jac = 27*y_short.

Thus this script extracts an explicit modular D12 section corresponding to the
second spinor, with the first spinor chosen as zero.
"""

import json
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"

q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next((
    p for p in q8_candidates
    if p.exists()
    and json.loads(p.read_text()).get("status") == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
), None)
if Q8 is None:
    raise SystemExit("missing exact q8/D13 child")

q8 = json.loads(Q8.read_text())
i9 = next(x for x in q8["child"]["finite_fibres"] if x["kodaira"] == "I9*")

RQ = PolynomialRing(QQ, "U")
fQ = RQ(str(i9["factor"]))
assert fQ.degree() == 1
alphaQ = -fQ[0] / fQ[1]

records = []
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig = json.loads(path.read_text())
        p = ZZ(sig["prime"])
    except Exception:
        continue
    if sig.get("status") != "PASS_Q32_MODP_SIGNATURE":
        continue
    records.append((p, sig, path))
records.sort(key=lambda z: int(z[0]))

if not records:
    raise SystemExit("no q32 signature artifacts")


def redq(q, p, F):
    q = QQ(q)
    den = ZZ(q.denominator()) % p
    if not den:
        raise ZeroDivisionError
    return F(ZZ(q.numerator()) % p) / F(den)


def rf_from_rec(rec, RV, K, F):
    n = RV([F(v) for v in rec["num"]])
    d = RV([F(v) for v in rec["den"]])
    return K(n) / K(d)


def poly_sqrt(P):
    R = P.parent()
    F = R.base_ring()
    if not P:
        return R.zero()
    fac = P.factor()
    unit = F(fac.unit())
    if not unit.is_square():
        return None
    out = R(unit.sqrt())
    for f, e in fac:
        e = int(e)
        if e % 2:
            return None
        out *= f ** (e // 2)
    assert out**2 == P
    return out


def rat_sqrt(value, RV, K):
    value = K(value)
    if not value:
        return K.zero()
    n = RV(value.numerator())
    d = RV(value.denominator())
    nr = poly_sqrt(n)
    dr = poly_sqrt(d)
    if nr is None or dr is None or not dr:
        return None
    out = K(nr) / K(dr)
    assert out**2 == value
    return out


def norm_rf(value, RV, K):
    value = K(value)
    n = RV(value.numerator())
    d = RV(value.denominator())
    lc = d.leading_coefficient()
    n /= lc
    d /= lc
    return {
        "num_degree": int(n.degree()),
        "den_degree": int(d.degree()),
        "num": [int(x) for x in n.list()],
        "den": [int(x) for x in d.list()],
    }


def general_weierstrass_holds(x, y, a1, a2, a3, a4, a6):
    return y**2 + a1*x*y + a3*y == x**3 + a2*x**2 + a4*x + a6


good = []
for p, sig, path in records:
    F = GF(p)
    RV = PolynomialRing(F, "V")
    V = RV.gen()
    K = RV.fraction_field()

    alpha = redq(alphaQ, p, F)

    coeffs = [
        rf_from_rec(rec, RV, K, F)
        for rec in sig["quartic_coefficients"]
    ]
    if len(coeffs) != 5:
        raise ArithmeticError(f"prime {p}: expected quartic degree 4")

    RT = PolynomialRing(K, "T")
    T = RT.gen()
    quartic = sum(RT(coeffs[i]) * T**i for i in range(5))
    shifted = quartic(T + K(alpha))

    r2 = K(shifted[0])
    r = rat_sqrt(r2, RV, K)
    if r is None or not r:
        raise ArithmeticError(f"prime {p}: spinor value is not a nonzero square")

    # w^2 = a*T^4 + b*T^3 + c*T^2 + d*T + r^2
    a = K(shifted[4])
    b = K(shifted[3])
    c = K(shifted[2])
    d = K(shifted[1])

    four = K(4)
    two = K(2)
    twelve = K(12)

    a1 = d / r
    a2 = c - d**2 / (four * r**2)
    a3 = two * r * b
    a4 = -four * r**2 * a
    a6 = a2 * a4

    b2 = a1**2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2

    c4p = b2**2 - 24*b4
    c6p = -b2**3 + 36*b2*b4 - 216*b6
    deltap = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6
    assert deltap

    jacA = rf_from_rec(sig["jacobian_A"], RV, K, F)
    jacB = rf_from_rec(sig["jacobian_B"], RV, K, F)

    # Compiler convention: A=-27 I, B=-27 J.
    # The pointed model has c4=16 I, c6=32 J, hence u=3 scaling.
    c4j = -48 * jacA
    c6j = -864 * jacB
    inv4 = (c4j == 81*c4p)
    inv6 = (c6j == 729*c6p)
    if not (inv4 and inv6):
        raise ArithmeticError(
            f"prime {p}: pointed/Jacobian invariant normalization mismatch"
        )

    # Opposite spinor relative to chosen spinor zero.
    xg = -a2
    yg = a1*a2 - a3
    assert general_weierstrass_holds(xg, yg, a1, a2, a3, a4, a6)

    # Complete square/cube:
    # Y = y + (a1*x+a3)/2
    # X = x + b2/12
    xs = xg + b2/twelve
    ys = yg + (a1*xg + a3)/two
    shortA = -c4p / 48
    shortB = -c6p / 864
    assert ys**2 == xs**3 + shortA*xs + shortB

    # Match compiler short Jacobian.
    xj = 9 * xs
    yj = 27 * ys
    assert jacA == 81*shortA
    assert jacB == 729*shortB
    assert yj**2 == xj**3 + jacA*xj + jacB

    # Choosing the other spinor as zero must negate the section.
    rn = -r
    na1 = d / rn
    na2 = c - d**2 / (four * rn**2)
    na3 = two * rn * b
    nxg = -na2
    nyg = na1*na2 - na3
    nb2 = na1**2 + 4*na2
    nxs = nxg + nb2/twelve
    nys = nyg + (na1*nxg + na3)/two
    nxj = 9*nxs
    nyj = 27*nys
    sign_ok = (nxj == xj and nyj == -yj)
    assert sign_ok

    xr = norm_rf(xj, RV, K)
    yr = norm_rf(yj, RV, K)
    rr = norm_rf(r, RV, K)

    print(
        "Q32POINTED|"
        f"prime={int(p)}|"
        f"c4scale={int(inv4)}|c6scale={int(inv6)}|"
        f"section_x_deg={xr['num_degree']}/{xr['den_degree']}|"
        f"section_y_deg={yr['num_degree']}/{yr['den_degree']}|"
        f"sign_swap={int(sign_ok)}|"
        "status=PASS_MARKED_SECTION",
        flush=True,
    )

    good.append({
        "prime": int(p),
        "alpha": int(alpha),
        "spinor_sqrt": rr,
        "marked_section_x": xr,
        "marked_section_y": yr,
        "sign_swap_negates_section": bool(sign_ok),
        "status": "PASS_MARKED_SECTION",
    })

status = (
    "PASS_POINTED_Q32_D12_SPINOR_MARKING"
    if len(good) == len(records)
    else "PARTIAL_POINTED_Q32_D12_SPINOR_MARKING"
)

out = LOCAL / "q32-pointed-spinor-weierstrass-anchor.json"
out.write_text(json.dumps({
    "schema": "elkies-k3.h3-q32-pointed-spinor-weierstrass-anchor.v1",
    "status": status,
    "old_I9star_root_QQ": str(alphaQ),
    "interpretation": (
        "Choose E10a as zero on the q32 quartic. The opposite spinor E10b "
        "maps to the displayed rational section on the compiler's short "
        "D12 Jacobian. Reversing the spinor choice fixes x and negates y."
    ),
    "normalization": {
        "pointed_to_compiler_x_scale": 9,
        "pointed_to_compiler_y_scale": 27,
        "compiler_c4_over_pointed_c4": 81,
        "compiler_c6_over_pointed_c6": 729,
    },
    "primes": good,
}, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{out}", flush=True)
print(
    "Q32POINTED_RESULT|"
    f"compatible={len(good)}/{len(records)}|"
    f"status={status}",
    flush=True,
)
