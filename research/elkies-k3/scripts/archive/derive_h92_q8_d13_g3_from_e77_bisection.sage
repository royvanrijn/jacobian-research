#!/usr/bin/env sage -python
"""
Reduce the explicit old_E7_7 q8 bisection to the missing primitive D13
Mordell-Weil generator G3.

Previously certified lattice facts in the anchored q8 frame:
    O  = II*_E8_1
    G1 = E6_1 = (1,0,0,0), height 3/4
    AJ(old_E7_7) = (1,0,1,0) = G1 + G3
    G3 = (0,0,1,0), height 11/4, P.O=1

The old_E7_7 curve has q8 degree two.  This script performs its exact
Abel-Jacobi reduction at equation level:

  * reconstruct the corrected q8 pencil from the certified QQ kernel;
  * restrict U=f1/f0 to the exact q6 section old_E7_7 and obtain the
    quadratic equation H_U(T)=0 for its two points over QQ(U);
  * recover the corresponding quartic W-coordinate from the actual old
    Weierstrass x-coordinate and the chord-discriminant square root;
  * pass to L=QQ(U)[theta]/H and add the two conjugate points on the
    branch-anchored elliptic curve;
  * prove the sum descends to QQ(U);
  * map it to the certified minimal I9*+9I1 D13 equation;
  * subtract both IV* branch points.  Exactly one result must have
    collision degree 1 (G3); the other has collision degree 2
    (2G1+G3).  This removes the old arbitrary quartic-sign assignment.

Run:
  sage -python ~/Downloads/derive_h92_q8_d13_g3_from_e77_bisection.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def rational(field, ring, data, numerator, denominator):
    return (
        field(polynomial(ring, data[numerator]))
        / field(polynomial(ring, data[denominator]))
    )


def monic_power_root(value, exponent):
    ring = value.parent()
    root = ring.one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def rational_map_degree(value, ring):
    value = value.parent()(value)
    return max(
        int(ring(value.numerator()).degree()),
        int(ring(value.denominator()).degree()),
    )


def squarefree_binary_quartic(radicand, old_base_ring):
    numerator = old_base_ring(radicand.numerator())
    denominator = old_base_ring(radicand.denominator())
    nf = numerator.factor()
    df = denominator.factor()
    odd_factors = tuple(
        factor
        for fac in (nf, df)
        for factor, exponent in fac
        if exponent % 2
    )
    quartic = old_base_ring(nf.unit() / df.unit())
    for factor in odd_factors:
        quartic *= factor
    assert quartic
    quotient = radicand / old_base_ring.fraction_field()(quartic)
    assert quotient.is_square()
    square_factor = quotient.sqrt()
    assert radicand == square_factor**2 * quartic
    return quartic, square_factor


def degree_or_minus_one(value):
    return -1 if not value else int(value.degree())


def infinity_order(value, ring):
    if not value:
        return 10**9
    return int(
        ring(value.denominator()).degree()
        - ring(value.numerator()).degree()
    )


def point_degrees(point, ring):
    if point.is_zero():
        return None
    x, y = point.xy()
    return {
        "x": [
            int(ring(x.numerator()).degree()),
            int(ring(x.denominator()).degree()),
        ],
        "y": [
            int(ring(y.numerator()).degree()),
            int(ring(y.denominator()).degree()),
        ],
    }


def collision_data(point, ring):
    if point.is_zero():
        return {"degree": -1, "h": None}
    x, y = point.xy()
    dx = ring(x.denominator())
    dy = ring(y.denominator())
    hx = monic_power_root(dx, 2)
    hy = monic_power_root(dy, 3)
    assert hx == hy
    return {"degree": int(hx.degree()), "h": hx}


def primitive_monic(poly):
    poly = poly.parent()(poly)
    assert poly
    return poly.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--q8-child", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"

CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
OLD_CURVES = LOCAL / "q8-explicit-old-curves.json"
BRANCH = LOCAL / "q8-d13-branch-anchor.json"

if args.q8_child:
    Q8_CHILD = args.q8_child.resolve()
else:
    candidates = [
        GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
        LOCAL / "q8-corrected2cover-qq-child.json",
    ]
    Q8_CHILD = next((p for p in candidates if p.exists()), candidates[0])

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-d13-g3-from-e77-bisection.json"
)

for path in (CHILD, ZERO, COMPONENTS, OLD_CURVES, BRANCH, Q8_CHILD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMPONENTS.read_text())
old_curves = json.loads(OLD_CURVES.read_text())
branch = json.loads(BRANCH.read_text())
q8child = json.loads(Q8_CHILD.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert old_curves["status"] == "PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE"
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert q8child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

e77_lattice = next(
    item for item in old_curves["curves"]
    if item["curve"] == "old_E7_7"
)
assert e77_lattice["q8_degree"] == 2
assert e77_lattice["mw_coordinates"] == [1, 0, 1, 0]
assert e77_lattice["index3_residue"] != 0

# ===========================================================================
# 1. Reconstruct q6 child and the corrected marked q8 section.
# ===========================================================================

R = PolynomialRing(QQ, "T")
T = R.gen()
Kold = R.fraction_field()

model = child["minimal_short_weierstrass"]
A = polynomial(R, model["A_coefficients_low_to_high"])
Bcurve = polynomial(R, model["B_coefficients_low_to_high"])
Eold = EllipticCurve(Kold, [0, 0, 0, Kold(A), Kold(Bcurve)])

zdata = zero["section"]
old_zero = Eold(
    rational(
        Kold, R, zdata,
        "x_numerator_coefficients_low_to_high",
        "x_denominator_coefficients_low_to_high",
    ),
    rational(
        Kold, R, zdata,
        "y_numerator_coefficients_low_to_high",
        "y_denominator_coefficients_low_to_high",
    ),
)

point_entries = {entry["sign"]: entry for entry in components["sections"]}
points = {
    sign: Eold(
        rational(
            Kold, R, entry,
            "x_numerator_coefficients_low_to_high",
            "x_denominator_coefficients_low_to_high",
        ),
        rational(
            Kold, R, entry,
            "y_numerator_coefficients_low_to_high",
            "y_denominator_coefficients_low_to_high",
        ),
    )
    for sign, entry in point_entries.items()
}

affine = points[components["source"]["affine_E7_sign"]]
e7_7 = points[components["source"]["E7_7_sign"]]

Pmap = e7_7 - old_zero
Qmap = e7_7 - affine
S = Pmap + Qmap
sx, sy = S.xy()

nx, dx = R(sx.numerator()), R(sx.denominator())
ny, dy = R(sy.numerator()), R(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

# Reconstruct the exact regular q-frame used by the QQ compiler.
ii = R(next(
    item["factor"] for item in child["finite_fibres"]
    if item["kodaira"] == "II*"
)).monic()
iv = R(next(
    item["factor"] for item in child["finite_fibres"]
    if item["kodaira"] == "IV*"
)).monic()
M = (ii**2 * iv**2).monic()

normalizer = (ny * dx * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
rho = (normalizer * nx.inverse_mod(M)).mod(M)

def pair_from_st(s, t):
    s, t = R(s), R(t)
    Bcoef = Kold(s) / Kold(h)
    Acoef = (
        -Kold(s)*p_fun/Kold(h)
        -Kold(s)*Kold(normalizer)/Kold(nx)
        +Kold(s*rho)
        +Kold(t*M)
    )
    return Acoef, Bcoef

kernel_polys = q8child["rr"]["kernel_polynomials"]
assert len(kernel_polys) == 2
pairs = [
    pair_from_st(R(entry["s"]), R(entry["t"]))
    for entry in kernel_polys
]
(A0, B0), (A1, B1) = pairs

# On old_E7_7 the chord through -S has this exact slope.
x7, y7 = e7_7.xy()
m7 = (y7 + sy) / (x7 - sx)
assert m7

U_on_e77 = (A1 + B1*m7) / (A0 + B0*m7)
assert rational_map_degree(U_on_e77, R) == 2

print(
    "Q8G3_BISECTION|old_curve=old_E7_7|q8_degree=2|"
    f"U_map_degree={rational_map_degree(U_on_e77,R)}|"
    "lattice_AJ=1,0,1,0|status=PASS",
    flush=True,
)

# ===========================================================================
# 2. Reconstruct the quartic over QQ(U), including its square factor.
# ===========================================================================

UR = PolynomialRing(QQ, "U")
U = UR.gen()
K = UR.fraction_field()
TR = PolynomialRing(K, "T")
TT = TR.gen()
KT = TR.fraction_field()

def lift_poly(value):
    value = R(value)
    return TR([K(c) for c in value.list()])

def lift_rat(value):
    value = Kold(value)
    return (
        KT(lift_poly(value.numerator()))
        / KT(lift_poly(value.denominator()))
    )

m_value = -(
    lift_rat(A1) - K(U)*lift_rat(A0)
) / (
    lift_rat(B1) - K(U)*lift_rat(B0)
)

sxU, syU = lift_rat(sx), lift_rat(sy)
AU, BU = lift_poly(A), lift_poly(Bcurve)

XR = PolynomialRing(KT, "x")
x = XR.gen()
y_line = XR(m_value)*(x-XR(sxU)) - XR(syU)
relation = y_line**2 - x**3 - XR(AU)*x - XR(BU)
quadratic, remainder = relation.quo_rem(x-XR(sxU))
assert not remainder and quadratic.degree() == 2

disc = KT(quadratic[1]**2 - 4*quadratic[2]*quadratic[0])
quartic, square_factor = squarefree_binary_quartic(disc, TR)
assert quartic.degree() == 4

# Check against the branch-anchor old-base points.
Tii = QQ(branch["zero"]["old_base_T"])
Tiv = QQ(branch["iv_sections"]["old_base_T"])
assert not K(quartic(K(Tii)))
qiv = K(quartic(K(Tiv)))
assert qiv and qiv.is_square()

# ===========================================================================
# 3. The E7_7 bisection equation H_U(T)=0 and its actual W coordinate.
# ===========================================================================

Ue77_lift = lift_rat(U_on_e77)
u_relation = Ue77_lift - KT(K(U))
H = primitive_monic(TR(u_relation.numerator()))
assert H.degree() == 2
assert H.gcd(TR(u_relation.denominator())) == 1

# Verify that the generic map really has no hidden factor.
assert H[2] == 1

x7U = lift_rat(x7)

# For q2*x^2+q1*x+q0=0, a chosen square root of the discriminant at
# the root x=x7 is 2*q2*x7+q1.  Divide by the compiler's certified
# square_factor to obtain W on W^2=quartic.
sqrt_disc_on_e77 = 2*KT(quadratic[2])*x7U + KT(quadratic[1])
W_e77 = sqrt_disc_on_e77 / KT(square_factor)

# It need only satisfy W^2=quartic modulo H, since x7 is a quadratic root
# precisely on the bisection equation H=0.
def remainder_mod_H(value):
    value = KT(value)
    num = TR(value.numerator())
    den = TR(value.denominator())
    assert den.gcd(H) == 1
    return (num * den.inverse_mod(H)).mod(H)

assert remainder_mod_H(W_e77**2 - KT(quartic)) == 0

print(
    "Q8G3_QUARTIC|H_degree=2|W_recovery=PASS|"
    f"H={H}",
    flush=True,
)

# ===========================================================================
# 4. Add the two conjugate bisection points on the branch-anchored curve.
# ===========================================================================

L = K.extension(H, names="theta")
theta = L.gen()
theta_bar = -L(H[1]) - theta

def eval_poly_L(poly, value):
    poly = TR(poly)
    return sum(L(poly[i]) * value**i for i in range(poly.degree()+1))

def eval_rat_L(value, at):
    value = KT(value)
    num = eval_poly_L(TR(value.numerator()), at)
    den = eval_poly_L(TR(value.denominator()), at)
    assert den
    return num / den

def descend_L(value):
    value = L(value)
    coeffs = list(value)
    if not coeffs:
        return K(0)
    assert all(not coefficient for coefficient in coeffs[1:]), (
        "quadratic trace did not descend", value
    )
    return K(coeffs[0])

assert H(theta) == 0
assert H(theta_bar) == 0

# Shifted quartic coefficients around the II* branch point.
rR = PolynomialRing(K, "r")
r = rR.gen()
shifted = rR(quartic(r + K(Tii)))
assert shifted[0] == 0 and shifted.degree() == 4
d = K(shifted[1])
c = K(shifted[2])
b = K(shifted[3])
a = K(shifted[4])
assert d

Eanchor = EllipticCurve(K, [0, c, 0, b*d, a*d**2])
EanchorL = Eanchor.change_ring(L)

def anchor_point(at):
    W = eval_rat_L(W_e77, at)
    Tval = at
    assert W**2 == eval_poly_L(quartic, Tval)
    delta = Tval - L(Tii)
    assert delta
    X = L(d) / delta
    Y = L(d) * W / delta**2
    P = EanchorL(X, Y)
    assert P in EanchorL
    return P

Ptheta = anchor_point(theta)
Pbar = anchor_point(theta_bar)
PsumL = Ptheta + Pbar
assert not PsumL.is_zero()

xsumL, ysumL = PsumL.xy()
xsum = descend_L(xsumL)
ysum = descend_L(ysumL)
Paj_anchor = Eanchor(xsum, ysum)
assert Paj_anchor in Eanchor

print(
    "Q8G3_TRACE|quadratic_points=2|sum_descends=1|"
    f"anchor_x_deg={point_degrees(Paj_anchor,UR)['x'][0]}/"
    f"{point_degrees(Paj_anchor,UR)['x'][1]}|"
    f"anchor_y_deg={point_degrees(Paj_anchor,UR)['y'][0]}/"
    f"{point_degrees(Paj_anchor,UR)['y'][1]}|status=PASS",
    flush=True,
)

# ===========================================================================
# 5. Move to the certified canonical D13 model.
# ===========================================================================

child_data = q8child["child"]
Amin = UR([QQ(v) for v in child_data["minimal_A_coefficients_low_to_high"]])
Bmin = UR([QQ(v) for v in child_data["minimal_B_coefficients_low_to_high"]])
Ecanon = EllipticCurve(K, [0, 0, 0, K(Amin), K(Bmin)])

assert Eanchor.is_isomorphic(Ecanon)
iso = Eanchor.isomorphism_to(Ecanon)
Paj = iso(Paj_anchor)
assert Paj in Ecanon

# The two old IV* fibre components are the two quartic points at T=Tiv.
Wiv = qiv.sqrt()
delta_iv = K(Tiv - Tii)
Xiv = d / delta_iv
Yiv = d * Wiv / delta_iv**2
Piv_plus = iso(Eanchor(Xiv, Yiv))
Piv_minus = iso(Eanchor(Xiv, -Yiv))
assert Piv_minus == -Piv_plus

candidate_plus = Paj - Piv_plus
candidate_minus = Paj - Piv_minus
cp = collision_data(candidate_plus, UR)
cm = collision_data(candidate_minus, UR)

degrees = sorted((cp["degree"], cm["degree"]))
assert degrees == [1, 2], (
    "expected G3 and 2G1+G3 collision degrees 1 and 2",
    cp["degree"], cm["degree"],
)

if cp["degree"] == 1:
    G1 = Piv_plus
    minus_G1 = Piv_minus
    G3 = candidate_plus
    other = candidate_minus
    e61_sign = "plus"
else:
    G1 = Piv_minus
    minus_G1 = Piv_plus
    G3 = candidate_minus
    other = candidate_plus
    e61_sign = "minus"

assert G1 + G3 == Paj
assert minus_G1 == -G1
assert collision_data(G1, UR)["degree"] == 0
assert collision_data(G3, UR)["degree"] == 1
assert collision_data(other, UR)["degree"] == 2

# Lattice-side identity and height/correction regression.
Hlat = matrix(QQ, [
    [QQ(v) for v in row]
    for row in old_curves["anchor"]["height_gram"]
])
zG1 = vector(ZZ, (1,0,0,0))
zAJ = vector(ZZ, e77_lattice["mw_coordinates"])
zG3 = zAJ - zG1
assert zG3 == vector(ZZ, (0,0,1,0))
assert zG3 * Hlat * zG3 == QQ(11)/4
assert QQ(4 + 2*1) - QQ(13)/4 == QQ(11)/4

g1deg = point_degrees(G1, UR)
g3deg = point_degrees(G3, UR)
ajdeg = point_degrees(Paj, UR)

print(
    "Q8G3_CANONICAL|"
    f"E6_1_quartic_sign={e61_sign}|"
    "AJ_mw=1,0,1,0|G1_mw=1,0,0,0|G3_mw=0,0,1,0|"
    "G3_height=11/4|G3_PdotO=1|G3_correction=13/4|"
    f"G1_xdeg={g1deg['x'][0]}/{g1deg['x'][1]}|"
    f"G3_xdeg={g3deg['x'][0]}/{g3deg['x'][1]}|"
    f"G3_ydeg={g3deg['y'][0]}/{g3deg['y'][1]}|"
    f"AJ_xdeg={ajdeg['x'][0]}/{ajdeg['x'][1]}|"
    "status=PASS",
    flush=True,
)

def point_payload(P):
    if P.is_zero():
        return {"zero": True}
    x, y = P.xy()
    return {
        "zero": False,
        "x": str(x),
        "y": str(y),
        "degrees": point_degrees(P, UR),
        "collision": {
            "degree": collision_data(P, UR)["degree"],
            "h": str(collision_data(P, UR)["h"]),
        },
    }

payload = {
    "schema": "elkies-k3.h92-q8-d13-g3-from-e77-bisection.v1",
    "status": "PASS_EXACT_D13_G3_FROM_E77_BISECTION",
    "source_bisection": {
        "curve": "old_E7_7",
        "q8_degree": 2,
        "anchored_D13_AJ_mw": [1,0,1,0],
        "U_on_curve": str(U_on_e77),
        "quadratic_H_over_QQ_U": str(H),
        "quartic_W_on_curve_mod_H": str(W_e77),
        "trace_descends_to_QQ_U": True,
    },
    "orientation": {
        "E6_1_quartic_sign": e61_sign,
        "E6_5_quartic_sign": "minus" if e61_sign == "plus" else "plus",
        "criterion": (
            "subtracting E6_1 from AJ(old_E7_7) gives collision degree 1; "
            "subtracting E6_5 gives collision degree 2"
        ),
    },
    "canonical_D13": {
        "AJ_old_E7_7": point_payload(Paj),
        "G1": point_payload(G1),
        "minus_G1": point_payload(minus_G1),
        "G3": point_payload(G3),
        "other_subtraction_2G1_plus_G3": point_payload(other),
    },
    "lattice_certificate": {
        "AJ_old_E7_7": [1,0,1,0],
        "G1": [1,0,0,0],
        "G3": [0,0,1,0],
        "G3_height": "11/4",
        "G3_P_dot_O": 1,
        "G3_D13_local_correction": "13/4",
    },
    "boundary": (
        "This gives explicit canonical-D13 rational coordinates for G3 and "
        "fixes the IV* quartic sign orientation. It does not yet construct "
        "G2, G4, the q24 marked point, or the D12 neighbour equation."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print("Q8G3_RESULT|status=PASS_EXACT_D13_G3_FROM_E77_BISECTION")
print(f"OUTPUT|{OUTPUT}")
