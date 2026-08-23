#!/usr/bin/env sage -python
"""
Audit the exact rational q6-child points used by the H3 q8 marking.

This stays entirely on the certified globally minimal E8+E6 elliptic K3.
It reconstructs

    P = E7_7 - old_zero
    Q = E7_7 - affine_E7
    R = P + Q
    S = 2P + 2Q

in the STANDARD Jacobian group law, then computes their geometric Shioda
heights from:

  * exact intersection multiplicity with the standard zero section O;
  * the two reducible fibres II* (E8) and IV* (E6).

For the globally minimal K3 model:
    <P,P> = 4 + 2(P.O) - sum contr_v(P).

At II* the component group is trivial and the unique multiplicity-one
component is the identity component, so the section correction is 0.

At IV* a section that specializes to a smooth point of the Weierstrass cubic
(including O) meets the identity component and contributes 0.  A section that
runs into the cusp (0,0) meets one of the two nonidentity multiplicity-one
components; both have E6 correction 4/3.

The script also checks the quadratic-height identity
    h(2R) = 4 h(R)
directly from the independently computed geometric data.  If this passes, the
point heights are internally certified without using the disputed MW bridge.

Finally it compares the actual Gram of (P,Q) to the claimed pinned MW
coordinates (-1,0,0),(0,-1,0) and searches the pinned rank-3 MW lattice for
integer vector pairs having the observed Gram.

Run:
  sage -python ~/Downloads/audit_h92_q6_child_q8_point_heights.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, pari, vector
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
    return ring([QQ(v) for v in values])


def rational(field, ring, data, nk, dk):
    return field(polynomial(ring, data[nk])) / field(polynomial(ring, data[dk]))


def monic_power_root(value, exponent):
    value = value.parent()(value)
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        if multiplicity % exponent:
            raise AssertionError(
                f"denominator is not an exact {exponent}-th power: "
                f"factor degree {factor.degree()}, multiplicity {multiplicity}"
            )
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def ord_at(value, factor):
    return int(
        value.numerator().valuation(factor)
        - value.denominator().valuation(factor)
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--search-cap", type=int, default=300)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"

CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
MARKING = GEN / "elkies-k3-h92-q6-child-q8-marking.json"

for path in (CHILD, ZERO, COMPONENTS, MARKING):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMPONENTS.read_text())
marking = json.loads(MARKING.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()

A = polynomial(R, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(R, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Delta = polynomial(R, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
assert A.degree() <= 8 and B.degree() <= 12 and Delta.degree() <= 24

E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

zdata = zero["section"]
P0 = E(
    rational(K, R, zdata,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(K, R, zdata,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)

point_by_sign = {}
for entry in components["sections"]:
    point_by_sign[entry["sign"]] = E(
        rational(K, R, entry,
                 "x_numerator_coefficients_low_to_high",
                 "x_denominator_coefficients_low_to_high"),
        rational(K, R, entry,
                 "y_numerator_coefficients_low_to_high",
                 "y_denominator_coefficients_low_to_high"),
    )

affine = point_by_sign[components["source"]["affine_E7_sign"]]
e77 = point_by_sign[components["source"]["E7_7_sign"]]

P = e77 - P0
Q = e77 - affine
Rsum = P + Q
S = 2 * Rsum

sdata = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
S_marked = E(
    rational(K, R, sdata,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(K, R, sdata,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
assert S == S_marked

fibres = {
    item["kodaira"]: R(item["factor"]).monic()
    for item in child["finite_fibres"]
    if item["kodaira"] in ("II*", "IV*")
}
assert set(fibres) == {"II*", "IV*"}
ii = fibres["II*"]
iv = fibres["IV*"]
assert ii.degree() == iv.degree() == 1

# Standard E6 section correction for either nonidentity component of IV*.
IVSTAR_NONIDENTITY_CORRECTION = QQ(4) / 3


def zero_intersection(point):
    """Return exact P.O from the reduced rational coordinates."""
    x, y = point.xy()
    dx = R(x.denominator())
    dy = R(y.denominator())
    hx = monic_power_root(dx, 2)
    hy = monic_power_root(dy, 3)
    assert hx == hy, "x/y denominator roots disagree"

    # Every finite denominator root is an O-intersection.  Verify with the
    # formal parameter z=-x/y: z/hroot is a unit at the finite divisor.
    z = -x / y
    z_over_h = z / K(hx)
    assert R(z_over_h.numerator()).gcd(hx) in QQ
    assert R(z_over_h.denominator()).gcd(hx) in QQ
    finite_degree = int(hx.degree())

    # In the minimal K3 scaling x=s^-4 X, y=s^-6 Y at infinity.
    xdiff = int(R(x.numerator()).degree() - R(x.denominator()).degree())
    ydiff = int(R(y.numerator()).degree() - R(y.denominator()).degree())
    x_extra = xdiff - 4
    y_extra = ydiff - 6
    assert x_extra % 2 == 0 and y_extra % 3 == 0
    ninf_x = max(0, x_extra // 2)
    ninf_y = max(0, y_extra // 3)
    assert ninf_x == ninf_y, (xdiff, ydiff, ninf_x, ninf_y)

    return finite_degree + ninf_x, {
        "finite_degree": finite_degree,
        "infinity_multiplicity": ninf_x,
        "x_degree_difference": xdiff,
        "y_degree_difference": ydiff,
        "finite_denominator_root_degree": finite_degree,
        "finite_denominator_gcd_delta_degree": int(hx.gcd(Delta).degree()),
    }


def fibre_status(point, kind, factor):
    x, y = point.xy()
    ox, oy = ord_at(x, factor), ord_at(y, factor)
    oz = ord_at(-x/y, factor)

    # If z has positive order, the section specializes to O.
    if oz > 0 and (ox < 0 or oy < 0):
        return "O", ox, oy, oz

    # If x,y are regular, inspect the affine specialization.
    if ox >= 0 and oy >= 0:
        base_point = -factor[0] / factor[1]
        x0, y0 = QQ(x(base_point)), QQ(y(base_point))
        if (x0, y0) == (0, 0):
            return "cusp", ox, oy, oz
        return "smooth_affine", ox, oy, oz

    return "unclassified", ox, oy, oz


def geometric_height(name, point):
    Oint, odetail = zero_intersection(point)
    ii_status, iix, iiy, iiz = fibre_status(point, "II*", ii)
    iv_status, ivx, ivy, ivz = fibre_status(point, "IV*", iv)

    # II*: E8 component group is trivial; section correction is zero.
    ii_contr = QQ(0)

    # IV*: standard zero / smooth affine points are identity component.
    # A cusp specialization reaches one of the two nonidentity mult.-1
    # components, both with correction 4/3.
    if iv_status == "cusp":
        iv_contr = IVSTAR_NONIDENTITY_CORRECTION
    elif iv_status in ("O", "smooth_affine"):
        iv_contr = QQ(0)
    else:
        raise RuntimeError(
            f"{name}: could not classify IV* specialization "
            f"(orders x={ivx}, y={ivy}, z={ivz})"
        )

    height = QQ(4 + 2*Oint) - ii_contr - iv_contr
    record = {
        "name": name,
        "O_intersection": Oint,
        "II_status": ii_status,
        "II_orders": (iix, iiy, iiz),
        "IV_status": iv_status,
        "IV_orders": (ivx, ivy, ivz),
        "II_correction": ii_contr,
        "IV_correction": iv_contr,
        "height": height,
        **odetail,
    }
    print(
        "Q8POINTHEIGHT|"
        f"name={name}|O={Oint}|"
        f"II={ii_status}|IIord={iix},{iiy},{iiz}|"
        f"IV={iv_status}|IVord={ivx},{ivy},{ivz}|"
        f"contr={ii_contr + iv_contr}|height={height}|"
        f"denroot={odetail['finite_denominator_root_degree']}|"
        f"denroot_gcd_delta={odetail['finite_denominator_gcd_delta_degree']}|"
        f"infO={odetail['infinity_multiplicity']}",
        flush=True,
    )
    return height, record


hP, recP = geometric_height("P=E7_7-old_zero", P)
hQ, recQ = geometric_height("Q=E7_7-affine", Q)
hR, recR = geometric_height("R=P+Q", Rsum)
hS, recS = geometric_height("S=2P+2Q", S)

quadratic_ok = (hS == 4*hR)
pair_actual = (hR - hP - hQ) / 2

H = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])
claimed_P = vector(QQ, (-1, 0, 0))
claimed_Q = vector(QQ, (0, -1, 0))
claimed_R = claimed_P + claimed_Q
claimed_S = 2*claimed_R

def h_of(v):
    return QQ(v * H * v)

claimed_pair = QQ(claimed_P * H * claimed_Q)

print(
    "Q8POINTGRAM|"
    f"hP={hP}|hQ={hQ}|pairPQ={pair_actual}|hR={hR}|hS={hS}|"
    f"quad_hS_eq_4hR={int(quadratic_ok)}",
    flush=True,
)
print(
    "Q8POINTCLAIM|"
    f"hP={h_of(claimed_P)}|hQ={h_of(claimed_Q)}|"
    f"pairPQ={claimed_pair}|hR={h_of(claimed_R)}|hS={h_of(claimed_S)}|"
    f"matches_actual={int((hP,hQ,pair_actual)==(h_of(claimed_P),h_of(claimed_Q),claimed_pair))}",
    flush=True,
)

# Search the pinned MW lattice for exact integer vectors of the observed
# lengths, then pairs with the observed inner product.  Work with G=3H to
# make the form integral.
G = (3*H).change_ring(ZZ)
targets = {}
for label, height in (("P", hP), ("Q", hQ)):
    target = ZZ(3*height)
    result = pari(G).qfminim(target)
    cols = matrix(ZZ, result[2]).columns()
    vectors = set()
    for col in cols:
        for sign in (1, -1):
            v = sign*vector(ZZ, col)
            if ZZ(v*G*v) == target:
                vectors.add(tuple(map(int, v)))
    targets[label] = sorted(vectors)
    print(
        f"Q8POINTSEARCH|name={label}|target3h={target}|"
        f"exact_vectors={len(vectors)}|sample={targets[label][:12]}",
        flush=True,
    )

pairs = []
target_pair3 = ZZ(3*pair_actual)
for pv in targets["P"]:
    vp = vector(ZZ, pv)
    for qv in targets["Q"]:
        vq = vector(ZZ, qv)
        if ZZ(vp*G*vq) == target_pair3:
            pairs.append((pv, qv))

print(
    "Q8POINTSEARCHPAIR|"
    f"target3pair={target_pair3}|pairs={len(pairs)}|sample={pairs[:20]}",
    flush=True,
)

status = (
    "PASS_ACTUAL_GRAM_MATCHES_CLAIM"
    if (hP, hQ, pair_actual) ==
       (h_of(claimed_P), h_of(claimed_Q), claimed_pair)
    else "FAIL_PRIMITIVE_POINT_TO_MW_BRIDGE"
)

print(
    "Q8POINTAUDIT|"
    f"quadratic_consistency={int(quadratic_ok)}|"
    f"status={status}",
    flush=True,
)

out = ROOT / "artifacts" / "local" / "elkies-k3" / "q8-point-mw-height-audit.json"
out.parent.mkdir(parents=True, exist_ok=True)
def jsonable(value):
    if hasattr(value, "parent") and value.parent() is QQ:
        return str(value)
    if isinstance(value, tuple):
        return [jsonable(v) for v in value]
    if isinstance(value, list):
        return [jsonable(v) for v in value]
    return value


out.write_text(json.dumps({
    "schema": "elkies-k3.h92-q6-child-q8-point-mw-height-audit.v1",
    "status": status,
    "points": {
        "P": {k: jsonable(v) for k,v in recP.items()},
        "Q": {k: jsonable(v) for k,v in recQ.items()},
        "R": {k: jsonable(v) for k,v in recR.items()},
        "S": {k: jsonable(v) for k,v in recS.items()},
    },
    "actual_gram": {
        "P_height": str(hP),
        "Q_height": str(hQ),
        "P_Q_pairing": str(pair_actual),
        "R_height": str(hR),
        "S_height": str(hS),
        "quadratic_S_equals_4R": quadratic_ok,
    },
    "claimed_gram": {
        "P_coordinates": [-1,0,0],
        "Q_coordinates": [0,-1,0],
        "P_height": str(h_of(claimed_P)),
        "Q_height": str(h_of(claimed_Q)),
        "P_Q_pairing": str(claimed_pair),
        "R_height": str(h_of(claimed_R)),
        "S_height": str(h_of(claimed_S)),
    },
    "pinned_lattice_search": {
        "P_vectors": [list(v) for v in targets["P"]],
        "Q_vectors": [list(v) for v in targets["Q"]],
        "matching_pairs": [[list(a),list(b)] for a,b in pairs],
    },
}, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{out}", flush=True)
