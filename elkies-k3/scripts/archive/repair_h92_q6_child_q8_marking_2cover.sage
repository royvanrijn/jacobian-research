#!/usr/bin/env sage -python
"""
Construct and certify the corrected H92 q6-child q8 marked section after
recognizing the binary-quartic covariant map as a 2-covering map.

The old marking code formed, in the standard child Jacobian,

    Pmap = phi(E7_7) - phi(old_O)
    Qmap = phi(E7_7) - phi(affine_E7)
    S_old = 2*Pmap + 2*Qmap.

But phi is the binary-quartic 2-covering map.  Hence these differences already
represent twice the corresponding primitive MW directions.  The lattice q8
target requires coordinates (-2,-2,0), so the correct standard-Jacobian point is

    S_new = Pmap + Qmap,

and the old point must satisfy S_old=2*S_new.

This script certifies from the actual minimal E8+E6 K3 equation that S_new:
  * is rational and exactly halves the old marked point;
  * has O-intersection 10;
  * is smooth/identity-component at II* and IV*;
  * has Shioda height 24, matching the pinned MW vector (-2,-2,0);
  * has a reduced smooth O-collision divisor of degree 10.

It writes a corrected marking payload under artifacts/local only.  It does not
overwrite the canonical generated marking.

Run:
  sage -python ~/Downloads/repair_h92_q6_child_q8_marking_2cover.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix, vector


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


def point_data(point, ring):
    x, y = point.xy()
    return {
        "x_numerator_coefficients_low_to_high":
            [str(v) for v in ring(x.numerator()).list()],
        "x_denominator_coefficients_low_to_high":
            [str(v) for v in ring(x.denominator()).list()],
        "y_numerator_coefficients_low_to_high":
            [str(v) for v in ring(y.numerator()).list()],
        "y_denominator_coefficients_low_to_high":
            [str(v) for v in ring(y.denominator()).list()],
        "coordinate_degrees": {
            "x": [int(ring(x.numerator()).degree()),
                  int(ring(x.denominator()).degree())],
            "y": [int(ring(y.numerator()).degree()),
                  int(ring(y.denominator()).degree())],
        },
    }


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0, (factor, multiplicity, exponent)
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def ord_at(value, factor):
    return int(value.numerator().valuation(factor)
               - value.denominator().valuation(factor))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument(
    "--output",
    type=Path,
    help="default: artifacts/local/elkies-k3/q8-marking-2cover-corrected.json",
)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
OUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts" / "local" / "elkies-k3" /
         "q8-marking-2cover-corrected.json"
)

CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
OLD_MARKING = GEN / "elkies-k3-h92-q6-child-q8-marking.json"

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMPONENTS.read_text())
old_marking = json.loads(OLD_MARKING.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert old_marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()

A = polynomial(R, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(R, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Delta = polynomial(
    R, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"]
)
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

points = {}
for entry in components["sections"]:
    points[entry["sign"]] = E(
        rational(K, R, entry,
                 "x_numerator_coefficients_low_to_high",
                 "x_denominator_coefficients_low_to_high"),
        rational(K, R, entry,
                 "y_numerator_coefficients_low_to_high",
                 "y_denominator_coefficients_low_to_high"),
    )
affine = points[components["source"]["affine_E7_sign"]]
e77 = points[components["source"]["E7_7_sign"]]

Pmap = e77 - P0
Qmap = e77 - affine
S_new = Pmap + Qmap
S_old = 2*Pmap + 2*Qmap

old_sdata = old_marking["selected_q8"][
    "relative_child_section_standard_jacobian_coordinates"
]
S_old_artifact = E(
    rational(K, R, old_sdata,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(K, R, old_sdata,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
assert S_old == S_old_artifact
assert S_old == 2*S_new
assert not S_new.is_zero()

# Collision divisor and exact O intersection.
sx, sy = S_new.xy()
hx = monic_power_root(R(sx.denominator()), 2)
hy = monic_power_root(R(sy.denominator()), 3)
assert hx == hy
h = hx
assert h.degree() == 10
assert h.gcd(Delta) in QQ

z = -sx/sy
z_over_h = z/K(h)
assert R(z_over_h.numerator()).gcd(h) in QQ
assert R(z_over_h.denominator()).gcd(h) in QQ

# No hidden O-intersection at infinity in the global minimal K3 scaling.
xdiff = int(R(sx.numerator()).degree() - R(sx.denominator()).degree())
ydiff = int(R(sy.numerator()).degree() - R(sy.denominator()).degree())
assert xdiff == 4 and ydiff == 6
O_intersection = int(h.degree())

fibres = {
    item["kodaira"]: R(item["factor"]).monic()
    for item in child["finite_fibres"]
    if item["kodaira"] in ("II*", "IV*")
}
assert set(fibres) == {"II*", "IV*"}

local = {}
for kind, factor in fibres.items():
    ox, oy = ord_at(sx, factor), ord_at(sy, factor)
    oz = ord_at(-sx/sy, factor)
    point = -factor[0]/factor[1]
    assert ox >= 0 and oy >= 0
    x0, y0 = QQ(sx(point)), QQ(sy(point))
    assert (x0, y0) != (0, 0)
    local[kind] = {
        "orders_x_y_z": [ox, oy, oz],
        "specialization": [str(x0), str(y0)],
        "component": "identity",
    }

# Shioda height: χ(K3)=2.  Both additive corrections vanish because S_new
# meets the identity component at II* and IV*.
height = QQ(4 + 2*O_intersection)
assert height == 24

H = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])
claimed = vector(QQ, old_marking["selected_q8"][
    "relative_child_section_MW_coordinates"
])
assert tuple(claimed) == (-2, -2, 0)
claimed_height = QQ(claimed*H*claimed)
assert claimed_height == height

# Old wrong section regression.
osx, osy = S_old.xy()
old_h = monic_power_root(R(osx.denominator()), 2)
assert old_h == monic_power_root(R(osy.denominator()), 3)
assert old_h.degree() == 46
old_height = QQ(4 + 2*old_h.degree())  # old section smooth at both additive fibres
assert old_height == 96 == 4*height

print(
    "Q8MARK2COVER|"
    f"old_equals_2new=1|old_collision_degree={old_h.degree()}|"
    f"new_collision_degree={h.degree()}|old_height={old_height}|"
    f"new_height={height}|claimed_height={claimed_height}",
    flush=True,
)
print(
    "Q8MARK2COVER_LOCAL|"
    f"II={local['II*']['component']}:{','.join(map(str,local['II*']['orders_x_y_z']))}|"
    f"IV={local['IV*']['component']}:{','.join(map(str,local['IV*']['orders_x_y_z']))}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-marking-2cover-corrected.v1",
    "status": "PASS_EXACT_Q8_MARKING_2COVER_CORRECTION",
    "correction": {
        "binary_quartic_covariant_role":
            "canonical degree-2 covering map to the Jacobian",
        "old_formula":
            "2*(phi(E7_7)-phi(old_O)) + 2*(phi(E7_7)-phi(affine_E7))",
        "correct_formula":
            "(phi(E7_7)-phi(old_O)) + (phi(E7_7)-phi(affine_E7))",
        "exact_group_check": "old_marked_point = 2 * corrected_marked_point",
    },
    "selected_q8": {
        "relative_child_section_MW_coordinates":
            [int(v) for v in claimed],
        "relative_child_section_standard_jacobian_coordinates":
            point_data(S_new, R),
        "height": str(height),
        "O_intersection": O_intersection,
        "collision_divisor": {
            "degree": int(h.degree()),
            "coefficients_low_to_high": [str(v) for v in h.list()],
            "squarefree": bool(h.gcd(h.derivative()) in QQ),
            "coprime_to_discriminant": bool(h.gcd(Delta) in QQ),
        },
        "additive_fibres": local,
    },
    "old_marking_regression": {
        "old_section_is_double": True,
        "old_height": str(old_height),
        "old_collision_degree": int(old_h.degree()),
    },
    "boundary": (
        "This repairs the rational point attached to the already-certified "
        "MW coordinate (-2,-2,0). It does not automatically reauthorize "
        "downstream q8 smooth/additive/global modules, all of which must be "
        "regenerated against the corrected degree-10 collision divisor."
    ),
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("Q8MARK2COVER_RESULT|status=PASS_EXACT_Q8_MARKING_2COVER_CORRECTION")
print(f"OUTPUT|{OUT}")
