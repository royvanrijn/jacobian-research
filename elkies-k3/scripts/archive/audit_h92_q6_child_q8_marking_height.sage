#!/usr/bin/env sage -python
"""
Audit the H92 q6-child q8 marking by comparing the actual rational section
against the MW height claimed for it.

The q8 marking constructs an exact rational section S on the globally minimal
E8+E6 elliptic K3 and assigns it MW coordinates (-2,-2,0) in the certified
height lattice

    [[8/3,1/3,-1],
     [1/3,8/3,1],
     [-1,1,46]].

Independently, the exact rational functions for S have x-denominator h^2 and
y-denominator h^3, with h squarefree of degree 46 and coprime to the
Weierstrass discriminant.  Thus S meets the standard zero section transversely
at 46 smooth fibres.

On an elliptic K3, Shioda's formula is

    <S,S> = 4 + 2(S.O) - sum_v contr_v(S).

The only reducible fibres are E8 and E6.  Whatever components S meets, the
local correction is at most the largest diagonal entry of the inverse Cartan
matrix at each fibre.  For the pinned E8 and E6 Cartan matrices these maxima
are 30 and 6.  Therefore the actual rational section must have height at least

    4 + 2*46 - (30+6) = 60.

If the claimed coordinate has height 24, the point-to-NS/MW identification is
inconsistent and all q8 child constructions that consume it must be treated as
untrusted until that bridge is repaired.

Run:
  sage -python ~/Downloads/audit_h92_q6_child_q8_marking_height.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, matrix, vector
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


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
MARKING = GEN / "elkies-k3-h92-q6-child-q8-marking.json"

for path in (CHILD, MARKING):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()

A = polynomial(R, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(R, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Delta = polynomial(R, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
assert A.degree() <= 8
assert B.degree() <= 12
assert Delta.degree() <= 24
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

sdata = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = K(polynomial(R, sdata["x_numerator_coefficients_low_to_high"])) / K(
    polynomial(R, sdata["x_denominator_coefficients_low_to_high"])
)
sy = K(polynomial(R, sdata["y_numerator_coefficients_low_to_high"])) / K(
    polynomial(R, sdata["y_denominator_coefficients_low_to_high"])
)
S = E(sx, sy)
assert S in E

dx = R(sx.denominator())
dy = R(sy.denominator())
h_x = monic_power_root(dx, 2)
h_y = monic_power_root(dy, 3)
assert h_x == h_y
h = h_x
assert h.degree() == 46
assert h.gcd(Delta) in QQ

# Check these really are O-intersections, not merely coordinate denominator
# artifacts: z=-x/y is the formal parameter at O and must have exactly one
# zero along each root of h.
z = -sx/sy
z_over_h = z/K(h)
assert R(z_over_h.numerator()).gcd(h) in QQ
assert R(z_over_h.denominator()).gcd(h) in QQ
smooth_O_intersection_degree = int(h.degree())

height_gram = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])
coords = vector(QQ, marking["selected_q8"]["relative_child_section_MW_coordinates"])
assert tuple(coords) == (-2, -2, 0)
claimed_height = QQ(coords * height_gram * coords)

E8 = matrix(QQ, [
    [2,0,-1,0,0,0,0,0],
    [0,2,0,-1,0,0,0,0],
    [-1,0,2,-1,0,0,0,0],
    [0,-1,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,0],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,-1],
    [0,0,0,0,0,0,-1,2],
])
E6 = matrix(QQ, [
    [2,-1,0,0,0,0],
    [-1,2,-1,0,0,0],
    [0,-1,2,-1,0,-1],
    [0,0,-1,2,-1,0],
    [0,0,0,-1,2,0],
    [0,0,-1,0,0,2],
])

max_E8 = max(E8.inverse().diagonal())
max_E6 = max(E6.inverse().diagonal())
max_total_correction = max_E8 + max_E6

# χ(O_K3)=2, so Shioda self-height begins at 2χ=4.
height_lower_bound = QQ(4 + 2*smooth_O_intersection_degree) - max_total_correction
required_correction_for_claim = QQ(4 + 2*smooth_O_intersection_degree) - claimed_height

print(
    "Q8MARKHEIGHT|"
    f"minimal_A_degree={A.degree()}|minimal_B_degree={B.degree()}|"
    f"delta_degree={Delta.degree()}|smooth_O_intersections={smooth_O_intersection_degree}|"
    f"claimed_coords={','.join(map(str,coords))}|claimed_height={claimed_height}",
    flush=True,
)
print(
    "Q8MARKHEIGHT_BOUND|"
    f"max_E8_correction={max_E8}|max_E6_correction={max_E6}|"
    f"max_total_correction={max_total_correction}|"
    f"required_correction_for_claim={required_correction_for_claim}|"
    f"height_lower_bound={height_lower_bound}",
    flush=True,
)

contradiction = claimed_height < height_lower_bound
print(
    "Q8MARKHEIGHT_RESULT|"
    f"contradiction={int(contradiction)}|"
    f"claimed_height={claimed_height}|lower_bound={height_lower_bound}|"
    f"gap={height_lower_bound-claimed_height}|"
    f"status={'FAIL_POINT_TO_MW_BRIDGE' if contradiction else 'PASS_NO_CONTRADICTION'}",
    flush=True,
)

# Optional direct Sage canonical-height diagnostic when the function-field
# backend exposes it.  This is not used in the proof above.
try:
    if hasattr(S, "height"):
        direct = S.height()
        print(f"Q8MARKHEIGHT_SAGE|height={direct}", flush=True)
except Exception as exc:
    print(f"Q8MARKHEIGHT_SAGE|unavailable={type(exc).__name__}:{exc}", flush=True)
