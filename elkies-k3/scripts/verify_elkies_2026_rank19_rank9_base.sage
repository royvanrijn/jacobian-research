#!/usr/bin/env sage-python
"""Verify the simplest discovered rank-at-least-nine paired base.

The pair consists of equation-level bisection orbit masks 42110 and 43109.
Both quadratic leading coefficients are squares, so their V4 fibre product is
a genus-one curve with a rational point at infinity.  The third quotient is
2-isogenous to that paired base.  This verifier recomputes its quartic
Jacobian and replays an exact nine-point finite-quotient independence
certificate.  It also parameterizes the first conic, constructs a pointed
quartic birational to the paired base, identifies its minimal Weierstrass
model, and transports the nine independent points through the degree-two
isogeny and the exact inverse birational map to nine explicit rational
``t``-values.  Finally it imports the complete paired-cover height theorem to
obtain generic surface rank at least 19.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOGUE = ROOT / "artifacts/generated-results/elkies-2026-immediate-point-pair-catalogue-full.json"
DEFAULT_LEDGER = ROOT / "artifacts/generated-results/elkies-2026-immediate-point-pair-rank-ledger.json"
DEFAULT_GEOMETRY = ROOT / "artifacts/generated-results/elkies-2026-bisection-pair-cover-geometry-full.json"
DEFAULT_DISJOINT_PAIRS = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority-disjoint-pairs-full.tsv"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank19-rank9-paired-base.json"
FINITE_QUOTIENT_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
SHORT_MODEL_HELPER = ROOT / "elliptic-curves/ecsearch/q12o5867_specialization.py"
PAIR_KEY = "42110:43109"
EXPECTED_MODEL = (
    1,
    0,
    0,
    -70087047578007713577216,
    3865770423647395544516350651140096,
)
EXPECTED_PAIRED_MODEL = (
    1,
    0,
    0,
    -60729194722297004073216,
    5758259762216167074332597509226496,
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def point_text(point) -> list[str]:
    if point.is_zero():
        return ["0"]
    return [rational_text(point[0]), rational_text(point[1])]


def polynomial_text(polynomial) -> str:
    return str(polynomial)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
parser.add_argument("--disjoint-pairs", type=Path, default=DEFAULT_DISJOINT_PAIRS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

catalogue = json.loads(args.catalogue.read_text())
ledger = json.loads(args.ledger.read_text())
geometry = json.loads(args.geometry.read_text())
if catalogue["status"] != "PASS_COMPLETE_5566_IMMEDIATE_POINT_PAIR_CATALOGUE":
    raise ArithmeticError("the complete immediate-point catalogue is not certified")
if ledger["status"] != "PASS_COMPLETE_EXACT_RANK_LOWER_BOUND_LEDGER":
    raise ArithmeticError("the complete rank lower-bound ledger is not certified")
if geometry["status"] != "PASS_COMPLETE_CONIC_AND_GENUS_ONE_PAIR_CLASSIFICATION":
    raise ArithmeticError("the complete pair geometry theorem is not certified")

row = next(record for record in catalogue["pairs"] if record["pair_key"] == PAIR_KEY)
rank_record = ledger["results"][PAIR_KEY]
if row["arithmetic_complexity_rank"] != 114:
    raise ArithmeticError("the selected pair's arithmetic complexity rank changed")
if rank_record["certified_rank_lower_bound"] != 9:
    raise ArithmeticError("the selected pair no longer has nine certified points")
if [cover["lattice_orbit_mask"] for cover in row["covers"]] != [42110, 43109]:
    raise ArithmeticError("the selected cover attachment changed")

R = PolynomialRing(QQ, "t")
t = R.gen()
q1, q2 = (
    R([QQ(value) for value in cover["q_coefficients_low_to_high"]])
    for cover in row["covers"]
)
if not q1.is_irreducible() or not q2.is_irreducible() or q1.gcd(q2) != 1:
    raise ArithmeticError("the pair does not have four distinct geometric branch points")
if not q1[2].is_square() or not q2[2].is_square():
    raise ArithmeticError("the displayed rational point at infinity disappeared")
infinity_coordinates = (QQ(q1[2]).sqrt(), QQ(q2[2]).sqrt())
if row["common_points"]["infinity"] != {
    "u_over_t": rational_text(infinity_coordinates[0]),
    "v_over_t": rational_text(infinity_coordinates[1]),
}:
    raise ArithmeticError("the stored point at infinity changed")

quartic = q1 * q2
e, d, c, b, a = (quartic[index] for index in range(5))
invariant_i = 12 * a * e - 3 * b * d + c**2
invariant_j = (
    72 * a * c * e
    + 9 * b * c * d
    - 27 * a * d**2
    - 27 * b**2 * e
    - 2 * c**3
)
raw_jacobian = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
minimal = raw_jacobian.global_minimal_model()
minimal_model = tuple(int(value) for value in minimal.a_invariants())
if minimal_model != EXPECTED_MODEL or list(minimal_model) != row["minimal_jacobian_a1_a2_a3_a4_a6"]:
    raise ArithmeticError("the minimal Jacobian changed")
if int(minimal.conductor()) != int(row["conductor"]) or int(minimal.root_number()) != -1:
    raise ArithmeticError("the conductor or global root number changed")

points = tuple(
    minimal(QQ(point[0]), QQ(point[1])) for point in rank_record["generators"]
)
if len(points) != 9:
    raise ArithmeticError("expected nine displayed Jacobian points")

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    source_point_to_target,
    verify_finite_quotient_certificate,
)

model_fractions = tuple(Fraction(value) for value in minimal_model)
short_model, short_change = short_certificate_model(model_fractions)
point_fractions = tuple(
    (Fraction(point[0]), Fraction(point[1])) for point in rank_record["generators"]
)
short_points = tuple(
    source_point_to_target(point, short_change) for point in point_fractions
)
independence = rank_record["independence_certificate"]
verify_finite_quotient_certificate(short_model, short_points, independence)
if not independence["certified_independent"] or independence["certified_rank_lower_bound"] != 9:
    raise ArithmeticError("the nine-point independence certificate failed")

# Parameterize q1 by writing u = 65*t+r.  The paired base then becomes a
# pointed quartic after substituting t=(289444-r^2)/(130*r-38636) into q2.
parameter_ring = PolynomialRing(QQ, "r")
r = parameter_ring.gen()
leading_root_1 = QQ(q1[2]).sqrt()
parameter_denominator = 2 * leading_root_1 * r - q1[1]
parameter_t = (q1[0] - r**2) / parameter_denominator
parameter_u = leading_root_1 * parameter_t + r
if parameter_u**2 != q1(parameter_t):
    raise ArithmeticError("the first-conic parameterization failed")
paired_quartic = parameter_ring(parameter_denominator**2 * q2(parameter_t))
expected_paired_quartic = parameter_ring(
    [
        46344697121074403584,
        168863136988245440,
        -309051947898044,
        -962433973020,
        1346816601,
    ]
)
if paired_quartic != expected_paired_quartic:
    raise ArithmeticError("the paired-base quartic changed")

# Put the rational point at infinity into the affine chart s=1/r,
# y=z/r^2.  The reciprocal quartic has base point (0,36699).
reciprocal_coefficients = tuple(
    QQ(paired_quartic[4 - index]) for index in range(5)
)
quartic_q = QQ(reciprocal_coefficients[0]).sqrt()
if quartic_q != infinity_coordinates[1]:
    raise ArithmeticError("the pointed reciprocal quartic lost its base point")
_, quartic_d, quartic_c, quartic_b, quartic_a = reciprocal_coefficients
pointed_model = EllipticCurve(
    QQ,
    [
        quartic_d / quartic_q,
        quartic_c - quartic_d**2 / (4 * quartic_q**2),
        2 * quartic_q * quartic_b,
        -4 * quartic_q**2 * quartic_a,
        quartic_a * (quartic_d**2 - 4 * quartic_q**2 * quartic_c),
    ],
)
paired_minimal = pointed_model.global_minimal_model()
if tuple(int(value) for value in paired_minimal.a_invariants()) != EXPECTED_PAIRED_MODEL:
    raise ArithmeticError("the paired-base minimal model changed")
if int(paired_minimal.conductor()) != int(minimal.conductor()):
    raise ArithmeticError("the two isogenous curves have different conductors")

candidate_isogenies = [
    morphism
    for morphism in minimal.isogenies_prime_degree(2)
    if morphism.codomain().is_isomorphic(paired_minimal)
]
if len(candidate_isogenies) != 1:
    raise ArithmeticError("the expected degree-two isogeny is not unique")
third_to_paired = candidate_isogenies[0]
if third_to_paired.codomain() != paired_minimal:
    codomain_to_paired = third_to_paired.codomain().isomorphism_to(paired_minimal)
else:
    codomain_to_paired = None
paired_to_pointed = paired_minimal.isomorphism_to(pointed_model)


def isogeny_image(point):
    image = third_to_paired(point)
    return codomain_to_paired(image) if codomain_to_paired is not None else image


def paired_point_to_cover(point):
    """Return (t,u,v,r,s,y) under the exact inverse pointed-quartic map."""

    if point.is_zero():
        raise ArithmeticError("the selected paired point maps to the quartic origin")
    pointed_point = paired_to_pointed(point)
    x_value, y_value = pointed_point[0], pointed_point[1]
    if y_value == 0:
        raise ArithmeticError("the inverse pointed-quartic map hit an exceptional point")
    s_value = (
        4 * quartic_q**2 * (x_value + quartic_c) - quartic_d**2
    ) / (2 * quartic_q * y_value)
    if s_value == 0:
        raise ArithmeticError("the inverse pointed-quartic map returned infinity")
    reciprocal_y = (
        (x_value * s_value**2 - quartic_d * s_value) / (2 * quartic_q)
        - quartic_q
    )
    reciprocal_value = sum(
        coefficient * s_value**index
        for index, coefficient in enumerate(reciprocal_coefficients)
    )
    if reciprocal_y**2 != reciprocal_value:
        raise ArithmeticError("the pointed-quartic inverse failed")
    r_value = 1 / s_value
    t_value = (q1[0] - r_value**2) / (2 * leading_root_1 * r_value - q1[1])
    u_value = leading_root_1 * t_value + r_value
    z_value = reciprocal_y * r_value**2
    v_value = z_value / (2 * leading_root_1 * r_value - q1[1])
    if u_value**2 != q1(t_value) or v_value**2 != q2(t_value):
        raise ArithmeticError("the explicit point missed the paired base")
    return t_value, u_value, v_value, r_value, s_value, reciprocal_y


paired_points = tuple(isogeny_image(point) for point in points)
if any(point.is_zero() for point in paired_points):
    raise ArithmeticError("an independent point fell into the isogeny kernel")

# Either sign remains an independent basis.  Store the sign yielding the
# shorter t-coordinate so that the certificate is usable by specialization
# searches without changing the exact argument.
explicit_cover_points = []
explicit_t_values = set()
for index, paired_point in enumerate(paired_points, start=1):
    alternatives = []
    for sign in (1, -1):
        signed_point = sign * paired_point
        cover_point = paired_point_to_cover(signed_point)
        t_value = cover_point[0]
        bit_score = (
            abs(t_value.numerator()).nbits() + t_value.denominator().nbits()
        )
        alternatives.append((bit_score, -sign, sign, signed_point, cover_point))
    _, _, sign, signed_point, cover_point = min(alternatives)
    t_value, u_value, v_value, r_value, s_value, reciprocal_y = cover_point
    if t_value in explicit_t_values:
        raise ArithmeticError("two transported generators produced the same t-value")
    explicit_t_values.add(t_value)
    explicit_cover_points.append(
        {
            "source_generator": index,
            "source_sign": sign,
            "paired_minimal_point": point_text(signed_point),
            "t": rational_text(t_value),
            "u": rational_text(u_value),
            "v": rational_text(v_value),
            "conic_parameter_r": rational_text(r_value),
            "reciprocal_parameter_s": rational_text(s_value),
            "reciprocal_quartic_y": rational_text(reciprocal_y),
        }
    )

target_masks = frozenset((42110, 43109))
in_disjoint_graph = False
pair_row_count = 0
with args.disjoint_pairs.open() as stream:
    for pair_row in csv.DictReader(stream, delimiter="\t"):
        pair_row_count += 1
        if frozenset(
            (int(pair_row["left_orbit_mask"]), int(pair_row["right_orbit_mask"]))
        ) == target_masks:
            in_disjoint_graph = True
if pair_row_count != 8895801 or in_disjoint_graph:
    raise ArithmeticError("the complete disjoint graph boundary changed")

pair_theorem = geometry["all_distinct_pairs"]
if pair_theorem["anti_invariant_height_matrix_on_pair_cover"] != [[24, 0], [0, 24]]:
    raise ArithmeticError("the imported paired-cover height matrix changed")
if pair_theorem["generic_mw_rank_lower_bound"] != 19:
    raise ArithmeticError("the imported surface rank consequence changed")

isogeny_x_map, isogeny_y_map = third_to_paired.rational_maps()
minimal_to_pointed_urst = tuple(QQ(value) for value in paired_to_pointed.tuple())

result = {
    "schema": "elkies-k3.elkies-2026-rank19-rank9-paired-base.v2",
    "status": "PASS_EXACT_RANK19_SURFACE_OVER_RANK_AT_LEAST_9_BASE",
    "inputs": {
        display_path(args.catalogue): digest(args.catalogue),
        display_path(args.ledger): digest(args.ledger),
        display_path(args.geometry): digest(args.geometry),
        display_path(args.disjoint_pairs): digest(args.disjoint_pairs),
        display_path(FINITE_QUOTIENT_HELPER): digest(FINITE_QUOTIENT_HELPER),
        display_path(SHORT_MODEL_HELPER): digest(SHORT_MODEL_HELPER),
    },
    "selection": {
        "pair_key": PAIR_KEY,
        "arithmetic_complexity_rank": row["arithmetic_complexity_rank"],
        "orbit_masks": [42110, 43109],
        "orbit_hex": ["0x0a47e", "0x0a865"],
        "outside_norm_four_disjoint_pair_graph": True,
    },
    "paired_base": {
        "equations": [f"u^2={q1}", f"v^2={q2}"],
        "rational_point_at_infinity": {
            "u_over_t": rational_text(infinity_coordinates[0]),
            "v_over_t": rational_text(infinity_coordinates[1]),
        },
        "genus": 1,
        "third_quotient_quartic_coefficients_low_to_high": [
            rational_text(quartic[index]) for index in range(5)
        ],
        "minimal_jacobian_a1_a2_a3_a4_a6": list(minimal_model),
        "conductor": int(minimal.conductor()),
        "global_root_number": int(minimal.root_number()),
        "independent_point_count": len(points),
        "points": rank_record["generators"],
        "independence_certificate": independence,
        "paired_base_rank_lower_bound": 9,
        "positive_rank_transfer": (
            "The fixed-point-free product involution gives a degree-two isogeny to the third "
            "quotient, so its nine independent Jacobian points prove paired-base rank at least 9."
        ),
        "explicit_parameter_map": {
            "first_conic_parameterization": {
                "parameter": "r",
                "formula_t": "(289444-r^2)/(130*r-38636)",
                "formula_u": "65*t+r",
            },
            "paired_quartic_equation": (
                "z^2=" + polynomial_text(paired_quartic)
            ),
            "paired_quartic_coefficients_low_to_high": [
                rational_text(paired_quartic[index]) for index in range(5)
            ],
            "reciprocal_quartic_coordinates": "s=1/r, y=z/r^2",
            "reciprocal_quartic_coefficients_low_to_high": [
                rational_text(value) for value in reciprocal_coefficients
            ],
            "reciprocal_quartic_origin": ["0", rational_text(quartic_q)],
            "pointed_weierstrass_a1_a2_a3_a4_a6": [
                rational_text(value) for value in pointed_model.a_invariants()
            ],
            "paired_minimal_a1_a2_a3_a4_a6": list(EXPECTED_PAIRED_MODEL),
            "minimal_to_pointed_urst": [
                rational_text(value) for value in minimal_to_pointed_urst
            ],
            "minimal_to_pointed_formula": (
                "X=(x-r0)/u0^2; Y=(y-s0*(x-r0)-t0)/u0^3 for "
                "(u0,r0,s0,t0)=minimal_to_pointed_urst"
            ),
            "pointed_to_quartic_formula": (
                "s=(4*q^2*(X+c)-d^2)/(2*q*Y); "
                "y=(X*s^2-d*s)/(2*q)-q, with q=36699 and "
                "(d,c) the linear and quadratic reciprocal-quartic coefficients"
            ),
            "quartic_to_cover_formula": (
                "r=1/s; t=(289444-r^2)/(130*r-38636); "
                "u=65*t+r; v=y*r^2/(130*r-38636)"
            ),
            "third_quotient_to_paired_isogeny": {
                "degree": int(third_to_paired.degree()),
                "kernel_polynomial": str(third_to_paired.kernel_polynomial()),
                "x_map": str(isogeny_x_map),
                "y_map": str(isogeny_y_map),
            },
            "transported_independent_points": explicit_cover_points,
            "distinct_explicit_t_value_count": len(explicit_t_values),
        },
    },
    "surface": {
        "anti_invariant_height_matrix": [[24, 0], [0, 24]],
        "generic_mw_rank_lower_bound": 19,
        "infinitely_many_rational_t_values": True,
    },
    "proof_boundary": (
        "The certificate proves the base rank is at least 9, not exactly 9, and proves generic "
        "surface rank at least 19. It asserts no exceptional specialization rank."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026R19R9BASE|pair=42110:43109|base_rank_lower_bound=9|"
    "explicit_t_values=9|surface_generic_rank_lower_bound=19|outside_disjoint_graph=true|"
    f"status={result['status']}|output={display_path(args.output)}"
)
