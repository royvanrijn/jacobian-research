#!/usr/bin/env sage-python
"""Verify a new low-complexity rank-19 paired bisection base.

Among the complete equation-ranked bisection catalogue, take the two cheapest
covers whose quadratic leading coefficients are rational squares.  Their V4
fibre product has a rational point over t=infinity.  This script verifies the
two exact branch equations, the genus-one quotient, its Jacobian, and three
independent rational points on that Jacobian.

The quotient map from the V4 curve to ``w^2=q1*q2`` is an unramified double
cover, hence a 2-isogeny after choosing the rational point at infinity.  The
rank-three lower bound therefore transfers to the paired base.  The two new
surface sections have distinct V4 characters and exact height matrix
diag(24,24), proving generic Mordell--Weil rank at least 19 over its function
field.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys

from sage.all import EllipticCurve, Infinity, PolynomialRing, QQ, ZZ, lcm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
DEFAULT_PAIRS = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority-disjoint-pairs-full.tsv"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank19-cheapest-infinity-pair.json"
FINITE_QUOTIENT_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
SHORT_MODEL_HELPER = ROOT / "elliptic-curves/ecsearch/q12o5867_specialization.py"
EXPECTED_MASKS = (28257, 81769)
EXPECTED_PRIMITIVE_Q = (
    (41627760409, 15206854416, 2278725696),
    (126480025, 108563070, 21650409),
)
EXPECTED_MINIMAL_MODEL = (
    1,
    1,
    1,
    -10283115414666818054869518594495,
    237013667219266831004461300259717721649810997157,
)
EXPECTED_GENERATORS = (
    (
        Fraction(-115272385602710865443, 19881),
        Fraction(894022058066911731326422170896, 2803221),
    ),
    (
        Fraction(-350052388948929430349333, 690060361),
        Fraction(-8919234683413439600875236866002684966, 18127195623109),
    ),
    (
        Fraction(6906836279496188575, 49),
        Fraction(-18147830468021807993312965366, 343),
    ),
)


def file_digest(path: Path) -> str:
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


def primitive_quadratic(coefficients):
    values = tuple(QQ(value) for value in coefficients)
    denominator = lcm([value.denominator() for value in values])
    integers = [ZZ(value * denominator) for value in values]
    content = reduce(gcd, (abs(int(value)) for value in integers))
    integers = tuple(value // content for value in integers)
    if integers[2] < 0:
        integers = tuple(-value for value in integers)
    return integers


def rational_square_root(value):
    value = QQ(value)
    if not value.is_square():
        raise ArithmeticError(f"{value} is not a rational square")
    return QQ(value.sqrt())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

batch = json.loads(args.input.read_text())
records = batch["bisections"]
by_mask = {int(record["lattice_orbit_mask"]): record for record in records}
if len(by_mask) != len(records):
    raise ArithmeticError("duplicate orbit masks in bisection batch")

leading_square_records = sorted(
    (
        record
        for record in records
        if QQ(record["residual_chord"]["q_coefficients"][2]).is_square()
    ),
    key=lambda record: int(record["priority_rank"]),
)
selected_masks = tuple(int(record["lattice_orbit_mask"]) for record in leading_square_records[:2])
if selected_masks != EXPECTED_MASKS:
    raise ArithmeticError(f"cheapest infinity pair changed: {selected_masks}")
selected = tuple(by_mask[mask] for mask in selected_masks)

primitive_q = tuple(
    primitive_quadratic(record["residual_chord"]["q_coefficients"])
    for record in selected
)
if primitive_q != EXPECTED_PRIMITIVE_Q:
    raise ArithmeticError(f"primitive branch equations changed: {primitive_q}")

square_rescalings = []
for record, coefficients in zip(selected, primitive_q):
    actual = tuple(QQ(value) for value in record["residual_chord"]["q_coefficients"])
    scale = actual[0] / coefficients[0]
    square_rescalings.append(rational_square_root(scale))
    if actual != tuple(scale * value for value in coefficients):
        raise ArithmeticError("branch normalization is not a scalar rescaling")

R = PolynomialRing(QQ, "t")
t = R.gen()
q1, q2 = (R(list(coefficients)) for coefficients in primitive_q)
if not q1.is_irreducible() or not q2.is_irreducible() or q1.gcd(q2) != 1:
    raise ArithmeticError("the selected covers do not have four distinct geometric branch points")
leading_roots = tuple(rational_square_root(q[2]) for q in (q1, q2))

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
jacobian = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
minimal = jacobian.global_minimal_model()
minimal_model = tuple(int(value) for value in minimal.a_invariants())
if minimal_model != EXPECTED_MINIMAL_MODEL:
    raise ArithmeticError(f"minimal Jacobian changed: {minimal_model}")

points = tuple(minimal(QQ(x), QQ(y)) for x, y in EXPECTED_GENERATORS)
if any(point.order() != Infinity for point in points):
    raise ArithmeticError("a displayed positive-rank point became torsion")

# Reuse the same independent finite-quotient verifier as the published paired
# cover certificate.  It works on a short integral model, so transport the
# points by the exact Weierstrass change first.
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    source_point_to_target,
    verify_finite_quotient_certificate,
)

minimal_fractions = tuple(Fraction(value) for value in minimal_model)
short_model, short_change = short_certificate_model(minimal_fractions)
short_points = tuple(
    source_point_to_target((Fraction(x), Fraction(y)), short_change)
    for x, y in EXPECTED_GENERATORS
)
independence = build_finite_quotient_certificate(
    short_model, short_points, relation_prime=3, prime_bound=500
)
verify_finite_quotient_certificate(short_model, short_points, independence)
if not independence["certified_independent"]:
    raise ArithmeticError("the three Jacobian points were not certified independent")

target_pair = frozenset(selected_masks)
in_disjoint_graph = False
pair_row_count = 0
with args.pairs.open() as stream:
    for row in csv.DictReader(stream, delimiter="\t"):
        pair_row_count += 1
        if frozenset((int(row["left_orbit_mask"]), int(row["right_orbit_mask"]))) == target_pair:
            in_disjoint_graph = True
if pair_row_count != 8895801:
    raise ArithmeticError(f"unexpected disjoint graph size: {pair_row_count}")
if in_disjoint_graph:
    raise ArithmeticError("the candidate unexpectedly lies in the norm-four disjointness graph")

height = int(batch["individual_base_change_certificate"]["anti_invariant_height"])
if int(batch["invariant_mw_rank"]) != 17 or height != 12:
    raise ArithmeticError("the imported rank/height certificate changed")

result = {
    "schema": "elkies-k3.elkies-2026-rank19-cheapest-infinity-pair.v1",
    "status": "PASS_EXACT_NEW_RANK19_GENUS_ONE_BASE_RANK_AT_LEAST_3",
    "inputs": {
        display_path(args.input): file_digest(args.input),
        display_path(args.pairs): file_digest(args.pairs),
        display_path(FINITE_QUOTIENT_HELPER): file_digest(FINITE_QUOTIENT_HELPER),
        display_path(SHORT_MODEL_HELPER): file_digest(SHORT_MODEL_HELPER),
    },
    "selection": {
        "rule": "two lowest priority ranks with square quadratic leading coefficient",
        "leading_square_cover_count": len(leading_square_records),
        "orbit_masks": list(selected_masks),
        "orbit_hex": [f"0x{mask:05x}" for mask in selected_masks],
        "priority_ranks": [int(record["priority_rank"]) for record in selected],
        "equation_ranks": [int(record["equation_rank"]) for record in selected],
        "outside_norm_four_disjoint_pair_graph": not in_disjoint_graph,
        "interpretation": (
            "The disjointness graph is a geometric priority heuristic, not a character-"
            "independence requirement; this exact candidate lies outside it."
        ),
    },
    "paired_base": {
        "equations": [
            "u^2=41627760409+15206854416*t+2278725696*t^2",
            "v^2=126480025+108563070*t+21650409*t^2",
        ],
        "square_rescalings_from_batch_q": [rational_text(value) for value in square_rescalings],
        "rational_point_at_infinity": {
            "u_over_t": rational_text(leading_roots[0]),
            "v_over_t": rational_text(leading_roots[1]),
            "uv_over_t_squared": rational_text(leading_roots[0] * leading_roots[1]),
        },
        "geometric_branch_points": 4,
        "genus": 1,
        "third_quotient": {
            "equation": f"w^2=({q1})*({q2})",
            "quartic_coefficients_low_to_high": [rational_text(quartic[index]) for index in range(5)],
            "binary_quartic_I": rational_text(invariant_i),
            "binary_quartic_J": rational_text(invariant_j),
        },
        "jacobian_minimal_model_a1_a2_a3_a4_a6": list(minimal_model),
        "jacobian_points": [
            [rational_text(x), rational_text(y)] for x, y in EXPECTED_GENERATORS
        ],
        "jacobian_independence": independence,
        "jacobian_rank_lower_bound": 3,
        "paired_base_rank_lower_bound": 3,
        "positive_rank_transfer": (
            "The fixed-point-free product involution gives a degree-two isogeny from the "
            "paired genus-one base to the third quotient; isogenous elliptic curves have "
            "the same rational rank."
        ),
    },
    "surface_rank_consequence": {
        "invariant_rank": 17,
        "anti_invariant_height_matrix": [[2 * height, 0], [0, 2 * height]],
        "height_matrix_rank": 2,
        "generic_mw_rank_lower_bound_over_paired_base": 19,
        "infinitely_many_rational_base_points": True,
        "infinitely_many_rational_t_values": True,
    },
    "proof_boundary": (
        "The finite-quotient certificate proves only Jacobian rank at least 3, not an upper "
        "bound. The surface conclusion is generic rank at least 19; no specialization rank "
        "beyond the inherited sections is asserted here."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026R19CHEAPPAIR|"
    f"masks={selected_masks}|pair_genus=1|base_rank_lower_bound=3|"
    "surface_generic_rank_lower_bound=19|outside_disjoint_graph=true|"
    f"status={result['status']}|output={display_path(args.output)}"
)
