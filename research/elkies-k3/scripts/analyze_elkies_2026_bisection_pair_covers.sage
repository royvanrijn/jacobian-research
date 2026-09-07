#!/usr/bin/env sage-python
"""Classify the rational and paired bases of the complete R17 bisection batch.

For every equation-level bisection the batch records a squarefree quadratic

    u^2 = q(t).

This script checks, over QQ and without a bounded point search, whether the
associated conic has a rational point.  It also compares the geometric branch
divisors of all records.  When every branch quadratic is irreducible and no
two are proportional, every pair has four distinct geometric branch points;
its connected V4 fibre product consequently has genus one.

The rank statement in the output uses only the already-certified inputs in
the batch: invariant rank 17 and anti-invariant height 12 on each quadratic
cover.  On a distinct paired cover the two sections have different V4
characters.  Their pullback height matrix is therefore diag(24,24).
"""

from __future__ import annotations

import argparse
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

from sage.all import Conic, PolynomialRing, QQ, ZZ, lcm


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-bisection-pair-cover-geometry-full.json"
EXPECTED_SCHEMA = "elkies-k3.bisection-extension-input.v1"
OUTPUT_SCHEMA = "elkies-k3.bisection-pair-cover-geometry.v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


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
    if len(values) != 3 or not values[2]:
        raise ValueError("each branch equation must be a genuine quadratic")
    denominator = lcm([value.denominator() for value in values])
    integers = [ZZ(value * denominator) for value in values]
    content = reduce(gcd, (abs(int(value)) for value in integers))
    integers = tuple(value // content for value in integers)
    if integers[2] < 0:
        integers = tuple(-value for value in integers)
    return integers


def integer_square(value) -> bool:
    value = ZZ(value)
    return value >= 0 and value.is_square()


def complexity_key(record):
    cost = record["equation_complexity"]
    return (
        int(cost["group_addition_upper_bound"]),
        int(cost["support_count"]),
        int(cost["dependency_count"]),
        int(cost["coordinate_input_bits"]),
        int(cost["maximum_absolute_coefficient"]),
        int(cost["coefficient_l1"]),
        tuple(int(value) for value in record["published_basis_w"]),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument(
    "--sample-size",
    type=int,
    default=100,
    help="number of cheapest exact conic points retained in the compact output",
)
args = parser.parse_args()

if args.sample_size < 0:
    raise ValueError("--sample-size must be nonnegative")

payload = json.loads(args.input.read_text())
if payload.get("schema") != EXPECTED_SCHEMA:
    raise ValueError(f"unexpected input schema: {payload.get('schema')!r}")
records = payload.get("bisections", [])
if not records:
    raise ValueError("input contains no bisections")

certificate = payload["individual_base_change_certificate"]
if int(payload["invariant_mw_rank"]) != 17:
    raise ValueError("paired rank calculation expects invariant rank 17")
if int(certificate["anti_invariant_height"]) != 12:
    raise ValueError("paired height calculation expects individual height 12")
if not certificate["branch_fibres_smooth_for_every_record"]:
    raise ValueError("the batch has not certified smooth branch fibres")

P = PolynomialRing(QQ, names=("T", "U", "Z"))
T, U, Z = P.gens()
branch_keys = {}
point_hash = sha256()
point_samples = []
rational_point_count = 0
affine_point_count = 0
infinity_point_count = 0
irreducible_count = 0
maximum_point_coordinate_bits = 0
square_constant_masks = set()
square_leading_masks = set()

for position, record in enumerate(sorted(records, key=complexity_key), start=1):
    q_coefficients = tuple(QQ(value) for value in record["residual_chord"]["q_coefficients"])
    coefficients = primitive_quadratic(q_coefficients)
    a, b, c = coefficients
    orbit_mask = int(record["lattice_orbit_mask"])
    if q_coefficients[0].is_square():
        square_constant_masks.add(orbit_mask)
    if q_coefficients[2].is_square():
        square_leading_masks.add(orbit_mask)
    discriminant = b * b - 4 * a * c
    irreducible = not integer_square(discriminant)
    irreducible_count += int(irreducible)
    branch_keys.setdefault(tuple(int(value) for value in coefficients), []).append(
        orbit_mask
    )

    q_a, q_b, q_c = q_coefficients
    conic = Conic(U**2 - (q_c * T**2 + q_b * T * Z + q_a * Z**2))
    soluble, point = conic.has_rational_point(point=True)
    if not soluble:
        continue
    rational_point_count += 1
    coordinates = tuple(QQ(value) for value in point)
    affine_point_count += int(bool(coordinates[2]))
    infinity_point_count += int(not coordinates[2])
    coordinate_text = tuple(rational_text(value) for value in coordinates)
    point_hash.update(
        (f"{orbit_mask}:" + ",".join(coordinate_text) + "\n").encode()
    )
    maximum_point_coordinate_bits = max(
        maximum_point_coordinate_bits,
        *(max(abs(int(value.numerator())), int(value.denominator())).bit_length() for value in coordinates),
    )
    if len(point_samples) < args.sample_size:
        point_samples.append(
            {
                "priority_rank": int(record["priority_rank"]),
                "equation_rank": int(record["equation_rank"]),
                "lattice_orbit_mask": orbit_mask,
                "orbit_hex": f"0x{orbit_mask:05x}",
                "primitive_q_coefficients_low_to_high": [int(value) for value in coefficients],
                "conic_point_T_U_Z": list(coordinate_text),
            }
        )

record_count = len(records)
duplicate_branch_groups = [
    {"primitive_q_coefficients_low_to_high": list(key), "orbit_masks": masks}
    for key, masks in branch_keys.items()
    if len(masks) > 1
]
all_pair_count = record_count * (record_count - 1) // 2
all_irreducible = irreducible_count == record_count
all_geometric_branches_distinct = len(branch_keys) == record_count
all_conics_rational = rational_point_count == record_count
pair_theorem_applies = all_irreducible and all_geometric_branches_distinct
square_both_masks = square_constant_masks & square_leading_masks


def pair_count(size):
    return size * (size - 1) // 2


easy_rational_pair_count = (
    pair_count(len(square_constant_masks))
    + pair_count(len(square_leading_masks))
    - pair_count(len(square_both_masks))
)

result = {
    "schema": OUTPUT_SCHEMA,
    "status": (
        "PASS_COMPLETE_CONIC_AND_GENUS_ONE_PAIR_CLASSIFICATION"
        if pair_theorem_applies
        else "INCOMPLETE_OR_COUNTEREXAMPLE"
    ),
    "input": {
        "path": display_path(args.input),
        "sha256": digest(args.input),
        "schema": payload["schema"],
    },
    "complete_conic_classification": {
        "record_count": record_count,
        "rational_point_count": rational_point_count,
        "anisotropic_conic_count": record_count - rational_point_count,
        "affine_solver_point_count": affine_point_count,
        "infinity_solver_point_count": infinity_point_count,
        "all_conics_Q_rational": all_conics_rational,
        "canonical_point_ledger_sha256": point_hash.hexdigest(),
        "maximum_solver_point_coordinate_bits": maximum_point_coordinate_bits,
        "proof_method": "Sage Conic.has_rational_point over QQ (exact Hasse-Minkowski algorithm)",
        "consequence": (
            f"Exactly {rational_point_count} quadratic bisection bases are Q-rational conics "
            "and can be parameterized, giving that many explicit P1-based "
            "generic-rank-at-least-18 families."
        ),
    },
    "complete_branch_classification": {
        "irreducible_quadratic_count": irreducible_count,
        "all_branch_quadratics_irreducible_over_Q": all_irreducible,
        "distinct_primitive_branch_quadratic_count": len(branch_keys),
        "all_geometric_branch_divisors_distinct": all_geometric_branches_distinct,
        "duplicate_branch_groups": duplicate_branch_groups,
    },
    "all_distinct_pairs": {
        "pair_count": all_pair_count,
        "v4_connected": pair_theorem_applies,
        "geometric_branch_point_count": 4 if pair_theorem_applies else None,
        "fiber_product_genus": 1 if pair_theorem_applies else None,
        "riemann_hurwitz": "2g-2=4*(-2)+4*(4-2)=0",
        "invariant_rank": 17,
        "new_character_count": 2,
        "anti_invariant_height_before_second_pullback": 12,
        "anti_invariant_height_matrix_on_pair_cover": [[24, 0], [0, 24]],
        "height_matrix_rank": 2,
        "generic_mw_rank_lower_bound": 19,
        "proof": (
            "Distinct quadratic squareclasses give a connected V4 cover. The two new sections "
            "belong to distinct nontrivial characters, so their cross-height is zero. A degree-two "
            "pullback doubles each certified height 12 to 24."
        ),
        "rational-base_boundary": (
            "No pair has a shared geometric branch point, so this catalogue contains no genus-zero "
            "paired base. Rational points and Mordell-Weil ranks of the genus-one pair bases require "
            "a separate arithmetic scan."
        ),
        "immediate_Q_point_subfamily": {
            "square_constant_cover_count": len(square_constant_masks),
            "square_leading_cover_count": len(square_leading_masks),
            "square_both_cover_count": len(square_both_masks),
            "distinct_pair_count": easy_rational_pair_count,
            "proof": (
                "Two square constants give a common fibre-product point over t=0; two square "
                "leading coefficients give a common point over t=infinity. Inclusion-exclusion "
                "removes pairs counted by both tests."
            ),
        },
    },
    "cheapest_conic_point_samples": point_samples,
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "ELKIES2026BISECTIONPAIRS|"
    f"records={record_count}|Q_rational_conics={rational_point_count}|"
    f"irreducible_quadratics={irreducible_count}|distinct_branches={len(branch_keys)}|"
    f"pairs={all_pair_count}|pair_genus={result['all_distinct_pairs']['fiber_product_genus']}|"
    f"pair_height_rank={result['all_distinct_pairs']['height_matrix_rank']}|"
    f"status={result['status']}|output={display_path(args.output)}"
)
