#!/usr/bin/env sage-python
"""Certify a bounded shortlist of rational alternate-Q80 V4 pair bases.

status: ACTIVE_PROOF
claim: exact rational points on a bounded shortlist of genus-one V4 bases
inputs: cheapest-1024 alternate bisections and the direct alternate equation
outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json
supersedes: none; gates the quartic product-twist P.O=0 search

A norm-ten rational bisection has class ``(2,2,w)`` in the direct
``U + MW(-1)`` splitting.  Hence two such curves have intersection

    (2,2,w_i).(2,2,w_j) = 8 - <w_i,w_j>.

When this number is one, the two QQ-rational curves have one geometric
intersection point counted with multiplicity.  That point is Galois fixed and
therefore gives an exact rational point on their V4 fibre product.  This
script recovers the point independently from the two explicit lifted-section
equations and verifies both square-cover equations coefficientwise.

The shortlist is selected only from the declared priority prefix.  No failure
outside the prefix, no Mordell--Weil rank upper bound for a base Jacobian, and
no section on a quartic product twist is claimed.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, prime_range, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
)
OUTPUT_SCHEMA = "elkies-k3.r17-norm12-11952-v4-pair-shortlist.v1"


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


def rational_bits(value) -> int:
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def polynomial(record, key, ring):
    return ring([QQ(value) for value in record["lifted_section"][f"{key}_coefficients"]])


def branch_polynomial(record, ring):
    numerator = ring([QQ(value) for value in record["branch"]["numerator_coefficients"]])
    denominator = ring([QQ(value) for value in record["branch"]["denominator_coefficients"]])
    if denominator.degree() != 0 or not denominator:
        raise ArithmeticError(f"{record['label']}: branch denominator is not constant")
    result = numerator / denominator[0]
    if result.degree() != 2 or not result.is_squarefree():
        raise ArithmeticError(f"{record['label']}: branch is not a squarefree quadratic")
    return ring(result)


def recover_unique_affine_intersection(left, right, ring):
    """Recover the unique intersection from two linear-in-radical lifts."""

    field = ring.fraction_field()
    q_left = branch_polynomial(left, ring)
    q_right = branch_polynomial(right, ring)
    x0_left, x1_left, y0_left, y1_left = (
        polynomial(left, key, ring) for key in ("x0", "x1", "y0", "y1")
    )
    x0_right, x1_right, y0_right, y1_right = (
        polynomial(right, key, ring) for key in ("x0", "x1", "y0", "y1")
    )
    delta_x = x0_right - x0_left
    delta_y = y0_right - y0_left
    determinant = -x1_left * y1_right + x1_right * y1_left
    if not determinant:
        raise ArithmeticError("the affine radical-recovery determinant vanishes identically")
    left_radical = field((-delta_x * y1_right + x1_right * delta_y) / determinant)
    right_radical = field((x1_left * delta_y - delta_x * y1_left) / determinant)
    left_residual = ring((left_radical**2 - q_left).numerator())
    right_residual = ring((right_radical**2 - q_right).numerator())
    common = left_residual.gcd(right_residual).monic()
    linear_factors = [factor for factor, _exponent in common.factor() if factor.degree() == 1]
    candidates = []
    for factor in linear_factors:
        parameter = QQ(-factor[0] / factor[1])
        if not left_radical.denominator()(parameter) or not right_radical.denominator()(parameter):
            continue
        left_root = QQ(left_radical(parameter))
        right_root = QQ(right_radical(parameter))
        if left_root**2 != q_left(parameter) or right_root**2 != q_right(parameter):
            continue
        x_left = QQ(x0_left(parameter) + x1_left(parameter) * left_root)
        y_left = QQ(y0_left(parameter) + y1_left(parameter) * left_root)
        x_right = QQ(x0_right(parameter) + x1_right(parameter) * right_root)
        y_right = QQ(y0_right(parameter) + y1_right(parameter) * right_root)
        if x_left != x_right or y_left != y_right:
            continue
        candidates.append((parameter, left_root, right_root, x_left, y_left))
    if len(candidates) != 1:
        raise ArithmeticError(
            f"expected one verified affine intersection, recovered {len(candidates)}"
        )
    parameter, left_root, right_root, x_value, y_value = candidates[0]
    return {
        "u": rational_text(parameter),
        "left_square_root": rational_text(left_root),
        "right_square_root": rational_text(right_root),
        "product_square_root": rational_text(left_root * right_root),
        "surface_x": rational_text(x_value),
        "surface_y": rational_text(y_value),
        "maximum_coordinate_bits": max(
            rational_bits(value)
            for value in (parameter, left_root, right_root, x_value, y_value)
        ),
        "recovery_gcd_degree": int(common.degree()),
        "both_cover_equations_verified": True,
        "two_lifted_surface_points_equal": True,
    }


def quartic_coefficients(left_q, right_q):
    product = left_q * right_q
    return product, tuple(QQ(product[index]) for index in range(5))


def finite_field_pair_data(left_q, right_q, model, prime):
    """Return exact good-reduction/local-density data, or ``None``."""

    field = GF(prime)
    ring = PolynomialRing(field, "u")
    try:
        left = ring([field(value) for value in left_q])
        right = ring([field(value) for value in right_q])
        coefficient_a = ring([field(QQ(value)) for value in model["A_coefficients_low_to_high"]])
        coefficient_b = ring([field(QQ(value)) for value in model["B_coefficients_low_to_high"]])
        discriminant = ring(
            [field(QQ(value)) for value in model["discriminant_coefficients_low_to_high"]]
        )
    except (ZeroDivisionError, TypeError, ValueError):
        return None
    if (
        left.degree() != 2
        or right.degree() != 2
        or discriminant.degree() != 24
        or not left.is_squarefree()
        or not right.is_squarefree()
        or not discriminant.is_squarefree()
        or left.gcd(right).degree()
        or (left * right).gcd(discriminant).degree()
    ):
        return None

    def square(value):
        return not value or value.is_square()

    simultaneous = sum(
        square(left(value)) and square(right(value)) for value in field
    )
    simultaneous += int(square(left[2]) and square(right[2]))

    twist = left * right
    twist_a = coefficient_a * twist**2
    twist_b = coefficient_b * twist**3

    def fibre_has_no_two_torsion(a_value, b_value):
        fibre_discriminant = -field(16) * (
            field(4) * a_value**3 + field(27) * b_value**2
        )
        return bool(fibre_discriminant) and not any(
            x**3 + a_value * x + b_value == 0 for x in field
        )

    chart_parameter = None
    for value in field:
        if fibre_has_no_two_torsion(twist_a(value), twist_b(value)):
            chart_parameter = int(value)
            break
    if chart_parameter is None:
        a_infinity = twist_a[16]
        b_infinity = twist_b[24]
        if fibre_has_no_two_torsion(a_infinity, b_infinity):
            chart_parameter = "infinity"
    return {
        "prime": int(prime),
        "simultaneous_split_points_on_P1_Fp": int(simultaneous),
        "density_numerator": int(simultaneous),
        "density_denominator": int(prime + 1),
        "complete_polynomial_section_chart_parameter": chart_parameter,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--direct", type=Path, default=DIRECT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--pool-size", type=int, default=1024)
    parser.add_argument("--shortlist-size", type=int, default=64)
    parser.add_argument("--candidate-window", type=int, default=512)
    parser.add_argument("--prime-bound", type=int, default=211)
    parser.add_argument("--recommended-primes", type=int, default=3)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.pool_size < 2 or args.shortlist_size < 1 or args.candidate_window < args.shortlist_size:
        parser.error("require pool-size >= 2 and candidate-window >= shortlist-size >= 1")
    if args.prime_bound < 7 or args.recommended_primes < 1:
        parser.error("prime-bound and recommended-primes are too small")

    batch = json.loads(args.bisections.read_text())
    direct = json.loads(args.direct.read_text())
    if batch.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("unexpected bisection schema")
    if batch.get("status") not in {
        "PASS_EXACT_ALTERNATE_BISECTION_EQUATION_CHUNK",
        "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
    }:
        raise ValueError("input is not an exact alternate-bisection compilation")
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("input is not the exact direct alternate-Q80 model")
    records = sorted(batch["bisections"], key=lambda row: int(row["priority_rank"]))
    if len(records) < args.pool_size:
        raise ValueError("the bisection artifact is shorter than --pool-size")
    records = records[: args.pool_size]
    if [int(row["priority_rank"]) for row in records] != list(range(1, args.pool_size + 1)):
        raise ArithmeticError("the selected input is not the complete declared priority prefix")

    gram = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    if gram.nrows() != 17 or gram.det() != 948:
        raise ArithmeticError("the direct alternate frame changed")
    vectors = [vector(ZZ, record["direct_alternate_w"]) for record in records]
    if any(value * gram * value != 10 for value in vectors):
        raise ArithmeticError("the selected prefix contains a non-norm-ten record")
    if any(
        int(record["lifted_section"]["anti_invariant_height"]) != 12
        or not record["lifted_section"]["constant_and_linear_identities_verified"]
        or not record["lifted_section"]["two_branches_verified"]
        for record in records
    ):
        raise ArithmeticError("one selected bisection lacks its exact lifted section")

    ring = PolynomialRing(QQ, "u")
    branches = [branch_polynomial(record, ring) for record in records]
    intersection_one = []
    intersection_histogram = {}
    for right_index in range(len(records)):
        for left_index in range(right_index):
            frame_pairing = int(vectors[left_index] * gram * vectors[right_index])
            intersection = 8 - frame_pairing
            intersection_histogram[str(intersection)] = (
                intersection_histogram.get(str(intersection), 0) + 1
            )
            if intersection != 1:
                continue
            product, coefficients = quartic_coefficients(
                branches[left_index], branches[right_index]
            )
            if product.degree() != 4 or not product.is_squarefree():
                raise ArithmeticError("an intersection-one pair has a singular product quartic")
            complexity = {
                "maximum_priority_rank": max(
                    int(records[left_index]["priority_rank"]),
                    int(records[right_index]["priority_rank"]),
                ),
                "priority_rank_sum": int(records[left_index]["priority_rank"])
                + int(records[right_index]["priority_rank"]),
                "quartic_maximum_coefficient_bits": max(
                    rational_bits(value) for value in coefficients
                ),
                "quartic_coefficient_bit_sum": sum(
                    rational_bits(value) for value in coefficients
                ),
            }
            selection_key = (
                complexity["maximum_priority_rank"],
                complexity["priority_rank_sum"],
                complexity["quartic_maximum_coefficient_bits"],
                complexity["quartic_coefficient_bit_sum"],
                records[left_index]["label"],
                records[right_index]["label"],
            )
            intersection_one.append(
                (selection_key, left_index, right_index, complexity)
            )
    intersection_one.sort(key=lambda item: item[0])
    candidate_window = intersection_one[: args.candidate_window]
    if len(candidate_window) < args.shortlist_size:
        raise ArithmeticError("too few intersection-one pairs in the candidate window")

    model = direct["weierstrass_model"]
    prime_candidates = list(prime_range(5, args.prime_bound + 1))
    shortlisted = []
    recovery_failures = []
    for _selection_key, left_index, right_index, complexity in candidate_window:
        left = records[left_index]
        right = records[right_index]
        try:
            point = recover_unique_affine_intersection(left, right, ring)
        except ArithmeticError as error:
            recovery_failures.append(
                {
                    "labels": [left["label"], right["label"]],
                    "reason": str(error),
                }
            )
            continue
        left_q = branches[left_index]
        right_q = branches[right_index]
        product, coefficients = quartic_coefficients(left_q, right_q)
        parameter = QQ(point["u"])
        if QQ(point["product_square_root"]) ** 2 != product(parameter):
            raise ArithmeticError("the recovered V4 point misses the product quartic")

        e, d, c, b, a = coefficients
        invariant_i = 12 * a * e - 3 * b * d + c**2
        invariant_j = (
            72 * a * c * e
            + 9 * b * c * d
            - 27 * a * d**2
            - 27 * b**2 * e
            - 2 * c**3
        )
        jacobian = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
        integral_jacobian = jacobian.global_integral_model()

        local_data = []
        for prime in prime_candidates:
            row = finite_field_pair_data(left_q, right_q, model, prime)
            if row is not None:
                local_data.append(row)
        recommended = [
            row["prime"]
            for row in local_data
            if row["complete_polynomial_section_chart_parameter"] is not None
        ][: args.recommended_primes]
        if len(recommended) < args.recommended_primes:
            continue

        pair_key = f"{left['label']}:{right['label']}"
        shortlisted.append(
            {
                "shortlist_rank": len(shortlisted) + 1,
                "pair_key": pair_key,
                "labels": [left["label"], right["label"]],
                "priority_ranks": [
                    int(left["priority_rank"]),
                    int(right["priority_rank"]),
                ],
                "lattice_orbit_masks": [
                    int(left["lattice_orbit_mask"]),
                    int(right["lattice_orbit_mask"]),
                ],
                "direct_frame_vectors": [
                    list(map(int, vectors[left_index])),
                    list(map(int, vectors[right_index])),
                ],
                "direct_frame_pairing": 7,
                "rational_bisection_intersection_number": 1,
                "v4_base_genus": 1,
                "v4_base_point": point,
                "branch_quadratics_coefficients_low_to_high": [
                    [rational_text(value) for value in left_q],
                    [rational_text(value) for value in right_q],
                ],
                "product_quartic_coefficients_low_to_high": [
                    rational_text(value) for value in coefficients
                ],
                "product_quartic_squarefree": True,
                "product_twist_chi": 4,
                "product_twist_P_dot_O_zero_degree_bounds": {
                    "X": 8,
                    "Y": 12,
                },
                "binary_quartic_invariants": {
                    "I": rational_text(invariant_i),
                    "J": rational_text(invariant_j),
                },
                "base_jacobian_integral_a1_a2_a3_a4_a6": [
                    rational_text(value) for value in integral_jacobian.a_invariants()
                ],
                "base_jacobian_rank_status": "UNKNOWN",
                "complexity": complexity,
                "good_small_prime_data": local_data,
                "recommended_complete_sieve_primes": recommended,
            }
        )
        print(
            "ALTV4PAIR|"
            f"rank={len(shortlisted)}|pair={pair_key}|"
            f"point_bits={point['maximum_coordinate_bits']}|"
            f"primes={','.join(map(str, recommended))}",
            flush=True,
        )
        if len(shortlisted) == args.shortlist_size:
            break
    if len(shortlisted) != args.shortlist_size:
        raise ArithmeticError(
            f"recovered only {len(shortlisted)} of {args.shortlist_size} requested pairs"
        )

    result = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_EXACT_BOUNDED_RATIONAL_V4_PAIR_SHORTLIST",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.bisections): digest(args.bisections),
            display_path(args.direct): digest(args.direct),
        },
        "scope": {
            "native_priority_prefix_size": args.pool_size,
            "all_unordered_pairs_in_prefix": args.pool_size * (args.pool_size - 1) // 2,
            "intersection_one_pair_count": len(intersection_one),
            "candidate_window_size": args.candidate_window,
            "shortlist_size": len(shortlisted),
            "small_prime_bound_exclusive": args.prime_bound + 1,
            "recommended_primes_per_pair": args.recommended_primes,
        },
        "exact_geometry": {
            "rational_bisection_class": "(2,2,w)",
            "intersection_formula": "C_i.C_j=8-<w_i,w_j>",
            "intersection_histogram": dict(
                sorted(intersection_histogram.items(), key=lambda item: int(item[0]))
            ),
            "rationality_argument": (
                "Intersection number one gives one geometric intersection point with "
                "multiplicity one. The two curves are defined over QQ, so the unique point "
                "is Galois fixed. The displayed coordinates independently replay this point."
            ),
        },
        "selection_order": [
            "maximum priority rank",
            "priority-rank sum",
            "quartic maximum coefficient bits",
            "quartic coefficient-bit sum",
            "labels",
        ],
        "affine_recovery_failures_in_candidate_window_before_completion": recovery_failures,
        "pairs": shortlisted,
        "software_assumptions": {"sage": SAGE_VERSION},
        "reproducing_command": (
            "sage -python "
            "elkies-k3/scripts/select_r17_norm12_11952_v4_pair_shortlist.sage"
        ),
        "proof_boundary": (
            "This exactly certifies the displayed rational points and genus-one V4 bases "
            "inside the declared cheapest-prefix/candidate-window search. Jacobian ranks "
            "remain UNKNOWN. The recommended primes only pass exact good-reduction and "
            "complete-chart gates; no modular system is solved here, and no quartic "
            "product-twist section or rank-20 base change is claimed."
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != serialized:
            raise ArithmeticError("stored V4 pair shortlist differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "ALTV4SHORTLIST|"
        f"prefix={args.pool_size}|intersection_one={len(intersection_one)}|"
        f"shortlist={len(shortlisted)}|status={result['status']}|"
        f"output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
