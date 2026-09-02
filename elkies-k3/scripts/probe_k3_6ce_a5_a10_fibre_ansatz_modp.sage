#!/usr/bin/env sage-python
"""Exhaust the normalized I6+I11 determinant-384 source chart modulo p."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-fibre-ansatz-mod5-v1.json"
SOURCE_ID = "K3-6ce16abb9de3c7c5-S0008"


def display_path(path):
    return str(path.resolve().relative_to(ROOT))


def local_square_root(unit_coefficients, root0):
    field = root0.parent()
    answer = [field.zero()] * len(unit_coefficients)
    answer[0] = root0
    for degree in range(1, len(answer)):
        known = sum(
            answer[left] * answer[degree - left] for left in range(1, degree)
        )
        answer[degree] = (unit_coefficients[degree] - known) / (2 * root0)
    return answer


def truncated_product(left, right, precision):
    field = left[0].parent()
    answer = [field.zero()] * precision
    for i, left_value in enumerate(left[:precision]):
        for j, right_value in enumerate(right[: precision - i]):
            answer[i + j] += left_value * right_value
    return answer


def multiplicative_branch(a_series, sign=1):
    field = a_series[0].parent()
    unit = [-value / field(3) for value in a_series]
    if not unit[0] or not unit[0].is_square():
        return None
    root0 = unit[0].sqrt()
    if sign == -1:
        root0 = -root0
    root = local_square_root(unit, root0)
    return [
        2 * value
        for value in truncated_product(
            truncated_product(root, root, len(root)), root, len(root)
        )
    ]


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCES)
    parser.add_argument("--source-id", default=SOURCE_ID)
    parser.add_argument("--prime", type=int, default=5)
    parser.add_argument("--examples", type=int, default=0)
    parser.add_argument("--max-samples", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.examples < 0 or arguments.max_samples < 0:
        parser.error("examples and max-samples must be nonnegative")

    source_path = arguments.source.resolve()
    source_payload = json.loads(source_path.read_text())
    matches = [
        row for row in source_payload["sources"] if row["source_id"] == arguments.source_id
    ]
    if len(matches) != 1:
        raise ValueError("expected one selected determinant-384 source")
    source = matches[0]["source"]
    if not (
        source["root_type"] == "A10+A5"
        and source["root_rank"] == 15
        and source["mw_rank_for_rho_19"] == 2
        and source["root_lattice_primitive"]
        and source["torsion"] == 1
    ):
        raise ArithmeticError("selected source invariants changed")

    order_zero, order_infinity = 6, 11
    field = GF(arguments.prime)
    if field.characteristic() in (2, 3):
        raise ValueError("prime must differ from 2 and 3")
    ring = PolynomialRing(field, "t")
    t = ring.gen()

    rows = []
    for jet in range(order_zero):
        rows.append([field(index == jet) for index in range(13)])
    for jet in range(order_infinity):
        rows.append([field(index == 12 - jet) for index in range(13)])
    hermite = matrix(field, rows)
    if hermite.ncols() != hermite.rank() or hermite.rank() != 13:
        raise ArithmeticError("two-support Hermite matrix lost full column rank")
    compatibility_matrix = hermite.left_kernel().basis_matrix()
    compatibility_equations = order_zero + order_infinity - 13
    if compatibility_matrix.nrows() != compatibility_equations:
        raise ArithmeticError("two-support compatibility rank changed")

    examples = []
    branch_eligible = compatible = exact_orders = squarefree = 0
    compatible_order_histogram = Counter()
    sample = 0
    exhaustive_total = arguments.prime**8
    for digits in itertools.product(range(arguments.prime), repeat=8):
        sample += 1
        if arguments.max_samples and sample > arguments.max_samples:
            sample -= 1
            break
        a_coefficients = [field(-3)] + [field(value) for value in digits]
        if not a_coefficients[8]:
            continue
        A = ring(a_coefficients)
        at_zero = a_coefficients[:order_zero]
        at_infinity = [
            a_coefficients[8 - jet] if 0 <= 8 - jet <= 8 else field.zero()
            for jet in range(order_infinity)
        ]
        zero_branch = multiplicative_branch(at_zero, 1)
        infinity_branch = multiplicative_branch(at_infinity, 1)
        if zero_branch is None or infinity_branch is None:
            continue
        for sign_infinity in (1, -1):
            branch_eligible += 1
            target = vector(
                field,
                zero_branch + [sign_infinity * value for value in infinity_branch],
            )
            if compatibility_matrix * target:
                continue
            b_coefficients = list(hermite.solve_right(target))
            compatible += 1
            B = ring(b_coefficients)
            discriminant_core = 4 * A**3 + 27 * B**2
            orders = (
                order_at(discriminant_core, field.zero()),
                24 - discriminant_core.degree(),
            )
            compatible_order_histogram[orders] += 1
            if orders != (order_zero, order_infinity):
                continue
            exact_orders += 1
            residual, remainder = discriminant_core.quo_rem(t**order_zero)
            residual_i1_count = 24 - order_zero - order_infinity
            if remainder or residual.degree() != residual_i1_count or residual(0) == 0:
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            squarefree += 1
            if not arguments.examples or len(examples) < arguments.examples:
                examples.append(
                    {
                        "sample_index": sample,
                        "branch_sign_at_infinity": sign_infinity,
                        "A_coefficients_low_to_high": [int(value) for value in a_coefficients],
                        "B_coefficients_low_to_high": [int(value) for value in b_coefficients],
                        "discriminant_orders": {
                            "0": order_zero,
                            "infinity": order_infinity,
                        },
                        "residual_discriminant_coefficients_low_to_high": [
                            int(value) for value in residual
                        ],
                        "residual_factorization": [
                            {"degree": int(factor.degree()), "multiplicity": int(power)}
                            for factor, power in residual.factor()
                        ],
                        "geometric_fibre_profile": "I6+I11+7I1",
                    }
                )

    exhausted = sample == exhaustive_total and not arguments.max_samples
    payload = {
        "schema": "elkies-k3.k3-6ce-a5-a10-mw2-fibre-ansatz-modp.v1",
        "status": (
            "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ"
            if exhausted and squarefree
            else (
                "PASS_BOUNDED_MODULAR_SOURCE_FIBRE_ANSATZ"
                if squarefree
                else "PASS_BOUNDED_NO_MODULAR_SOURCE_FIBRE_ANSATZ"
            )
        ),
        "prime": arguments.prime,
        "scan": {
            "normalized_A_polynomials": exhaustive_total,
            "samples_consumed": sample,
            "exhausted": exhausted,
        },
        "accounting": {
            "branch_eligible_with_signs": branch_eligible,
            "hermite_compatible_with_signs": compatible,
            "exact_prescribed_orders": exact_orders,
            "squarefree_examples_with_signs": squarefree,
            "stored_examples": len(examples),
            "compatible_discriminant_order_histogram": {
                ",".join("null" if value is None else str(value) for value in key): count
                for key, count in sorted(
                    compatible_order_histogram.items(),
                    key=lambda item: tuple(
                        -1 if value is None else value for value in item[0]
                    ),
                )
            },
        },
        "ansatz": {
            "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
            "degree_bounds": {"A": 8, "B": 12},
            "normalization": "A(0)=-3; reducible supports at 0 and infinity",
            "normalized_reducible_supports": ["0:I6", "infinity:I11"],
            "hermite_conditions": 17,
            "B_coefficient_rank": 13,
            "compatibility_equations_on_A": compatibility_equations,
            "expected_fibre_stratum_dimension": 8 - compatibility_equations,
            "expected_K3_6ce_MW2_locus_dimension": 1,
            "expected_MW_conditions_still_missing": 2,
            "section_marking": {
                "status": "NOT_IMPOSED_AT_FIBRE_GATE",
                "basis_pole_profile": [0, 1],
                "height_gram": source["mw_height_gram"],
            },
        },
        "examples": examples,
        "source": {
            "artifact": display_path(source_path),
            "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "source_id": arguments.source_id,
            "source_gram_sha256": source["gram_sha256"],
            "surface_id": "K3-6ce16abb9de3c7c5",
            "determinant": 384,
            "root_type": source["root_type"],
            "mw_height_gram": source["mw_height_gram"],
        },
        "proof_boundary": {
            "proved": (
                "Every stored example is an exact short-Weierstrass K3 model over "
                "the displayed finite field with fibre profile I6+I11+7I1."
            ),
            "not_proved": (
                "The two MW sections, full lattice marking, rational parameterization, "
                "characteristic-zero lifting, arithmetic moduli curve, and neighbour "
                "route are not proved."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/probe_k3_6ce_a5_a10_fibre_ansatz_modp.sage "
            f"--source-id {arguments.source_id} --prime {arguments.prime} "
            f"--examples {arguments.examples} --output {display_path(arguments.output)}"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "K36CEA5A10FIBRE|"
        f"p={arguments.prime}|samples={sample}|compatible={compatible}|"
        f"squarefree={squarefree}|stored={len(examples)}|exhausted={int(exhausted)}|"
        f"status={'PASS' if squarefree else 'BOUNDED_NEGATIVE'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
