#!/usr/bin/env python3
"""Select and parameterize the cheapest bisections spanning ICARM 394.

The optimization is exact and finite.  Among minimum-cardinality subsets of
the 25 covers split at ``t=3/8``, it minimizes maximum equation rank, sum of
equation ranks, maximum priority rank, sum of priority ranks, then mask tuple.
Each selected conic is parameterized through its certified control point.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import shlex
import sys


ROOT = Path(__file__).resolve().parents[2]
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
CONTROLS = ROOT / "artifacts/generated-results/elliptic-curves/elkies_2026_bisection_specialization_controls_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-control-spanning-bisections.json"
T0 = Fraction(3, 8)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational(value: object) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(str(value))


def rational_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def binary_rank(rows) -> int:
    pivots: list[int] = []
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        for pivot in pivots:
            value = min(value, value ^ pivot)
        if value:
            pivots.append(value)
            pivots.sort(reverse=True)
    return len(pivots)


def evaluate(coefficients, value: Fraction) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + rational(coefficient)
    return answer


def multiply(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return result


def add(left, right, scale=Fraction(1)):
    result = [Fraction(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += scale * value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def compose_quadratic(q, numerator, denominator):
    return add(
        add(
            [q[0] * value for value in multiply(denominator, denominator)],
            [q[1] * value for value in multiply(numerator, denominator)],
        ),
        [q[2] * value for value in multiply(numerator, numerator)],
    )


def conic_parameterization(q, root: Fraction) -> dict[str, object]:
    a, b, _c = q[2], q[1], q[0]
    derivative = 2 * a * T0 + b
    denominator = [-a, Fraction(0), Fraction(1)]
    t_numerator = [-a * T0 + derivative, -2 * root, T0]
    u_numerator = [-a * root, derivative, -root]
    identity = add(
        multiply(u_numerator, u_numerator),
        compose_quadratic(q, t_numerator, denominator),
        scale=Fraction(-1),
    )
    if any(identity):
        raise ArithmeticError("the conic parameterization identity failed")
    return {
        "parameter": "s",
        "t_numerator_low_to_high": [rational_text(value) for value in t_numerator],
        "u_numerator_low_to_high": [rational_text(value) for value in u_numerator],
        "common_denominator_low_to_high": [rational_text(value) for value in denominator],
        "formula": (
            "t(s)=T/den and u(s)=U/den for the displayed coefficient vectors; "
            "u(s)^2=q(t(s)) holds coefficientwise"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--controls", type=Path, default=CONTROLS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    sys.set_int_max_str_digits(0)

    batch = json.loads(args.bisections.read_text())
    controls = json.loads(args.controls.read_text())
    target = next(fibre for fibre in controls["fibres"] if fibre["parameter"] == "3/8")
    if len(target["hits"]) != 25 or target["split_bisection_count"] != 25:
        raise ArithmeticError("the ICARM 394 split set changed")
    by_mask = {int(record["lattice_orbit_mask"]): record for record in batch["bisections"]}
    hits = {int(hit["lattice_orbit_mask"]): hit for hit in target["hits"]}
    masks = sorted(hits)
    quotient_rows = {
        mask: tuple(
            int(value)
            for value in hits[mask]["finite_quotient_class_modulo_generic_17"][
                "coordinates_over_f2"
            ]
        )
        for mask in masks
    }
    dimension = binary_rank(quotient_rows.values())
    spanning = [
        subset
        for subset in combinations(masks, dimension)
        if binary_rank(quotient_rows[mask] for mask in subset) == dimension
    ]
    if dimension != 4 or not spanning:
        raise ArithmeticError("the four-dimensional exceptional span changed")

    def complexity(subset):
        equation = [int(by_mask[mask]["equation_rank"]) for mask in subset]
        priority = [int(by_mask[mask]["priority_rank"]) for mask in subset]
        return max(equation), sum(equation), max(priority), sum(priority), tuple(subset)

    selected = min(spanning, key=complexity)
    records = []
    for mask in selected:
        record = by_mask[mask]
        hit = hits[mask]
        q = tuple(rational(value) for value in record["residual_chord"]["q_coefficients"])
        root = rational(hit["canonical_positive_square_root"])
        if root**2 != evaluate(q, T0):
            raise ArithmeticError("a selected control root changed")
        records.append(
            {
                "lattice_orbit_mask": mask,
                "orbit_hex": f"0x{mask:05x}",
                "priority_rank": int(record["priority_rank"]),
                "equation_rank": int(record["equation_rank"]),
                "exceptional_coordinates_over_f2": list(quotient_rows[mask]),
                "q_coefficients_low_to_high": [rational_text(value) for value in q],
                "control_point": {"t": "3/8", "u": rational_text(root)},
                "conic_parameterization": conic_parameterization(q, root),
                "lifted_section": record["lifted_section"],
                "generic_rank_statement": (
                    "After this quadratic base change, the pulled generic R17 subgroup plus the "
                    "displayed anti-invariant lifted section has rank at least 18."
                ),
            }
        )

    result = {
        "schema": "elkies-k3.elkies-2026-control-spanning-bisections.v1",
        "status": "PASS_EXACT_MINIMUM_COMPLEXITY_SPANNING_PACKET",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.bisections): digest(args.bisections),
            display_path(args.controls): digest(args.controls),
        },
        "control": {
            "label": target["label"],
            "t": "3/8",
            "split_bisection_count": 25,
            "exceptional_dimension": dimension,
        },
        "selection": {
            "minimum_cardinality": dimension,
            "spanning_subset_count_at_minimum_cardinality": len(spanning),
            "optimization_order": [
                "maximum equation rank",
                "sum of equation ranks",
                "maximum priority rank",
                "sum of priority ranks",
                "orbit masks",
            ],
            "complexity_key": list(complexity(selected)[:-1]),
            "orbit_masks": list(selected),
            "orbit_hex": [f"0x{mask:05x}" for mask in selected],
        },
        "rank18_families": records,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Minimum cardinality and the declared complexity optimum are exhaustive finite statements. "
            "The conic identities and control roots are exact. Each rank-at-least-18 statement uses the "
            "already certified bisection theorem; no specialization rank or rank upper bound is inferred."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLSPANNING|"
        f"masks={','.join(map(str, selected))}|spanning_subsets={len(spanning)}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
