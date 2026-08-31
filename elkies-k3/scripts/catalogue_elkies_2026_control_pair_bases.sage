#!/usr/bin/env sage-python
"""Catalogue the 300 pair bases selected by the 25-split t=3/8 fibre.

status: ACTIVE_SEARCH
claim: exact control-selected pair-base Jacobian catalogue
inputs: complete bisection batch and exact control-specialization certificate
outputs: artifacts/generated-results/elkies-2026-control-pair-base-catalogue.json
supersedes: none; complements the t=0/infinity immediate-point catalogue

Every pair among S(3/8) has an exact rational point on its V4 base.  The
associated third quotient is the pointed quartic w^2=q_i(t)q_j(t); its
Jacobian is 2-isogenous to the pair-base Jacobian and has the same rank.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import shlex
import sys
from time import perf_counter

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
CONTROLS = ROOT / "artifacts/generated-results/elliptic-curves/elkies_2026_bisection_specialization_controls_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-control-pair-base-catalogue.json"
OUTPUT_SCHEMA = "elkies-k3.elkies-2026-control-pair-base-catalogue.v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def evaluate(coefficients, value):
    answer = QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + QQ(coefficient)
    return answer


def rational_bits(value) -> int:
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def binary_rank(rows) -> int:
    """Return the rank of short binary rows without a CAS-dependent matrix."""

    pivots = []
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        for pivot in pivots:
            value = min(value, value ^ pivot)
        if value:
            pivots.append(value)
            pivots.sort(reverse=True)
    return len(pivots)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--controls", type=Path, default=CONTROLS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    started = perf_counter()

    batch = json.loads(args.bisections.read_text())
    controls = json.loads(args.controls.read_text())
    if batch.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("unexpected bisection schema")
    if controls.get("schema") != "elliptic-curves.elkies-2026-bisection-specialization-controls.v1":
        raise ValueError(f"unexpected control schema: {controls.get('schema')!r}")
    target = next(fibre for fibre in controls["fibres"] if fibre["parameter"] == "3/8")
    if target["split_bisection_count"] != 25 or len(target["hits"]) != 25:
        raise ArithmeticError("the t=3/8 split set changed")
    by_mask = {int(record["lattice_orbit_mask"]): record for record in batch["bisections"]}
    masks = sorted(int(hit["lattice_orbit_mask"]) for hit in target["hits"])
    hit_by_mask = {int(hit["lattice_orbit_mask"]): hit for hit in target["hits"]}
    incidence = {
        fibre["label"]: {int(hit["lattice_orbit_mask"]): hit for hit in fibre["hits"]}
        for fibre in controls["fibres"]
    }
    control_parameters = {fibre["label"]: QQ(fibre["parameter"]) for fibre in controls["fibres"]}

    quotient_rows = {
        mask: tuple(
            int(value)
            for value in hit_by_mask[mask]["finite_quotient_class_modulo_generic_17"][
                "coordinates_over_f2"
            ]
        )
        for mask in masks
    }
    quotient_dimension = binary_rank(quotient_rows.values())
    spanning_subsets = [
        subset
        for subset in combinations(masks, quotient_dimension)
        if binary_rank(quotient_rows[mask] for mask in subset) == quotient_dimension
    ]
    if quotient_dimension != 4 or not spanning_subsets:
        raise ArithmeticError("the four-dimensional t=3/8 quotient span changed")

    def spanning_complexity(subset):
        equation_ranks = [int(by_mask[mask]["equation_rank"]) for mask in subset]
        priority_ranks = [int(by_mask[mask]["priority_rank"]) for mask in subset]
        return (
            max(equation_ranks),
            sum(equation_ranks),
            max(priority_ranks),
            sum(priority_ranks),
            tuple(subset),
        )

    minimum_spanning_subset = min(spanning_subsets, key=spanning_complexity)

    ring = PolynomialRing(QQ, "t")
    rows = []
    root_numbers = Counter()
    for left_mask, right_mask in combinations(masks, 2):
        left = by_mask[left_mask]
        right = by_mask[right_mask]
        left_q = ring([QQ(value) for value in left["residual_chord"]["q_coefficients"]])
        right_q = ring([QQ(value) for value in right["residual_chord"]["q_coefficients"]])
        if left_q.degree() != 2 or right_q.degree() != 2 or left_q.gcd(right_q) != 1:
            raise ArithmeticError(f"pair {left_mask}:{right_mask} is not a smooth four-branch quartic")
        quartic = left_q * right_q
        e, d, c, b, a = (QQ(quartic[index]) for index in range(5))
        invariant_i = 12 * a * e - 3 * b * d + c**2
        invariant_j = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
        raw = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
        minimal = raw.global_minimal_model()
        minimal_model = tuple(ZZ(value) for value in minimal.a_invariants())
        root_number = int(minimal.root_number())
        root_numbers[root_number] += 1

        common_controls = []
        for label, hits in incidence.items():
            if left_mask not in hits or right_mask not in hits:
                continue
            parameter = control_parameters[label]
            left_root = QQ(hits[left_mask]["canonical_positive_square_root"])
            right_root = QQ(hits[right_mask]["canonical_positive_square_root"])
            if left_root**2 != left_q(parameter) or right_root**2 != right_q(parameter):
                raise ArithmeticError("control square roots no longer match the batch equations")
            common_controls.append(
                {
                    "label": label,
                    "t": rational_text(parameter),
                    "u": rational_text(left_root),
                    "v": rational_text(right_root),
                    "third_quotient_w": rational_text(left_root * right_root),
                }
            )
        if not any(point["t"] == "3/8" for point in common_controls):
            raise ArithmeticError("a selected pair lost its defining t=3/8 point")

        boundary_points = []
        for label, index in (("zero", 0), ("infinity", 2)):
            left_value = QQ(left_q[index])
            right_value = QQ(right_q[index])
            if left_value.is_square() and right_value.is_square():
                boundary_points.append(
                    {
                        "label": label,
                        "left_root": rational_text(left_value.sqrt()),
                        "right_root": rational_text(right_value.sqrt()),
                    }
                )
        row = {
            "pair_key": f"{left_mask}:{right_mask}",
            "orbit_masks": [left_mask, right_mask],
            "orbit_hex": [f"0x{left_mask:05x}", f"0x{right_mask:05x}"],
            "priority_ranks": sorted((int(left["priority_rank"]), int(right["priority_rank"]))),
            "equation_ranks": sorted((int(left["equation_rank"]), int(right["equation_rank"]))),
            "rational_point_sources": [point["label"] for point in common_controls] + [
                point["label"] for point in boundary_points
            ],
            "common_control_points": common_controls,
            "boundary_points": boundary_points,
            "q_coefficients_low_to_high": [
                [rational_text(left_q[index]) for index in range(3)],
                [rational_text(right_q[index]) for index in range(3)],
            ],
            "quartic_coefficients_low_to_high": [rational_text(quartic[index]) for index in range(5)],
            "binary_quartic_I": rational_text(invariant_i),
            "binary_quartic_J": rational_text(invariant_j),
            "minimal_jacobian_a1_a2_a3_a4_a6": [int(value) for value in minimal_model],
            "conductor": int(minimal.conductor()),
            "global_root_number": root_number,
            "complexity": {
                "minimal_model_max_coefficient_bits": max(abs(value).nbits() for value in minimal_model),
                "quartic_max_rational_bits": max(rational_bits(quartic[index]) for index in range(5)),
                "priority_rank_max": max(int(left["priority_rank"]), int(right["priority_rank"])),
                "priority_rank_sum": int(left["priority_rank"]) + int(right["priority_rank"]),
            },
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            row["complexity"]["minimal_model_max_coefficient_bits"],
            row["complexity"]["quartic_max_rational_bits"],
            row["complexity"]["priority_rank_max"],
            row["complexity"]["priority_rank_sum"],
            row["orbit_masks"],
        )
    )
    for rank, row in enumerate(rows, start=1):
        row["arithmetic_complexity_rank"] = rank
    if len(rows) != 300:
        raise ArithmeticError(f"expected 300 pairs, got {len(rows)}")

    repeated = [row for row in rows if len(row["rational_point_sources"]) > 1]
    result = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_EXACT_CONTROL_SELECTED_PAIR_BASE_CATALOGUE",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.bisections): digest(args.bisections),
            display_path(args.controls): digest(args.controls),
        },
        "selection": {
            "control_label": target["label"],
            "parameter": target["parameter"],
            "split_mask_count": len(masks),
            "split_masks": masks,
            "pair_count": len(rows),
            "minimum_complexity_exceptional_spanning_subset": {
                "optimization_order": [
                    "maximum equation rank",
                    "sum of equation ranks",
                    "maximum priority rank",
                    "sum of priority ranks",
                    "orbit masks",
                ],
                "minimum_cardinality": quotient_dimension,
                "full_quotient_dimension": quotient_dimension,
                "spanning_subset_count_at_minimum_cardinality": len(spanning_subsets),
                "orbit_masks": list(minimum_spanning_subset),
                "orbit_hex": [f"0x{mask:05x}" for mask in minimum_spanning_subset],
                "quotient_coordinates_over_f2": [
                    list(quotient_rows[mask]) for mask in minimum_spanning_subset
                ],
                "equation_ranks": [
                    int(by_mask[mask]["equation_rank"]) for mask in minimum_spanning_subset
                ],
                "priority_ranks": [
                    int(by_mask[mask]["priority_rank"]) for mask in minimum_spanning_subset
                ],
                "complexity_key": list(spanning_complexity(minimum_spanning_subset)[:-1]),
            },
        },
        "summary": {
            "root_number_counts": {str(key): value for key, value in sorted(root_numbers.items())},
            "pairs_with_an_additional_known_base_point": len(repeated),
            "pairs_with_multiple_point_sources": [row["pair_key"] for row in repeated],
        },
        "pairs": rows,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The 300 pair bases, their displayed rational points, quartic invariants, minimal Jacobians, "
            "conductors, and root numbers are exact. Positive Mordell-Weil rank is not inferred here; it "
            "requires a separately certified non-torsion point on the Jacobian."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLPAIRCATALOGUE|"
        f"pairs={len(rows)}|repeated_point_sources={len(repeated)}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
