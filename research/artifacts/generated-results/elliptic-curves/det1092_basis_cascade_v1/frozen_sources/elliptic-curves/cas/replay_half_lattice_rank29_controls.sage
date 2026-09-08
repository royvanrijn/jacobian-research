#!/usr/bin/env sage -python
"""Blind half-lattice replay on three exact rank-29 family controls."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
INPUT = ROOT / "elliptic-curves/data/half_lattice_rank29_control_inputs_v1.json"
ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_rank29_controls_blind_v1.json"
sys.path.insert(0, str(CAS))

from mod2_reduction_independence import combined_mod2_rank, find_mod2_reduction_certificate
from search_nagao_u135_alternate_covers import relation_proposals


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def binary_rank(values) -> int:
    pivots: dict[int, int] = {}
    for value in values:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def restricted_rank(signatures, columns) -> int:
    values = []
    for signature in signatures:
        for row in signature.rows:
            values.append(sum((int(row[column]) & 1) << offset for offset, column in enumerate(columns)))
    return binary_rank(values)


def point_record(point):
    return {"x": str(point[0]), "y": str(point[1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()

    frozen = json.loads(args.input.read_text())
    if frozen.get("status") != "FROZEN_GENERIC_SUBGROUP_ONLY_NO_EXCEPTIONAL_COORDINATES":
        raise ValueError("rank29 generic-only input is not frozen")
    if frozen["boundary"]["output_contains_exceptional_point_coordinates"] is not False:
        raise ValueError("rank29 input leaked an exceptional point")
    engine = SourceFileLoader("half_lattice_rank29_engine", str(ENGINE)).load_module()

    results = []
    for case in frozen["cases"]:
        label = case["label"]
        model = tuple(Fraction(value) for value in case["short_model"])
        generic_points = tuple(
            (Fraction(point[0]), Fraction(point[1])) for point in case["generic_points"]
        )
        generic_gram = tuple(tuple(int(value) for value in row) for row in case["generic_height_gram"])
        signatures0 = find_mod2_reduction_certificate(model, generic_points, prime_bound=1000)
        if combined_mod2_rank(signatures0, 17) != 17:
            raise ArithmeticError(f"{label}: transported generic subgroup lost independence")

        generic_oracle = engine.CosetOracle(generic_gram)
        generic_rows = []
        for mask in range(1 << 17):
            norm, representative, unused_error = generic_oracle.solve(mask)
            generic_rows.append((norm, mask, representative))
        generic_rows.sort(key=lambda row: (-row[0], row[1]))
        generic_top = generic_rows[:43]
        if len(generic_top) != 43 or {row[0] for row in generic_top} != {12}:
            raise ArithmeticError(f"{label}: generic deepest-hole census changed")

        specialized_gram = engine.canonical_height_gram(model, generic_points)
        specialized_runs = {}
        for scale in (100_000, 1_000_000):
            rounded = tuple(
                tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
                for row in specialized_gram
            )
            oracle = engine.CosetOracle(rounded)
            rows = []
            for mask in range(1 << 17):
                unused_norm, representative, unused_error = oracle.solve(mask)
                rows.append((engine.quadratic_decimal(specialized_gram, representative) / 4, mask, representative))
            rows.sort(key=lambda row: (-row[0], row[1]))
            specialized_runs[scale] = rows[:43]
        specialized_top = specialized_runs[1_000_000]
        specialized_top_masks = {row[1] for row in specialized_top}
        generic_top_masks = {row[1] for row in generic_top}
        scale_stable = [row[1] for row in specialized_runs[100_000]] == [row[1] for row in specialized_top]

        specialized_representative = {mask: representative for unused_depth, mask, representative in specialized_top}
        # Recompute representatives for generic-only masks in the specialized metric.
        rounded6 = tuple(
            tuple(int((value * Decimal(1_000_000)).to_integral_value()) for value in row)
            for row in specialized_gram
        )
        specialized_oracle = engine.CosetOracle(rounded6)
        selected_masks = sorted(generic_top_masks | specialized_top_masks)
        discoveries = {}
        cover_records = []
        for position, mask in enumerate(selected_masks, 1):
            representative = specialized_representative.get(mask)
            if representative is None:
                unused_norm, representative, unused_error = specialized_oracle.solve(mask)
            depth = engine.quadratic_decimal(specialized_gram, representative) / 4
            outcome = engine.run_quartic_search(
                mask=mask,
                representative=representative,
                short_model=model,
                generic_points=generic_points,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            for point in outcome.curve_points:
                discoveries.setdefault(point, set()).add(mask)
            record = outcome.record
            cover_records.append(
                {
                    "position_in_union": position,
                    "mask": mask,
                    "hex": f"0x{mask:05x}",
                    "generic_deepest": mask in generic_top_masks,
                    "specialized_top43": mask in specialized_top_masks,
                    "specialized_representative": list(representative),
                    "specialized_depth": str(depth),
                    "status": record["status"],
                    "finite_curve_point_count": len(outcome.curve_points),
                    "integral_coefficient_bits": record["integral_model_maximum_coefficient_bits"],
                    "reduced_coefficient_bits": (
                        record["reduced_model"]["maximum_coefficient_bits"]
                        if record["status"] == "bounded_search_complete" else None
                    ),
                    "modular_density_product": (
                        record["local_stage"]["joint_independent_density_product"]
                        if record["status"] == "bounded_search_complete" else None
                    ),
                    "search_milliseconds": record.get("search_milliseconds"),
                    "wall_seconds": record["wall_seconds"],
                }
            )
            print(
                f"RANK29HALF|case={label}|cover={position}/{len(selected_masks)}|mask={mask:#x}|"
                f"status={record['status']}|points={len(outcome.curve_points)}",
                flush=True,
            )

        basis_signs = {signed for point in generic_points for signed in (point, (point[0], -point[1]))}
        candidates = tuple(
            sorted(
                (point for point in discoveries if point not in basis_signs),
                key=lambda point: (
                    max(abs(point[0].numerator).bit_length(), point[0].denominator.bit_length()),
                    max(abs(point[1].numerator).bit_length(), point[1].denominator.bit_length()),
                    point,
                ),
            )
        )
        proposals = relation_proposals(
            model, generic_points, candidates, timeout=120.0, stack_bytes=args.stack_bytes
        ) if candidates else ()
        unexplained = tuple(point for point, (unused_relation, exact) in zip(candidates, proposals) if not exact)
        signatures = find_mod2_reduction_certificate(
            model, generic_points + unexplained, prime_bound=1800
        ) if unexplained else signatures0
        combined_rank = combined_mod2_rank(signatures, 17 + len(unexplained))
        baseline_rank = restricted_rank(signatures, range(17))
        certificate_valid = baseline_rank == 17 and combined_rank >= 17
        columns = list(range(17))
        current_rank = baseline_rank
        selected_offsets = []
        if certificate_valid:
            for offset in range(len(unexplained)):
                trial = columns + [17 + offset]
                trial_rank = restricted_rank(signatures, trial)
                if trial_rank > current_rank:
                    columns = trial
                    current_rank = trial_rank
                    selected_offsets.append(offset)
        unexplained_index = {point: index for index, point in enumerate(unexplained)}
        candidate_rows = []
        for point, (unused_relation, exact) in zip(candidates, proposals):
            if exact:
                continue
            offset = unexplained_index[point]
            candidate_rows.append(
                {
                    "point": point_record(point),
                    "source_masks": sorted(discoveries[point]),
                    "selected_for_independent_quotient_basis": certificate_valid and offset in selected_offsets,
                }
            )
        results.append(
            {
                "label": label,
                "lineage": case["lineage"],
                "generic_deepest_count": 43,
                "specialized_top_count": 43,
                "generic_specialized_intersection_count": len(generic_top_masks & specialized_top_masks),
                "selected_union_count": len(selected_masks),
                "specialized_top43_order_stable_scale_1e5_vs_1e6": scale_stable,
                "cover_records": cover_records,
                "blind_result": {
                    "distinct_nonbasis_candidates": len(candidates),
                    "unexplained_candidate_count": len(unexplained),
                    "finite_reduction_certificate_valid": certificate_valid,
                    "finite_mod2_rank_lower_bound": combined_rank if certificate_valid else None,
                    "finite_mod2_quotient_gain": combined_rank - 17 if certificate_valid else None,
                    "selected_independent_candidate_count": len(selected_offsets) if certificate_valid else None,
                    "candidate_points": candidate_rows,
                },
            }
        )

    payload = {
        "schema": "elliptic-curves.half-lattice-rank29-controls-blind.v1",
        "status": "PASS_BOUNDED_BLIND_RANK29_CONTROL_SEARCH",
        "blindness_boundary": {
            "search_input": str(args.input.relative_to(ROOT)),
            "search_loaded_exceptional_coordinates": False,
            "verification_fixture": "separate verifier",
        },
        "declared_budget": {
            "generic_deepest_classes": 43,
            "specialized_top_classes": 43,
            "height_bound_each": args.height_bound,
            "timeout_seconds_each": args.timeout_seconds,
            "single_pass_no_retry": True,
        },
        "results": results,
        "claim_boundary": [
            "Generic CVPs use exact integral height forms; specialized CVPs use rounded canonical-height forms.",
            "Every returned curve point is checked exactly, and finite-reduction gains are exact lower bounds when marked valid.",
            "Search misses are bounded and no exact rank upper bound is asserted.",
        ],
        "input_hashes": {
            str(args.input.relative_to(ROOT)): digest(args.input),
            str(ENGINE.relative_to(ROOT)): digest(ENGINE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "reproducing_command": "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_half_lattice_rank29_controls.sage",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "RANK29HALF|status=PASS|gains="
        + ",".join(str(row["blind_result"]["finite_mod2_quotient_gain"]) for row in results)
        + f"|output={args.output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
