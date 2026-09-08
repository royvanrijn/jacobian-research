#!/usr/bin/env sage -python
"""Blind fixed-rule replay of the 43 exact deepest R17 half-lattice holes.

The parameter list is frozen below.  This search imports no public exceptional
points and uses the same rule on four known high-rank parameters and four
censored ordinary controls: enumerate the exact norm-12 generic holes, choose
the shortest specialized representative in each class, minimize/reduce the
associated pointed quartic, and search to reduced-coordinate height 10^5.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BLIND_ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_r17_matrix_blind_v1.json"
)

# The four positives are calibration parameters, not features.  The four small
# parameters are fixed censored controls; no exact-rank assertion is made for them.
FIBRES = (
    ("r17-control-a", -2, 377, "known_positive_label_withheld_from_search"),
    ("r17-control-b", -308, 251, "known_positive_label_withheld_from_search"),
    ("r17-control-c", 2456, 135, "known_positive_label_withheld_from_search"),
    ("r17-control-d", -9529, 5471, "known_positive_label_withheld_from_search"),
    ("r17-censored-t4", 4, 1, "censored_ordinary_control"),
    ("r17-censored-t6", 6, 1, "censored_ordinary_control"),
    ("r17-censored-t7", 7, 1, "censored_ordinary_control"),
    ("r17-censored-t8", 8, 1, "censored_ordinary_control"),
)

sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def binary_rank(rows) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = int(row)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def restricted_signature_rank(signatures, columns: Sequence[int]) -> int:
    packed = []
    for signature in signatures:
        for row in signature.rows:
            packed.append(
                sum((int(row[column]) & 1) << offset for offset, column in enumerate(columns))
            )
    return binary_rank(packed)


def point_record(point) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid bounded-search budget")

    engine = SourceFileLoader("half_lattice_blind_engine", str(BLIND_ENGINE)).load_module()
    generic_oracle = engine.CosetOracle(engine.GENERIC_GRAM)
    deepest = []
    for mask in range(1 << engine.DIMENSION):
        norm, representative, error = generic_oracle.solve(mask)
        if norm == 12:
            deepest.append((mask, representative))
    if len(deepest) != 43:
        raise ArithmeticError("the exact R17 deepest-hole count changed")

    family = load_q12o5867_data(MODEL, SECTIONS)
    fibre_rows = []
    for fibre_index, (identifier, numerator, denominator, role) in enumerate(FIBRES, 1):
        specialization = evaluate_projective_specialization(family, numerator, denominator)
        minimal_model, minimal_change, unused_metadata = global_minimal_model_with_change(
            specialization.model
        )
        minimal_points = tuple(
            source_point_to_target(point, minimal_change) for point in specialization.points
        )
        short_model, short_change = short_certificate_model(minimal_model)
        generic_points = tuple(
            source_point_to_target(point, short_change) for point in minimal_points
        )
        generic_signatures = find_mod2_reduction_certificate(
            short_model, generic_points, prime_bound=500
        )
        if combined_mod2_rank(generic_signatures, 17) != 17:
            raise ArithmeticError(f"{identifier}: specialized generic subgroup lost independence")
        height_gram = engine.canonical_height_gram(short_model, generic_points)
        rounded_1e5 = tuple(
            tuple(int((value * Decimal(100_000)).to_integral_value()) for value in row)
            for row in height_gram
        )
        rounded_1e6 = tuple(
            tuple(int((value * Decimal(1_000_000)).to_integral_value()) for value in row)
            for row in height_gram
        )
        oracle_1e5 = engine.CosetOracle(rounded_1e5)
        oracle_1e6 = engine.CosetOracle(rounded_1e6)

        cover_records = []
        discovered: dict[tuple[Fraction, Fraction], set[int]] = {}
        representative_disagreements = 0
        for cover_index, (mask, unused_generic_representative) in enumerate(deepest, 1):
            unused_norm5, representative5, unused_error5 = oracle_1e5.solve(mask)
            unused_norm6, representative, unused_error6 = oracle_1e6.solve(mask)
            representative_disagreements += representative5 != representative
            depth = engine.quadratic_decimal(height_gram, representative) / 4
            outcome = engine.run_quartic_search(
                mask=mask,
                representative=representative,
                short_model=short_model,
                generic_points=generic_points,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            for point in outcome.curve_points:
                discovered.setdefault(point, set()).add(mask)
            record = outcome.record
            compact = {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "specialized_representative": list(representative),
                "specialized_depth": str(depth),
                "status": record["status"],
                "wall_seconds": record["wall_seconds"],
                "finite_curve_point_count": len(outcome.curve_points),
                "raw_rational_coefficient_maximum_bits": record[
                    "raw_rational_coefficient_maximum_bits"
                ],
                "integral_model_maximum_coefficient_bits": record[
                    "integral_model_maximum_coefficient_bits"
                ],
            }
            if record["status"] == "bounded_search_complete":
                compact.update(
                    {
                        "reduced_model_maximum_coefficient_bits": record["reduced_model"][
                            "maximum_coefficient_bits"
                        ],
                        "reduced_model_discriminant": record["reduced_model"]["discriminant"],
                        "independent_modular_density_product": record["local_stage"][
                            "joint_independent_density_product"
                        ],
                        "search_milliseconds": record["search_milliseconds"],
                    }
                )
            else:
                compact["error"] = record.get("error")
            cover_records.append(compact)
            print(
                f"R17DEEPMATRIX|fibre={fibre_index}/{len(FIBRES)}|id={identifier}|"
                f"cover={cover_index}/43|mask={mask:#07x}|status={record['status']}|"
                f"points={len(outcome.curve_points)}",
                flush=True,
            )

        basis_with_signs = {
            signed for point in generic_points for signed in (point, (point[0], -point[1]))
        }
        candidates = tuple(
            sorted(
                (point for point in discovered if point not in basis_with_signs),
                key=lambda point: (
                    max(abs(point[0].numerator).bit_length(), point[0].denominator.bit_length()),
                    max(abs(point[1].numerator).bit_length(), point[1].denominator.bit_length()),
                    point,
                ),
            )
        )
        proposals = (
            relation_proposals(
                short_model,
                generic_points,
                candidates,
                timeout=60.0,
                stack_bytes=args.stack_bytes,
            )
            if candidates
            else ()
        )
        unexplained = tuple(
            point for point, (unused_relation, exact) in zip(candidates, proposals) if not exact
        )
        signatures = (
            find_mod2_reduction_certificate(
                short_model, generic_points + unexplained, prime_bound=1200
            )
            if unexplained
            else generic_signatures
        )
        combined_rank = combined_mod2_rank(signatures, 17 + len(unexplained))
        selected_offsets = []
        selected_columns = list(range(17))
        current_rank = restricted_signature_rank(signatures, selected_columns)
        for offset in range(len(unexplained)):
            trial = selected_columns + [17 + offset]
            trial_rank = restricted_signature_rank(signatures, trial)
            if trial_rank > current_rank:
                selected_offsets.append(offset)
                selected_columns = trial
                current_rank = trial_rank
        if current_rank != combined_rank:
            raise ArithmeticError(f"{identifier}: greedy quotient basis lost finite-code rank")

        unexplained_index = {point: index for index, point in enumerate(unexplained)}
        candidate_rows = []
        for point, (relation, exact) in zip(candidates, proposals):
            if exact:
                continue
            offset = unexplained_index[point]
            candidate_rows.append(
                {
                    "point": point_record(point),
                    "source_masks": sorted(discovered[point]),
                    "source_hex": [f"0x{mask:05x}" for mask in sorted(discovered[point])],
                    "selected_for_independent_quotient_basis": offset in selected_offsets,
                }
            )
        fibre_rows.append(
            {
                "id": identifier,
                "role": role,
                "parameter": f"{numerator}/{denominator}",
                "projective_parameter": [numerator, denominator],
                "short_model": [str(value) for value in short_model],
                "specialized_generic_mod2_rank": 17,
                "specialized_generic_certificate_primes": [
                    signature.prime for signature in generic_signatures
                ],
                "representative_disagreements_scale_1e5_vs_1e6": representative_disagreements,
                "cover_records": cover_records,
                "bounded_search_result": {
                    "distinct_nonbasis_candidates": len(candidates),
                    "candidates_unexplained_by_exact_generic_group_law": len(unexplained),
                    "finite_mod2_certified_rank_lower_bound": combined_rank,
                    "finite_mod2_certified_quotient_gain": combined_rank - 17,
                    "selected_independent_candidate_count": len(selected_offsets),
                    "candidate_points": candidate_rows,
                },
            }
        )

    payload: dict[str, Any] = {
        "schema": "elliptic-curves.half-lattice-fake-descent-r17-blind-matrix.v1",
        "status": "PASS_FIXED_GENERIC_DEEPEST43_BLIND_MATRIX",
        "blindness_boundary": {
            "loaded_public_exceptional_points": False,
            "loaded_positive_control_artifact": False,
            "parameter_roles_not_used_by_search": True,
            "fixed_rule": "all 43 exact generic R17 half-lattice classes of norm 12",
        },
        "input_hashes": {
            str(MODEL.relative_to(ROOT)): digest(MODEL),
            str(SECTIONS.relative_to(ROOT)): digest(SECTIONS),
            str(BLIND_ENGINE.relative_to(ROOT)): digest(BLIND_ENGINE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "declared_budget": {
            "class_count_per_fibre": 43,
            "height_bound_each": args.height_bound,
            "timeout_seconds_each": args.timeout_seconds,
            "stack_bytes_each": args.stack_bytes,
            "single_pass_no_retry": True,
        },
        "fibre_count": len(fibre_rows),
        "fibres": fibre_rows,
        "claim_boundary": [
            "Curve membership, group-law exclusion, and finite-reduction independence are exact.",
            "The generic hole list is the exact complete norm-12 R17 CVP stratum.",
            "Specialized representatives use rounded high-precision canonical heights and are numerical lattice evidence.",
            "Every miss is bounded by the declared reduced-coordinate height and timeout.",
            "The ordinary controls are censored search controls, not certified low-rank fibres.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/replay_r17_generic_deep_holes_matrix.sage"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "R17DEEPMATRIX|status=PASS|gains="
        + ",".join(
            str(row["bounded_search_result"]["finite_mod2_certified_quotient_gain"])
            for row in fibre_rows
        )
        + f"|output={args.output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
