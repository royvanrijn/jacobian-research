#!/usr/bin/env sage -python
"""Verification-only labeling of the frozen blind rank-28 replay.

This program is intentionally separate from ``half_lattice_fake_descent_replay.sage``.
Only this second pass opens the public exceptional-point fixture.  It proves
that every blind non-generic point lies in the certified public rank-28
subgroup, labels its quotient coordinates, and checks that the eleven points
selected blindly span the full eleven-dimensional public complement.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
BLIND_SCRIPT = CAS / "half_lattice_fake_descent_replay.sage"
BLIND_RESULT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_rank28_blind_v1.json"
)
PUBLIC_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
KUMMER_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_kummer_quotients_controls_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_rank28_verification_v1.json"
)

sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from build_elkies_2026_rank28_relative_descent_magma import (  # noqa: E402
    load_relative_input,
)
from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def point_from_record(row: dict[str, str]) -> tuple[Fraction, Fraction]:
    return Fraction(row["x"]), Fraction(row["y"])


def binary_rank(values: Iterable[int]) -> int:
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


def bit_mask(values: Sequence[int]) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(values))


def xor_rows(rows: Sequence[Sequence[int]], mask: int) -> list[int]:
    if not rows:
        return []
    answer = [0] * len(rows[0])
    for index, row in enumerate(rows):
        if (mask >> index) & 1:
            answer = [left ^ int(right) for left, right in zip(answer, row)]
    return answer


def median(values: Sequence[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) & 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def load_kummer_basis() -> tuple[list[list[int]], dict[str, Any]]:
    artifact = json.loads(KUMMER_CONTROLS.read_text())
    matches = [row for row in artifact["runs"] if row["parameter"] == "-9529/5471"]
    if len(matches) != 1:
        raise ArithmeticError("the Kummer controls have no unique rank-28 row")
    run = matches[0]
    points = [
        row for row in run["points"] if row["exceptional_quotient_coordinates"] is not None
    ]
    if len(points) != 11:
        raise ArithmeticError("the public Kummer complement no longer has eleven rows")
    for index, row in enumerate(points):
        expected = [0] * 11
        expected[index] = 1
        if row["label"] != f"Q{index + 1}" or row["exceptional_quotient_coordinates"] != expected:
            raise ArithmeticError("the Kummer rows are not in the public Q1,...,Q11 order")
    basis = [list(map(int, row["local_squareclass_row"])) for row in points]
    return basis, {
        "fingerprint_dimension": run["fingerprint_dimension"],
        "selected_auxiliary_primes": run["selected_auxiliary_primes"],
        "known_exceptional_quotient_dimension": run["known_exceptional_quotient_dimension"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind-result", type=Path, default=BLIND_RESULT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()

    blind = json.loads(args.blind_result.read_text())
    if blind.get("schema") != "elliptic-curves.half-lattice-fake-descent-rank28-blind.v1":
        raise ValueError("unexpected blind replay schema")
    if blind.get("status") != "PASS_BLIND_BOUNDED_RANK28_HALF_LATTICE_REPLAY":
        raise ValueError("the blind replay is not a completed search")
    if blind["blindness_boundary"]["search_loaded_exceptional_point_fixture"] is not False:
        raise ValueError("the input artifact does not assert a clean blindness boundary")
    if blind["input_hashes"][str(BLIND_SCRIPT.relative_to(ROOT))] != digest(BLIND_SCRIPT):
        raise ValueError("the blind search script changed after the artifact was frozen")

    # This is the first point at which public exceptional data are loaded.
    source = load_relative_input(PUBLIC_CONTROLS)
    short_model, short_change = short_certificate_model(source.model)
    public_short = tuple(
        source_point_to_target(point, short_change) for point in source.public_complement
    )
    generic_short = tuple(
        point_from_record(row) for row in blind["fibre"]["generic_points"]
    )
    if generic_short != tuple(
        source_point_to_target(point, short_change) for point in source.generic_points
    ):
        raise ArithmeticError("the blind and verification generic bases differ")
    if [str(value) for value in short_model] != blind["fibre"]["short_model"]:
        raise ArithmeticError("the blind and verification short models differ")

    candidate_source_rows = [
        row
        for row in blind["blind_results"]["candidate_points"]
        if not row["exact_relation_in_generic_subgroup"]
    ]
    candidate_points = tuple(point_from_record(row["point"]) for row in candidate_source_rows)
    relations = relation_proposals(
        short_model,
        generic_short + public_short,
        candidate_points,
        timeout=args.timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if not all(exact for unused_relation, exact in relations):
        raise ArithmeticError("a blind candidate is outside the certified public rank-28 subgroup")

    kummer_basis, kummer_metadata = load_kummer_basis()
    labeled = []
    cover_quotients: dict[int, set[int]] = defaultdict(set)
    selected_quotients = []
    for source_row, point, (relation, exact) in zip(
        candidate_source_rows, candidate_points, relations
    ):
        quotient_coordinates = relation[17:]
        quotient_mask = bit_mask(quotient_coordinates)
        if quotient_mask == 0:
            raise ArithmeticError("a non-generic blind point has zero public quotient parity")
        barcode = xor_rows(kummer_basis, quotient_mask)
        for center_mask in source_row["source_masks"]:
            cover_quotients[int(center_mask)].add(quotient_mask)
        if source_row["selected_for_blind_independent_quotient_basis"]:
            selected_quotients.append(quotient_mask)
        labeled.append(
            {
                "point": source_row["point"],
                "source_half_lattice_masks": source_row["source_masks"],
                "source_half_lattice_hex": source_row["source_hex"],
                "exact_relation_in_generic17_plus_public11": list(relation),
                "exceptional_quotient_coordinates": list(quotient_coordinates),
                "exceptional_quotient_mod2_mask": quotient_mask,
                "exceptional_quotient_mod2_hex": f"0x{quotient_mask:03x}",
                "auxiliary_local_kummer_barcode": barcode,
                "selected_blindly_for_independent_basis": source_row[
                    "selected_for_blind_independent_quotient_basis"
                ],
            }
        )

    all_quotients = [row["exceptional_quotient_mod2_mask"] for row in labeled]
    if binary_rank(all_quotients) != 11 or binary_rank(selected_quotients) != 11:
        raise ArithmeticError("the blind quotient does not equal the public 11-space")

    class_rows = {row["mask"]: row for row in blind["ranking"]["selected_classes"]}
    search_rows = {row["mask"]: row for row in blind["search_records"]}
    ordered_masks = blind["ranking"]["selected_masks"]
    prefix_vectors: list[int] = []
    prefix_rank = 0
    per_cover = []
    first_full_rank_at = None
    for mask in ordered_masks:
        vectors = sorted(cover_quotients.get(mask, ()))
        before = prefix_rank
        prefix_vectors.extend(vectors)
        prefix_rank = binary_rank(prefix_vectors)
        if vectors:
            row = class_rows[mask]
            per_cover.append(
                {
                    "half_lattice_mask": mask,
                    "half_lattice_hex": f"0x{mask:05x}",
                    "generic_depth": row["generic_depth"],
                    "generic_rank_with_mask_tiebreak": row[
                        "generic_rank_with_mask_tiebreak"
                    ],
                    "specialized_depth": row["specialized_depth"],
                    "specialized_rank": row["specialized_rank"],
                    "recovered_nonzero_quotient_masks": vectors,
                    "recovered_nonzero_quotient_hex": [
                        f"0x{value:03x}" for value in vectors
                    ],
                    "quotient_span_dimension_within_cover": binary_rank(vectors),
                    "prefix_quotient_rank_before": before,
                    "prefix_quotient_rank_after": prefix_rank,
                    "incremental_quotient_gain": prefix_rank - before,
                    "integral_model_maximum_coefficient_bits": search_rows[mask][
                        "integral_model_maximum_coefficient_bits"
                    ],
                    "reduced_model_maximum_coefficient_bits": search_rows[mask][
                        "reduced_model"
                    ]["maximum_coefficient_bits"],
                    "search_milliseconds": search_rows[mask]["search_milliseconds"],
                }
            )
        if prefix_rank == 11 and first_full_rank_at is None:
            first_full_rank_at = {
                "selected_union_position": ordered_masks.index(mask) + 1,
                "specialized_rank": class_rows[mask]["specialized_rank"],
                "half_lattice_hex": f"0x{mask:05x}",
            }

    productive_masks = sorted(cover_quotients)
    blind_module = SourceFileLoader("blind_half_lattice_replay", str(BLIND_SCRIPT)).load_module()
    generic_oracle = blind_module.CosetOracle(blind_module.GENERIC_GRAM)
    xor_depth_histogram: Counter[str] = Counter()
    xor_masks = set()
    for left_index, left in enumerate(productive_masks):
        for right in productive_masks[left_index + 1 :]:
            difference = left ^ right
            xor_masks.add(difference)
            norm, unused_representative, unused_error = generic_oracle.solve(difference)
            xor_depth_histogram[str(Fraction(norm, 4))] += 1

    productive_set = set(productive_masks)
    groups: dict[str, list[dict[str, Any]]] = {
        "productive_exceptional": [],
        "searched_not_exceptional": [],
    }
    for mask in ordered_masks:
        groups[
            "productive_exceptional" if mask in productive_set else "searched_not_exceptional"
        ].append(search_rows[mask])
    contrasts = {}
    for label, rows in groups.items():
        reduced_bits = [row["reduced_model"]["maximum_coefficient_bits"] for row in rows]
        density = [
            float(row["local_stage"]["joint_independent_density_product"])
            for row in rows
        ]
        contrasts[label] = {
            "class_count": len(rows),
            "generic_exact_deepest_count": sum(row["generic_depth"] == "3" for row in rows),
            "median_reduced_model_maximum_coefficient_bits": median(reduced_bits),
            "median_independent_modular_density_product": median(density),
        }

    payload = {
        "schema": "elliptic-curves.half-lattice-fake-descent-rank28-verification.v1",
        "status": "PASS_POSTHOC_PUBLIC_FIXTURE_VERIFICATION_OF_BLIND_GAIN_11",
        "phase_boundary": {
            "blind_search_artifact_frozen_before_fixture_load": True,
            "blind_search_loaded_public_exceptional_points": False,
            "verification_loaded_public_exceptional_points": True,
            "verification_role": "post-hoc labeling and completeness check only",
        },
        "input_hashes": {
            str(args.blind_result.relative_to(ROOT)): digest(args.blind_result),
            str(PUBLIC_CONTROLS.relative_to(ROOT)): digest(PUBLIC_CONTROLS),
            str(KUMMER_CONTROLS.relative_to(ROOT)): digest(KUMMER_CONTROLS),
            str((CAS / "elkies_rank28.py").relative_to(ROOT)): digest(CAS / "elkies_rank28.py"),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "exact_result": {
            "blind_non_generic_candidate_count": len(labeled),
            "all_candidates_exactly_in_public_rank28_group": True,
            "all_recovered_quotient_span_dimension": binary_rank(all_quotients),
            "blindly_selected_basis_count": len(selected_quotients),
            "blindly_selected_basis_quotient_dimension": binary_rank(selected_quotients),
            "public_exceptional_quotient_dimension": 11,
            "blind_search_recovers_full_public_exceptional_quotient": True,
            "first_full_quotient_rank_in_specialized_order": first_full_rank_at,
        },
        "half_lattice_class_summary": {
            "productive_center_count": len(productive_masks),
            "productive_centers": per_cover,
            "productive_center_f2_span_dimension": binary_rank(productive_masks),
            "distinct_pairwise_xor_count": len(xor_masks),
            "pairwise_xor_generic_depth_histogram": dict(sorted(xor_depth_histogram.items())),
        },
        "kummer_fingerprint": {
            **kummer_metadata,
            "interpretation": (
                "exact auxiliary-prime squareclass fingerprints inherited linearly from the "
                "certified public Q1,...,Q11 basis; not a complete Selmer computation"
            ),
        },
        "labeled_blind_points": labeled,
        "bounded_search_contrasts_within_selected_union": contrasts,
        "interpretation_boundary": [
            "The exact 11-dimensional quotient recovery and every displayed group-law relation are theorem-level certificates.",
            "The half-lattice center depths use the exact generic form and a high-precision rounded specialized height form as stated in the blind artifact.",
            "The hit rates and complexity contrasts are bounded-search observations conditioned on the predeclared 64-class union.",
            "The fake-descent quartics are pointed models birational to E; their local solubility is automatic and cannot explain rank-jump simultaneity.",
            "Distinct productive half-lattice centers are search charts, not distinct Selmer torsors unless a separate genuine 2-descent identifies them as such.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/verify_half_lattice_fake_descent_replay.sage"
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing output: {args.output}")
        print(f"HALFLATTICEVERIFY|status=PASS_CHECK|output={args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        "HALFLATTICEVERIFY|status=PASS|blind_gain=11|public_gain=11|"
        f"productive_centers={len(productive_masks)}|output={args.output.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
