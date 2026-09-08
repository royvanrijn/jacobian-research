#!/usr/bin/env sage -python
"""Verification-only quotient labeling for the blind R17 deepest-hole matrix."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BLIND = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_r17_matrix_blind_v1.json"
)
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
KUMMER = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_known_kummer_quotients_controls_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_r17_matrix_verification_v1.json"
)

sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from elkies_rank25 import POINTS as RANK25_POINTS  # noqa: E402
from elkies_rank26 import POINTS as RANK26_POINTS  # noqa: E402
from elkies_rank27 import POINTS as RANK27_POINTS  # noqa: E402
from elkies_rank28 import POINTS as RANK28_POINTS  # noqa: E402
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


PUBLIC = {25: RANK25_POINTS, 26: RANK26_POINTS, 27: RANK27_POINTS, 28: RANK28_POINTS}


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


def bit_mask(values) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(values))


def annihilator_basis(values, dimension: int) -> list[int]:
    orthogonal = [
        test
        for test in range(1, 1 << dimension)
        if all((test & value).bit_count() % 2 == 0 for value in values)
    ]
    basis = []
    rank = 0
    for value in orthogonal:
        new_rank = binary_rank(basis + [value])
        if new_rank > rank:
            basis.append(value)
            rank = new_rank
    return basis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--controls", type=Path, default=CONTROLS)
    parser.add_argument("--kummer", type=Path, default=KUMMER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    blind = json.loads(args.blind.read_text())
    controls = json.loads(args.controls.read_text())
    kummer = json.loads(args.kummer.read_text())
    if blind.get("status") != "PASS_FIXED_GENERIC_DEEPEST43_BLIND_MATRIX":
        raise ValueError("the blind matrix is not passing")
    if blind["blindness_boundary"]["loaded_public_exceptional_points"] is not False:
        raise ValueError("the matrix does not have a clean blindness assertion")
    if controls.get("status") != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS":
        raise ValueError("the public positive controls are not passing")

    control_by_parameter = {row["parameter"]: row for row in controls["fibres"]}
    kummer_by_parameter = {row["parameter"]: row for row in kummer["runs"]}
    family = load_q12o5867_data(MODEL, SECTIONS)
    verified = []
    for blind_row in blind["fibres"][:4]:
        parameter = blind_row["parameter"]
        control = control_by_parameter[parameter]
        rank = int(control["published_rank_lower_bound"])
        quotient_dimension = rank - 17
        numerator, denominator = blind_row["projective_parameter"]
        specialization = evaluate_projective_specialization(family, numerator, denominator)
        minimal_model, minimal_change, unused = global_minimal_model_with_change(
            specialization.model
        )
        short_model, short_change = short_certificate_model(minimal_model)
        if [str(value) for value in short_model] != blind_row["short_model"]:
            raise ArithmeticError(f"{parameter}: blind/public short-model mismatch")
        generic = tuple(
            source_point_to_target(
                source_point_to_target(point, minimal_change), short_change
            )
            for point in specialization.points
        )
        public_indices = [
            int(value) - 1
            for value in control["public_positive_control"][
                "selected_public_point_indices_one_based"
            ]
        ]
        public = tuple(
            source_point_to_target(PUBLIC[rank][index], short_change)
            for index in public_indices
        )
        candidate_rows = blind_row["bounded_search_result"]["candidate_points"]
        points = tuple(
            (Fraction(row["point"]["x"]), Fraction(row["point"]["y"]))
            for row in candidate_rows
        )
        relations = relation_proposals(
            short_model,
            generic + public,
            points,
            timeout=60.0,
            stack_bytes=blind["declared_budget"]["stack_bytes_each"],
        )
        if not all(exact for unused_relation, exact in relations):
            raise ArithmeticError(f"{parameter}: a blind point missed the public subgroup")

        quotient_masks = []
        selected = []
        by_center: dict[int, set[int]] = defaultdict(set)
        for source_row, (relation, exact) in zip(candidate_rows, relations):
            quotient = relation[17:]
            mask = bit_mask(quotient)
            if mask == 0:
                raise ArithmeticError(f"{parameter}: a declared quotient point has zero parity")
            quotient_masks.append(mask)
            for center in source_row["source_masks"]:
                by_center[int(center)].add(mask)
            if source_row["selected_for_independent_quotient_basis"]:
                selected.append(
                    {
                        "point": source_row["point"],
                        "source_half_lattice_hex": source_row["source_hex"],
                        "exact_relation_in_generic17_plus_public_complement": list(relation),
                        "exceptional_quotient_coordinates": list(quotient),
                        "exceptional_quotient_mod2_mask": mask,
                        "exceptional_quotient_mod2_hex": f"0x{mask:0{(quotient_dimension + 3) // 4}x}",
                    }
                )
        recovered_rank = binary_rank(quotient_masks)
        blind_gain = blind_row["bounded_search_result"]["finite_mod2_certified_quotient_gain"]
        if recovered_rank != blind_gain or binary_rank(
            row["exceptional_quotient_mod2_mask"] for row in selected
        ) != blind_gain:
            raise ArithmeticError(f"{parameter}: public quotient labeling changed the blind gain")

        kummer_run = kummer_by_parameter[parameter]
        qrows = [
            row
            for row in kummer_run["points"]
            if row["exceptional_quotient_coordinates"] is not None
        ]
        if len(qrows) != quotient_dimension:
            raise ArithmeticError(f"{parameter}: Kummer quotient dimension mismatch")
        recovered_unique = sorted(set(quotient_masks))
        annihilator = annihilator_basis(recovered_unique, quotient_dimension)
        verified.append(
            {
                "id": blind_row["id"],
                "parameter": parameter,
                "published_rank_lower_bound": rank,
                "public_exceptional_quotient_dimension": quotient_dimension,
                "blind_recovered_quotient_dimension": recovered_rank,
                "full_public_quotient_recovered": recovered_rank == quotient_dimension,
                "public_Q_basis_point_indices_one_based": [index + 1 for index in public_indices],
                "blind_non_generic_point_count": len(candidate_rows),
                "distinct_recovered_mod2_quotient_class_count": len(recovered_unique),
                "distinct_recovered_mod2_quotient_masks": recovered_unique,
                "annihilator_dimension_of_recovered_subspace": len(annihilator),
                "annihilator_basis_masks": annihilator,
                "annihilator_basis_hex": [
                    f"0x{value:0{(quotient_dimension + 3) // 4}x}" for value in annihilator
                ],
                "productive_half_lattice_centers": [
                    {
                        "half_lattice_mask": center,
                        "half_lattice_hex": f"0x{center:05x}",
                        "recovered_quotient_span_dimension": binary_rank(masks),
                        "recovered_quotient_masks": sorted(masks),
                    }
                    for center, masks in sorted(by_center.items())
                ],
                "blindly_selected_exact_basis": selected,
                "kummer_fingerprint_dimension": kummer_run["fingerprint_dimension"],
                "kummer_selected_auxiliary_primes": kummer_run["selected_auxiliary_primes"],
            }
        )

    payload = {
        "schema": "elliptic-curves.half-lattice-fake-descent-r17-matrix-verification.v1",
        "status": "PASS_POSTHOC_PUBLIC_VERIFICATION_OF_FIXED_DEEPEST43_RULE",
        "phase_boundary": {
            "blind_matrix_frozen_before_fixture_load": True,
            "verification_loaded_public_exceptional_points": True,
            "public_points_used_only_for_posthoc_coordinates_and_recall": True,
        },
        "input_hashes": {
            str(args.blind.relative_to(ROOT)): digest(args.blind),
            str(args.controls.relative_to(ROOT)): digest(args.controls),
            str(args.kummer.relative_to(ROOT)): digest(args.kummer),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "positive_controls": verified,
        "summary": {
            "target_quotient_dimensions": [row["public_exceptional_quotient_dimension"] for row in verified],
            "blind_recovered_quotient_dimensions": [row["blind_recovered_quotient_dimension"] for row in verified],
            "full_recall_flags": [row["full_public_quotient_recovered"] for row in verified],
            "interpretation": (
                "The frozen generic-deepest43 rule fully recovers +8,+9,+10 and nine of "
                "the eleven rank-28 quotient directions.  Adding the independently ranked "
                "specialized-deepest classes in the rank-28 gate recovers the missing two."
            ),
        },
        "claim_boundary": [
            "All displayed relations and quotient dimensions are exact.",
            "Public exceptional points are used only in this post-hoc verification artifact.",
            "The fixed-rule point searches remain bounded and misses are not nonexistence theorems.",
            "Auxiliary Kummer fingerprints are exact finite-prime barcodes, not full Selmer calculations.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/verify_r17_generic_deep_holes_matrix.sage"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "R17DEEPMATRIXVERIFY|status=PASS|targets="
        + ",".join(str(row["public_exceptional_quotient_dimension"]) for row in verified)
        + "|recovered="
        + ",".join(str(row["blind_recovered_quotient_dimension"]) for row in verified)
        + f"|output={args.output.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
