#!/usr/bin/env sage-python
"""Freeze the rank-blind small-field R17 class-quotient cohort.

The candidate universe is every reduced finite parameter ``a/b`` of projective
height at most 24 in the published R17 coordinate.  Selection uses only exact
structural data: nonsingularity, irreducibility of the completed-square
2-division cubic, distinct rational j-invariant, and the absolute discriminant
of that explicit cubic order.  No point search, rank label, Nagao score, or
exceptional-point coordinate is loaded.

The output deliberately contains neither class-group features nor detector
outcomes.  Those are separate, ordered phases of the laboratory.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from typing import Any

from sage.all import PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
MW16 = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_curve398_hidden_a1_mw16_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-small-field-class-quotient-cohort-v1.json"
)

SCHEMA = "elkies-k3.r17-small-field-class-quotient-cohort.v1"
STATUS = "FROZEN_RANK_BLIND_PRE_CLASS_GROUP_COHORT"
SALT = "r17-small-field-class-quotient-cohort-v1"
HEIGHT_BOUND = 24
COHORT_SIZE = 100


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def homogeneous_value(coefficients: list[int], numerator: int, denominator: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(coefficients)
    )


def selection_record(row: dict[str, Any]) -> dict[str, Any]:
    return {
        key: row[key]
        for key in (
            "sample_id",
            "family",
            "parameter",
            "projective_pair",
            "projective_height",
            "cubic_order_discriminant",
            "absolute_cubic_order_discriminant_bits",
            "j_invariant",
        )
    }


def build(*, height_bound: int = HEIGHT_BOUND, cohort_size: int = COHORT_SIZE) -> dict[str, Any]:
    if height_bound < 1:
        raise ValueError("the projective-height bound must be positive")
    if not 100 <= cohort_size <= 500:
        raise ValueError("the committed laboratory size must lie between 100 and 500")

    model = json.loads(MODEL.read_text())
    sections = json.loads(SECTIONS.read_text())
    mw16 = json.loads(MW16.read_text())
    if model.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
        raise ArithmeticError("the published R17 model is not certified")
    if sections.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS":
        raise ArithmeticError("the published R17 section input is not certified")
    if len(sections.get("sections", [])) != 17:
        raise ArithmeticError("the published generic basis stopped having rank 17")
    if mw16.get("status") != "PASS_EXACT_HIDDEN_A1_MW16_FIBRATION_PARAMETER_AND_PUBLIC_SUBGROUP":
        raise ArithmeticError("the recovered MW16 expansion input is not certified")

    coefficients_a = [int(value) for value in model["A_coefficients_low_to_high"]]
    coefficients_b = [int(value) for value in model["B_coefficients_low_to_high"]]
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    candidates = []
    structural_rejections = {
        "singular_specialization": 0,
        "reducible_two_division_cubic": 0,
    }
    universe_count = 0
    for denominator in range(1, height_bound + 1):
        for numerator in range(-height_bound, height_bound + 1):
            if max(abs(numerator), denominator) > height_bound:
                continue
            if gcd(abs(numerator), denominator) != 1:
                continue
            universe_count += 1
            coefficient_a = homogeneous_value(coefficients_a, numerator, denominator)
            coefficient_b = homogeneous_value(coefficients_b, numerator, denominator)
            discriminant = -256 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
            if discriminant == 0:
                structural_rejections["singular_specialization"] += 1
                continue
            polynomial = z**3 + 16 * coefficient_a * z + 64 * coefficient_b
            if not polynomial.is_irreducible():
                structural_rejections["reducible_two_division_cubic"] += 1
                continue
            j_invariant = Fraction(
                6912 * coefficient_a**3,
                4 * coefficient_a**3 + 27 * coefficient_b**2,
            )
            candidates.append(
                {
                    "numerator": numerator,
                    "denominator": denominator,
                    "projective_height": max(abs(numerator), denominator),
                    "coefficient_a": coefficient_a,
                    "coefficient_b": coefficient_b,
                    "cubic_order_discriminant": discriminant,
                    "absolute_cubic_order_discriminant_bits": abs(discriminant).bit_length(),
                    "j_invariant": j_invariant,
                }
            )

    candidates.sort(
        key=lambda row: (
            abs(row["cubic_order_discriminant"]),
            row["projective_height"],
            row["denominator"],
            row["numerator"],
        )
    )
    selected = []
    seen_j: set[Fraction] = set()
    duplicate_j_rejections = 0
    for candidate in candidates:
        if candidate["j_invariant"] in seen_j:
            duplicate_j_rejections += 1
            continue
        seen_j.add(candidate["j_invariant"])
        selected.append(candidate)
        if len(selected) == cohort_size:
            break
    if len(selected) != cohort_size:
        raise ArithmeticError("the frozen universe does not contain enough admissible fibres")

    rows = []
    for selection_rank, candidate in enumerate(selected, start=1):
        numerator = candidate["numerator"]
        denominator = candidate["denominator"]
        parameter = str(Fraction(numerator, denominator))
        sample_key = {
            "salt": SALT,
            "family": "published-R17-MW17",
            "projective_pair": [numerator, denominator],
        }
        row = {
            "sample_id": canonical_hash(sample_key)[:24],
            "selection_rank_one_based": selection_rank,
            "family": "published-R17-MW17",
            "generic_rank": 17,
            "parameter": parameter,
            "projective_pair": [numerator, denominator],
            "projective_height": candidate["projective_height"],
            "cubic_order_discriminant": str(candidate["cubic_order_discriminant"]),
            "absolute_cubic_order_discriminant_bits": candidate[
                "absolute_cubic_order_discriminant_bits"
            ],
            "j_invariant": (
                str(candidate["j_invariant"].numerator)
                if candidate["j_invariant"].denominator == 1
                else f"{candidate['j_invariant'].numerator}/{candidate['j_invariant'].denominator}"
            ),
            "selection_status": "FROZEN_BEFORE_CLASS_GROUP_AND_POINT_SEARCH",
            "class_quotient_features": None,
            "feature_status": "NOT_OPENED",
            "detector_outcome": None,
            "outcome_status": "SEALED_UNTIL_ALL_FEATURES_FREEZE",
        }
        rows.append(row)

    candidate_commitment = [selection_record(row) for row in rows]
    return {
        "schema": SCHEMA,
        "status": STATUS,
        "laboratory_id": "r17-small-field-class-quotient-v1",
        "question": "Does pre-search dim Q_t predict later certified Mordell-Weil quotient gain?",
        "commitment": {
            "candidate_count": len(rows),
            "candidate_list_sha256": canonical_hash(candidate_commitment),
            "frozen_before_class_group_computation": True,
            "frozen_before_detector_search": True,
            "selection_used_rank_or_point_search_labels": False,
            "selection_used_exceptional_point_coordinates": False,
            "selection_used_nagao_or_selmer_scores": False,
        },
        "selection_protocol": {
            "family": "published R17 in the exact published t-coordinate",
            "parameter_universe": (
                "all reduced finite a/b with b>0 and max(|a|,b)<=24"
            ),
            "projective_height_bound": height_bound,
            "raw_parameter_count": universe_count,
            "structural_rejections": {
                **structural_rejections,
                "duplicate_exact_j_before_cohort_filled": duplicate_j_rejections,
            },
            "ordering": (
                "absolute discriminant of z^3+16*A_h(a,b)*z+64*B_h(a,b), "
                "then projective height, denominator, numerator"
            ),
            "cohort_size": cohort_size,
            "ordinary_means": (
                "nonsingular specialization with irreducible 2-division cubic; it does "
                "not assert absence from every historical computation"
            ),
            "field_discriminant_boundary": (
                "The selection metric is the explicit cubic-order discriminant.  The "
                "maximal-order field discriminant is computed and frozen in Phase 1."
            ),
        },
        "family_expansion_gate": {
            "recovered_A1_MW16_input": relative(MW16),
            "included_in_v1": False,
            "reason": (
                "The current exact lambda gauge did not pass the small-field feasibility "
                "probe: maximal-order reduction at lambda=0 exceeded 30 seconds.  This is "
                "an operational observation, not a mathematical exclusion."
            ),
            "admission_rule": (
                "Freeze a rank-blind rational reparameterization and demonstrate complete "
                "certified BNF plus generic-MW16 S-class images under the same per-row "
                "resource envelope before extending the cohort."
            ),
        },
        "phase_order": [
            "0_freeze_rank_blind_cohort",
            "1_compute_and_freeze_complete_unconditional_class_quotient_features",
            "2_freeze_uniform_blind_half_lattice_protocol_bound_to_phase_1_hash",
            "3_run_detector_without_loading_feature_values",
            "4_join_features_and_certified_outcomes_once",
        ],
        "rows": rows,
        "inputs": {
            relative(path): digest(path) for path in (MODEL, SECTIONS, MW16)
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                "sage -python elkies-k3/scripts/"
                "build_r17_small_field_class_quotient_cohort.sage"
            ),
        },
        "claim_boundary": [
            "This artifact freezes a rank-blind finite cohort, not any class-group result.",
            "Small explicit order discriminant does not by itself prove that BNF certification will finish.",
            "A later bounded detector miss will not prove rank 17 or absence of Mordell-Weil escape.",
            "Any conclusion is scoped to this discriminant-selected R17 population.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=HEIGHT_BOUND)
    parser.add_argument("--cohort-size", type=int, default=COHORT_SIZE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build(height_bound=args.height_bound, cohort_size=args.cohort_size)
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored small-field cohort differs from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17SMALLFIELDCOHORT"
        f"|rows={len(document['rows'])}"
        f"|hash={document['commitment']['candidate_list_sha256']}"
        f"|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
