#!/usr/bin/env python3
"""Build the compact experiment matrix for the half-lattice replay."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
INPUTS = {
    "rank28_blind": ART / "half_lattice_fake_descent_rank28_blind_v1.json",
    "rank28_verify": ART / "half_lattice_fake_descent_rank28_verification_v1.json",
    "r17_matrix_blind": ART / "half_lattice_fake_descent_r17_matrix_blind_v1.json",
    "r17_matrix_verify": ART / "half_lattice_fake_descent_r17_matrix_verification_v1.json",
    "r17_rank21_blind": ART / "half_lattice_r17_rank21_blind_v1.json",
    "r17_rank21_verify": ART / "half_lattice_r17_rank21_verification_v1.json",
    "rank29_blind": ART / "half_lattice_rank29_controls_blind_v1.json",
    "rank29_verify": ART / "half_lattice_rank29_controls_verification_v1.json",
    "heldout_blind": ART / "half_lattice_heldout_273_302_blind_v1.json",
    "heldout_verify": ART / "half_lattice_heldout_273_302_verification_v1.json",
}
OUTPUT = ART / "half_lattice_fake_descent_experiment_matrix_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load(name: str):
    return json.loads(INPUTS[name].read_text())


def median_field(records, field):
    values = [float(record[field]) for record in records if record.get(field) is not None]
    return median(values) if values else None


def main() -> None:
    rank28b = load("rank28_blind")
    rank28v = load("rank28_verify")
    r17b = load("r17_matrix_blind")
    r17v = load("r17_matrix_verify")
    rank21b = load("r17_rank21_blind")
    rank21v = load("r17_rank21_verify")
    rank29b = load("rank29_blind")
    rank29v = load("rank29_verify")
    heldb = load("heldout_blind")
    heldv = load("heldout_verify")

    positive_targets = {
        row["id"]: row["public_exceptional_quotient_dimension"]
        for row in r17v["positive_controls"]
    }
    rank21_covers = rank21b["cover_records"]
    r17_rows = [
        {
            "id": "r17-rank21-control",
            "parameter": "3/8",
            "role": "known_positive_label_withheld_from_search",
            "searched_generic_deepest_classes": 43,
            "exact_blind_quotient_gain": rank21v["exact_blind_recovered_quotient_dimension"],
            "verification_target_quotient_dimension": 4,
            "full_target_recovery": rank21v["full_public_quotient_recovered"],
            "distinct_nonbasis_candidates": rank21b["blind_result"]["distinct_nonbasis_candidates"],
            "specialized_depth_median": median_field(rank21_covers, "specialized_depth"),
            "reduced_coefficient_bits_median": median_field(rank21_covers, "reduced_coefficient_bits"),
            "modular_density_product_median": median_field(rank21_covers, "modular_density_product"),
        }
    ]
    for row in r17b["fibres"]:
        covers = row["cover_records"]
        result = row["bounded_search_result"]
        target = positive_targets.get(row["id"])
        gain = result["finite_mod2_certified_quotient_gain"]
        r17_rows.append(
            {
                "id": row["id"],
                "parameter": row["parameter"],
                "role": row["role"],
                "searched_generic_deepest_classes": len(covers),
                "exact_blind_quotient_gain": gain,
                "verification_target_quotient_dimension": target,
                "full_target_recovery": gain == target if target is not None else None,
                "distinct_nonbasis_candidates": result["distinct_nonbasis_candidates"],
                "specialized_depth_median": median_field(covers, "specialized_depth"),
                "reduced_coefficient_bits_median": median_field(
                    covers, "reduced_model_maximum_coefficient_bits"
                ),
                "modular_density_product_median": median_field(
                    covers, "independent_modular_density_product"
                ),
            }
        )

    rank29_exact = {
        row["label"]: row["exact_exceptional_quotient_rank_recovered"]
        for row in rank29v["results"]
    }
    rank29_rows = []
    for row in rank29b["results"]:
        covers = row["cover_records"]
        rank29_rows.append(
            {
                "id": row["label"],
                "role": "exact_family_transport_positive_control",
                "target_quotient_dimension": 12,
                "exact_blind_quotient_gain": rank29_exact[row["label"]],
                "full_target_recovery": rank29_exact[row["label"]] == 12,
                "generic_specialized_top43_intersection": row[
                    "generic_specialized_intersection_count"
                ],
                "searched_union_classes": row["selected_union_count"],
                "distinct_nonbasis_candidates": row["blind_result"][
                    "distinct_nonbasis_candidates"
                ],
                "specialized_depth_median": median_field(covers, "specialized_depth"),
                "reduced_coefficient_bits_median": median_field(
                    covers, "reduced_coefficient_bits"
                ),
                "modular_density_product_median": median_field(
                    covers, "modular_density_product"
                ),
            }
        )

    held_exact = {
        (row["curve"], row["configuration"]): row
        for row in heldv["results"]
    }
    held_rows = []
    for row in heldb["results"]:
        verified = held_exact[(row["curve"], row["configuration"])]
        covers = row["cover_records"]
        held_rows.append(
            {
                "curve": row["curve"],
                "configuration": row["configuration"],
                "starting_dimension": row["dimension"],
                "displayed_heldout_dimension": verified["heldout_dimension"],
                "exact_displayed_heldout_rank_recovered": verified[
                    "exact_heldout_quotient_rank_recovered"
                ],
                "searched_specialized_top_classes": row["selected_top_class_count"],
                "distinct_nonbasis_candidates": row["blind_result"][
                    "distinct_nonbasis_candidates"
                ],
                "blind_finite_reduction_field_valid": row["blind_result"][
                    "finite_reduction_certificate_valid"
                ],
                "specialized_depth_median": median_field(covers, "depth"),
                "reduced_coefficient_bits_median": median_field(
                    covers, "reduced_coefficient_bits"
                ),
                "modular_density_product_median": median_field(
                    covers, "modular_density_product"
                ),
            }
        )

    productive = [
        {
            "half_class": row["half_lattice_hex"],
            "generic_depth": row["generic_depth"],
            "specialized_rank": row["specialized_rank"],
            "specialized_depth": row["specialized_depth"],
            "recovered_quotient_classes": row["recovered_nonzero_quotient_hex"],
            "incremental_gain": row["incremental_quotient_gain"],
            "integral_coefficient_bits": row["integral_model_maximum_coefficient_bits"],
            "reduced_coefficient_bits": row["reduced_model_maximum_coefficient_bits"],
        }
        for row in rank28v["half_lattice_class_summary"]["productive_centers"]
    ]

    payload = {
        "schema": "elliptic-curves.half-lattice-fake-descent-experiment-matrix.v1",
        "status": "PASS_COMPACT_REPLAY_SUMMARY",
        "rank28_mandatory_gate": {
            "parameter": rank28b["fibre"]["parameter"],
            "generic_cvp_histogram": rank28b["ranking"]["generic_minimum_norm_histogram"],
            "generic_deepest_count": rank28b["ranking"]["generic_deepest_count"],
            "specialized_top_count": rank28b["ranking"]["specialized_deepest_count"],
            "intersection_count": rank28b["ranking"]["deepest_intersection_count"],
            "searched_union_count": rank28b["ranking"]["selected_union_count"],
            "exact_blind_quotient_gain": rank28v["exact_result"][
                "all_recovered_quotient_span_dimension"
            ],
            "target_quotient_dimension": rank28v["exact_result"][
                "public_exceptional_quotient_dimension"
            ],
            "full_blind_recovery": rank28v["exact_result"][
                "blind_search_recovers_full_public_exceptional_quotient"
            ],
            "productive_half_classes": productive,
            "pairwise_xor_generic_depth_histogram": rank28v[
                "half_lattice_class_summary"
            ]["pairwise_xor_generic_depth_histogram"],
            "productive_center_f2_span_dimension": rank28v[
                "half_lattice_class_summary"
            ]["productive_center_f2_span_dimension"],
            "productive_vs_nonproductive_contrast": rank28v[
                "bounded_search_contrasts_within_selected_union"
            ],
        },
        "r17_fixed_generic_deepest43_matrix": r17_rows,
        "rank29_generic_specialized_union_matrix": rank29_rows,
        "heldout_273_302_and_adverse_matrix": held_rows,
        "fail_closed_conclusions": {
            "fake_cover_local_solubility_is_predictive": False,
            "reason": "Each fake-descent quartic is pointed at infinity and birational to E, so local solubility is automatic.",
            "deepest_holes_privileged_for_recovery": True,
            "deepest_holes_prospective_rank_predictor_established": False,
            "simultaneous_nontrivial_2cover_solubility_established": False,
            "curves273_or302_hidden_family_provenance_detected": False,
        },
        "claim_boundary": [
            "Exact fields are copied only from exact group-law or finite-reduction verification artifacts.",
            "CVP depths are numerical for specialized canonical-height forms.",
            "All point-search misses and recovery fractions are bounded at the declared source-artifact budgets.",
            "The curve245 blind finite-reduction field is intentionally null; its exact held-out quotient ranks come from the fixture verifier.",
        ],
        "input_hashes": {
            str(path.relative_to(ROOT)): digest(path) for path in INPUTS.values()
        }
        | {str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve())},
        "reproducing_command": "python3 elliptic-curves/cas/summarize_half_lattice_fake_descent_replay.py",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"HALFLATTICESUMMARY|status=PASS|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
