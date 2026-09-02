#!/usr/bin/env sage-python
"""Certify the equation-first promotion of the determinant-500 MW17 K3.

The expanded rank-seven catalogue first exposed this surface through a
rootless MW17 frame.  The ordered same-genus Niemeier search subsequently
found three reduced Gram representatives with root system A3+A4+A9 and MW
rank one.  This script proves that those representatives are one integral
isometry class and records the exact source invariants, the target's current
multisection audit, and the remaining arithmetic/equation gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, matrix, pari


ROOT = Path(__file__).resolve().parents[2]
TARGET = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-source-search-target-v1.json"
)
SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
)
MULTISECTIONS = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-multisection-spectrum-v1.json"
)
DEGREE_THREE = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-degree3-complete-v1.json"
)
FIBRES_MOD7 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod7-v1.json"
)
MARKING_MOD7 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-v1.json"
)
MARKING_MOD7_NONSQUARE = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-nonsquare-v1.json"
)
HENSEL_MOD7 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-marked-gf7-hensel-v1.json"
)
FORMAL_MOD7 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-formal-smoothness-v1.json"
)
QQ_REJECTION = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-source-qq-rejection-v1.json"
)
RATIONAL_SCAN = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-rational-parameter-scan-v1.json"
)
MW2_ISOMETRIES = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-semistable-mw2-source-isometries-v1.json"
)
MW2_FIBRES_MOD5 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod5-v1.json"
)
MW2_FIBRES_MOD7 = (
    ROOT / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod7-v1.json"
)
MW2_MARKINGS = (
    ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod5-v1.json",
    ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod5-nonsquare-v1.json",
    ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod7-v1.json",
    ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod7-nonsquare-v1.json",
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-k3-04b86146cc6b284b-equation-first-promotion-v1.json"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def integral_isometry(left, right):
    """Return Q with Q*left*Q^t=right, normalizing PARI's convention."""

    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw)
    for value in (candidate, candidate.transpose()):
        if value * left * value.transpose() == right:
            return value
        if value * right * value.transpose() == left:
            inverse = value.inverse()
            if inverse.denominator() == 1:
                inverse = inverse.change_ring(ZZ)
                if inverse * left * inverse.transpose() == right:
                    return inverse
    raise ArithmeticError("unrecognized PARI qfisom orientation")


def build(arguments):
    target_payload = json.loads(arguments.target.read_text())
    source_payload = json.loads(arguments.sources.read_text())
    multisection_payload = json.loads(arguments.multisections.read_text())
    degree_three_payload = json.loads(arguments.degree_three.read_text())
    fibres_mod7 = json.loads(arguments.fibres_mod7.read_text())
    marking_mod7 = json.loads(arguments.marking_mod7.read_text())
    marking_mod7_nonsquare = json.loads(arguments.marking_mod7_nonsquare.read_text())
    hensel_mod7 = json.loads(arguments.hensel_mod7.read_text())
    formal_mod7 = json.loads(arguments.formal_mod7.read_text())
    qq_rejection = json.loads(arguments.qq_rejection.read_text())
    rational_scan = json.loads(arguments.rational_scan.read_text())
    mw2_isometries = json.loads(arguments.mw2_isometries.read_text())
    mw2_fibres_mod5 = json.loads(arguments.mw2_fibres_mod5.read_text())
    mw2_fibres_mod7 = json.loads(arguments.mw2_fibres_mod7.read_text())
    mw2_markings = [json.loads(path.read_text()) for path in arguments.mw2_markings]

    assert target_payload["surface_id"] == "K3-04b86146cc6b284b"
    assert source_payload["status"] == "PASS_SUCCESS_CONDITION_HIT"
    target = target_payload["frame"]
    target_gram = matrix(ZZ, target["gram"])
    assert target_gram.det() == 500
    assert int(pari(target_gram).qfminim(2)[0]) == 0
    assert target["root_rank"] == 0 and target["mw_rank_for_rho_19"] == 17

    candidates = [
        row
        for row in source_payload["sources"]
        if row["source"]["mw_rank_for_rho_19"] == 1
        and row["source"]["root_type"] == "A3+A4+A9"
        and row["source"]["support_count"] == 3
        and row["source"]["torsion"] == 1
        and row["source"]["pole_audit"][
            "minimum_nonzero_section_pole_order"
        ]
        == 1
    ]
    assert len(candidates) == 3
    candidates.sort(key=lambda row: row["source_id"])
    representative = candidates[0]
    representative_gram = matrix(ZZ, representative["source"]["gram"])
    isometries = []
    for candidate in candidates:
        candidate_gram = matrix(ZZ, candidate["source"]["gram"])
        isometry = integral_isometry(candidate_gram, representative_gram)
        assert isometry is not None and abs(isometry.det()) == 1
        isometries.append(
            {
                "source_id": candidate["source_id"],
                "to_representative": rows(isometry),
                "determinant": int(isometry.det()),
            }
        )

    source = representative["source"]
    assert source["determinant"] == 500
    assert source["root_rank"] == 16
    assert source["root_determinant"] == 200
    assert QQ(source["mw_regulator"]) == QQ(5) / 2
    assert QQ(source["root_determinant"]) * QQ(source["mw_regulator"]) == 500
    assert source["pole_audit"]["basis_with_all_poles_at_most_two"] is True

    audited_target = next(
        row
        for row in multisection_payload["targets"]
        if row["frame_id"] == target["frame_id"]
    )
    richness = audited_target["richness_coordinates"]
    degree_three = next(
        row
        for row in degree_three_payload["spectra"]
        if row["frame_id"] == target["frame_id"]
    )
    assert fibres_mod7["scan"]["exhausted"] is True
    assert fibres_mod7["accounting"]["squarefree_examples_with_signs"] == 6
    assert marking_mod7["accounting"]["marked_mw1_sections"] == 2
    assert marking_mod7_nonsquare["accounting"]["marked_mw1_sections"] == 4
    assert hensel_mod7["jacobian_certificate"]["rank_mod_prime"] == 39
    assert hensel_mod7["jacobian_certificate"]["tangent_dimension"] == 1
    assert hensel_mod7["finite_precision_lift"]["achieved_precision_exponent"] == 8
    assert formal_mod7["status"] == (
        "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY"
    )
    assert qq_rejection["status"] == (
        "PASS_EXACT_QQ_POINT_REJECTED_PRIMITIVE_CLOSURE_DET20_FIFTH_ROOT"
    )
    assert rational_scan["search"]["candidate_count"] == 87
    assert rational_scan["search"]["status_counts"] == {
        "EXACT_QQ_POINT_REJECTED_PRIMITIVE_DET20": 1,
        "NO_FULL_RR": 86,
    }
    assert mw2_isometries["accounting"] == {
        "class_sizes": [9],
        "integral_isometry_classes": 1,
        "reduced_gram_rows": 9,
        "selected_physical_basis_profile_sizes": [3, 6],
        "selected_physical_basis_profiles": 2,
    }
    assert mw2_fibres_mod5["accounting"]["squarefree_examples_with_signs"] == 30
    assert mw2_fibres_mod7["accounting"]["squarefree_examples_with_signs"] == 114
    assert all(row["accounting"]["marked_ordered_basis_pairs"] == 0 for row in mw2_markings)
    assert mw2_markings[-1]["accounting"]["marked_generator_sections"] == [20, 32]
    assert mw2_markings[-1]["accounting"]["component_matched_pair_candidates"] == 28
    assert mw2_markings[-1]["accounting"]["pairs_meeting_singular_fibres"] == 28

    return {
        "schema": "elkies-k3.rank7-equation-first-promotion.v1",
        "status": "PASS_FORMALLY_SMOOTH_Z7_PROMOTION_FIRST_QQ_POINT_REJECTED_DET20",
        "surface_id": target_payload["surface_id"],
        "transcendental_lattice": {
            "gram": [[0, 0, 5], [0, 20, 0], [5, 0, 0]],
            "decomposition": "U(5) + <20>",
            "determinant": -500,
            "signature": [2, 1],
        },
        "target": {
            "frame_id": target["frame_id"],
            "determinant": 500,
            "root_type": "",
            "mw_rank_for_rho_19": 17,
            "rootless": True,
            "gram_sha256": target["gram_sha256"],
            "multisection_richness": richness,
            "complete_degree_three": {
                "status": degree_three_payload["status"],
                "translation_cosets": degree_three["translation_cosets"],
                "rational_trisection_translation_cosets": degree_three[
                    "rational_trisection_translation_cosets"
                ],
                "genus_one_trisection_translation_cosets": degree_three[
                    "genus_one_trisection_translation_cosets"
                ],
                "maximum_coset_minimum_norm": degree_three[
                    "maximum_coset_minimum_norm"
                ],
            },
        },
        "source": {
            "promoted_source_id": representative["source_id"],
            "integral_isometry_class_representatives_merged": len(candidates),
            "root_type": "A3+A4+A9",
            "multiplicative_reducible_fibres": ["I4", "I5", "I10"],
            "reducible_fibre_support_count": 3,
            "mw_rank_for_rho_19": 1,
            "torsion_order": 1,
            "mw_height_gram": [["5/2"]],
            "mw_regulator": "5/2",
            "minimum_nonzero_section_pole_order": 1,
            "complete_mw_basis_pole_profile": [1],
            "frame_gram": source["gram"],
            "frame_gram_sha256": source["gram_sha256"],
            "root_determinant_times_regulator_over_torsion_squared": "500",
            "ambient_provenance": sorted(
                {
                    ambient
                    for candidate in candidates
                    for ambient in candidate["ambient_provenance"]
                }
            ),
        },
        "integral_isometry_classification": {
            "candidate_source_ids": [row["source_id"] for row in candidates],
            "classes": 1,
            "explicit_isometries": isometries,
        },
        "equation_first_score": {
            "target_mw_rank": 17,
            "source_mw_rank": 1,
            "reducible_fibre_support_count": 3,
            "semistable_configuration_compatible": True,
            "expected_ns_locus_dimension": 1,
            "expected_fibre_stratum_dimension": 2,
            "expected_additional_section_conditions": 1,
            "minimum_section_pole_order": 1,
            "rational_source_marking": None,
            "source_marking_galois_orbit_size": None,
            "rational_parameterization": None,
            "certified_neighbor_route": None,
        },
        "equation_gate": {
            "prime": 7,
            "normalized_fibre_polynomials_exhausted": 7**8,
            "squarefree_I4_I5_I10_5I1_models": 6,
            "square_twist_marked_sections": 2,
            "nonsquare_twist_marked_sections": 4,
            "marked_system_variables": 40,
            "marked_system_equations": 53,
            "jacobian_rank": 39,
            "tangent_dimension": 1,
            "unit_minor_mod_7": 6,
            "finite_lift_precision_exponent": 8,
            "formally_smooth_relative_dimension": 1,
            "first_rational_reconstruction": {
                "free_parameter_integer": -20,
                "result": "REJECTED_DISPLAYED_SECTION_IS_FIVE_DIVISIBLE",
                "primitive_closure_determinant": 20,
            },
            "bounded_integral_parameter_scan": {
                "smooth_GF7_disks": 3,
                "integer_parameters": 87,
                "no_full_rational_reconstruction": 86,
                "exact_QQ_points": 1,
                "intended_determinant_500_points": 0,
            },
            "status": "PASS_EXACT_EXHAUSTIVE_GF7_MARKING_AND_FORMAL_SMOOTHNESS",
        },
        "same_surface_mw2_fallback": {
            "source_id_tested": "K3-04b86146cc6b284b-S2021",
            "root_type": "A3+A4+A8",
            "multiplicative_reducible_fibres": ["I4", "I5", "I9"],
            "mw_rank_for_rho_19": 2,
            "basis_pole_profile": [0, 0],
            "reduced_gram_rows": 9,
            "integral_isometry_classes": 1,
            "selected_basis_profiles": 2,
            "common_component_depth_profiles": [[1, 0, 1], [1, 0, 3]],
            "required_smooth_pair_intersection": 1,
            "GF5_squarefree_fibre_models": 30,
            "GF7_squarefree_fibre_models": 114,
            "GF7_nonsquare_individual_generator_sections": [20, 32],
            "GF7_nonsquare_component_matched_pair_candidates": 28,
            "GF7_nonsquare_pairs_meeting_reducible_fibres": 28,
            "marked_basis_pairs_over_GF5_or_GF7": 0,
            "decision": "RETAIN_AS_EQUATION_RICH_BUT_MARKING_POOR_FALLBACK",
        },
        "promotion_decision": {
            "decision": "PROMOTE_TO_FORMAL_SMOOTHNESS_AND_RATIONAL_PARAMETRIZATION_GATE",
            "reason": (
                "The same determinant-500 Picard-19 lattice family contains a "
                "rootless MW17 target and a single exact MW1 A3+A4+A9 source class "
                "with three supports, trivial torsion, and a pole-one generator."
            ),
            "next_gates": [
                "derive or identify a rational parameterization of the marked I4-I5-I10 source",
                "enforce nondivisibility of the height-5/2 marking on every rational candidate",
                "determine the source marking field and Galois orbit",
                "certify a low-pole neighbour corridor to the rootless target",
                "widen beyond the single semistable pole-[0,0] MW2 isometry class only if the MW1 rational gate remains empty",
                "compare other equation-friendly MW15--17 surfaces on the same exact multisection coordinates",
            ],
        },
        "proof_boundary": {
            "proved": (
                "The target is an exact rootless determinant-500 rank-17 frame. "
                "The promoted source is an exact determinant-500 rank-17 frame "
                "with primitive A3+A4+A9 roots, MW rank one, trivial torsion, "
                "height 5/2, and an exact pole-one MW basis. The three retained "
                "reduced Gram representatives are explicitly integrally isometric."
            ),
            "inferred": (
                "The lattice-polarized moduli dimension is one. A-type roots are "
                "compatible with I4, I5, and I10, but the remaining singular fibres "
                "and full semistability require an equation."
            ),
            "not_proved": (
                "No QQ equation, rational source marking, marking orbit, section "
                "specialization, or elliptic-neighbour route is yet certified."
            ),
        },
        "inputs": {
            relative(arguments.target): digest(arguments.target),
            relative(arguments.sources): digest(arguments.sources),
            relative(arguments.multisections): digest(arguments.multisections),
            relative(arguments.degree_three): digest(arguments.degree_three),
            relative(arguments.fibres_mod7): digest(arguments.fibres_mod7),
            relative(arguments.marking_mod7): digest(arguments.marking_mod7),
            relative(arguments.marking_mod7_nonsquare): digest(arguments.marking_mod7_nonsquare),
            relative(arguments.hensel_mod7): digest(arguments.hensel_mod7),
            relative(arguments.formal_mod7): digest(arguments.formal_mod7),
            relative(arguments.qq_rejection): digest(arguments.qq_rejection),
            relative(arguments.rational_scan): digest(arguments.rational_scan),
            relative(arguments.mw2_isometries): digest(arguments.mw2_isometries),
            relative(arguments.mw2_fibres_mod5): digest(arguments.mw2_fibres_mod5),
            relative(arguments.mw2_fibres_mod7): digest(arguments.mw2_fibres_mod7),
            **{relative(path): digest(path) for path in arguments.mw2_markings},
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_k3_04b_equation_first_promotion.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=TARGET)
    parser.add_argument("--sources", type=Path, default=SOURCES)
    parser.add_argument("--multisections", type=Path, default=MULTISECTIONS)
    parser.add_argument("--degree-three", type=Path, default=DEGREE_THREE)
    parser.add_argument("--fibres-mod7", type=Path, default=FIBRES_MOD7)
    parser.add_argument("--marking-mod7", type=Path, default=MARKING_MOD7)
    parser.add_argument(
        "--marking-mod7-nonsquare", type=Path, default=MARKING_MOD7_NONSQUARE
    )
    parser.add_argument("--hensel-mod7", type=Path, default=HENSEL_MOD7)
    parser.add_argument("--formal-mod7", type=Path, default=FORMAL_MOD7)
    parser.add_argument("--qq-rejection", type=Path, default=QQ_REJECTION)
    parser.add_argument("--rational-scan", type=Path, default=RATIONAL_SCAN)
    parser.add_argument("--mw2-isometries", type=Path, default=MW2_ISOMETRIES)
    parser.add_argument("--mw2-fibres-mod5", type=Path, default=MW2_FIBRES_MOD5)
    parser.add_argument("--mw2-fibres-mod7", type=Path, default=MW2_FIBRES_MOD7)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    arguments.target = arguments.target.resolve()
    arguments.sources = arguments.sources.resolve()
    arguments.multisections = arguments.multisections.resolve()
    arguments.degree_three = arguments.degree_three.resolve()
    arguments.fibres_mod7 = arguments.fibres_mod7.resolve()
    arguments.marking_mod7 = arguments.marking_mod7.resolve()
    arguments.marking_mod7_nonsquare = arguments.marking_mod7_nonsquare.resolve()
    arguments.hensel_mod7 = arguments.hensel_mod7.resolve()
    arguments.formal_mod7 = arguments.formal_mod7.resolve()
    arguments.qq_rejection = arguments.qq_rejection.resolve()
    arguments.rational_scan = arguments.rational_scan.resolve()
    arguments.mw2_isometries = arguments.mw2_isometries.resolve()
    arguments.mw2_fibres_mod5 = arguments.mw2_fibres_mod5.resolve()
    arguments.mw2_fibres_mod7 = arguments.mw2_fibres_mod7.resolve()
    arguments.mw2_markings = tuple(path.resolve() for path in MW2_MARKINGS)
    arguments.output = arguments.output.resolve()

    result = build(arguments)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != serialized:
            raise SystemExit("determinant-500 equation-first promotion is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "K304BEPROMOTE|source=A3+A4+A9/MW1|target=rootless/MW17|"
        "pole=1|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
