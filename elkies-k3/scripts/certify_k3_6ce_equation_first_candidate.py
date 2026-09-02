#!/usr/bin/env python3
"""Certify the equation-first status of the determinant-384 MW2-to-MW15 K3."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
ADAPTER = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-source-search-target-partner1-lattice-only-v1.json"
SOURCES = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
CLASSIFIER = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-source-isometries-v1.json"
FIBRES = {
    5: GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-fibre-ansatz-mod5-v1.json",
    7: GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-fibre-ansatz-mod7-v1.json",
}
MARKINGS = {
    (5, "square"): GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod5-square-v1.json",
    (5, "nonsquare"): GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod5-nonsquare-v1.json",
    (7, "square"): GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod7-square-v1.json",
    (7, "nonsquare"): GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod7-nonsquare-v1.json",
}
CORRIDOR = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-same-ns-compiler-routes-rankfirst-cap2000-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-equation-first-candidate-v1.json"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    adapter = load(ADAPTER)
    sources = load(SOURCES)
    classifier = load(CLASSIFIER)
    fibres = {prime: load(path) for prime, path in FIBRES.items()}
    markings = {key: load(path) for key, path in MARKINGS.items()}
    corridor = load(CORRIDOR)

    if adapter["status"] != "PASS_EXACT_CATALOGUE_LATTICE_ONLY_SOURCE_SEARCH_TARGET_EXTRACTION":
        raise ArithmeticError("lattice-only adapter status changed")
    if adapter.get("equation_work_authorized") is not False:
        raise ArithmeticError("unresolved T-arithmetic gate was promoted")
    if adapter["frame"]["mw_rank_for_rho_19"] != 15 or adapter["determinant"] != 384:
        raise ArithmeticError("target invariants changed")
    if sources["status"] != "PASS_SUCCESS_CONDITION_HIT" or len(sources["sources"]) != 144:
        raise ArithmeticError("low-MW source inventory changed")
    root_histogram = Counter(row["source"]["root_type"] for row in sources["sources"])
    if root_histogram != Counter({"A10+A5": 69, "A15": 45, "A2+A5+A8": 30}):
        raise ArithmeticError("source root histogram changed")
    if {row["source"]["mw_rank_for_rho_19"] for row in sources["sources"]} != {2}:
        raise ArithmeticError("source inventory is no longer purely MW2")
    good_basis = [
        row
        for row in sources["sources"]
        if row["source"]["pole_audit"]["basis_with_all_poles_at_most_two"]
    ]
    if len(good_basis) != 50:
        raise ArithmeticError("physical MW-basis count changed")
    if classifier["accounting"] != {
        "class_sizes": [24],
        "integral_isometry_classes": 1,
        "reduced_gram_rows": 24,
        "selected_physical_basis_profiles": 1,
    }:
        raise ArithmeticError("best source classification changed")

    expected_fibres = {5: 152, 7: 1032}
    for prime, payload in fibres.items():
        if payload["status"] != "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ":
            raise ArithmeticError(f"fibre status changed at p={prime}")
        if payload["accounting"]["squarefree_examples_with_signs"] != expected_fibres[prime]:
            raise ArithmeticError(f"fibre count changed at p={prime}")
    for (prime, square_class), payload in markings.items():
        if payload["status"] != "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_BASIS_LOCUS":
            raise ArithmeticError(f"marking status changed at p={prime}, {square_class}")
        if payload["accounting"]["marked_ordered_basis_pairs"]:
            raise ArithmeticError("empty marking chart contains a basis")
    if corridor["status"] != "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_EMPTY":
        raise ArithmeticError("corridor status changed")
    route = corridor["results"]
    if len(route) != 1 or route[0]["case"] != "k36ce" or route[0]["best_routes_by_target"]:
        raise ArithmeticError("corridor result changed")

    marking_summary = []
    for key in sorted(markings):
        payload = markings[key]
        accounting = payload["accounting"]
        marking_summary.append(
            {
                "prime": key[0],
                "twist_square_class": key[1],
                "fibre_models": payload["scope"]["fibre_models"],
                "marked_generator_sections": accounting["marked_generator_sections"],
                "models_with_both_generator_section_classes": accounting[
                    "models_with_both_generator_section_classes"
                ],
                "component_matched_pair_candidates": accounting[
                    "component_matched_pair_candidates"
                ],
                "marked_basis_pairs": accounting["marked_ordered_basis_pairs"],
            }
        )

    paths = [ADAPTER, SOURCES, CLASSIFIER, *FIBRES.values(), *MARKINGS.values(), CORRIDOR]
    payload = {
        "schema": "elkies-k3.k3-6ce-equation-first-candidate.v1",
        "status": "PASS_LATTICE_SOURCE_PROMOTION_NORMALIZED_MARKING_GATES_EMPTY",
        "surface": {
            "surface_id": "K3-6ce16abb9de3c7c5",
            "determinant": 384,
            "complex_moduli_dimension": 1,
            "t_arithmetic_status": adapter["t_arithmetic_pre_solver_gate"]
            ["pre_solver_gate"]["status"],
            "equation_work_authorized_from_t_arithmetic": False,
        },
        "target": {
            "frame_id": "K3-6ce16abb9de3c7c5-F001",
            "root_type": "A2",
            "mw_rank": 15,
        },
        "source_inventory": {
            "semistable_mw2_sources": 144,
            "root_type_histogram": dict(sorted(root_histogram.items())),
            "sources_with_physical_mw_basis_through_pole_two": len(good_basis),
        },
        "best_source": {
            "representative_source_id": "K3-6ce16abb9de3c7c5-S0008",
            "root_type": "A10+A5",
            "fibre_profile": "I6+I11+7I1",
            "mw_rank": 2,
            "support_count": 2,
            "basis_pole_profile": [0, 1],
            "integral_frame_classes": 1,
            "marked_basis_profiles": 1,
            "normalized_equation_conditions": {
                "component_depths_at_I6_I11": [[0, 4], [0, 2]],
                "required_smooth_pair_intersection": 1,
                "fibre_Hermite_compatibility_equations_on_A": 4,
                "residual_base_scaling_dimension": 1,
            },
        },
        "equation_gate": {
            "squarefree_fibre_models": expected_fibres,
            "marking_charts": marking_summary,
            "decision": (
                "Do not launch characteristic-zero lifting from the A10+A5 chart. "
                "Test the other determinant-384 MW2 marking profiles first."
            ),
        },
        "corridor_gate": {
            "status": corridor["status"],
            "best_root_rank_by_depth": [
                [row["depth"], row["best_root_rank"]] for row in route[0]["accounting"]
            ],
            "decision": "No exact MW2-to-MW15 target route in the declared depth-eight beam.",
        },
        "comparison": (
            "This is the strongest new lattice-source candidate and validates the MW-rank "
            "tradeoff, but determinant 500 remains the leading equation branch because its "
            "MW1 chart has a formally smooth Z7 marked family."
        ),
        "proof_boundary": {
            "proved": (
                "The same determinant-384 Picard-19 lattice has the displayed MW15 target "
                "and 144 exact semistable MW2 source Grams in the declared complete cut. "
                "The best source class and all four normalized finite-field marking charts "
                "are exact, as is every retained edge in the bounded corridor beam."
            ),
            "not_proved": (
                "No rational source marking, characteristic-zero equation, modular-curve "
                "identification, multisection spectrum for the rootful target, or complete "
                "neighbour-graph obstruction is asserted."
            ),
        },
        "inputs": {relative(path): digest(path) for path in paths},
        "reproduce": "python3 elkies-k3/scripts/certify_k3_6ce_equation_first_candidate.py",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print("K36CEEQUATIONFIRST|sources=144|best=A10+A5/MW2|marked=0|status=PASS")


if __name__ == "__main__":
    main()
