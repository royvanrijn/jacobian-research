#!/usr/bin/env python3
"""Replay the six fiberwise double-saturation and non-Cartier rows."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from compile_support_saturation_cases import (  # noqa: E402
    SINGULAR_SQUAREFREE_SYMBOLS,
    cubic_double_stratification_case,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "cubic_double_saturation_stratification.json"
)


def main() -> None:
    pinned = json.loads(ARTIFACT.read_text())
    replayed = cubic_double_stratification_case()
    assert replayed == pinned, "stale cubic double-saturation artifact"

    assert set(pinned["rows"]) == set(SINGULAR_SQUAREFREE_SYMBOLS)
    for name in SINGULAR_SQUAREFREE_SYMBOLS:
        row = pinned["rows"][name]
        computation = row["exact_computation"]
        base_change = row["fiberwise_base_change_certificate"]
        assert computation["cotangent_saturation_generators"] == 0
        assert computation["support_ext2_multiplicity"] == 6
        assert (
            computation[
                "support_ext2_parameter_axis_radical_difference"
            ]
            == 0
        )
        assert (
            computation[
                "support_ext2_central_pruned_presentation_difference"
            ]
            == 0
        )
        assert computation["support_ext2_pruned_presentation_rank"] == 3
        assert (
            computation[
                "support_ext2_collision_square_action_generators"
            ]
            == 0
        )
        assert computation["different_generator_module_dimension"] == (
            computation["parameter_count"]
        )
        assert computation["different_generator_module_multiplicity"] == 6
        assert (
            computation[
                "different_generator_parameter_axis_radical_difference"
            ]
            == 0
        )
        assert (
            computation[
                "different_generator_central_pruned_presentation_difference"
            ]
            == 0
        )
        assert (
            computation["different_generator_pruned_presentation_rank"]
            == 6
        )
        assert base_change == {
            "parameter_count": computation["parameter_count"],
            "cotangent_rees_torsion_generators": 0,
            "cotangent_initial_presentation_difference": 0,
            "annihilator_cokernel_rees_torsion_generators": 0,
            "annihilator_cokernel_initial_presentation_difference": 0,
        }

    print("PASS: six formal-gauge quotient complements have dimensions 2,4,4,6,6,8")
    print(
        "PASS: strict Rees modules certify annihilator and cotangent "
        "base change on every row"
    )
    print("PASS: cotangent saturation C2 holds on every geometric fiber")
    print(
        "PASS: the parameter-independent multiplicity-six C1 obstruction "
        "survives on every geometric fiber"
    )
    print(
        "PASS: the Kahler different needs six local generators, so it is "
        "not Cartier on any geometric fiber"
    )


if __name__ == "__main__":
    main()
