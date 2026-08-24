#!/usr/bin/env python3
"""Classify the rank-seven bases in the multi-prime max-200 Mestre replay.

This is an exact *base-point* classification only.  It does not infer that a
particular recovered pair extends to either known component, nor does it make
a Mordell--Weil independence claim.  Its purpose is to distinguish genuinely
new moduli points from points already lying on a known rank-13 base family.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from icarm_curve245_mestre import fermigier_roots
from screen_mestre_two_section_transverse_seeds import normalized_moduli
from verify_mestre_transverse_two_section_component import component_coordinates


Q = Fraction
DEFAULT_INPUT = Path(
    "artifacts/generated-results/elliptic_mestre_two_section_seed_screen_max200_all_primes.json"
)
FERMIGIER_PARAMETERS = (Q(-3), -Q(8, 3))
FERMIGIER_ROOTS = (0, 8, 58, 77, 85, 102)


def on_rational_component(roots: tuple[int, ...]) -> bool:
    """Test the exact moduli identity for the previously audited base curve."""

    moduli = normalized_moduli(roots)
    z = moduli[0] + 35
    return component_coordinates(z)[:4] == moduli


def fermigier_normalization_matches(roots: tuple[int, ...]) -> bool:
    """Compare the candidate root set with alpha_1 -> 0, alpha_2 -> 1."""

    source = fermigier_roots(*FERMIGIER_PARAMETERS)
    normalized = tuple(sorted((root - source[0]) / (source[1] - source[0]) for root in source))
    target = tuple(sorted(Q(root, roots[1]) for root in roots))
    return normalized == target


def replay(input_path: Path = DEFAULT_INPUT) -> dict[str, object]:
    source = json.loads(input_path.read_text())
    rank_seven = [
        record
        for record in source["records"]
        if record["transverse_rank_seven_pairs"]
    ]
    classifications = []
    for record in rank_seven:
        roots = tuple(record["roots"])
        if on_rational_component(roots):
            kind = "previous rational two-section component base curve"
        elif roots == FERMIGIER_ROOTS and fermigier_normalization_matches(roots):
            kind = "Fermigier two-parameter base family at u=-3, v=-8/3"
        else:
            raise AssertionError(f"unclassified rank-seven base point: {roots}")
        classifications.append(
            {
                "roots": list(roots),
                "reconstructed_section_count": len(record["reconstructed_sections"]),
                "rank_seven_pair_count": len(record["transverse_rank_seven_pairs"]),
                "base_classification": kind,
            }
        )
    if len(rank_seven) != 8 or sum(
        item["rank_seven_pair_count"] for item in classifications
    ) != 38:
        raise AssertionError("the pinned all-prime rank-seven population changed")
    return {
        "status": "exact base-point classification of multi-prime rank-seven seeds completed",
        "input": str(input_path),
        "input_scope": source["scope"],
        "rank_seven_base_point_count": len(classifications),
        "rank_seven_pair_count": sum(item["rank_seven_pair_count"] for item in classifications),
        "classifications": classifications,
        "Fermigier_normalization": {
            "parameters": [str(value) for value in FERMIGIER_PARAMETERS],
            "source_affine_normalization": "alpha_1 -> 0, alpha_2 -> 1",
            "target_roots": list(FERMIGIER_ROOTS),
            "exact_root_set_match": True,
        },
        "conclusion": (
            "every rank-seven base point recovered in this bounded multi-prime replay "
            "lies on a previously known rank-13 base family"
        ),
        "not_established": [
            "that every recovered labelled pair has a generic visible-subgroup relation",
            "pair intersections, Shioda data, saturation, or independence at the Fermigier point",
            "any absence statement outside the bounded reconstruction screen",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(args.input), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
