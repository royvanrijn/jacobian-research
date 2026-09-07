#!/usr/bin/env python3
"""Exact known-branch audit for the rank-seven points in diameter 211--215.

The two recovered rank-seven bases lie on the previously verified rational
two-section base curve.  This audit compares each exact incidence tangent with
the derivative of that curve.  It is not a generic relation proof for every
specialized affine line, and it makes no Mordell--Weil claim.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from probe_mestre_two_section_local_continuation import Field, residuals, solve_square_over_q
from screen_mestre_two_section_transverse_seeds import normalized_moduli
from verify_mestre_transverse_two_section_component import component_coordinates


Q = Fraction
DEFAULT_INPUT = Path(
    "artifacts/generated-results/elliptic_mestre_two_section_transverse_band_211_215_all_primes.json"
)


def component_tangent(z: Fraction) -> tuple[Fraction, ...]:
    """Derivative of ``component_coordinates(z)[:4]``, with dc1/dz = 1."""

    return (
        Q(1),
        -Q(77, 3) + Q(13, 18) * z,
        Q(652, 3) - Q(217, 18) * z + z**2 / 6,
        -608 + Q(449, 9) * z - Q(49, 36) * z**2 + z**3 / 81,
    )


def incidence_tangent(moduli: tuple[Fraction, ...], sections: tuple[tuple[Fraction, ...], ...]) -> tuple[Fraction, ...]:
    """Normalize the rank-seven kernel vector by dc1=1 exactly."""

    jacobian = [list(value.gradient) for value in residuals((*moduli, *sections[0], *sections[1]), Field())]
    columns = range(1, 8)
    minor = [[row[index] for index in columns] for row in jacobian]
    correction = solve_square_over_q(minor, [-row[0] for row in jacobian])
    tangent = (Q(1), *correction)
    if any(sum(row[index] * tangent[index] for index in range(8)) for row in jacobian):
        raise AssertionError("computed vector is not in the incidence tangent kernel")
    return tangent


def replay(input_path: Path = DEFAULT_INPUT) -> dict[str, object]:
    source = json.loads(input_path.read_text())
    hits = source["rank_seven_hits"]
    audits = []
    for hit in hits:
        roots = tuple(hit["roots"])
        moduli = normalized_moduli(roots)
        z = moduli[0] + 35
        if component_coordinates(z)[:4] != moduli:
            raise AssertionError(f"base point is not on the rational component: {roots}")
        expected = component_tangent(z)
        pairs = []
        for record in hit["transverse_pairs"]:
            sections = tuple(tuple(Q(value) for value in section) for section in record["sections"])
            tangent = incidence_tangent(moduli, sections)
            if tangent[:4] != expected:
                raise AssertionError(f"incidence tangent diverged from the known branch at {roots}")
            pairs.append(
                {
                    "sections": [[str(value) for value in section] for section in sections],
                    "normalized_moduli_tangent": [str(value) for value in tangent[:4]],
                    "matches_known_component_tangent": True,
                }
            )
        audits.append({"roots": list(roots), "z": str(z), "pairs": pairs})
    if len(audits) != 2 or sum(len(item["pairs"]) for item in audits) != 6:
        raise AssertionError("the pinned 211--215 rank-seven population changed")
    return {
        "status": "exact tangent audit of diameter-211--215 rank-seven pairs completed",
        "input": str(input_path),
        "rank_seven_base_point_count": len(audits),
        "rank_seven_pair_count": sum(len(item["pairs"]) for item in audits),
        "audits": audits,
        "conclusion": (
            "every recovered rank-seven pair has the tangent direction of the "
            "previously verified rational two-section base curve"
        ),
        "not_established": [
            "a generic visible-subgroup relation for every specialized extra line",
            "pair intersections, Shioda data, saturation, or a rank claim",
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
