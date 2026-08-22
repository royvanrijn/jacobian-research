#!/usr/bin/env python3
"""Exact visible-subgroup audit at two conjugate-slope Mestre seeds.

The local conjugate-slope germ has two distant exact seeds.  This checker
maps their visible and affine quartic points to the short Jacobian at ``T=1``
and verifies explicit three-term visible-subgroup relations.  It also records
the exact mod-3 quotient rank before and after the pair.  This is deliberately
seed-level rank-neutrality evidence, not a generic relation on the unresolved
formal surface.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from alternate_quartic_covers import short_add
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from screen_mestre_fermigier_two_section_escape import square_root
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
T = Q(1)
PRIME_BOUND = 251
SEEDS = (
    {
        "roots": (0, 7, 127, 128, 225, 233),
        "lines": ((Q(233, 113), -Q(97, 113)), (Q(233, 113), Q(97, 113))),
        "relations": (((6, -1), (7, 1), (11, -1)), ((6, -1), (7, 1), (10, 1))),
    },
    {
        "roots": (0, 21, 151, 169, 200, 239),
        "lines": ((Q(239, 109), -Q(31, 109)), (Q(239, 109), Q(31, 109))),
        "relations": (((6, -1), (7, 1), (11, -1)), ((6, 1), (7, -1), (10, -1))),
    },
)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def negative(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def visible_combination(
    coefficients: tuple[Fraction, ...],
    visible: tuple[tuple[Fraction, Fraction], ...],
    relation: tuple[tuple[int, int], ...],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for index, sign in relation:
        point = visible[index] if sign > 0 else negative(visible[index])
        answer = short_add(coefficients, answer, point)
    return answer


def relation_text(relation: tuple[tuple[int, int], ...]) -> str:
    return "".join(
        ("+" if sign > 0 else "-") + f"V{index + 1}"
        for index, sign in relation
    ).lstrip("+")


def replay() -> dict[str, object]:
    records = []
    for seed in SEEDS:
        integer_roots = seed["roots"]
        construction = SixRootMestreConstruction(
            tuple(Q(value, integer_roots[1]) for value in integer_roots)
        )
        quartic = primitive_quartic_coefficients(construction, T)
        coefficients = short_jacobian_coefficients(construction, T)
        visible = tuple(
            quartic_point_to_short_jacobian(construction, T, point)
            for point in primitive_visible_points(construction, T)
        )
        affine = tuple(
            quartic_point_to_short_jacobian(
                construction,
                T,
                (intercept + slope * T, square_root(quartic_value(quartic, intercept + slope * T))),
            )
            for intercept, slope in seed["lines"]
        )
        relations = seed["relations"]
        if tuple(
            visible_combination(coefficients, visible, relation) for relation in relations
        ) != affine:
            raise AssertionError("the conjugate seed visible-subgroup relation changed")
        baseline = mod3_independence_certificate(
            coefficients, visible, prime_bound=PRIME_BOUND
        )
        augmented = mod3_independence_certificate(
            coefficients, (*visible, *affine), prime_bound=PRIME_BOUND
        )
        if (
            baseline["combined_exact_rank_over_F3"] != 9
            or augmented["combined_exact_rank_over_F3"] != 9
        ):
            raise AssertionError("the conjugate seed quotient rank changed")
        records.append(
            {
                "integer_roots": list(integer_roots),
                "lines": [[rational_text(value) for value in line] for line in seed["lines"]],
                "exact_affine_relations": [f"P{index + 1}={relation_text(relation)}" for index, relation in enumerate(relations)],
                "visible_mod3_rank": baseline["combined_exact_rank_over_F3"],
                "augmented_mod3_rank": augmented["combined_exact_rank_over_F3"],
            }
        )
    return {
        "status": "conjugate-slope seed visible-subgroup audit completed",
        "Mestre_parameter": str(T),
        "prime_bound": PRIME_BOUND,
        "records": records,
        "conclusion": "both affine points lie in the visible subgroup at each of the two audited conjugate-slope seeds",
        "not_established": [
            "a generic relation on the local conjugate-slope germ",
            "a global component identity, saturation, heights, intersections, or a Shioda Gram matrix",
            "generic rank at least 14",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
