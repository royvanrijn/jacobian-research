#!/usr/bin/env python3
"""Bounded search for one affine target line carrying all quintic groups.

The pinned anchors are the exceptional transitive-group output of the
exhaustive projective-height search through height 30.  For every line
through anchors of two different groups, this script scans reduced rational
parameters and sends only exact-screen survivors to PARI/GP.  It is a
discovery/minimality experiment, not an oracle-free group certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from pathlib import Path

from search_universal_quintic_calculator import (
    SCREEN_PRIMES,
    Target,
    discriminant_square_class,
    is_square,
    local_f20_allowed,
    pari_classify,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "universal_quintic_target_line_search.json"
)
EXCEPTIONAL_GROUPS = ("A5", "C5", "D5", "F20")

# Exhaustive output at projective height 30, modulo
# (B,C) -> (-B,-C).  The involution is expanded below.
ANCHORS_HEIGHT_30 = (
    ("A5", (5, 5, 0, 2)),
    ("A5", (4, 8, 16, -1)),
    ("A5", (3, 3, 4, -18)),
    ("A5", (5, 2, 21, -5)),
    ("C5", (10, 10, 0, 7)),
    ("D5", (10, 4, 21, -20)),
    ("D5", (1, -1, 27, 28)),
    ("F20", (10, 5, 15, 4)),
)

Point = tuple[Fraction, Fraction, Fraction]
LineKey = tuple[tuple[int, int, int], Point]


def affine_point(projective: tuple[int, int, int, int]) -> Point:
    w, pi, b, c = projective
    return Fraction(pi, w), Fraction(b, w), Fraction(c, w)


def target_from_point(point: Point) -> Target:
    denominator = math.lcm(*(coordinate.denominator for coordinate in point))
    values = [denominator, *(int(coordinate * denominator) for coordinate in point)]
    content = reduce(math.gcd, (abs(value) for value in values))
    values = [value // content for value in values]
    if values[0] < 0:
        values = [-value for value in values]
    return Target(*values)


def rational_parameters(bound: int) -> list[Fraction]:
    return [
        Fraction(numerator, denominator)
        for denominator in range(1, bound + 1)
        for numerator in range(-bound, bound + 1)
        if math.gcd(abs(numerator), denominator) == 1
    ]


def primitive_direction(direction: Point) -> tuple[int, int, int]:
    denominator = math.lcm(*(coordinate.denominator for coordinate in direction))
    values = [int(coordinate * denominator) for coordinate in direction]
    content = reduce(math.gcd, (abs(value) for value in values))
    values = [value // content for value in values]
    if next(value for value in values if value) < 0:
        values = [-value for value in values]
    return values[0], values[1], values[2]


def line_key(point: Point, direction: Point) -> LineKey:
    primitive = primitive_direction(direction)
    moment = (
        point[1] * primitive[2] - point[2] * primitive[1],
        point[2] * primitive[0] - point[0] * primitive[2],
        point[0] * primitive[1] - point[1] * primitive[0],
    )
    return primitive, moment


def parameter_on_line(point: Point, origin: Point, direction: Point) -> Fraction:
    index = next(index for index, value in enumerate(direction) if value)
    parameter = (point[index] - origin[index]) / direction[index]
    assert all(
        coordinate == base + parameter * step
        for coordinate, base, step in zip(point, origin, direction, strict=True)
    )
    return parameter


def expanded_anchors() -> list[tuple[str, Point]]:
    result: set[tuple[str, Point]] = set()
    for group, projective in ANCHORS_HEIGHT_30:
        point = affine_point(projective)
        result.add((group, point))
        result.add((group, (point[0], -point[1], -point[2])))
    return sorted(result)


def exceptional_candidate(target: Target, screen_prime_count: int) -> bool:
    square_class = discriminant_square_class(target)
    if square_class <= 0:
        return False
    if is_square(square_class):
        return True
    return all(
        local_f20_allowed(target.primitive_inverse_coefficients, prime)
        for prime in SCREEN_PRIMES[:screen_prime_count]
    )


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def point_json(point: Point) -> list[str]:
    return [fraction_text(value) for value in point]


def run_search(parameter_bound: int, screen_prime_count: int, gp: str) -> dict:
    anchors = expanded_anchors()

    # Recheck the pinned discovery catalogue before using it.
    anchor_targets = {target_from_point(point) for _, point in anchors}
    anchor_classifications = pari_classify(anchor_targets, gp)
    for expected_group, point in anchors:
        assert anchor_classifications[target_from_point(point)] == expected_group

    lines: dict[LineKey, tuple[Point, Point]] = {}
    for (left_group, left), (right_group, right) in itertools.combinations(
        anchors, 2
    ):
        if left_group == right_group:
            continue
        direction = tuple(
            right_coordinate - left_coordinate
            for left_coordinate, right_coordinate in zip(left, right, strict=True)
        )
        lines.setdefault(line_key(left, direction), (left, direction))

    parameters = rational_parameters(parameter_bound)
    occurrences: dict[Target, set[tuple[LineKey, Fraction]]] = defaultdict(set)
    candidates: set[Target] = set()
    for key, (origin, direction) in lines.items():
        for parameter in parameters:
            point = tuple(
                base + parameter * step
                for base, step in zip(origin, direction, strict=True)
            )
            if point[0] == 0:
                continue
            target = target_from_point(point)
            if exceptional_candidate(target, screen_prime_count):
                candidates.add(target)
                occurrences[target].add((key, parameter))

    classifications = pari_classify(candidates, gp)
    coverage: dict[LineKey, dict[str, list[tuple[Fraction, Target]]]] = {
        key: defaultdict(list) for key in lines
    }
    for target, group in classifications.items():
        for key, parameter in occurrences[target]:
            coverage[key][group].append((parameter, target))

    # Include every catalogue anchor geometrically on the line, regardless of
    # its parameter in this particular anchor normalization.
    for key, (origin, direction) in lines.items():
        primitive, _ = key
        key_direction = tuple(Fraction(value) for value in primitive)
        for group, point in anchors:
            if line_key(point, key_direction) != key:
                continue
            parameter = parameter_on_line(point, origin, direction)
            coverage[key][group].append((parameter, target_from_point(point)))

    def exceptional_count(key: LineKey) -> int:
        return len(set(coverage[key]).intersection(EXCEPTIONAL_GROUPS))

    def line_score(key: LineKey) -> tuple[int, int, int, LineKey]:
        return (
            exceptional_count(key),
            len(coverage[key]),
            sum(len(values) for values in coverage[key].values()),
            key,
        )

    best_key = max(lines, key=line_score)
    best_groups = sorted(coverage[best_key])
    best_origin, best_direction = lines[best_key]
    best_rows = []
    for group in best_groups:
        parameter, target = min(
            set(coverage[best_key][group]),
            key=lambda item: (
                max(
                    abs(item[0].numerator).bit_length(),
                    item[0].denominator.bit_length(),
                ),
                abs(item[0]),
                item[1].height,
            ),
        )
        best_rows.append(
            {
                "group": group,
                "t": fraction_text(parameter),
                "projective_target": [target.w, target.pi, target.b, target.c],
            }
        )

    histogram: dict[str, int] = defaultdict(int)
    for key in coverage:
        histogram[str(exceptional_count(key))] += 1

    maximum_exceptional = max(exceptional_count(key) for key in lines)
    found = maximum_exceptional == len(EXCEPTIONAL_GROUPS)
    return {
        "status": "found" if found else "not_found",
        "claim_kind": "bounded_computation",
        "fixed_map_inverse_polynomial": (
            "Pi^5*S^5 - 5*Pi*S^3 - 2*B*S^2 + 4*S - 2*C"
        ),
        "search": {
            "anchor_projective_height_bound": 30,
            "anchor_count_after_sign_involution": len(anchors),
            "line_count": len(lines),
            "parameter_height_bound": parameter_bound,
            "parameter_count": len(parameters),
            "f20_screen_primes": list(SCREEN_PRIMES[:screen_prime_count]),
            "pari_role": "classification of exact-screen survivors only",
        },
        "anchor_group_counts": {
            group: sum(anchor_group == group for anchor_group, _ in anchors)
            for group in EXCEPTIONAL_GROUPS
        },
        "screen_survivor_count": len(candidates),
        "classified_transitive_count": len(classifications),
        "line_exceptional_group_count_histogram": dict(sorted(histogram.items())),
        "maximum_exceptional_groups_on_one_line": maximum_exceptional,
        "best_line": {
            "u": point_json(best_origin),
            "v": point_json(best_direction),
            "groups": best_groups,
            "rows": best_rows,
        },
        "interpretation": (
            "No line in this finite anchor/parameter search realizes all four "
            "exceptional groups; S5 cannot repair two missing exceptional groups. "
            "This is not a global nonexistence result."
            if not found
            else (
                "Candidate only: add S5 and replay all five rows with an "
                "oracle-free certificate before promotion."
            )
        ),
        "future_witness_cost": {
            "rational_encoding": (
                "reduced signed numerator/positive denominator; normalize two "
                "witness parameters to 0 and 1 and minimize over all ten pairs"
            ),
            "coefficient_bit_length": (
                "lexicographic (maximum bit length, total bit length) over "
                "u, v, and the five normalized t values"
            ),
            "largest_witness_prime": (
                "minimum possible maximum over the exact modular witnesses"
            ),
            "total_resolvent_degree": (
                "sum of displayed resolvent degrees; degree-6 Dummit certificates "
                "for D5 and F20 give target total 12"
            ),
            "certificate_byte_count": (
                "UTF-8 bytes of canonical minified JSON with sorted keys, "
                "including exact factors/automorphisms but excluding prose"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parameter-bound", type=int, default=20)
    parser.add_argument("--screen-primes", type=int, default=12)
    parser.add_argument("--gp", default="gp")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    assert args.parameter_bound >= 1
    assert 1 <= args.screen_primes <= len(SCREEN_PRIMES)

    result = run_search(args.parameter_bound, args.screen_primes, args.gp)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write:
        args.artifact.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
