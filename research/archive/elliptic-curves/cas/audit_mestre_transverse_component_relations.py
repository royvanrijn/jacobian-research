#!/usr/bin/env python3
"""Exact rank-neutrality audit at r=2,T=1 on the rational component.

This is not a generic relation proof.  It records two exact group-law
relations between the affine points and labelled visible Mestre points at one
clean specialization, then combines an exact finite-reduction lower bound
with optional PARI/GP height and ellrank diagnostics.  It demonstrates that
this specialization cannot witness the desired rank increase.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import shutil
import subprocess

from alternate_quartic_covers import short_add
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_dsquare_four import rational_square_root
from search_mestre_root_tuple_scale import (
    height_matrix_replay,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from verify_mestre_transverse_two_section_component import (
    component_coordinates,
    split_roots,
)


Q = Fraction
R = Q(2)
T = Q(1)
PRIME_BOUND = 499
EXPECTED_COEFFICIENTS = (Q(0), Q(0), Q(0), Q(-436481168886243), Q(1186790811178179337758))
EXPECTED_PIVOTS = [1, 2, 3, 4, 5, 7, 8, 9, 10]


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def negate(point: tuple[Fraction, Fraction] | None) -> tuple[Fraction, Fraction] | None:
    return None if point is None else (point[0], -point[1])


def add_terms(
    coefficients: tuple[Fraction, ...],
    terms: tuple[tuple[tuple[Fraction, Fraction], int], ...],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point, multiplier in terms:
        for _ in range(abs(multiplier)):
            answer = short_add(coefficients, answer, point if multiplier > 0 else negate(point))
    return answer


def labelled_visible_point(
    construction: SixRootMestreConstruction,
    root: Fraction,
    sign: int,
) -> tuple[Fraction, Fraction]:
    """Map the visible point at the named root, independently of root sorting."""

    approximant = construction.square_approximant_coefficients(T)
    x_value = root + sign * T
    y_value = sum(coefficient * x_value**degree for degree, coefficient in enumerate(approximant))
    y_value /= T * construction.quartic_square_scale
    return quartic_point_to_jacobian(construction, T, (x_value, y_value))


def pari_diagnostics(
    coefficients: tuple[Fraction, ...],
    points: tuple[tuple[Fraction, Fraction], ...],
) -> dict[str, object]:
    if shutil.which("gp") is None:
        return {"status": "skipped: gp unavailable"}
    heights = height_matrix_replay(
        coefficients,
        points,
        precisions=(64, 128),
        timeout=25,
        stack_bytes=256_000_000,
    )
    a_value, b_value = coefficients[3:]
    program = (
        f"E=ellinit([0,0,0,{a_value},{b_value}]);"
        'print("ELLRANK_BEGIN");print(ellrank(E));print("ELLRANK_END");quit\n'
    )
    completed = subprocess.run(
        ["gp", "-q"], input=program, text=True, capture_output=True, check=True, timeout=25
    )
    match = re.search(r"ELLRANK_BEGIN\s*\n(.*?)\nELLRANK_END", completed.stdout, re.S)
    if match is None:
        raise AssertionError("PARI/GP omitted the ellrank result")
    text = match.group(1).strip()
    if not text.startswith("[9, 9,"):
        raise AssertionError(f"unexpected PARI/GP ellrank bounds: {text[:200]}")
    return {
        "status": "PARI/GP height and ellrank diagnostics completed",
        "height_matrix": heights,
        "ellrank_output": text,
        "ellrank_bounds": [9, 9],
    }


def replay(*, include_pari: bool = True) -> dict[str, object]:
    roots = split_roots(R)
    construction = SixRootMestreConstruction(roots)
    z = 12 * (R + 3) / (1 - R**2)
    coordinates = component_coordinates(z)
    coefficients = construction.primitive_jacobian_coefficients(T)
    if coefficients != EXPECTED_COEFFICIENTS:
        raise AssertionError("the specialization short model changed")
    quartic = construction.primitive_quartic_coefficients(T)
    visible = tuple(
        quartic_point_to_jacobian(construction, T, point)
        for point in primitive_visible_points(construction, T)
    )
    affine = []
    for intercept, slope in ((coordinates[4], coordinates[5]), (coordinates[6], coordinates[7])):
        x_value = intercept + slope * T
        y_value = rational_square_root(quartic_value(quartic, x_value))
        if y_value in (None, 0):
            raise AssertionError("an affine component point lost its ordinate")
        affine.append(quartic_point_to_jacobian(construction, T, (x_value, y_value)))

    # The labels use the unsorted split-root formula: roots[0]=0, roots[1]=1,
    # and roots[3]=(42-z)/6.  They make the relations stable independently of
    # the sorting convention in SixRootMestreConstruction.
    v0_minus = labelled_visible_point(construction, roots[0], -1)
    v0_plus = labelled_visible_point(construction, roots[0], 1)
    v1_plus = labelled_visible_point(construction, roots[1], 1)
    vb_minus = labelled_visible_point(construction, roots[3], -1)
    first_relation = add_terms(coefficients, ((v0_minus, 1), (v0_plus, -1), (v1_plus, 1)))
    second_relation = add_terms(coefficients, ((v0_minus, -1), (v1_plus, -1), (vb_minus, -1)))
    if first_relation != affine[0] or second_relation != affine[1]:
        raise AssertionError("the exact component relations changed")

    visible_certificate = mod3_independence_certificate(
        coefficients, visible, prime_bound=PRIME_BOUND
    )
    augmented_certificate = mod3_independence_certificate(
        coefficients, (*visible, *affine), prime_bound=PRIME_BOUND
    )
    if (
        visible_certificate["combined_exact_rank_over_F3"] != 9
        or augmented_certificate["combined_exact_rank_over_F3"] != 9
        or visible_certificate["independent_subset_indices_one_based"] != EXPECTED_PIVOTS
        or augmented_certificate["independent_subset_indices_one_based"] != EXPECTED_PIVOTS
    ):
        raise AssertionError("the finite-reduction non-promotion certificate changed")
    return {
        "status": "exact specialization relation and non-promotion audit completed",
        "component_parameter": "r=2",
        "mestre_parameter": "T=1",
        "short_weierstrass_coefficients": [rational_text(value) for value in coefficients],
        "exact_affine_relations": [
            "P1=V(0,-)-V(0,+)+V(1,+)",
            "P2=-V(0,-)-V(1,+)-V((42-z)/6,-)",
        ],
        "visible_mod3_rank": 9,
        "augmented_mod3_rank": 9,
        "finite_reduction_pivots": EXPECTED_PIVOTS,
        "pari_diagnostics": pari_diagnostics(coefficients, (*visible, *affine)) if include_pari else {"status": "skipped by caller"},
        "conclusion": (
            "at this specialization both affine points are exactly in the visible subgroup; "
            "it cannot witness a new Mordell-Weil direction"
        ),
        "not_established": [
            "a generic relation on the rational component",
            "saturation of the visible rank-nine sublattice",
            "a Shioda Gram matrix or generic rank statement",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pari", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(include_pari=not args.skip_pari), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
