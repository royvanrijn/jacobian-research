#!/usr/bin/env python3
"""Exact two-anchor exceptional-direction transport in the Fermigier family.

This is a structural calculation, not a parameter score sweep.  It identifies
the exceptional quotients at the rank-22 and rank-20 anchors, classifies every
cross-anchor affine interpolant, applies discriminant collision tests to the
one-parameter quadratic and Mobius interpolant families, and constructs the
actual biquadratic fiber products for every pair of independently exceptional
affine directions.

The polynomial coordinate is the legacy literal shift ``T=s=2u``.  Anchor
metadata always records the canonical adapter coordinate ``u`` as well.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
from itertools import combinations
import json
import multiprocessing as mp
from pathlib import Path
import platform
import sys
from typing import Any, Iterable, Sequence

import sympy as sp

from ecsearch.fermigier_near_miss import canonical_ratpoints_output
from ecsearch.fermigier_rank import (
    parse_ratpoints_output,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (
    IndependenceCertificate,
    add_rational_points,
    negate_rational_point,
    select_independent_subset,
    verify_independence_certificate,
)
from search_fermigier_rank22_accidental_slices import (
    FERMIGIER_BIVARIATE_COEFFICIENTS,
    published_accidental_points,
    select_reconstruction_convention,
)


Q = Fraction
E22_U = Q(19754, 39)
E22_T = 2 * E22_U
RANK20_U = Q(28917, 20)
RANK20_T = 2 * RANK20_U
RANK20_EXCEPTIONAL_X = (
    Q(-8545),
    Q(23004, 5),
    Q(-8817, 10),
    Q(8183, 10),
    Q(76563, 10),
    Q(-69561, 20),
    Q(-431673, 70),
    Q(-408943, 110),
)
RANK20_SELECTED_INDICES = tuple(range(12)) + (12, 14, 28, 30, 40, 44, 57, 59)
E22_INDEPENDENT_EXCEPTIONAL_LABELS = tuple(f"P{index}" for index in range(13, 23))
E22_RELATION_COEFFICIENTS = (
    -1,
    -1,
    -1,
    -1,
    2,
    2,
    -1,
    -1,
    2,
    2,
    -1,
    -3,
    3,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
)
RANK20_ARTIFACT_RELATIVE = Path(
    "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"
)
RANK20_ARTIFACT_SHA256 = "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1"
RANK22_RECONSTRUCTION_RELATIVE = Path(
    "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_accidental_slices.json"
)
RANK22_RECONSTRUCTION_SHA256 = "3794f23d37685edcc1ad5c8279d48fcd247b1d31b7c88aa4939ca5efd40f79a7"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_text(value: Fraction | sp.Rational) -> str:
    value = Q(int(value.p), int(value.q)) if isinstance(value, sp.Rational) else Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def multiply_point(
    coefficients: Sequence[Fraction], point: tuple[Fraction, Fraction], scalar: int
) -> tuple[Fraction, Fraction] | None:
    if scalar < 0:
        negated = negate_rational_point(coefficients, point)
        assert negated is not None
        return multiply_point(coefficients, negated, -scalar)
    answer = None
    addend: tuple[Fraction, Fraction] | None = point
    while scalar:
        if scalar & 1:
            answer = add_rational_points(coefficients, answer, addend)
        addend = add_rational_points(coefficients, addend, addend)
        scalar //= 2
    return answer


def exact_linear_combination(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    scalars: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    if len(points) != len(scalars):
        raise ValueError("point/scalar lengths differ")
    answer = None
    for point, scalar in zip(points, scalars, strict=True):
        answer = add_rational_points(
            coefficients, answer, multiply_point(coefficients, point, scalar)
        )
    return answer


def e22_exceptional_quotient() -> tuple[dict[str, Any], tuple[tuple[str, Fraction], ...]]:
    specialization = specialize_fermigier_rank_sections(E22_U)
    _, _, reconstruction, _ = select_reconstruction_convention()
    accidentals = published_accidental_points(reconstruction)
    labels_and_x = tuple((label, point[0]) for label, point in accidentals)
    positive_points = parse_ratpoints_output(
        specialization.quartic_model,
        canonical_ratpoints_output(tuple(value for _, value in labels_and_x)),
    )[::2]
    cloud = section_and_point_cloud_differences(specialization, positive_points)
    if len(cloud) != 23:
        raise AssertionError("the E22 generic-plus-accidental cloud changed")
    selected, certificate = select_independent_subset(
        specialization.canonical_model,
        cloud,
        relation_prime=5,
        maximum_reduction_prime=2_000,
    )
    expected = tuple(range(12)) + tuple(range(13, 23))
    if selected != expected:
        raise AssertionError("the E22 exceptional quotient selection changed")
    if exact_linear_combination(
        specialization.canonical_model, cloud, E22_RELATION_COEFFICIENTS
    ) is not None:
        raise AssertionError("the exact P6/generic relation failed")
    return (
        {
            "generic_rank_lower_bound": 12,
            "accidental_source_labels": [label for label, _ in labels_and_x],
            "accidental_source_x": {
                label: rational_text(value) for label, value in labels_and_x
            },
            "independent_exceptional_labels_modulo_generic": list(
                E22_INDEPENDENT_EXCEPTIONAL_LABELS
            ),
            "exceptional_quotient_rank_lower_bound": 10,
            "independent_union_indices_zero_based": list(selected),
            "independence_certificate": certificate.to_json_object(),
            "P6_exact_relation": {
                "ordered_points": [
                    *[f"G{index}" for index in range(1, 13)],
                    *[label for label, _ in labels_and_x],
                ],
                "coefficients": list(E22_RELATION_COEFFICIENTS),
                "human_readable": (
                    "3*P6-(G1+G2+G3+G4+G7+G8+G11)"
                    "+2*(G5+G6+G9+G10)-3*G12=O"
                ),
                "verified_by_exact_group_law": True,
            },
        },
        labels_and_x,
    )


def rank20_exceptional_quotient(root: Path) -> tuple[dict[str, Any], tuple[tuple[str, Fraction], ...]]:
    path = root / RANK20_ARTIFACT_RELATIVE
    if sha256_file(path) != RANK20_ARTIFACT_SHA256:
        raise AssertionError("the pinned rank-20 artifact changed")
    artifact = json.loads(path.read_text())
    specialization = specialize_fermigier_rank_sections(RANK20_U)
    abscissas = tuple(Q(value) for value in artifact["bounded_search"]["abscissas"])
    searched = parse_ratpoints_output(
        specialization.quartic_model, canonical_ratpoints_output(abscissas)
    )
    cloud = section_and_point_cloud_differences(specialization, searched)
    selected_indices = tuple(artifact["point_cloud"]["selected_indices"])
    if selected_indices != RANK20_SELECTED_INDICES:
        raise AssertionError("the pinned rank-20 selected subset changed")
    selected_points = tuple(cloud[index] for index in selected_indices)
    certificate = IndependenceCertificate.from_json_object(
        artifact["point_cloud"]["certificate"]
    )
    verify_independence_certificate(
        specialization.canonical_model, selected_points, certificate
    )

    point_to_x: dict[tuple[Fraction, Fraction], Fraction] = {}
    for x_value in abscissas:
        pair = parse_ratpoints_output(
            specialization.quartic_model,
            canonical_ratpoints_output((x_value,)),
        )
        one_x_cloud = section_and_point_cloud_differences(specialization, pair)
        for point in one_x_cloud[12:]:
            point_to_x.setdefault(point, x_value)
    exceptional_x = tuple(point_to_x[cloud[index]] for index in selected_indices[12:])
    if exceptional_x != RANK20_EXCEPTIONAL_X:
        raise AssertionError("the rank-20 exceptional preimage list changed")
    labels_and_x = tuple(
        (f"R20E{index}", value)
        for index, value in enumerate(exceptional_x, start=1)
    )
    return (
        {
            "generic_rank_lower_bound": 12,
            "exceptional_quotient_rank_lower_bound": 8,
            "independent_union_indices_zero_based": list(selected_indices),
            "independent_exceptional_preimages": {
                label: rational_text(value) for label, value in labels_and_x
            },
            "independence_certificate": certificate.to_json_object(),
            "source_artifact": str(RANK20_ARTIFACT_RELATIVE),
            "source_artifact_sha256": RANK20_ARTIFACT_SHA256,
        },
        labels_and_x,
    )


def family_expression(T: sp.Symbol, X: sp.Symbol) -> sp.Expr:
    answer = sp.Integer(0)
    for x_degree, coefficients in enumerate(FERMIGIER_BIVARIATE_COEFFICIENTS):
        coefficient = sum(
            sp.Rational(value.numerator, value.denominator) * T**t_degree
            for t_degree, value in enumerate(coefficients)
        )
        answer += coefficient * X**x_degree
    return sp.expand(answer)


def primitive_poly(expression: sp.Expr, variable: sp.Symbol) -> sp.Poly:
    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(numerator, variable, domain=sp.QQ).primitive()[1]
    return -polynomial if polynomial.LC() < 0 else polynomial


def polynomial_line(polynomial: sp.Poly) -> str:
    return ",".join(rational_text(value) for value in reversed(polynomial.all_coeffs()))


def factor_signature(polynomial: sp.Poly) -> tuple[tuple[int, int], ...]:
    return tuple(
        (factor.degree(), exponent)
        for factor, exponent in sp.factor_list(polynomial)[1]
    )


def hyperelliptic_genus(squarefree_degree: int) -> int:
    return max(0, (squarefree_degree - 1) // 2)


def squareclass_kernel_degree(polynomial: sp.Poly) -> int:
    """Return the degree left after removing every even-power factor."""

    return sum(
        factor.degree()
        for factor, exponent in sp.factor_list(polynomial)[1]
        if exponent % 2
    )


def _mobius_pair_record(
    row: tuple[str, Fraction, str, Fraction],
) -> dict[str, Any]:
    """Classify every rational degeneration in one Mobius interpolation pencil.

    The finite chart is

        x(T) = (a(c) T + b(c)) / (c T + 1),

    where ``a(c),b(c)`` are uniquely determined by the two anchor values.  The
    missing projective value ``c=infinity`` is checked separately using the
    denominator ``T``.  Factoring the exact discriminant in ``QQ[c]`` is a
    complete rational branch-collision test for this pencil.
    """

    left_label, left_x_q, right_label, right_x_q = row
    T, X, c = sp.symbols("T X c")
    family = family_expression(T, X)
    first_anchor = sp.Rational(E22_T.numerator, E22_T.denominator)
    second_anchor = sp.Rational(RANK20_T.numerator, RANK20_T.denominator)
    left_x = sp.Rational(left_x_q.numerator, left_x_q.denominator)
    right_x = sp.Rational(right_x_q.numerator, right_x_q.denominator)

    numerator_slope = sp.cancel(
        (
            left_x * (c * first_anchor + 1)
            - right_x * (c * second_anchor + 1)
        )
        / (first_anchor - second_anchor)
    )
    numerator_intercept = sp.cancel(
        left_x * (c * first_anchor + 1) - numerator_slope * first_anchor
    )
    denominator = c * T + 1
    mobius_x = sp.cancel((numerator_slope * T + numerator_intercept) / denominator)
    cleared = sp.together(family.subs(X, mobius_x) * denominator**4).as_numer_denom()[0]
    generic = sp.Poly(cleared, T, domain=sp.QQ.frac_field(c))
    discriminant = sp.Poly(
        sp.discriminant(generic.as_expr(), T), c, domain=sp.QQ
    ).primitive()[1]
    discriminant_factors = sp.factor_list(discriminant)[1]
    signature = tuple(
        (factor.degree(), exponent) for factor, exponent in discriminant_factors
    )

    degenerations = []
    nonlinear_factors = []
    for factor, exponent in discriminant_factors:
        if factor.degree() != 1:
            nonlinear_factors.append(factor)
            continue
        value = sp.cancel(-factor.nth(0) / factor.nth(1))
        specialized = primitive_poly(generic.as_expr().subs(c, value), T)
        kernel_degree = squareclass_kernel_degree(specialized)
        if value == 0:
            status = "valid affine limit"
        elif value == -1 / first_anchor:
            status = "invalid: denominator and numerator vanish at E22 anchor"
        elif value == -1 / second_anchor:
            status = "invalid: denominator and numerator vanish at rank20 anchor"
        else:
            status = "unexpected rational discriminant root"
        degenerations.append(
            {
                "c": rational_text(value),
                "discriminant_multiplicity": exponent,
                "status": status,
                "cleared_polynomial_degree": specialized.degree(),
                "factor_signature": [list(item) for item in factor_signature(specialized)],
                "squareclass_kernel_degree": kernel_degree,
                "squareclass_genus": hyperelliptic_genus(kernel_degree),
            }
        )

    # The d=0 projective chart: x=(a*T+b)/T.  It is the c=infinity member
    # omitted by the normalization c*T+1.
    infinity_slope = sp.cancel(
        (left_x * first_anchor - right_x * second_anchor)
        / (first_anchor - second_anchor)
    )
    infinity_intercept = sp.cancel(
        left_x * first_anchor - infinity_slope * first_anchor
    )
    infinity_x = sp.cancel((infinity_slope * T + infinity_intercept) / T)
    infinity_polynomial = primitive_poly(
        sp.together(family.subs(X, infinity_x) * T**4).as_numer_denom()[0], T
    )
    infinity_kernel_degree = squareclass_kernel_degree(infinity_polynomial)

    return {
        "left": left_label,
        "right": right_label,
        "generic_degree": generic.degree(),
        "generic_genus": hyperelliptic_genus(generic.degree()),
        "discriminant_degree_in_c": discriminant.degree(),
        "discriminant_factor_signature": [list(item) for item in signature],
        "rational_degenerations": sorted(degenerations, key=lambda item: Q(item["c"])),
        "nonlinear_discriminant_factor_sha256": sha256_lines(
            polynomial_line(factor) for factor in nonlinear_factors
        ),
        "infinity_chart": {
            "denominator": "T",
            "degree": infinity_polynomial.degree(),
            "factor_signature": [
                list(item) for item in factor_signature(infinity_polynomial)
            ],
            "squareclass_kernel_degree": infinity_kernel_degree,
            "genus": hyperelliptic_genus(infinity_kernel_degree),
            "polynomial_sha256": sha256_lines((polynomial_line(infinity_polynomial),)),
        },
    }


def mobius_transport(
    e22: Sequence[tuple[str, Fraction]],
    rank20: Sequence[tuple[str, Fraction]],
    *,
    workers: int,
) -> dict[str, Any]:
    """Exhaust the projective one-parameter Mobius pencil for all 88 pairs."""

    if workers <= 0:
        raise ValueError("workers must be positive")
    rows = [
        (left_label, left_x, right_label, right_x)
        for left_label, left_x in e22
        for right_label, right_x in rank20
    ]
    if workers == 1:
        records = [_mobius_pair_record(row) for row in rows]
    else:
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
        with mp.get_context(start_method).Pool(workers) as pool:
            records = pool.map(_mobius_pair_record, rows)
    records.sort(key=lambda item: (item["left"], item["right"]))

    histogram = Counter(
        (
            record["generic_degree"],
            tuple(map(tuple, record["discriminant_factor_signature"])),
            tuple(
                (item["c"], item["discriminant_multiplicity"], item["status"])
                for item in record["rational_degenerations"]
            ),
            record["infinity_chart"]["degree"],
            tuple(map(tuple, record["infinity_chart"]["factor_signature"])),
        )
        for record in records
    )
    unexpected = [
        {"left": record["left"], "right": record["right"], **item}
        for record in records
        for item in record["rational_degenerations"]
        if item["status"] == "unexpected rational discriminant root"
    ]
    low_genus = [
        {"left": record["left"], "right": record["right"], **item}
        for record in records
        for item in record["rational_degenerations"]
        if item["status"].startswith("valid") and item["squareclass_genus"] <= 1
    ]
    manifest_lines = (
        f"{record['left']}|{record['right']}|{record['generic_degree']}|"
        f"{record['discriminant_degree_in_c']}|{record['discriminant_factor_signature']}|"
        f"{record['rational_degenerations']}|{record['nonlinear_discriminant_factor_sha256']}|"
        f"{record['infinity_chart']}"
        for record in records
    )
    return {
        "definition": (
            "the complete projective pencil x(T)=(a*T+b)/(c*T+d) through each "
            "anchor pair; d=1 is the finite chart and d=0 is checked separately"
        ),
        "pair_count": len(records),
        "histogram": [
            {
                "generic_degree": key[0],
                "discriminant_factor_signature": [list(item) for item in key[1]],
                "rational_degenerations": [
                    {"c": item[0], "multiplicity": item[1], "status": item[2]}
                    for item in key[2]
                ],
                "infinity_degree": key[3],
                "infinity_factor_signature": [list(item) for item in key[4]],
                "count": count,
            }
            for key, count in sorted(histogram.items())
        ],
        "manifest_sha256": sha256_lines(manifest_lines),
        "unexpected_rational_degenerations": unexpected,
        "low_genus_candidates": low_genus,
        "records": records,
        "interpretation": (
            "Every finite non-affine rational member is squarefree degree 10 (genus 4). "
            "The only valid rational degeneration is c=0, the already classified affine "
            "genus-2 member.  The other two rational discriminant roots put a pole at an "
            "anchor, and the d=0 chart is irreducible degree 10.  Thus no genus-0 or "
            "genus-1 Mobius transport exists in the declared pencils."
        ),
    }


def transport_polynomials(
    e22: Sequence[tuple[str, Fraction]],
    rank20: Sequence[tuple[str, Fraction]],
) -> tuple[dict[str, Any], dict[tuple[str, str], sp.Poly]]:
    T, X, curvature = sp.symbols("T X curvature")
    family = family_expression(T, X)
    first_anchor = sp.Rational(E22_T.numerator, E22_T.denominator)
    second_anchor = sp.Rational(RANK20_T.numerator, RANK20_T.denominator)
    normalized = (T - first_anchor) / (second_anchor - first_anchor)
    affine_records: list[dict[str, Any]] = []
    quadratic_records: list[dict[str, Any]] = []
    affine_polynomials: dict[tuple[str, str], sp.Poly] = {}

    for left_label, left_x_q in e22:
        for right_label, right_x_q in rank20:
            left_x = sp.Rational(left_x_q.numerator, left_x_q.denominator)
            right_x = sp.Rational(right_x_q.numerator, right_x_q.denominator)
            affine_x = left_x + (right_x - left_x) * normalized
            affine = primitive_poly(family.subs(X, affine_x), T)
            affine_factors = factor_signature(affine)
            affine_gcd_degree = sp.gcd(affine, affine.diff()).degree()
            affine_records.append(
                {
                    "left": left_label,
                    "right": right_label,
                    "degree": affine.degree(),
                    "factor_signature": [list(value) for value in affine_factors],
                    "squarefree_gcd_degree": affine_gcd_degree,
                    "squarefree_degree": affine.degree() - affine_gcd_degree,
                    "genus": hyperelliptic_genus(affine.degree()),
                    "polynomial_sha256": sha256_lines((polynomial_line(affine),)),
                }
            )
            affine_polynomials[left_label, right_label] = affine

            quadratic_x = affine_x + curvature * normalized * (normalized - 1)
            numerator = sp.together(family.subs(X, quadratic_x)).as_numer_denom()[0]
            generic_polynomial = sp.Poly(
                numerator, T, domain=sp.QQ.frac_field(curvature)
            )
            discriminant = sp.Poly(
                sp.discriminant(generic_polynomial.as_expr(), T),
                curvature,
                domain=sp.QQ,
            ).primitive()[1]
            discriminant_factors = sp.factor_list(discriminant)[1]
            roots: list[sp.Rational] = []
            signature = []
            for factor, exponent in discriminant_factors:
                signature.append((factor.degree(), exponent))
                if factor.degree() == 1:
                    roots.append(-factor.nth(0) / factor.nth(1))
            root_records = []
            for root in roots:
                specialized = primitive_poly(numerator.subs(curvature, root), T)
                specialized_signature = factor_signature(specialized)
                squarefree_part = sp.Poly(1, T, domain=sp.QQ)
                for factor, _ in sp.factor_list(specialized)[1]:
                    squarefree_part *= factor
                root_records.append(
                    {
                        "curvature": rational_text(root),
                        "degree": specialized.degree(),
                        "factor_signature": [list(value) for value in specialized_signature],
                        "squarefree_degree": squarefree_part.degree(),
                        "genus": hyperelliptic_genus(squarefree_part.degree()),
                    }
                )
            quadratic_records.append(
                {
                    "left": left_label,
                    "right": right_label,
                    "generic_degree": generic_polynomial.degree(),
                    "generic_genus": hyperelliptic_genus(generic_polynomial.degree()),
                    "discriminant_degree_in_curvature": discriminant.degree(),
                    "discriminant_factor_signature": [list(value) for value in signature],
                    "rational_collision_parameters": root_records,
                    "nonlinear_discriminant_factor_sha256": sha256_lines(
                        polynomial_line(factor)
                        for factor, _ in discriminant_factors
                        if factor.degree() > 1
                    ),
                }
            )

    affine_histogram = Counter(
        (row["degree"], tuple(map(tuple, row["factor_signature"])), row["genus"])
        for row in affine_records
    )
    quadratic_histogram = Counter(
        (
            row["generic_degree"],
            tuple(map(tuple, row["discriminant_factor_signature"])),
            tuple(item["curvature"] for item in row["rational_collision_parameters"]),
        )
        for row in quadratic_records
    )
    affine_lines = (
        f"{row['left']}|{row['right']}|{row['degree']}|{row['factor_signature']}|"
        f"{row['squarefree_gcd_degree']}|{row['polynomial_sha256']}"
        for row in affine_records
    )
    quadratic_lines = (
        f"{row['left']}|{row['right']}|{row['generic_degree']}|"
        f"{row['discriminant_degree_in_curvature']}|{row['discriminant_factor_signature']}|"
        f"{row['rational_collision_parameters']}|{row['nonlinear_discriminant_factor_sha256']}"
        for row in quadratic_records
    )
    return (
        {
            "cross_anchor_pair_count": len(affine_records),
            "affine": {
                "definition": "the unique x(T)=a*T+b through each anchor pair",
                "histogram": [
                    {
                        "degree": key[0],
                        "factor_signature": [list(value) for value in key[1]],
                        "genus": key[2],
                        "count": count,
                    }
                    for key, count in sorted(affine_histogram.items())
                ],
                "manifest_sha256": sha256_lines(affine_lines),
                "low_genus_candidates": [
                    row for row in affine_records if row["genus"] <= 1
                ],
            },
            "quadratic": {
                "definition": (
                    "x(T)=L(T)+k*z*(z-1), z=(T-T_E22)/(T_rank20-T_E22); "
                    "all rational branch-collision k are roots of Disc_T"
                ),
                "histogram": [
                    {
                        "generic_degree": key[0],
                        "discriminant_factor_signature": [list(value) for value in key[1]],
                        "rational_collision_parameters": list(key[2]),
                        "count": count,
                    }
                    for key, count in sorted(quadratic_histogram.items())
                ],
                "manifest_sha256": sha256_lines(quadratic_lines),
                "low_genus_candidates": [
                    {
                        "left": row["left"],
                        "right": row["right"],
                        **collision,
                    }
                    for row in quadratic_records
                    for collision in row["rational_collision_parameters"]
                    if collision["genus"] <= 1
                ],
                "interpretation": (
                    "The only rational collision is k=0, which is exactly the affine genus-2 case."
                ),
            },
        },
        affine_polynomials,
    )


def fiber_products(
    affine_polynomials: dict[tuple[str, str], sp.Poly],
    independent_e22_labels: set[str],
) -> dict[str, Any]:
    keys = sorted(
        key for key in affine_polynomials if key[0] in independent_e22_labels
    )
    records = []
    for first, second in combinations(keys, 2):
        left = affine_polynomials[first]
        right = affine_polynomials[second]
        common = sp.gcd(left, right)
        third = (left * right).exquo(common * common)
        third_squarefree = sp.gcd(third, third.diff()).degree() == 0
        third_genus = hyperelliptic_genus(third.degree())
        connected = common.degree() == 0 and left != right
        fiber_genus = 2 + 2 + third_genus if connected else None
        records.append(
            {
                "first": list(first),
                "second": list(second),
                "endpoint_stratum": (
                    "shared-E22-endpoint"
                    if first[0] == second[0]
                    else (
                        "shared-rank20-endpoint"
                        if first[1] == second[1]
                        else "distinct-at-both-anchors"
                    )
                ),
                "common_branch_gcd_degree": common.degree(),
                "third_quotient_squarefree_degree": third.degree(),
                "third_quotient_is_squarefree": third_squarefree,
                "third_quotient_genus": third_genus,
                "fiber_product_connected": connected,
                "fiber_product_genus": fiber_genus,
            }
        )
    histogram = Counter(
        (
            row["common_branch_gcd_degree"],
            row["third_quotient_squarefree_degree"],
            row["third_quotient_genus"],
            row["fiber_product_genus"],
        )
        for row in records
    )
    lines = (
        f"{row['first']}|{row['second']}|{row['endpoint_stratum']}|{row['common_branch_gcd_degree']}|"
        f"{row['third_quotient_squarefree_degree']}|{row['third_quotient_genus']}|"
        f"{row['fiber_product_genus']}"
        for row in records
    )
    return {
        "construction": (
            "actual biquadratic cover y1^2=f1(T), y2^2=f2(T); the third character "
            "quotient is y3^2=squarefree(f1*f2), never a product-square surrogate"
        ),
        "pair_scope": (
            "all unordered pairs of the 80 affine transports from P13..P22 to the eight rank20 directions"
        ),
        "pair_count": len(records),
        "endpoint_strata": dict(
            sorted(Counter(row["endpoint_stratum"] for row in records).items())
        ),
        "histogram": [
            {
                "common_branch_gcd_degree": key[0],
                "third_quotient_squarefree_degree": key[1],
                "third_quotient_genus": key[2],
                "fiber_product_genus": key[3],
                "count": count,
            }
            for key, count in sorted(histogram.items())
        ],
        "manifest_sha256": sha256_lines(lines),
        "low_genus_third_quotients": [
            row for row in records if row["third_quotient_genus"] <= 1
        ],
        "interpretation": (
            "Every pair, including pairs sharing one anchor endpoint, has disjoint six-point "
            "branch loci; the third quotient has genus 5 and the connected fiber product has genus 9."
        ),
    }


def result_digest(artifact: dict[str, Any]) -> str:
    stable = {
        "anchors": artifact["anchors"],
        "exceptional_quotients": artifact["exceptional_quotients"],
        "transport": artifact["transport"],
        "fiber_products": artifact["fiber_products"],
        "outcome": artifact["outcome"],
        "sources": artifact["sources"],
    }
    return hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def run(root: Path, *, workers: int = 4) -> dict[str, Any]:
    reconstruction_path = root / RANK22_RECONSTRUCTION_RELATIVE
    if sha256_file(reconstruction_path) != RANK22_RECONSTRUCTION_SHA256:
        raise AssertionError("the pinned rank-22 reconstruction artifact changed")
    e22_record, e22_all = e22_exceptional_quotient()
    rank20_record, rank20 = rank20_exceptional_quotient(root)
    transport, affine_polynomials = transport_polynomials(e22_all, rank20)
    transport["mobius"] = mobius_transport(e22_all, rank20, workers=workers)
    products = fiber_products(
        affine_polynomials, set(E22_INDEPENDENT_EXCEPTIONAL_LABELS)
    )
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-exceptional-transport.v1",
        "status": "complete exact bounded structural classification",
        "claim_level": "exact computation; no new section or specialization found",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "anchors": {
            "E22": {
                "canonical_adapter_u": rational_text(E22_U),
                "literal_shift_T": rational_text(E22_T),
                "certified_rank_lower_bound": 22,
            },
            "rank20": {
                "canonical_adapter_u": rational_text(RANK20_U),
                "literal_shift_T": rational_text(RANK20_T),
                "certified_rank_lower_bound": 20,
            },
        },
        "exceptional_quotients": {"E22": e22_record, "rank20": rank20_record},
        "transport": transport,
        "fiber_products": products,
        "outcome": {
            "independent_E22_exceptional_directions": 10,
            "independent_rank20_exceptional_directions": 8,
            "affine_low_genus_candidates": len(
                transport["affine"]["low_genus_candidates"]
            ),
            "quadratic_low_genus_candidates": len(
                transport["quadratic"]["low_genus_candidates"]
            ),
            "mobius_low_genus_candidates": len(
                transport["mobius"]["low_genus_candidates"]
            ),
            "fiber_product_low_genus_quotients": len(
                products["low_genus_third_quotients"]
            ),
            "new_sections": 0,
            "new_specializations": 0,
            "target_met": False,
        },
        "sources": {
            "rank22_reconstruction_artifact": str(RANK22_RECONSTRUCTION_RELATIVE),
            "rank22_reconstruction_artifact_sha256": RANK22_RECONSTRUCTION_SHA256,
            "rank20_artifact": str(RANK20_ARTIFACT_RELATIVE),
            "rank20_artifact_sha256": RANK20_ARTIFACT_SHA256,
            "script_sha256": sha256_file(Path(__file__)),
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/analyze_fermigier_exceptional_transport.py"
        ),
    }
    artifact["result_sha256"] = result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts/generated-results/elliptic-curves/elliptic_fermigier_exceptional_transport.json",
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    artifact = run(root, workers=args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
