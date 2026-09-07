#!/usr/bin/env python3
"""Independently verify the six-root rank-at-least-13 construction.

The family formulas live in :mod:`mestre_rank13_02393128133175`.  This
verifier does not use its pre-expanded quartic as an input: it reconstructs
Mestre's square approximant from the six roots over ``QQ[T]``, removes the
fixed square content, and checks the displayed sections symbolically.

At ``u=1`` it then maps the affine sections and one split point at infinity
to the short Jacobian.  Exact reductions modulo good primes certify thirteen
independent specialized points.  Specialization therefore proves generic
rank at least thirteen.  The discriminant calculation is conductor geometry,
not an exact conductor formula for every specialization.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import isqrt
from typing import Any, Sequence

import sympy as sp

import mestre_rank13_02393128133175 as formulas
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points, quartic_point_to_short_jacobian
from search_mestre_root_tuple_scale import point_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
ROOTS = (0, 23, 93, 128, 133, 175)
CENTER = Q(92)
CENTERED_ROOTS = (-92, -69, 1, 36, 41, 83)
BASE_CHANGE_CONSTANT = 14_406
SPECIALIZATION_U = Q(1)
SPECIALIZATION_T = Q(14_405, 2)
QUARTIC_SQUARE_CONTENT = 5_760_000
QUARTIC_SQUARE_SCALE = 2_400


# These are the six-section diagnostic inventory supplied by the exact linear-
# section calculation.  Only the first section is used in the rank theorem;
# this verifier checks their identities but does not classify every possible
# linear-abscissa section.
# Each tuple is (slope numerator over 21, intercept, numerator of y, y-denom).
COMPANION_DATA = (
    (
        31,
        Q(619, 3),
        (520, 148_862, 17_234_623, 577_212_405),
        147,
    ),
    (
        -31,
        Q(619, 3),
        (520, -148_862, 17_234_623, -577_212_405),
        147,
    ),
    (19, Q(304, 3), (80, -7_448, -55_468, -52_005_660), 147),
    (-19, Q(304, 3), (80, 7_448, -55_468, 52_005_660), 147),
    (1, Q(115, 3), (440, 2_254, -3_750_607, 6_792_429), 147),
    (-1, Q(115, 3), (440, -2_254, -3_750_607, -6_792_429), 147),
)


def rational_text(value: Fraction | sp.Rational | int) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def coefficient_digest(coefficients: Sequence[sp.Expr]) -> str:
    payload = "".join(f"{sp.Rational(value)}\n" for value in coefficients)
    return hashlib.sha256(payload.encode()).hexdigest()


def _primitive_positive(poly: sp.Poly) -> tuple[sp.Rational, sp.Poly]:
    content, primitive = poly.primitive()
    if primitive.LC() < 0:
        content = -content
        primitive = -primitive
    return sp.Rational(content), primitive


@dataclass(frozen=True)
class SymbolicConstruction:
    x: sp.Symbol
    t: sp.Symbol
    u: sp.Symbol
    root_polynomial: sp.Expr
    shifted_product: sp.Expr
    square_approximant: sp.Expr
    primitive_quartic: sp.Expr
    quartic_coefficients: tuple[sp.Expr, ...]
    discriminant_content: sp.Rational
    discriminant_core: sp.Poly
    base_changed_content: sp.Rational
    base_changed_core: sp.Poly
    infinity_x: sp.Expr
    infinity_y: sp.Expr


def reconstruct_symbolically() -> SymbolicConstruction:
    """Rebuild the quartic, sections, split infinity, and discriminant."""

    x, t, u = sp.symbols("x t u")
    field_t = sp.QQ.frac_field(t)
    root_polynomial = sp.prod(x - root for root in ROOTS)
    shifted_product = sp.expand(
        root_polynomial.subs(x, x - t) * root_polynomial.subs(x, x + t)
    )
    product_poly = sp.Poly(shifted_product, x, domain=field_t)

    approximant = x**6
    for lower_degree in range(5, -1, -1):
        current = sp.Poly(approximant, x, domain=field_t)
        correction = sp.cancel(
            (
                product_poly.nth(6 + lower_degree)
                - (current * current).nth(6 + lower_degree)
            )
            / 2
        )
        approximant = sp.expand(approximant + correction * x**lower_degree)

    remainder_over_t2 = sp.cancel(
        (approximant**2 - shifted_product) / t**2
    )
    remainder_poly = sp.Poly(remainder_over_t2, x, t, domain=sp.QQ)
    square_content, primitive_multivariate = remainder_poly.primitive()
    if square_content != QUARTIC_SQUARE_CONTENT:
        raise AssertionError("the independently derived quartic content changed")
    if isqrt(int(square_content)) != QUARTIC_SQUARE_SCALE:
        raise AssertionError("the fixed quartic content ceased to be a square")
    primitive_quartic = primitive_multivariate.as_expr()
    quartic_poly = sp.Poly(primitive_quartic, x, domain=field_t)
    if quartic_poly.degree() != 4:
        raise AssertionError("the Mestre obstruction did not leave a quartic")
    quartic_coefficients = tuple(
        sp.factor(quartic_poly.nth(index)) for index in range(5)
    )

    expected_coefficients = (
        9 * t**6
        - 253_406 * t**4
        + 1_434_086_185 * t**2
        + 7_050_150_764_944,
        12 * (276 * t**4 - 2_861_579 * t**2 - 31_373_984_992),
        -3 * (6 * t**4 - 158_072 * t**2 - 2_245_309_213),
        -36 * (92 * t**2 + 1_377_831),
        9 * (t**2 + BASE_CHANGE_CONSTANT),
    )
    if any(
        sp.expand(actual - expected) != 0
        for actual, expected in zip(quartic_coefficients, expected_coefficients)
    ):
        raise AssertionError("the reconstructed quartic missed the formula module")

    # All twelve paired-root sections follow directly from a root of one
    # shifted factor.  Check the actual primitive ordinates, not just values at
    # a finite set of parameters.
    for root in ROOTS:
        for sign in (-1, 1):
            abscissa = root + sign * t
            ordinate = sp.cancel(
                approximant.subs(x, abscissa) / (t * QUARTIC_SQUARE_SCALE)
            )
            if sp.cancel(primitive_quartic.subs(x, abscissa) - ordinate**2) != 0:
                raise AssertionError("a generic paired-root section failed")

    # Verify the full six-section diagnostic inventory.  These identities do
    # not assert six new independent Mordell-Weil directions.
    for slope_numerator, intercept, y_numerators, y_denominator in COMPANION_DATA:
        abscissa = sp.Rational(slope_numerator, 21) * t + sp.Rational(intercept)
        y3, y2, y1, y0 = y_numerators
        ordinate = sp.Rational(1, y_denominator) * (
            y3 * t**3 + y2 * t**2 + y1 * t + y0
        )
        if sp.cancel(primitive_quartic.subs(x, abscissa) - ordinate**2) != 0:
            raise AssertionError("a nonvisible companion identity failed")

    # The base change parametrizes a^2 = 9(T^2+C).  Evaluate the homogeneous
    # quartic covariants at [x:z]=[1:0] and verify their image on the Jacobian.
    base_parameter = (BASE_CHANGE_CONSTANT - u**2) / (2 * u)
    leading_square = 3 * (BASE_CHANGE_CONSTANT + u**2) / (2 * u)
    base_changed_coefficients = tuple(
        sp.cancel(value.subs(t, base_parameter)) for value in quartic_coefficients
    )
    e, d, c, b, a = base_changed_coefficients
    if sp.cancel(a - leading_square**2) != 0:
        raise AssertionError("the base change did not split quartic infinity")
    invariant_i = sp.cancel(12 * a * e - 3 * b * d + c**2)
    invariant_j = sp.cancel(
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    infinity_x = sp.cancel(36 * g0 / a)
    infinity_y = sp.cancel(
        54 * (a * g1 - b * g0) / leading_square**3
    )
    if sp.cancel(
        infinity_y**2
        - infinity_x**3
        + 27 * invariant_i * infinity_x
        + 27 * invariant_j
    ) != 0:
        raise AssertionError("the generic split-infinity image missed the Jacobian")

    # The nonconstant discriminant factor is enough for exact generic
    # singular-fiber geometry.  It is not, by itself, a specialized conductor.
    e0, d0, c0, b0, a0 = quartic_coefficients
    invariant_i_t = sp.expand(12 * a0 * e0 - 3 * b0 * d0 + c0**2)
    invariant_j_t = sp.expand(
        72 * a0 * c0 * e0
        + 9 * b0 * c0 * d0
        - 27 * a0 * d0**2
        - 27 * b0**2 * e0
        - 2 * c0**3
    )
    discriminant = sp.Poly(
        sp.expand((4 * invariant_i_t**3 - invariant_j_t**2) / 27),
        t,
        domain=sp.QQ,
    )
    discriminant_content, discriminant_core = _primitive_positive(discriminant)
    if (
        discriminant_core.degree() != 20
        or any(discriminant_core.nth(index) for index in range(1, 21, 2))
        or sp.gcd(discriminant_core, discriminant_core.diff()).degree() != 0
        or not discriminant_core.is_irreducible
    ):
        raise AssertionError("the original degree-20 discriminant geometry changed")

    pulled_back = sp.Poly(
        sp.cancel(
            (2 * u) ** 20
            * discriminant_core.as_expr().subs(t, base_parameter)
        ),
        u,
        domain=sp.QQ,
    )
    base_changed_content, base_changed_core = _primitive_positive(pulled_back)
    if (
        base_changed_core.degree() != 40
        or any(base_changed_core.nth(index) for index in range(1, 41, 2))
        or sp.gcd(base_changed_core, base_changed_core.diff()).degree() != 0
        or not base_changed_core.is_irreducible
    ):
        raise AssertionError("the pulled-back degree-40 frontier changed")
    for degree in range(41):
        if (
            base_changed_core.nth(40 - degree)
            * BASE_CHANGE_CONSTANT ** (40 - degree)
            != BASE_CHANGE_CONSTANT**20 * base_changed_core.nth(degree)
        ):
            raise AssertionError("the base-change reciprocity identity failed")

    return SymbolicConstruction(
        x=x,
        t=t,
        u=u,
        root_polynomial=root_polynomial,
        shifted_product=shifted_product,
        square_approximant=approximant,
        primitive_quartic=primitive_quartic,
        quartic_coefficients=quartic_coefficients,
        discriminant_content=discriminant_content,
        discriminant_core=discriminant_core,
        base_changed_content=base_changed_content,
        base_changed_core=base_changed_core,
        infinity_x=infinity_x,
        infinity_y=infinity_y,
    )


def companion_point(
    parameter_t: Fraction, companion_index: int = 0
) -> tuple[Fraction, Fraction]:
    slope_numerator, intercept, y_numerators, y_denominator = COMPANION_DATA[
        companion_index
    ]
    parameter_t = Q(parameter_t)
    y3, y2, y1, y0 = y_numerators
    return (
        Q(slope_numerator, 21) * parameter_t + Q(intercept),
        Q(
            y3 * parameter_t**3
            + y2 * parameter_t**2
            + y1 * parameter_t
            + y0,
            y_denominator,
        ),
    )


def split_infinity_point_at_one(
    construction: SixRootMestreConstruction,
) -> tuple[Fraction, Fraction]:
    coefficients = construction.primitive_quartic_coefficients(SPECIALIZATION_T)
    _, d, c, b, a = coefficients
    square_root = Q(3 * (BASE_CHANGE_CONSTANT + 1), 2)
    if square_root**2 != a:
        raise AssertionError("u=1 did not split the leading coefficient")
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    return (
        36 * g0 / a,
        54 * (a * g1 - b * g0) / square_root**3,
    )


def specialization_certificate() -> dict[str, Any]:
    """Return exact mod-3 certificates at the good specialization ``u=1``."""

    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    coefficients = construction.primitive_jacobian_coefficients(SPECIALIZATION_T)
    visible_quartic = primitive_visible_points(construction, SPECIALIZATION_T)
    chosen_companion = companion_point(SPECIALIZATION_T)
    affine_quartic = visible_quartic + (chosen_companion,)
    affine_jacobian = tuple(
        quartic_point_to_short_jacobian(construction, SPECIALIZATION_T, point)
        for point in affine_quartic
    )
    infinity = split_infinity_point_at_one(construction)
    all_points = affine_jacobian + (infinity,)
    all_companion_quartic = tuple(
        companion_point(SPECIALIZATION_T, index)
        for index in range(len(COMPANION_DATA))
    )
    all_companion_jacobian = tuple(
        quartic_point_to_short_jacobian(construction, SPECIALIZATION_T, point)
        for point in visible_quartic + all_companion_quartic
    )

    # Cross-check, rather than assume, the separate formula module.
    if formulas.ROOTS != ROOTS or formulas.BASE_CHANGE_CONSTANT != BASE_CHANGE_CONSTANT:
        raise AssertionError("the formula-module family identifier changed")
    if formulas.base_parameter(SPECIALIZATION_U) != SPECIALIZATION_T:
        raise AssertionError("the formula-module base change changed")
    if formulas.linear_extra_point(SPECIALIZATION_T) != chosen_companion:
        raise AssertionError("the formula-module chosen companion changed")
    if formulas.split_infinity_jacobian_point(SPECIALIZATION_U) != infinity:
        raise AssertionError("the formula-module split infinity changed")
    if formulas.known_jacobian_points(SPECIALIZATION_U) != all_points:
        raise AssertionError("the formula-module point order changed")

    visible_certificate = mod3_independence_certificate(
        coefficients, affine_jacobian[:12], prime_bound=251
    )
    affine_certificate = mod3_independence_certificate(
        coefficients, affine_jacobian, prime_bound=251
    )
    all_companion_certificate = mod3_independence_certificate(
        coefficients, all_companion_jacobian, prime_bound=251
    )
    full_certificate = mod3_independence_certificate(
        coefficients, all_points, prime_bound=251
    )
    if visible_certificate["combined_exact_rank_over_F3"] != 11:
        raise AssertionError("the visible mod-3 image dimension changed")
    if affine_certificate["combined_exact_rank_over_F3"] != 12:
        raise AssertionError("the affine mod-3 image dimension changed")
    if all_companion_certificate["combined_exact_rank_over_F3"] != 12:
        raise AssertionError("the six-companion mod-3 image dimension changed")
    if full_certificate["combined_exact_rank_over_F3"] != 13:
        raise AssertionError("the split-infinity rank certificate failed")
    if 14 not in full_certificate["independent_subset_indices_one_based"]:
        raise AssertionError("split infinity ceased to be a certificate pivot")

    return {
        "parameter_u": "1",
        "parameter_T": rational_text(SPECIALIZATION_T),
        "short_jacobian_coefficients": [rational_text(value) for value in coefficients],
        "quartic_discriminant_nonzero": (
            construction.primitive_quartic_discriminant(SPECIALIZATION_T) != 0
        ),
        "displayed_point_count": len(all_points),
        "point_sha256": point_digest(all_points),
        "visible_mod3": visible_certificate,
        "affine_with_one_companion_mod3": affine_certificate,
        "affine_with_all_six_companions_mod3": all_companion_certificate,
        "with_split_infinity_mod3": full_certificate,
    }


def build_verification(*, full_polynomials: bool = False) -> dict[str, Any]:
    symbolic = reconstruct_symbolically()
    specialization = specialization_certificate()
    centered_polynomial = sp.Poly(
        sp.prod(symbolic.x - root for root in CENTERED_ROOTS),
        symbolic.x,
        domain=sp.QQ,
    )
    c4 = centered_polynomial.nth(4)
    c3 = centered_polynomial.nth(3)
    c1 = centered_polynomial.nth(1)
    if 2 * c1 != c3 * c4:
        raise AssertionError("the intrinsic centered Mestre relation failed")

    discriminant_coefficients = tuple(
        symbolic.discriminant_core.nth(index) for index in range(21)
    )
    base_changed_coefficients = tuple(
        symbolic.base_changed_core.nth(index) for index in range(41)
    )
    result: dict[str, Any] = {
        "status": "exact generic algebraic rank lower bound",
        "claim": {
            "generic_rank_lower_bound": 13,
            "reason": (
                "thirteen selected sections specialize at u=1 to points whose "
                "exact good-reduction images are independent modulo 3"
            ),
            "not_claimed": [
                "generic rank equals 13",
                "the six companion sections give six independent directions",
                "degree-40 discriminant geometry determines every exact conductor",
            ],
        },
        "family": {
            "roots": list(ROOTS),
            "mean": rational_text(CENTER),
            "centered_roots": list(CENTERED_ROOTS),
            "centered_coefficients": {
                "c4": rational_text(c4),
                "c3": rational_text(c3),
                "c1": rational_text(c1),
                "identity_2c1_equals_c3c4": True,
            },
            "quartic_fixed_square_content": str(QUARTIC_SQUARE_CONTENT),
            "quartic_fixed_square_scale": str(QUARTIC_SQUARE_SCALE),
            "primitive_quartic_coefficients": [
                str(value) for value in symbolic.quartic_coefficients
            ],
        },
        "sections": {
            "paired_root_section_count": 12,
            "all_paired_root_identities_symbolic": True,
            "displayed_nonvisible_linear_companion_count": len(COMPANION_DATA),
            "companion_identities_symbolic": True,
            "companion_inventory": [
                {
                    "x": (
                        f"({slope}/21)*T+({rational_text(intercept)})"
                    ),
                    "y_numerator_coefficients_descending": list(y_numerators),
                    "y_denominator": y_denominator,
                }
                for slope, intercept, y_numerators, y_denominator in COMPANION_DATA
            ],
            "rank_certificate_uses_companion_index_one_based": 1,
            "split_infinity_identity_symbolic": True,
        },
        "base_change": {
            "constant": BASE_CHANGE_CONSTANT,
            "formula": "T=(14406-u^2)/(2u)",
            "leading_coefficient": "9*(T^2+14406)",
            "leading_square_root": "3*(14406+u^2)/(2u)",
        },
        "specialization_certificate": specialization,
        "discriminant_geometry": {
            "original_discriminant_content": rational_text(
                symbolic.discriminant_content
            ),
            "original_core_degree_in_T": symbolic.discriminant_core.degree(),
            "original_core_even": True,
            "original_core_squarefree": True,
            "original_core_irreducible_over_Q": True,
            "original_core_coefficients_ascending_sha256": coefficient_digest(
                discriminant_coefficients
            ),
            "base_changed_clear_factor": "(2u)^20",
            "base_changed_content_removed": rational_text(
                symbolic.base_changed_content
            ),
            "base_changed_core_degree_in_u": symbolic.base_changed_core.degree(),
            "base_changed_core_even": True,
            "base_changed_core_squarefree": True,
            "base_changed_core_irreducible_over_Q": True,
            "base_changed_core_coefficients_ascending_sha256": coefficient_digest(
                base_changed_coefficients
            ),
            "reciprocity": "u^40*P(14406/u)=14406^20*P(u)",
            "conductor_caveat": (
                "This is the nonconstant discriminant frontier. Exact minimal "
                "conductors also depend on local minimalization, denominators, "
                "and fixed small primes."
            ),
        },
        "kihara_comparison": {
            "this_family_generic_rank_lower_bound": 13,
            "this_family_squarefree_frontier_degree": 40,
            "kihara_generic_rank_lower_bound": 14,
            "kihara_squarefree_frontier_degree": 398,
            "frontier_degree_per_certified_generic_rank": {
                "this_family": "40/13",
                "kihara": "199/7",
            },
            "interpretation": (
                "The degree comparison is an asymptotic conductor-growth proxy, "
                "not a conductor theorem for individual fibers."
            ),
        },
    }
    if full_polynomials:
        result["discriminant_geometry"][
            "original_core_coefficients_ascending"
        ] = [str(value) for value in discriminant_coefficients]
        result["discriminant_geometry"][
            "base_changed_core_coefficients_ascending"
        ] = [str(value) for value in base_changed_coefficients]
        result["base_change"]["split_infinity_jacobian_x"] = str(
            symbolic.infinity_x
        )
        result["base_change"]["split_infinity_jacobian_y"] = str(
            symbolic.infinity_y
        )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-polynomials",
        action="store_true",
        help="include degree-20/40 coefficient arrays and generic infinity coordinates",
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    print(
        json.dumps(
            build_verification(full_polynomials=arguments.full_polynomials),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
