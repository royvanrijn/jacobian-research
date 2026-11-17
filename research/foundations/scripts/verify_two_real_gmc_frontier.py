#!/usr/bin/env python3
"""Exact regressions for the two-real-variable GMC frontier.

The theorem-level arguments are written in
``extended-geometry/TWO_REAL_GMC_FRONTIER.md``.  This script checks their
coefficient identities and performs three finite Groebner exclusions for the
most direct attempt to replace Long's third real Gaussian by the existing
circular pair.

Every expectation is computed from

    E(Z^a W^b) = delta_(a,b) a!.

The Groebner conclusions concern only the displayed finite supports.  They
are not advertised as a proof of GMC(2).
"""

from __future__ import annotations

import json
from itertools import combinations, product
from math import gcd
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "two_real_gmc_frontier.json"

Z, W, T = sp.symbols("Z W T")


def expectation(expression: sp.Expr) -> sp.Expr:
    """Apply circular Gaussian contraction exactly."""

    polynomial = sp.Poly(sp.expand(expression), Z, W)
    return sp.expand(
        sum(
            coefficient * sp.factorial(z_degree)
            for (z_degree, w_degree), coefficient in polynomial.terms()
            if z_degree == w_degree
        )
    )


def truncated_exponential_moments(
    polynomial: sp.Expr,
    multiplier: sp.Expr,
    order: int,
) -> sp.Expr:
    power = sp.Integer(1)
    answer = sp.Integer(0)
    for exponent in range(order + 1):
        if exponent:
            power = sp.expand(power * polynomial)
        answer += expectation(multiplier * power) * T**exponent / sp.factorial(
            exponent
        )
    return sp.expand(answer)


def verify_two_weight_identity() -> None:
    """Check the exact binomial/factorial reduction on several weight pairs."""

    U = Z * W
    A = 1 + 2 * U + U**2
    B = 2 - U + 3 * U**2
    for positive_weight, negative_weight in ((1, 1), (2, 3), (4, 6)):
        divisor = gcd(positive_weight, negative_weight)
        a0 = positive_weight // divisor
        b0 = negative_weight // divisor
        period = a0 + b0
        radial_exponent = divisor * a0 * b0
        polynomial = (
            Z**positive_weight * A + W**negative_weight * B
        )
        power = sp.Integer(1)
        for moment in range(1, 2 * period + 1):
            power = sp.expand(power * polynomial)
            actual = expectation(power)
            if moment % period:
                expected = 0
            else:
                r = moment // period
                radial = sp.expand(
                    (U**radial_exponent * A**b0 * B**a0) ** r
                )
                expected = (
                    sp.binomial(moment, b0 * r) * expectation(radial)
                )
            assert sp.expand(actual - expected) == 0


def verify_affine_source_identity() -> None:
    """Regress the one-source Lagrange identity through a finite order."""

    z = sp.symbols("z")
    h = 1 + z + 2 * z**2
    v = z - z**2
    polynomial = W * h.subs(z, Z) + v.subs(z, Z)
    order = 7

    # Solve g=t*h(g) coefficient by coefficient by fixed-point iteration.
    g = sp.Integer(0)
    for _ in range(order + 1):
        g = sp.series(T * h.subs(z, g), T, 0, order + 1).removeO().expand()

    right = sp.exp(T * v.subs(z, g)) / (
        1 - T * sp.diff(h, z).subs(z, g)
    )
    right = sp.series(right, T, 0, order + 1).removeO().expand()
    left = truncated_exponential_moments(polynomial, sp.Integer(1), order)
    assert sp.expand(left - right) == 0

    right_z = sp.series(g * right, T, 0, order + 1).removeO().expand()
    left_z = truncated_exponential_moments(polynomial, Z, order)
    assert sp.expand(left_z - right_z) == 0


def verify_quadratic_tilt_identity() -> None:
    """Check the quadratic Gaussian determinant and tilted-moment formulas."""

    t, s1, s2 = sp.symbols("t s1 s2")
    x1, x2 = sp.symbols("x1 x2")
    # A nonzero symmetric nilpotent matrix for the standard complex bilinear
    # form.  The resulting quadratic is a scalar multiple of Z^2.
    A = sp.Matrix([[1, sp.I], [sp.I, -1]])
    assert A**2 == sp.zeros(2)
    b = sp.Matrix([1, sp.I])
    assert (b.T * b)[0] == 0
    assert (b.T * A * b)[0] == 0

    identity = sp.eye(2)
    resolvent = identity + 2 * t * A
    assert sp.simplify((identity - 2 * t * A) * resolvent) == identity
    assert sp.expand((identity - 2 * t * A).det()) == 1

    source = sp.Matrix([s1, s2])
    tilted = sp.exp(
        sp.Rational(1, 2) * (source.T * resolvent * source)[0]
        + t * (source.T * resolvent * b)[0]
    )

    polynomial = (
        (sp.Matrix([x1, x2]).T * A * sp.Matrix([x1, x2]))[0]
        + (b.T * sp.Matrix([x1, x2]))[0]
    )
    # Convert the real-coordinate test polynomial to circular coordinates.
    circular = sp.expand(
        polynomial.subs(
            {
                x1: (Z + W) / sp.sqrt(2),
                x2: (Z - W) / (sp.I * sp.sqrt(2)),
            },
            simultaneous=True,
        )
    )
    for q_real, q_source in (
        (sp.Integer(1), sp.Integer(1)),
        (x1, sp.diff(tilted, s1).subs({s1: 0, s2: 0})),
        (
            x1 * x2,
            sp.diff(tilted, s1, s2).subs({s1: 0, s2: 0}),
        ),
    ):
        circular_q = sp.expand(
            q_real.subs(
                {
                    x1: (Z + W) / sp.sqrt(2),
                    x2: (Z - W) / (sp.I * sp.sqrt(2)),
                },
                simultaneous=True,
            )
        )
        left = truncated_exponential_moments(circular, circular_q, 6)
        right = sp.series(q_source, t, 0, 7).removeO().subs(t, T)
        assert sp.simplify(left - right) == 0


def long_like_exclusion(max_degree: int, moment_order: int) -> dict[str, object]:
    """Exclude ``W(1+Z)+W^2 k(Z)`` for one bounded polynomial support."""

    coefficients = sp.symbols(f"c0:{max_degree + 1}")
    k = sum(coefficients[j] * Z**j for j in range(max_degree + 1))
    polynomial = W * (1 + Z) + W**2 * k
    equations: list[sp.Expr] = []
    power = sp.Integer(1)
    for _ in range(moment_order):
        power = sp.expand(power * polynomial)
        equations.append(expectation(power))
    basis = sp.groebner(equations, *coefficients, order="grevlex")
    assert list(basis) == [1]
    return {
        "ansatz": f"P=W(1+Z)+W^2*sum_(j=0)^{max_degree} c_j Z^j",
        "coefficient_count": max_degree + 1,
        "moments_used": moment_order,
        "groebner_basis": ["1"],
        "conclusion": "no complex coefficients solve the displayed finite moment system",
    }


def cubic_weight_component(weight: int) -> list[tuple[sp.Symbol, int, int]]:
    """Return the complete total-degree-three component of one weight."""

    radial_degree = (3 - abs(weight)) // 2
    prefix = (
        f"negative_{-weight}"
        if weight < 0
        else f"positive_{weight}"
        if weight > 0
        else "zero"
    )
    terms = []
    for radial_power in range(radial_degree + 1):
        coefficient = sp.Symbol(f"{prefix}_{radial_power}")
        if weight >= 0:
            terms.append(
                (
                    coefficient,
                    weight + radial_power,
                    radial_power,
                )
            )
        else:
            terms.append(
                (
                    coefficient,
                    radial_power,
                    -weight + radial_power,
                )
            )
    return terms


def cubic_three_weight_exclusions() -> dict[str, object]:
    """Exclude every mixed-sign cubic with exactly three nonzero weights.

    For a support ``k1 < k2 < k3``, multiplication of ``P`` by a nonzero
    scalar and the torus change ``(Z,W) -> (r Z,r^-1 W)`` normalize one
    selected nonzero coefficient in each of the first two components to one.
    The charts choose a nonzero coefficient in every component; a localization
    variable keeps the selected coefficient of the third component nonzero.
    """

    weights = tuple(range(-3, 4))
    records = []
    total_charts = 0
    for support in combinations(weights, 3):
        if min(support) >= 0 or max(support) <= 0:
            continue
        components = {
            weight: cubic_weight_component(weight) for weight in support
        }
        chart_orders = []
        for chart in product(
            *(range(len(components[weight])) for weight in support)
        ):
            selected = [
                components[weight][choice][0]
                for weight, choice in zip(support, chart)
            ]
            normalization = {selected[0]: 1, selected[1]: 1}
            polynomial = sp.expand(
                sum(
                    coefficient * Z**z_degree * W**w_degree
                    for weight in support
                    for coefficient, z_degree, w_degree in components[weight]
                ).subs(normalization)
            )
            localizer = sp.Symbol("localizer")
            equations = [localizer * selected[2] - 1]
            variables = sorted(
                polynomial.free_symbols - {Z, W}, key=str
            ) + [localizer]
            power = sp.Integer(1)
            certificate_order = None
            for moment in range(1, 9):
                power = sp.expand(power * polynomial)
                equation = expectation(power)
                if equation != 0:
                    equations.append(equation)
                if moment not in (3, 4, 5, 6, 8):
                    continue
                basis = sp.groebner(
                    equations,
                    *variables,
                    order="grevlex",
                )
                if list(basis) == [1]:
                    certificate_order = moment
                    break
            assert certificate_order is not None
            chart_orders.append(certificate_order)
            total_charts += 1
        records.append(
            {
                "rotational_weight_support": list(support),
                "nonvanishing_charts": len(chart_orders),
                "maximum_moment_order": max(chart_orders),
                "chart_order_histogram": {
                    str(order): chart_orders.count(order)
                    for order in sorted(set(chart_orders))
                },
                "groebner_basis_on_every_chart": ["1"],
            }
        )

    assert len(records) == 27
    assert total_charts == 72
    assert max(record["maximum_moment_order"] for record in records) == 8
    return {
        "total_degree_bound": 3,
        "mixed_sign_three_weight_supports": len(records),
        "nonvanishing_coefficient_charts": total_charts,
        "moments_used_at_most": 8,
        "normalization": (
            "P -> lambda*P(r*Z,r^-1*W) sets one selected coefficient "
            "in each of the first two distinct weight components to one"
        ),
        "conclusion": (
            "no total-degree-at-most-three polynomial with exactly three "
            "nonzero rotational weights of both signs has all pure moments zero"
        ),
        "supports": records,
    }


def negate_weight_support(support: tuple[int, ...]) -> tuple[int, ...]:
    """Swap ``Z`` and ``W`` on a rotational-weight support."""

    return tuple(sorted(-weight for weight in support))


def cubic_four_weight_reduction() -> dict[str, object]:
    """Exclude every four-weight cubic outside three sign-symmetry orbits."""

    exceptional_representatives = {
        (-3, -1, 1, 3),
        (-2, -1, 0, 1),
        (-2, -1, 1, 2),
    }
    all_supports = [
        support
        for support in combinations(range(-3, 4), 4)
        if min(support) < 0 < max(support)
    ]
    representatives = sorted(
        {
            min(support, negate_weight_support(support))
            for support in all_supports
        }
    )
    records = []
    excluded_support_count = 0
    excluded_chart_count = 0
    for support in representatives:
        if support in exceptional_representatives:
            continue
        components = {
            weight: cubic_weight_component(weight) for weight in support
        }
        chart_orders = []
        for chart in product(
            *(range(len(components[weight])) for weight in support)
        ):
            selected = [
                components[weight][choice][0]
                for weight, choice in zip(support, chart)
            ]
            normalization = {selected[0]: 1, selected[1]: 1}
            polynomial = sp.expand(
                sum(
                    coefficient * Z**z_degree * W**w_degree
                    for weight in support
                    for coefficient, z_degree, w_degree in components[weight]
                ).subs(normalization)
            )
            nonzero_product = sp.prod(selected[2:])
            localizer = sp.Symbol("localizer")
            equations = [localizer * nonzero_product - 1]
            variables = sorted(
                polynomial.free_symbols - {Z, W}, key=str
            ) + [localizer]
            power = sp.Integer(1)
            certificate_order = None
            for moment in range(1, 7):
                power = sp.expand(power * polynomial)
                equation = expectation(power)
                if equation != 0:
                    equations.append(equation)
                if moment not in (3, 4, 5, 6):
                    continue
                basis = sp.groebner(
                    equations,
                    *variables,
                    order="grevlex",
                )
                if list(basis) == [1]:
                    certificate_order = moment
                    break
            assert certificate_order is not None
            chart_orders.append(certificate_order)

        reflected = negate_weight_support(support)
        orbit_size = 1 if reflected == support else 2
        excluded_support_count += orbit_size
        excluded_chart_count += orbit_size * len(chart_orders)
        records.append(
            {
                "sign_orbit_representative": list(support),
                "reflected_support": list(reflected),
                "orbit_size": orbit_size,
                "charts_checked_on_representative": len(chart_orders),
                "maximum_moment_order": max(chart_orders),
                "chart_order_histogram": {
                    str(order): chart_orders.count(order)
                    for order in sorted(set(chart_orders))
                },
                "groebner_basis_on_every_chart": ["1"],
            }
        )

    exceptional_supports = sorted(
        {
            support
            for support in all_supports
            if min(support, negate_weight_support(support))
            in exceptional_representatives
        }
    )
    exceptional_chart_count = sum(
        2
        ** sum(weight in (-1, 0, 1) for weight in support)
        for support in exceptional_supports
    )
    assert len(all_supports) == 33
    assert len(representatives) == 18
    assert len(records) == 15
    assert excluded_support_count == 29
    assert excluded_chart_count == 97
    assert exceptional_supports == [
        (-3, -1, 1, 3),
        (-2, -1, 0, 1),
        (-2, -1, 1, 2),
        (-1, 0, 1, 2),
    ]
    assert exceptional_chart_count == 24
    assert max(record["maximum_moment_order"] for record in records) == 6
    return {
        "total_degree_bound": 3,
        "mixed_sign_four_weight_supports": len(all_supports),
        "sign_symmetry_orbits": len(representatives),
        "excluded_sign_orbits": len(records),
        "excluded_supports": excluded_support_count,
        "excluded_nonvanishing_charts": excluded_chart_count,
        "moments_used_at_most": 6,
        "exceptional_sign_orbit_representatives": [
            list(support) for support in sorted(exceptional_representatives)
        ],
        "remaining_supports": [
            list(support) for support in exceptional_supports
        ],
        "remaining_nonvanishing_charts": exceptional_chart_count,
        "conclusion": (
            "a total-degree-three four-weight counterexample must lie on one "
            "of four explicit supports, comprising 24 nonvanishing charts"
        ),
        "excluded_orbits": records,
    }


def main() -> None:
    verify_two_weight_identity()
    verify_affine_source_identity()
    verify_quadratic_tilt_identity()
    cubic_three_weights = cubic_three_weight_exclusions()
    cubic_four_weights = cubic_four_weight_reduction()
    exclusions = [
        long_like_exclusion(2, 2),
        long_like_exclusion(3, 5),
        long_like_exclusion(4, 6),
    ]
    payload = {
        "format": "two-real-gmc-frontier-v2",
        "theorem_level_results": {
            "quadratic_GMC": (
                "GMC holds for every polynomial P of total degree at most two "
                "in every Gaussian dimension"
            ),
            "two_rotational_weights": (
                "in two real variables, a pure-moment-zero polynomial with at "
                "most two rotational weights reduces to one-sided weight support"
            ),
            "affine_circular_source": (
                "in two real variables, P=W*h(Z)+v(Z) with all pure moments "
                "zero is one-sided after the exact classification"
            ),
            "necessary_shape_of_GMC2_counterexample": (
                "degree at least three, at least three rotational weights, "
                "and degree at least two in each of Z and W"
            ),
            "cubic_three_weight_exclusion": (
                "a total-degree-three counterexample must have at least four "
                "nonzero rotational weights"
            ),
            "cubic_four_weight_reduction": (
                "among 33 mixed-sign four-weight cubic supports, 29 are "
                "excluded and only four explicit supports remain"
            ),
        },
        "exact_cubic_three_weight_exclusion": cubic_three_weights,
        "exact_cubic_four_weight_reduction": cubic_four_weights,
        "exact_finite_support_exclusions": exclusions,
        "scope_warning": (
            "the Groebner exclusions are bounded ansatz results; GMC(2) remains open"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS two-real GMC: two-weight contraction identity")
    print("PASS two-real GMC: affine-source Lagrange identity through order 7")
    print("PASS quadratic GMC: nilpotent quadratic tilt regression")
    print(
        "PASS cubic GMC(2):",
        cubic_three_weights["mixed_sign_three_weight_supports"],
        "supports and",
        cubic_three_weights["nonvanishing_coefficient_charts"],
        "charts excluded through moment order",
        cubic_three_weights["moments_used_at_most"],
    )
    print(
        "PASS cubic four-weight GMC(2):",
        cubic_four_weights["excluded_supports"],
        "of",
        cubic_four_weights["mixed_sign_four_weight_supports"],
        "supports excluded; remaining charts",
        cubic_four_weights["remaining_nonvanishing_charts"],
    )
    for exclusion in exclusions:
        print(
            "PASS finite exclusion:",
            exclusion["coefficient_count"],
            "coefficients,",
            exclusion["moments_used"],
            "moments, Groebner basis [1]",
        )
    print(f"PASS two-real GMC: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
