#!/usr/bin/env python3
"""Exact low-degree stable-target stabilizer audit for the base quintic."""

from __future__ import annotations

from math import gcd
import subprocess
import sys

from flint import fmpz_mat
import sympy as sp


S, P, B, C = sp.symbols("S P B C")
TARGET_VARIABLES = (P, B, C)
LAMBDA, MU = sp.symbols("lambda mu")

SEED_4 = sp.Rational(-61712472, 10440125)
SEED_5 = sp.Rational(5636405776, 4437053125)
SELECTED_TARGET = {
    P: sp.Rational(85, 274),
    B: sp.Rational(225, 137),
    C: sp.Rational(120, 137),
}


inverse_relation = (
    SEED_5 * P**5 * S**5
    + SEED_4 * P**4 * S**4
    + P * S**3
    - B * S**2 / 2
    + S
    - C / 2
)
discriminant = sp.factor(sp.discriminant(inverse_relation, S))
discriminant_polynomial = sp.Poly(
    discriminant,
    *TARGET_VARIABLES,
    domain=sp.QQ,
)
minimum_p_order = min(
    monomial[0]
    for monomial, _coefficient in discriminant_polynomial.terms()
)
assert minimum_p_order == 8

H = sp.factor(discriminant / P**minimum_p_order)
H_polynomial = sp.Poly(H, *TARGET_VARIABLES, domain=sp.QQ)
assert H_polynomial.is_irreducible
assert H_polynomial.total_degree() == 16
assert len(H_polynomial.terms()) == 59
top_terms = [
    (monomial, coefficient)
    for monomial, coefficient in H_polynomial.terms()
    if sum(monomial) == H_polynomial.total_degree()
]
assert len(top_terms) == 1
assert top_terms[0][0] == (12, 0, 4)
assert H.subs(SELECTED_TARGET) != 0

# The positive upper Newton hull is an all-degree statement.  Write
#
#   A=(12,0,4), D=(2,5,1).
#
# Every support exponent u is coordinatewise dominated by a point of the
# segment [D,A].  Hence a strictly positive weight exposes A, D, or their
# common edge, and no other face.  The interval below is the exact set of
# segment parameters t for which D+t(A-D) dominates u.
NEWTON_VERTEX_A = (12, 0, 4)
NEWTON_VERTEX_D = (2, 5, 1)
for monomial, _coefficient in H_polynomial.terms():
    first, second, third = monomial
    lower_bound = max(
        sp.Rational(0),
        sp.Rational(first - 2, 10),
        sp.Rational(third - 1, 3),
    )
    upper_bound = min(
        sp.Rational(1),
        sp.Rational(5 - second, 5),
    )
    assert lower_bound <= upper_bound
    if monomial not in (NEWTON_VERTEX_A, NEWTON_VERTEX_D):
        segment_parameter = sp.Rational(5 - second, 5)
        lies_on_segment = (
            sp.Rational(first) == 2 + 10 * segment_parameter
            and sp.Rational(third) == 1 + 3 * segment_parameter
        )
        assert not lies_on_segment

coefficient_A = H_polynomial.coeff_monomial(NEWTON_VERTEX_A)
coefficient_D = H_polynomial.coeff_monomial(NEWTON_VERTEX_D)
assert coefficient_A != 0
assert coefficient_D != 0

# D is an intruder in the sense of Derksen--Hadas--Makar-Limanov: every
# exponent is positive, and the positive weight (1,3,1) exposes it uniquely.
# Their coordinate-polynomial theorem therefore turns T_i=x_i+H*V_i into an
# all-degree unstabilized obstruction whenever V_i is nonzero.
intruder_weight = (1, 3, 1)
intruder_degree = sum(
    exponent * weight
    for exponent, weight in zip(
        NEWTON_VERTEX_D,
        intruder_weight,
        strict=True,
    )
)
assert all(exponent > 0 for exponent in NEWTON_VERTEX_D)
for monomial, _coefficient in H_polynomial.terms():
    if monomial == NEWTON_VERTEX_D:
        continue
    assert sum(
        exponent * weight
        for exponent, weight in zip(
            monomial,
            intruder_weight,
            strict=True,
        )
    ) < intruder_degree

# The two vertices really can tie, so the preceding Newton reduction is not
# a monomial-avoidance theorem.  An explicit logarithmic field V=H*W gives a
# cancellable leading wall.  Take
#
#   W=(1, k*P^39*C^14, P*C^2),  k^5=-coefficient_A^17/coefficient_D.
#
# Then T=x+H*V=x+H^2*W has component degrees (32,85,35).  The A and D terms
# have the same exponent after substitution and their coefficients cancel.
wall_weights = (32, 85, 35)
assert sum(
    exponent * weight
    for exponent, weight in zip(
        NEWTON_VERTEX_A,
        wall_weights,
        strict=True,
    )
) == sum(
    exponent * weight
    for exponent, weight in zip(
        NEWTON_VERTEX_D,
        wall_weights,
        strict=True,
    )
)
leading_P = (24, 0, 8)
leading_B = (63, 0, 22)
leading_C = (25, 0, 10)
substituted_A = tuple(
    12 * leading_P[index] + 4 * leading_C[index]
    for index in range(3)
)
substituted_D = tuple(
    2 * leading_P[index]
    + 5 * leading_B[index]
    + leading_C[index]
    for index in range(3)
)
assert substituted_A == substituted_D == (388, 0, 136)
k_fifth_power = -(coefficient_A**17) / coefficient_D
assert (
    coefficient_A**33
    + coefficient_D * k_fifth_power * coefficient_A**16
) == 0

# A smaller cancellable wall already occurs in the Koszul submodule.  The
# P-zero two-generator ladder is
#
#   V_B=L*H_C,
#   V_C=-L*H_B-k*G*H,
#   Q=-k*G*H_C.
#
# Since deg(H_B)=13 and deg(H_C)=15, put deg(L)=1+3*n and
# deg(G)=18+5*n.  Its component weights are
#
#   (w_P,w_B,w_C)=(1,32+3*n,50+5*n),
#
# and every rung lies on the A-D wall.
gradient_total_degrees = tuple(
    sp.Poly(
        sp.diff(H, variable),
        *TARGET_VARIABLES,
        domain=sp.QQ,
    ).total_degree()
    for variable in TARGET_VARIABLES
)
assert gradient_total_degrees == (15, 13, 15)
for ladder_index in range(5):
    ladder_weights = (
        1,
        32 + 3 * ladder_index,
        50 + 5 * ladder_index,
    )
    assert (
        10 * ladder_weights[0]
        - 5 * ladder_weights[1]
        + 3 * ladder_weights[2]
    ) == 0

# At n=0, leading cancellation would force R^3=P^7*C*L with L linear.
# If L is P or C the displayed exponent pairs are not divisible by three;
# every other linear L contributes a new irreducible factor to exponent one.
assert any(exponent % 3 for exponent in (8, 1))
assert any(exponent % 3 for exponent in (7, 2))

# At n=1 the obstruction disappears:
#
#   L=P^2*C^2, R=P^3*C, G=P*C^2*R^5=P^16*C^7.
#
# The leading B- and C-components of T=x+H*V are respectively
# 4*h_A^2*P^26*C^9 and -k*h_A^2*P^40*C^15.  The two Newton terms become the
# same monomial P^172*C^60 and cancel when
# k^3=4^5*h_D*h_A^3.
koszul_L = (2, 0, 2)
koszul_R = (3, 0, 1)
koszul_G = tuple(
    unit + 5 * exponent
    for unit, exponent in zip((1, 0, 2), koszul_R, strict=True)
)
assert koszul_G == (16, 0, 7)
koszul_leading_B = tuple(
    NEWTON_VERTEX_A[index]
    + koszul_L[index]
    + (12, 0, 3)[index]
    for index in range(3)
)
koszul_leading_C = tuple(
    2 * NEWTON_VERTEX_A[index] + koszul_G[index]
    for index in range(3)
)
assert koszul_leading_B == (26, 0, 9)
assert koszul_leading_C == (40, 0, 15)
koszul_substituted_A = tuple(
    12 * (1, 0, 0)[index] + 4 * koszul_leading_C[index]
    for index in range(3)
)
koszul_substituted_D = tuple(
    2 * (1, 0, 0)[index]
    + 5 * koszul_leading_B[index]
    + koszul_leading_C[index]
    for index in range(3)
)
assert koszul_substituted_A == koszul_substituted_D == (172, 0, 60)
koszul_k_cubed = 4**5 * coefficient_D * coefficient_A**3
assert (
    koszul_k_cubed * coefficient_A**9
    - 4**5 * coefficient_D * coefficient_A**12
) == 0


def sparse_top_product(
    left: dict[tuple[int, int, int], sp.Expr],
    left_maximum_degree: int,
    right: dict[tuple[int, int, int], sp.Expr],
    right_maximum_degree: int,
    depth: int,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Multiply sparse polynomials, retaining the top ``depth`` layers."""

    maximum_degree = left_maximum_degree + right_maximum_degree
    result: dict[tuple[int, int, int], sp.Expr] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_entry + right_entry
                for left_entry, right_entry in zip(
                    left_monomial,
                    right_monomial,
                    strict=True,
                )
            )
            if maximum_degree - sum(monomial) > depth:
                continue
            result[monomial] = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
    return {
        monomial: sp.expand(coefficient)
        for monomial, coefficient in result.items()
        if coefficient != 0
    }


def sparse_top_power(
    polynomial: dict[tuple[int, int, int], sp.Expr],
    maximum_degree: int,
    exponent: int,
    depth: int,
) -> tuple[dict[tuple[int, int, int], sp.Expr], int]:
    """Top layers of a small nonnegative power."""

    result = {(0, 0, 0): sp.Integer(1)}
    result_maximum_degree = 0
    for _index in range(exponent):
        result = sparse_top_product(
            result,
            result_maximum_degree,
            polynomial,
            maximum_degree,
            depth,
        )
        result_maximum_degree += maximum_degree
    return result, result_maximum_degree


def research_koszul_wall_layers(depth: int = 6) -> None:
    """Expand exact top layers of the target-degree-55 Koszul wall."""

    eta = sp.symbols("eta")
    lower_l, lower_g = sp.symbols("lower_l lower_g")
    h_sparse = {
        monomial: coefficient
        for monomial, coefficient in H_polynomial.terms()
    }
    hc_sparse = {
        (first, second, third - 1): third * coefficient
        for (first, second, third), coefficient in H_polynomial.terms()
        if third
    }

    # T_B=B+H*P^2*C^2*H_C.  The identity B is thirty-four layers below
    # the leading correction and is irrelevant at the requested depth.
    h_times_hc = sparse_top_product(
        h_sparse,
        16,
        hc_sparse,
        15,
        depth,
    )
    leading_b = {
        tuple(
            entry + shift
            for entry, shift in zip(monomial, koszul_L, strict=True)
        ): coefficient
        for monomial, coefficient in h_times_hc.items()
        if 35 - sum(
            entry + shift
            for entry, shift in zip(monomial, koszul_L, strict=True)
        )
        <= depth
    }
    lower_l_monomial = (2, 0, 1)
    for monomial, coefficient in h_times_hc.items():
        shifted = tuple(
            entry + shift
            for entry, shift in zip(
                monomial,
                lower_l_monomial,
                strict=True,
            )
        )
        if 35 - sum(shifted) <= depth:
            leading_b[shifted] = (
                leading_b.get(shifted, 0)
                + lower_l * coefficient
            )

    # The high part of T_C is -eta*P^16*C^7*H^2.  The other Koszul term
    # H*L*H_B begins in degree 33, twenty-two layers lower.
    h_squared = sparse_top_product(
        h_sparse,
        16,
        h_sparse,
        16,
        depth,
    )
    leading_c = {
        tuple(
            entry + shift
            for entry, shift in zip(monomial, koszul_G, strict=True)
        ): -eta * coefficient
        for monomial, coefficient in h_squared.items()
        if 55 - sum(
            entry + shift
            for entry, shift in zip(monomial, koszul_G, strict=True)
        )
        <= depth
    }
    lower_g_monomial = (16, 0, 6)
    for monomial, coefficient in h_squared.items():
        shifted = tuple(
            entry + shift
            for entry, shift in zip(
                monomial,
                lower_g_monomial,
                strict=True,
            )
        )
        if 55 - sum(shifted) <= depth:
            leading_c[shifted] = (
                leading_c.get(shifted, 0)
                - eta * lower_g * coefficient
            )

    maximum_weight = 232
    composed: dict[tuple[int, int, int], sp.Expr] = {}
    for outer_monomial, outer_coefficient in H_polynomial.terms():
        first, second, third = outer_monomial
        outer_weight = first + 35 * second + 55 * third
        outer_gap = maximum_weight - outer_weight
        if outer_gap < 0 or outer_gap > depth:
            continue
        local_depth = depth - outer_gap
        b_power, b_maximum = sparse_top_power(
            leading_b,
            35,
            second,
            local_depth,
        )
        c_power, c_maximum = sparse_top_power(
            leading_c,
            55,
            third,
            local_depth,
        )
        product = sparse_top_product(
            b_power,
            b_maximum,
            c_power,
            c_maximum,
            local_depth,
        )
        for monomial, coefficient in product.items():
            shifted = (
                monomial[0] + first,
                monomial[1],
                monomial[2],
            )
            if maximum_weight - sum(shifted) > depth:
                continue
            composed[shifted] = (
                composed.get(shifted, 0)
                + outer_coefficient * coefficient
            )

    eta_relation = sp.Poly(
        eta**3 - koszul_k_cubed,
        eta,
        domain=sp.EX,
    )
    reduced = {}
    for monomial, coefficient in composed.items():
        remainder = sp.Poly(
            sp.expand(coefficient),
            eta,
            domain=sp.EX,
        ).rem(eta_relation).as_expr()
        if remainder != 0:
            reduced[monomial] = remainder

    by_degree: dict[int, dict[tuple[int, int, int], sp.Expr]] = {}
    for monomial, coefficient in reduced.items():
        by_degree.setdefault(sum(monomial), {})[monomial] = coefficient
    assert maximum_weight not in by_degree

    print("R5KOSZUL_WALL_TOP_DEGREE=232")
    print("R5KOSZUL_WALL_TOP_CANCELS=true")
    first_layer = by_degree.get(maximum_weight - 1, {})
    assert len(first_layer) == 1
    first_equation = next(iter(first_layer.values()))
    lower_l_solutions = sp.solve(
        sp.together(first_equation),
        lower_l,
        dict=True,
    )
    assert len(lower_l_solutions) == 1
    lower_l_solution = lower_l_solutions[0][lower_l]
    print("R5KOSZUL_WALL_DEGREE_231_TERMS=1")
    print("R5KOSZUL_WALL_DEGREE_231_CORRECTABLE=true")
    print(
        "R5KOSZUL_WALL_DEGREE_231_FREE_PARAMETERS="
        f"{int(lower_l_solution.has(lower_g))}"
    )

    corrected_by_degree = {
        degree: {
            monomial: sp.factor(
                coefficient.subs(lower_l, lower_l_solution)
            )
            for monomial, coefficient in layer.items()
            if coefficient.subs(lower_l, lower_l_solution) != 0
        }
        for degree, layer in by_degree.items()
    }
    for degree in range(maximum_weight - 2, maximum_weight - depth - 1, -1):
        layer = corrected_by_degree.get(degree, {})
        layer = {
            monomial: coefficient
            for monomial, coefficient in layer.items()
            if coefficient != 0
        }
        if layer:
            print(
                "R5KOSZUL_WALL_FIRST_SURVIVOR_AFTER_DEGREE_231="
                f"{degree}:{len(layer)}"
            )
            print(
                "R5KOSZUL_WALL_FIRST_SURVIVOR_MONOMIALS="
                f"{sorted(layer)}"
            )
            parameter_derivatives = sp.Matrix(
                [
                    [
                        sp.diff(coefficient, lower_g)
                        for coefficient in layer.values()
                    ]
                ]
            )
            print(
                "R5KOSZUL_WALL_FIRST_SURVIVOR_PARAMETER_RANK="
                f"{parameter_derivatives.rank()}"
            )
            break


if "--research-koszul-wall-layers" in sys.argv:
    research_koszul_wall_layers()
    raise SystemExit(0)


def research_koszul_hensel(maximum_depth: int) -> None:
    """Lift the degree-55 P-zero Koszul wall one homogeneous layer at a time."""

    assert 1 <= maximum_depth < 22
    eta = sp.symbols("eta")
    lower_l_monomials = [
        (first, 0, third)
        for degree in range(4)
        for first in range(degree + 1)
        for third in [degree - first]
    ]
    lower_l_symbols = sp.symbols(
        " ".join(
            f"l_{first}_{third}"
            for first, _second, third in lower_l_monomials
        )
    )
    lower_l_parameters = dict(
        zip(lower_l_monomials, lower_l_symbols, strict=True)
    )
    forced_zero_arguments = [
        argument.removeprefix("--research-zero-l=")
        for argument in sys.argv
        if argument.startswith("--research-zero-l=")
    ]
    assert len(forced_zero_arguments) <= 1
    forced_zero_l_names = (
        set(forced_zero_arguments[0].split(","))
        if forced_zero_arguments
        else set()
    )
    assert forced_zero_l_names <= {
        str(parameter)
        for parameter in lower_l_symbols
    }
    continue_constraints = "--research-continue-constraints" in sys.argv
    active_lower_l_symbols = tuple(
        parameter
        for parameter in lower_l_symbols
        if str(parameter) not in forced_zero_l_names
    )
    h_sparse = {
        monomial: coefficient
        for monomial, coefficient in H_polynomial.terms()
    }
    hc_sparse = {
        (first, second, third - 1): third * coefficient
        for (first, second, third), coefficient in H_polynomial.terms()
        if third
    }
    h_times_hc = sparse_top_product(
        h_sparse,
        16,
        hc_sparse,
        15,
        maximum_depth,
    )
    target_b = {
        tuple(
            entry + shift
            for entry, shift in zip(monomial, koszul_L, strict=True)
        ): coefficient
        for monomial, coefficient in h_times_hc.items()
        if 35
        - sum(
            entry + shift
            for entry, shift in zip(monomial, koszul_L, strict=True)
        )
        <= maximum_depth
    }
    # Keep every lower B-free homogeneous layer of L as a parameter.
    for l_monomial, l_parameter in lower_l_parameters.items():
        if str(l_parameter) in forced_zero_l_names:
            continue
        for monomial, coefficient in h_times_hc.items():
            shifted = tuple(
                entry + shift
                for entry, shift in zip(
                    monomial,
                    l_monomial,
                    strict=True,
                )
            )
            if 35 - sum(shifted) <= maximum_depth:
                target_b[shifted] = (
                    target_b.get(shifted, 0)
                    + l_parameter * coefficient
                )
    h_squared = sparse_top_product(
        h_sparse,
        16,
        h_sparse,
        16,
        maximum_depth,
    )
    eta_relation = sp.Poly(
        eta**3 - koszul_k_cubed,
        eta,
        domain=sp.EX,
    )
    constraint_equations: list[sp.Expr] = []
    constraint_basis = None

    def reduce_eta(expression: sp.Expr) -> sp.Expr:
        eta_remainder = sp.Poly(
            sp.expand(expression),
            eta,
            domain=sp.EX,
        ).rem(eta_relation).as_expr()
        if constraint_basis is None:
            return eta_remainder

        reduced_expression = 0
        eta_polynomial = sp.Poly(
            eta_remainder,
            eta,
            domain=sp.EX,
        )
        for (eta_power,), coefficient in eta_polynomial.terms():
            numerator, denominator = sp.together(
                coefficient
            ).as_numer_denom()
            assert not any(
                denominator.has(parameter)
                for parameter in active_lower_l_symbols
            )
            reduced_numerator = constraint_basis.reduce(
                sp.Poly(
                    numerator,
                    *active_lower_l_symbols,
                    domain=sp.QQ,
                ).as_expr()
            )[1]
            reduced_expression += (
                reduced_numerator * eta**eta_power / denominator
            )
        return sp.expand(reduced_expression)

    def compose(
        g_polynomial: dict[tuple[int, int, int], sp.Expr],
        depth: int,
    ) -> dict[tuple[int, int, int], sp.Expr]:
        target_c: dict[tuple[int, int, int], sp.Expr] = {}
        for h_monomial, h_coefficient in h_squared.items():
            for g_monomial, g_coefficient in g_polynomial.items():
                shifted = tuple(
                    entry + shift
                    for entry, shift in zip(
                        h_monomial,
                        g_monomial,
                        strict=True,
                    )
                )
                if 55 - sum(shifted) > depth:
                    continue
                target_c[shifted] = (
                    target_c.get(shifted, 0)
                    - eta * h_coefficient * g_coefficient
                )

        maximum_weight = 232
        result: dict[tuple[int, int, int], sp.Expr] = {}
        for outer_monomial, outer_coefficient in H_polynomial.terms():
            first, second, third = outer_monomial
            outer_weight = first + 35 * second + 55 * third
            outer_gap = maximum_weight - outer_weight
            if outer_gap < 0 or outer_gap > depth:
                continue
            local_depth = depth - outer_gap
            b_power, b_maximum = sparse_top_power(
                target_b,
                35,
                second,
                local_depth,
            )
            c_power, c_maximum = sparse_top_power(
                target_c,
                55,
                third,
                local_depth,
            )
            product = sparse_top_product(
                b_power,
                b_maximum,
                c_power,
                c_maximum,
                local_depth,
            )
            for monomial, coefficient in product.items():
                shifted = (
                    monomial[0] + first,
                    monomial[1],
                    monomial[2],
                )
                if maximum_weight - sum(shifted) > depth:
                    continue
                result[shifted] = (
                    result.get(shifted, 0)
                    + outer_coefficient * coefficient
                )

        return {
            monomial: remainder
            for monomial, coefficient in result.items()
            if (
                remainder := reduce_eta(coefficient)
            )
            != 0
        }

    # G_23=P^16*C^7.  At depth r, the new homogeneous piece G_(23-r)
    # enters the boundary equation by multiplication with
    # 3*4^5*h_D*h_A^12*eta*P^156*C^53.
    g_polynomial: dict[tuple[int, int, int], sp.Expr] = {
        koszul_G: sp.Integer(1)
    }
    response_scalar_without_eta = (
        3 * 4**5 * coefficient_D * coefficient_A**12
    )
    response_monomial = (156, 0, 53)

    print(f"R5KOSZUL_HENSEL_REQUESTED_DEPTH={maximum_depth}")
    for depth in range(1, maximum_depth + 1):
        composed = compose(g_polynomial, depth)
        layer_degree = 232 - depth
        layer = {
            monomial: coefficient
            for monomial, coefficient in composed.items()
            if sum(monomial) == layer_degree
        }
        correction: dict[tuple[int, int, int], sp.Expr] = {}
        divisible = True
        obstruction_monomial = None
        obstruction_equations: list[sp.Expr] = []
        for monomial, coefficient in layer.items():
            quotient_monomial = tuple(
                exponent - base
                for exponent, base in zip(
                    monomial,
                    response_monomial,
                    strict=True,
                )
            )
            if min(quotient_monomial) < 0:
                divisible = False
                if obstruction_monomial is None:
                    obstruction_monomial = monomial
                obstruction_equations.append(coefficient)
                continue
            if sum(quotient_monomial) != 23 - depth:
                divisible = False
                if obstruction_monomial is None:
                    obstruction_monomial = monomial
                obstruction_equations.append(coefficient)
                continue
            # eta^(-1)=eta^2/koszul_k_cubed.
            quotient_coefficient = reduce_eta(
                -coefficient
                * eta**2
                / (
                    response_scalar_without_eta
                    * koszul_k_cubed
                )
            )
            correction[quotient_monomial] = (
                correction.get(quotient_monomial, 0)
                + quotient_coefficient
            )

        print(
            f"R5KOSZUL_HENSEL_DEPTH_{depth}="
            f"terms:{len(layer)},divisible:{str(divisible).lower()}"
        )
        if not divisible:
            print(
                "R5KOSZUL_HENSEL_FIRST_OBSTRUCTION_DEPTH="
                f"{depth}"
            )
            print(
                "R5KOSZUL_HENSEL_FIRST_OBSTRUCTION_MONOMIAL="
                f"{obstruction_monomial}"
            )
            scalar_equations = []
            for equation in obstruction_equations:
                eta_polynomial = sp.Poly(
                    reduce_eta(equation),
                    eta,
                    domain=sp.EX,
                )
                scalar_equations.extend(
                    sp.together(coefficient)
                    for coefficient in eta_polynomial.all_coeffs()
                    if coefficient != 0
                )
            active_l_parameters = sorted(
                {
                    parameter
                    for equation in scalar_equations
                    for parameter in lower_l_symbols
                    if str(parameter) not in forced_zero_l_names
                    if equation.has(parameter)
                },
                key=str,
            )
            print(
                "R5KOSZUL_HENSEL_OBSTRUCTION_EQUATIONS="
                f"{len(scalar_equations)}"
            )
            print(
                "R5KOSZUL_HENSEL_ACTIVE_LOWER_L_PARAMETERS="
                f"{len(active_l_parameters)}"
            )
            print(
                "R5KOSZUL_HENSEL_ACTIVE_LOWER_L_PARAMETER_NAMES="
                f"{active_l_parameters}"
            )
            if len(active_l_parameters) <= 2:
                normalized_obstructions = [
                    sp.Poly(
                        equation.as_numer_denom()[0],
                        *active_l_parameters,
                        domain=sp.QQ,
                    ).primitive()[1].as_expr()
                    for equation in scalar_equations
                ]
                print(
                    "R5KOSZUL_HENSEL_NORMALIZED_OBSTRUCTIONS="
                    f"{[sp.factor(equation) for equation in normalized_obstructions]}"
                )
            if len(active_l_parameters) == 1:
                active_parameter = active_l_parameters[0]
                active_polynomials = [
                    sp.Poly(
                        equation.as_numer_denom()[0],
                        active_parameter,
                        domain=sp.QQ.frac_field(
                            *[
                                parameter
                                for parameter in lower_l_symbols
                                if parameter != active_parameter
                            ]
                        ),
                    )
                    for equation in scalar_equations
                ]
                common_polynomial = active_polynomials[0]
                for polynomial in active_polynomials[1:]:
                    common_polynomial = sp.gcd(
                        common_polynomial,
                        polynomial,
                    )
                print(
                    "R5KOSZUL_HENSEL_SINGLE_PARAMETER_DEGREE="
                    f"{common_polynomial.degree()}"
                )
                print(
                    "R5KOSZUL_HENSEL_SINGLE_PARAMETER_POLYNOMIAL="
                    f"{sp.factor(common_polynomial.monic().as_expr())}"
                )
            if not continue_constraints:
                return
            assert active_lower_l_symbols
            constraint_equations.extend(
                sp.Poly(
                    equation.as_numer_denom()[0],
                    *active_lower_l_symbols,
                    domain=sp.QQ,
                ).primitive()[1].as_expr()
                for equation in scalar_equations
            )
            constraint_basis = sp.groebner(
                constraint_equations,
                *active_lower_l_symbols,
                order="grevlex",
                domain=sp.QQ,
            )
            print(
                "R5KOSZUL_HENSEL_CONSTRAINT_BASIS_SIZE="
                f"{len(constraint_basis.polys)}"
            )
            if constraint_basis.reduce(sp.Integer(1))[1] == 0:
                print(
                    "R5KOSZUL_HENSEL_CONSTRAINT_IDEAL=unit"
                )
                return
            print(
                "R5KOSZUL_HENSEL_CONSTRAINT_ZERO_DIMENSIONAL="
                f"{str(constraint_basis.is_zero_dimensional).lower()}"
            )
        for monomial, coefficient in correction.items():
            g_polynomial[monomial] = (
                g_polynomial.get(monomial, 0)
                + coefficient
            )

        verification = compose(g_polynomial, depth)
        assert not any(
            sum(monomial) >= layer_degree
            for monomial in verification
        )

    print(
        "R5KOSZUL_HENSEL_NO_OBSTRUCTION_THROUGH_DEPTH="
        f"{maximum_depth}"
    )
    print(f"R5KOSZUL_HENSEL_G_TERMS={len(g_polynomial)}")


if "--research-koszul-hensel" in sys.argv:
    depth_arguments = [
        argument.removeprefix("--research-depth=")
        for argument in sys.argv
        if argument.startswith("--research-depth=")
    ]
    assert len(depth_arguments) <= 1
    requested_depth = int(depth_arguments[0]) if depth_arguments else 8
    research_koszul_hensel(requested_depth)
    raise SystemExit(0)

# A fixed rank-five coefficient point has alpha^5=alpha^6=1, hence the
# displayed coefficient-torus stabilizer is trivial.
assert gcd(5, 6) == 1


_denominator, primitive_H = H_polynomial.clear_denoms()
primitive_H = sp.Poly(
    primitive_H,
    *TARGET_VARIABLES,
    domain=sp.ZZ,
)
primitive_H_expression = primitive_H.as_expr()
primitive_gradient = [
    sp.Poly(
        sp.diff(primitive_H_expression, variable),
        *TARGET_VARIABLES,
        domain=sp.ZZ,
    )
    for variable in TARGET_VARIABLES
]

# The logarithmic multiplier is one in every degree.  An exact algebraic
# triple-root point supplies a singular point of H=0, so evaluating
# V(H)-QH-kappa there forces kappa=0.
triple_root = sp.symbols("triple_root")
triple_root_relation = (
    15 * SEED_5 * triple_root**4
    + 8 * SEED_4 * triple_root**3
    + 3 * triple_root**2
    - 1
)
triple_root_B = (
    20 * SEED_5 * triple_root**3
    + 12 * SEED_4 * triple_root**2
    + 6 * triple_root
)
triple_root_C = 2 * (
    SEED_5 * triple_root**5
    + SEED_4 * triple_root**4
    + triple_root**3
    - triple_root_B * triple_root**2 / 2
    + triple_root
)
triple_root_substitution = {
    P: 1,
    B: triple_root_B,
    C: triple_root_C,
}
for discriminant_equation in [
    primitive_H_expression,
    *[derivative.as_expr() for derivative in primitive_gradient],
]:
    numerator = sp.together(
        discriminant_equation.subs(triple_root_substitution)
    ).as_numer_denom()[0]
    assert sp.rem(
        sp.Poly(numerator, triple_root, domain=sp.QQ),
        sp.Poly(triple_root_relation, triple_root, domain=sp.QQ),
    ).is_zero


def macaulay2_homogenized_h() -> str:
    """Return the degree-sixteen homogenization in Macaulay2 syntax."""

    Z = sp.symbols("Z")
    homogenized_H = sum(
        int(coefficient)
        * P**first
        * B**second
        * C**third
        * Z ** (16 - first - second - third)
        for (first, second, third), coefficient in primitive_H.terms()
    )
    return str(sp.expand(homogenized_H)).replace("**", "^")


def research_ring() -> str:
    """Coefficient ring selected for optional long research calculations."""

    prefix = "--research-characteristic="
    values = [
        argument.removeprefix(prefix)
        for argument in sys.argv
        if argument.startswith(prefix)
    ]
    assert len(values) <= 1
    if not values:
        return "QQ"
    characteristic = int(values[0])
    assert characteristic > 0
    return f"ZZ/{characteristic}"


def run_macaulay2_research(task: str) -> None:
    """Run an optional long module/singular-locus calculation."""

    coefficient_ring = research_ring()
    macaulay2_H = macaulay2_homogenized_h()
    common_program = f"""
R={coefficient_ring}[P,B,C,Z];
H={macaulay2_H};
HP=diff(P,H); HB=diff(B,H); HC=diff(C,H);
"""
    if task == "koszul-homology":
        task_program = """
F=R^{3:-15,1:-16};
phi=map(R^1,F,matrix{{HP,HB,HC,-H}});
Q=prune HH_1(koszul phi);
print("R5WALL_TASK=koszul-homology");
print("R5WALL_RING=" | toString coefficientRing R);
print("R5WALL_DIM=" | toString dim Q);
print("R5WALL_DEGREE=" | toString degree Q);
print("R5WALL_BETTI=" | toString betti res Q);
print("R5WALL_HILBERT=" | toString hilbertSeries Q);
"""
    elif task == "singular-primes":
        task_program = """
J=ideal(H,HP,HB,HC);
primeList=minimalPrimes J;
print("R5WALL_TASK=singular-primes");
print("R5WALL_RING=" | toString coefficientRing R);
print("R5WALL_SINGULAR_DIM=" | toString dim J);
print("R5WALL_SINGULAR_DEGREE=" | toString degree J);
print("R5WALL_PRIME_COUNT=" | toString(#primeList));
scan(0..((#primeList)-1), index -> (
    componentIdeal=primeList#index;
    print(
        "R5WALL_PRIME_" | toString index
        | "_DIM=" | toString dim componentIdeal
        | "_DEGREE=" | toString degree componentIdeal
        | "_GENERATOR_DEGREES=" | toString degrees source gens componentIdeal
    );
));
"""
    elif task == "singular-boundary":
        task_program = """
J=ideal(H,HP,HB,HC);
affinePZero=saturate(J+ideal(P),ideal(Z));
infinityIdeal=trim(J+ideal(Z));
infinityRadical=radical infinityIdeal;
expectedInfinityRadical=ideal(Z,P*C);
affinePZeroIsEmpty=isSubset(ideal(1_R),affinePZero);
infinityRadicalIsExpected=(
    isSubset(infinityRadical,expectedInfinityRadical)
    and isSubset(expectedInfinityRadical,infinityRadical)
);
assert affinePZeroIsEmpty;
assert infinityRadicalIsExpected;
print("R5WALL_TASK=singular-boundary");
print("R5WALL_RING=" | toString coefficientRing R);
print("R5WALL_SINGULAR_DIM=" | toString dim J);
print("R5WALL_SINGULAR_DEGREE=" | toString degree J);
print(
    "R5WALL_P_ZERO_AFFINE_EMPTY="
    | toString affinePZeroIsEmpty
);
print(
    "R5WALL_INFINITY_RADICAL_EXPECTED="
    | toString infinityRadicalIsExpected
);
print("R5WALL_INFINITY_COMPONENTS=(Z,P),(Z,C)");
"""
    else:
        raise AssertionError(task)

    capture_result = task == "singular-boundary"
    result = subprocess.run(
        ["M2", "--no-readline", "--silent"],
        input=common_program + task_program + "\nexit 0;\n",
        text=True,
        capture_output=capture_result,
        check=False,
    )
    if capture_result:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        assert "error:" not in result.stdout + result.stderr
        assert "R5WALL_P_ZERO_AFFINE_EMPTY=true" in result.stdout
        assert "R5WALL_INFINITY_RADICAL_EXPECTED=true" in result.stdout
    raise SystemExit(result.returncode)


def run_root_partition_research(task: str) -> None:
    """Compute a targeted affine singular-locus component by elimination."""

    coefficient_ring = research_ring()
    seed_4 = str(SEED_4)
    seed_5 = str(SEED_5)
    macaulay2_affine_h = str(primitive_H_expression).replace("**", "^")
    if task == "triple-root":
        program = f"""
R={coefficient_ring}[s,P,B,C,MonomialOrder=>Eliminate 1];
u4={seed_4}; u5={seed_5};
H={macaulay2_affine_h};
tripleEquation=15*u5*P^5*s^4+8*u4*P^4*s^3+3*P*s^2-1;
rootB=20*u5*P^5*s^3+12*u4*P^4*s^2+6*P*s;
rootC=2*(u5*P^5*s^5+u4*P^4*s^4+P*s^3-rootB*s^2/2+s);
parameterIdeal=ideal(tripleEquation,B-rootB,C-rootC);
targetIdeal=trim eliminate({{s}},parameterIdeal);
singularIdeal=ideal(H,diff(P,H),diff(B,H),diff(C,H));
parameterIsPrime=isPrime parameterIdeal;
assert parameterIsPrime;
print("R5ROOT_TASK=triple-root");
print("R5ROOT_RING=" | toString coefficientRing R);
print("R5ROOT_TARGET_DIM=" | toString(dim targetIdeal-1));
print("R5ROOT_DEGREE=" | toString degree targetIdeal);
print(
    "R5ROOT_PARAMETER_IS_PRIME="
    | toString parameterIsPrime
);
print("R5ROOT_IS_PRIME_BY_CONTRACTION=true");
print(
    "R5ROOT_CONTAINS_SINGULAR_IDEAL="
    | toString isSubset(singularIdeal,targetIdeal)
);
print("R5ROOT_GENERATORS=" | toString flatten entries gens targetIdeal);
"""
    elif task == "two-double-root":
        program = f"""
R={coefficient_ring}[x,y,u,P,B,C,MonomialOrder=>Eliminate 3];
u4={seed_4}; u5={seed_5};
H={macaulay2_affine_h};
leading=u5*P^5;
coefficient4=leading*(u+2*x)+u4*P^4;
coefficient3=leading*(2*x*u+x^2+2*y)-P;
coefficient1=leading*(2*x*y*u+y^2)-1;
rootB=2*leading*((x^2+2*y)*u+2*x*y);
rootC=2*leading*y^2*u;
parameterIdeal=ideal(
    coefficient4,
    coefficient3,
    coefficient1,
    B-rootB,
    C-rootC
);
targetIdeal=trim eliminate({{x,y,u}},parameterIdeal);
singularIdeal=ideal(H,diff(P,H),diff(B,H),diff(C,H));
parameterIsPrime=isPrime parameterIdeal;
assert parameterIsPrime;
print("R5ROOT_TASK=two-double-root");
print("R5ROOT_RING=" | toString coefficientRing R);
print("R5ROOT_TARGET_DIM=" | toString(dim targetIdeal-3));
print("R5ROOT_DEGREE=" | toString degree targetIdeal);
print(
    "R5ROOT_PARAMETER_IS_PRIME="
    | toString parameterIsPrime
);
print("R5ROOT_IS_PRIME_BY_CONTRACTION=true");
print(
    "R5ROOT_CONTAINS_SINGULAR_IDEAL="
    | toString isSubset(singularIdeal,targetIdeal)
);
print("R5ROOT_GENERATORS=" | toString flatten entries gens targetIdeal);
"""
    else:
        raise AssertionError(task)

    result = subprocess.run(
        ["M2", "--no-readline", "--silent"],
        input=program + "\nexit 0;\n",
        text=True,
        capture_output=True,
        check=False,
    )
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    assert "error:" not in result.stdout + result.stderr
    assert "R5ROOT_TARGET_DIM=1" in result.stdout
    assert "R5ROOT_PARAMETER_IS_PRIME=true" in result.stdout
    assert "R5ROOT_IS_PRIME_BY_CONTRACTION=true" in result.stdout
    assert "R5ROOT_CONTAINS_SINGULAR_IDEAL=true" in result.stdout
    raise SystemExit(result.returncode)


if "--research-koszul-homology" in sys.argv:
    run_macaulay2_research("koszul-homology")

if "--research-singular-primes" in sys.argv:
    run_macaulay2_research("singular-primes")

if "--research-singular-boundary" in sys.argv:
    run_macaulay2_research("singular-boundary")

if "--research-triple-root-prime" in sys.argv:
    run_root_partition_research("triple-root")

if "--research-two-double-root-prime" in sys.argv:
    run_root_partition_research("two-double-root")


def verify_generic_fibre_newton_topology() -> None:
    """Certify Newton nondegeneracy and the generic-fibre Euler characteristic."""

    from fractions import Fraction
    from itertools import combinations

    origin = (0, 0, 0)
    support_coefficients = {
        monomial: sp.Rational(coefficient)
        for monomial, coefficient in H_polynomial.terms()
    }
    support_coefficients[origin] = -sp.symbols("h")

    def subtract(
        left: tuple[int, ...],
        right: tuple[int, ...],
    ) -> tuple[int, ...]:
        return tuple(
            left_entry - right_entry
            for left_entry, right_entry in zip(left, right, strict=True)
        )

    def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
        return sum(
            left_entry * right_entry
            for left_entry, right_entry in zip(left, right, strict=True)
        )

    def cross(
        left: tuple[int, int, int],
        right: tuple[int, int, int],
    ) -> tuple[int, int, int]:
        return (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )

    def primitive_hyperplane(
        normal: tuple[int, int, int],
        constant: int,
    ) -> tuple[tuple[int, int, int], int]:
        divisor = gcd(
            gcd(abs(normal[0]), abs(normal[1])),
            gcd(abs(normal[2]), abs(constant)),
        )
        divisor = divisor or 1
        return tuple(entry // divisor for entry in normal), constant // divisor

    def facets_three(
        points: tuple[tuple[int, int, int], ...],
    ) -> dict[
        tuple[tuple[int, int, int], int],
        frozenset[tuple[int, int, int]],
    ]:
        facets = {}
        for first, second, third in combinations(points, 3):
            normal = cross(
                subtract(second, first),
                subtract(third, first),
            )
            if normal == origin:
                continue
            constant = dot(normal, first)
            values = [dot(normal, point) - constant for point in points]
            if min(values) < 0 < max(values):
                continue
            if max(values) > 0:
                normal = tuple(-entry for entry in normal)
                constant = -constant
            normal, constant = primitive_hyperplane(normal, constant)
            facets[(normal, constant)] = frozenset(
                point
                for point in points
                if dot(normal, point) == constant
            )
        return facets

    def matrix_rank(rows: list[tuple[int, int, int]]) -> int:
        return int(sp.Matrix(rows).rank())

    def three_dimensional_data(
        points: tuple[tuple[int, int, int], ...],
    ) -> tuple[
        set[frozenset[tuple[int, int, int]]],
        Fraction,
        tuple[tuple[int, int, int], ...],
    ]:
        facets = facets_three(points)
        vertices = tuple(
            point
            for point in points
            if matrix_rank(
                [
                    normal
                    for normal, constant in facets
                    if dot(normal, point) == constant
                ]
            )
            == 3
        )
        faces = {frozenset(points)}
        facet_point_sets = tuple(facets.values())
        for count in range(1, len(facet_point_sets) + 1):
            for selected in combinations(facet_point_sets, count):
                intersection = frozenset.intersection(*selected)
                if intersection:
                    faces.add(intersection)

        centre = tuple(
            Fraction(sum(point[index] for point in vertices), len(vertices))
            for index in range(3)
        )
        normalized_volume = Fraction(0)
        for (normal, constant), facet_points in facets.items():
            facet_vertices = [
                point for point in vertices if point in facet_points
            ]
            dropped_coordinate = max(
                range(3),
                key=lambda index: abs(normal[index]),
            )
            retained = [
                index for index in range(3) if index != dropped_coordinate
            ]
            projected_centre = tuple(
                sum(point[index] for point in facet_vertices)
                / len(facet_vertices)
                for index in retained
            )
            facet_vertices.sort(
                key=lambda point: float(
                    sp.atan2(
                        point[retained[1]] - projected_centre[1],
                        point[retained[0]] - projected_centre[0],
                    )
                )
            )
            for index in range(1, len(facet_vertices) - 1):
                tetrahedron = sp.Matrix(
                    [
                        [
                            sp.Rational(point[coordinate])
                            - sp.Rational(
                                centre[coordinate].numerator,
                                centre[coordinate].denominator,
                            )
                            for coordinate in range(3)
                        ]
                        for point in (
                            facet_vertices[0],
                            facet_vertices[index],
                            facet_vertices[index + 1],
                        )
                    ]
                )
                determinant = tetrahedron.det()
                normalized_volume += abs(
                    Fraction(
                        int(sp.numer(determinant)),
                        int(sp.denom(determinant)),
                    )
                )
        assert normalized_volume.denominator == 1
        return faces, normalized_volume, vertices

    def convex_hull_two(
        points: tuple[tuple[int, int], ...],
    ) -> tuple[tuple[int, int], ...]:
        unique_points = sorted(set(points))

        def orientation(
            first: tuple[int, int],
            second: tuple[int, int],
            third: tuple[int, int],
        ) -> int:
            return (
                (second[0] - first[0]) * (third[1] - first[1])
                - (second[1] - first[1]) * (third[0] - first[0])
            )

        lower: list[tuple[int, int]] = []
        for point in unique_points:
            while (
                len(lower) >= 2
                and orientation(lower[-2], lower[-1], point) <= 0
            ):
                lower.pop()
            lower.append(point)
        upper: list[tuple[int, int]] = []
        for point in reversed(unique_points):
            while (
                len(upper) >= 2
                and orientation(upper[-2], upper[-1], point) <= 0
            ):
                upper.pop()
            upper.append(point)
        return tuple(lower[:-1] + upper[:-1])

    def two_dimensional_data(
        points: tuple[tuple[int, int], ...],
    ) -> tuple[set[frozenset[tuple[int, int]]], int]:
        vertices = convex_hull_two(points)
        faces = {frozenset(points)}
        for index, vertex in enumerate(vertices):
            next_vertex = vertices[(index + 1) % len(vertices)]
            direction = subtract(next_vertex, vertex)
            faces.add(
                frozenset(
                    point
                    for point in points
                    if (
                        direction[0] * (point[1] - vertex[1])
                        - direction[1] * (point[0] - vertex[0])
                    )
                    == 0
                    and min(vertex[0], next_vertex[0])
                    <= point[0]
                    <= max(vertex[0], next_vertex[0])
                    and min(vertex[1], next_vertex[1])
                    <= point[1]
                    <= max(vertex[1], next_vertex[1])
                )
            )
            faces.add(frozenset({vertex}))
        normalized_area = abs(
            sum(
                vertices[index][0]
                * vertices[(index + 1) % len(vertices)][1]
                - vertices[index][1]
                * vertices[(index + 1) % len(vertices)][0]
                for index in range(len(vertices))
            )
        )
        return faces, normalized_area

    active_coordinate_sets = tuple(
        selected
        for count in range(1, 4)
        for selected in combinations(range(3), count)
    )
    all_face_cases: list[
        tuple[tuple[int, ...], frozenset[tuple[int, ...]]]
    ] = []
    normalized_volumes: dict[tuple[int, ...], int] = {}
    full_vertices: tuple[tuple[int, int, int], ...] | None = None
    for active in active_coordinate_sets:
        projected_points = tuple(
            sorted(
                {
                    tuple(monomial[index] for index in active)
                    for monomial in support_coefficients
                    if all(
                        monomial[index] == 0
                        for index in range(3)
                        if index not in active
                    )
                }
            )
        )
        dimension = len(active)
        if dimension == 1:
            minimum = min(point[0] for point in projected_points)
            maximum = max(point[0] for point in projected_points)
            faces = {
                frozenset(projected_points),
                frozenset({(minimum,)}),
                frozenset({(maximum,)}),
            }
            normalized_volume = maximum - minimum
        elif dimension == 2:
            faces, normalized_volume = two_dimensional_data(
                projected_points
            )
        else:
            faces, volume_fraction, full_vertices = three_dimensional_data(
                projected_points
            )
            normalized_volume = int(volume_fraction)
        normalized_volumes[active] = normalized_volume
        all_face_cases.extend((active, face) for face in faces)

    assert full_vertices == (
        (0, 0, 0),
        (0, 2, 0),
        (0, 3, 1),
        (2, 0, 2),
        (2, 4, 0),
        (2, 5, 1),
        (4, 3, 0),
        (4, 4, 1),
        (8, 0, 0),
        (12, 0, 3),
        (12, 0, 4),
    )
    assert normalized_volumes == {
        (0,): 8,
        (1,): 2,
        (2,): 0,
        (0, 1): 38,
        (0, 2): 52,
        (1, 2): 2,
        (0, 1, 2): 328,
    }

    h = sp.symbols("h")
    coordinate_names = ("P", "B", "C")
    singular_checks = [
        "ring r=(0,h),(P,B,C,z),dp;",
        "option(redSB);",
    ]
    nontrivial_face_count = 0
    for case_index, (active, face) in enumerate(all_face_cases):
        if len(face) == 1:
            continue
        face_expression = sp.Integer(0)
        for projected_monomial in face:
            full_monomial = [0, 0, 0]
            for coordinate, exponent in zip(
                active,
                projected_monomial,
                strict=True,
            ):
                full_monomial[coordinate] = exponent
            full_monomial_tuple = tuple(full_monomial)
            coefficient = (
                -h
                if full_monomial_tuple == origin
                else support_coefficients[full_monomial_tuple]
            )
            face_expression += coefficient * sp.prod(
                variable**exponent
                for variable, exponent in zip(
                    TARGET_VARIABLES,
                    full_monomial_tuple,
                    strict=True,
                )
            )
        polynomial_name = f"f{case_index}"
        ideal_name = f"I{case_index}"
        basis_name = f"G{case_index}"
        singular_expression = str(sp.expand(face_expression)).replace(
            "**",
            "^",
        )
        derivatives = ",".join(
            f"diff({polynomial_name},{coordinate_names[index]})"
            for index in active
        )
        torus_product = "*".join(
            coordinate_names[index] for index in active
        )
        singular_checks.extend(
            [
                f"poly {polynomial_name}={singular_expression};",
                (
                    f"ideal {ideal_name}={polynomial_name},{derivatives},"
                    f"z*{torus_product}-1;"
                ),
                f"ideal {basis_name}=std({ideal_name});",
                (
                    f'if (reduce(1,{basis_name})!=0) '
                    f'{{ print("R5NEWTON_DEGENERATE_CASE={case_index}"); '
                    "exit(1); }"
                ),
            ]
        )
        nontrivial_face_count += 1
    singular_checks.append(
        f'print("R5NEWTON_NONTRIVIAL_FACES={nontrivial_face_count}");'
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(singular_checks) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "R5NEWTON_DEGENERATE_CASE" not in result.stdout
    assert (
        f"R5NEWTON_NONTRIVIAL_FACES={nontrivial_face_count}"
        in result.stdout
    )

    euler_contributions = {
        active: (-1) ** (len(active) - 1) * volume
        for active, volume in normalized_volumes.items()
    }
    assert euler_contributions == {
        (0,): 8,
        (1,): 2,
        (2,): 0,
        (0, 1): -38,
        (0, 2): -52,
        (1, 2): -2,
        (0, 1, 2): 328,
    }
    assert sum(euler_contributions.values()) == 246
    print(
        "PASS: all "
        f"{nontrivial_face_count} nontrivial Newton faces are torus nondegenerate"
    )
    print(
        "PASS: coordinate-stratum normalized volumes are "
        "8,2,0;38,52,2;328"
    )
    print("PASS: the generic fibre has Euler characteristic 246")
    print("SCOPE: the H^2=0 stable-rigidity criterion is unavailable")


if "--generic-fibre-newton" in sys.argv:
    verify_generic_fibre_newton_topology()
    raise SystemExit(0)


def verify_graded_logarithmic_module() -> None:
    """Verify the homogenized logarithmic module resolution in Macaulay2."""

    macaulay2_H = macaulay2_homogenized_h()
    macaulay2_program = f"""
R=QQ[P,B,C,Z];
H={macaulay2_H};
HP=diff(P,H); HB=diff(B,H); HC=diff(C,H);
F=R^{{3:-15,1:-16}};
phi=map(R^1,F,matrix{{{{HP,HB,HC,-H}}}});
M=kernel phi;
Cplx=res M;
print("R5LOG_G0=" | toString degrees source gens M);
print("R5LOG_G1=" | toString degrees source Cplx.dd_1);
print("R5LOG_G2=" | toString degrees source Cplx.dd_2);
print("R5LOG_LENGTH=" | toString length Cplx);
exit 0;
"""
    result = subprocess.run(
        ["M2", "--no-readline", "--silent"],
        input=macaulay2_program,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    compact_output = "".join(result.stdout.split())
    expected_generators = (
        "R5LOG_G0={{22},{22},"
        + ",".join("{23}" for _index in range(13))
        + "}"
    )
    expected_relations = (
        "R5LOG_G1={"
        + ",".join("{24}" for _index in range(18))
        + "}"
    )
    expected_second_relations = (
        "R5LOG_G2={"
        + ",".join("{25}" for _index in range(6))
        + "}"
    )
    assert expected_generators in compact_output, result.stdout
    assert expected_relations in compact_output, result.stdout
    assert expected_second_relations in compact_output, result.stdout
    assert "R5LOG_LENGTH=2" in compact_output, result.stdout
    print(
        "PASS: the graded logarithmic module has "
        "2 degree-7 and 13 degree-8 generators"
    )
    print(
        "PASS: its minimal relations have shifted Betti counts "
        "15, 18, 6"
    )
    print(
        "PASS: the filtered Hilbert numerator is "
        "2*t^7+13*t^8-18*t^9+6*t^10"
    )


if "--module-resolution" in sys.argv:
    verify_graded_logarithmic_module()
    raise SystemExit(0)


def monomial_exponents(bound: int) -> list[tuple[int, int, int]]:
    """Exponent triples of total degree at most ``bound``."""

    if bound < 0:
        return []
    return [
        (first, second, third)
        for first in range(bound + 1)
        for second in range(bound + 1 - first)
        for third in range(bound + 1 - first - second)
    ]


def add_exponents(
    left: tuple[int, int, int],
    right: tuple[int, int, int],
) -> tuple[int, int, int]:
    return tuple(
        left_entry + right_entry
        for left_entry, right_entry in zip(left, right, strict=True)
    )


def logarithmic_matrix(
    bound: int,
) -> tuple[
    fmpz_mat,
    list[tuple[int, int, int]],
    list[tuple[int, int, int]],
]:
    """Integer matrix for V(H)-QH-kappa=0."""

    vector_monomials = monomial_exponents(bound)
    quotient_monomials = monomial_exponents(bound - 1)
    vector_block_size = len(vector_monomials)
    quotient_offset = 3 * vector_block_size
    column_count = quotient_offset + len(quotient_monomials) + 1
    rows: dict[
        tuple[int, int, int],
        dict[int, int],
    ] = {}

    def add_coefficient(
        monomial: tuple[int, int, int],
        column: int,
        coefficient: int,
    ) -> None:
        row = rows.setdefault(monomial, {})
        row[column] = row.get(column, 0) + int(coefficient)
        if row[column] == 0:
            del row[column]

    for component, derivative in enumerate(primitive_gradient):
        for vector_index, vector_monomial in enumerate(vector_monomials):
            column = component * vector_block_size + vector_index
            for derivative_monomial, coefficient in derivative.terms():
                add_coefficient(
                    add_exponents(vector_monomial, derivative_monomial),
                    column,
                    coefficient,
                )

    for quotient_index, quotient_monomial in enumerate(quotient_monomials):
        column = quotient_offset + quotient_index
        for h_monomial, coefficient in primitive_H.terms():
            add_coefficient(
                add_exponents(quotient_monomial, h_monomial),
                column,
                -coefficient,
            )

    add_coefficient((0, 0, 0), column_count - 1, -1)
    ordered_rows = sorted(rows)
    entries: list[int] = []
    for monomial in ordered_rows:
        row = rows[monomial]
        entries.extend(row.get(column, 0) for column in range(column_count))

    return (
        fmpz_mat(len(ordered_rows), column_count, entries),
        vector_monomials,
        quotient_monomials,
    )


expected_nullities = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    2,
    21,
    54,
    104,
    174,
    267,
)
degree_seven_data = None
degree_eight_data = None
newton_data = {}
for degree_bound, expected_nullity in enumerate(expected_nullities):
    matrix, vector_monomials, quotient_monomials = logarithmic_matrix(
        degree_bound
    )
    nullspace, nullity = matrix.nullspace()
    assert nullity == expected_nullity
    if degree_bound == 7:
        degree_seven_data = (
            nullspace,
            vector_monomials,
            quotient_monomials,
        )
    if degree_bound == 8:
        degree_eight_data = (
            nullspace,
            vector_monomials,
            quotient_monomials,
        )
    if degree_bound in (8, 9, 10, 11, 12):
        newton_data[degree_bound] = (
            nullspace,
            vector_monomials,
            quotient_monomials,
        )

assert degree_seven_data is not None
assert degree_eight_data is not None
(
    degree_seven_nullspace,
    degree_seven_vector_monomials,
    degree_seven_quotient_monomials,
) = degree_seven_data
vector_block_size = len(degree_seven_vector_monomials)
quotient_offset = 3 * vector_block_size
column_count = (
    quotient_offset
    + len(degree_seven_quotient_monomials)
    + 1
)

vector_basis: list[sp.Matrix] = []
for basis_index in range(2):
    coefficients = [
        int(degree_seven_nullspace[row, basis_index])
        for row in range(column_count)
    ]
    content = sp.gcd_list(
        [coefficient for coefficient in coefficients if coefficient]
    )
    coefficients = [
        coefficient // content
        for coefficient in coefficients
    ]
    assert coefficients[-1] == 0

    components = []
    for component in range(3):
        offset = component * vector_block_size
        components.append(
            sum(
                coefficients[offset + monomial_index]
                * P**first
                * B**second
                * C**third
                for monomial_index, (first, second, third) in enumerate(
                    degree_seven_vector_monomials
                )
            )
        )
    vector_basis.append(sp.Matrix(components))

selected_evaluation = sp.Matrix.hstack(
    *[
        vector.subs(SELECTED_TARGET)
        for vector in vector_basis
    ]
)
assert selected_evaluation.rank() == 2

# Exact preservation H(T)=H eliminates this apparent stable frontier before
# any determinant condition is imposed.  The unique top term of H is
# P^12*C^4, and no nonzero combination of the two fields kills the
# degree-seven leading part of either its P- or C-component.
for component in (0, 2):
    leading_monomials = sorted(
        {
            monomial
            for vector in vector_basis
            for monomial, _coefficient in sp.Poly(
                vector[component],
                *TARGET_VARIABLES,
            ).terms()
            if sum(monomial) == 7
        }
    )
    leading_matrix = sp.Matrix(
        [
            [
                sp.Poly(
                    vector[component],
                    *TARGET_VARIABLES,
                ).coeff_monomial(monomial)
                for vector in vector_basis
            ]
            for monomial in leading_monomials
        ]
    )
    assert leading_matrix.rank() == 2


# The next logarithmic layer is already much larger.  It evaluates
# surjectively at the selected target, but the unique top term of H forces
# every exact-boundary candidate into the kernel of at least one of the
# leading P- and C-component maps.
(
    degree_eight_nullspace,
    degree_eight_vector_monomials,
    degree_eight_quotient_monomials,
) = degree_eight_data
degree_eight_block_size = len(degree_eight_vector_monomials)
degree_eight_quotient_offset = 3 * degree_eight_block_size
degree_eight_column_count = (
    degree_eight_quotient_offset
    + len(degree_eight_quotient_monomials)
    + 1
)
assert all(
    int(degree_eight_nullspace[degree_eight_column_count - 1, index]) == 0
    for index in range(21)
)

degree_eight_evaluations: list[sp.Matrix] = []
for basis_index in range(21):
    components = []
    for component in range(3):
        offset = component * degree_eight_block_size
        components.append(
            sum(
                int(
                    degree_eight_nullspace[
                        offset + monomial_index,
                        basis_index,
                    ]
                )
                * P**first
                * B**second
                * C**third
                for monomial_index, (first, second, third) in enumerate(
                    degree_eight_vector_monomials
                )
            )
        )
    degree_eight_evaluations.append(
        sp.Matrix(components).subs(SELECTED_TARGET)
    )
assert sp.Matrix.hstack(*degree_eight_evaluations).rank() == 3

degree_eight_index = {
    monomial: index
    for index, monomial in enumerate(degree_eight_vector_monomials)
}
degree_eight_top_monomials = [
    monomial
    for monomial in degree_eight_vector_monomials
    if sum(monomial) == 8
]
leading_blocks = []
for component in (0, 2):
    block = sp.Matrix(
        [
            [
                int(
                    degree_eight_nullspace[
                        component * degree_eight_block_size
                        + degree_eight_index[monomial],
                        basis_index,
                    ]
                )
                for basis_index in range(21)
            ]
            for monomial in degree_eight_top_monomials
        ]
    )
    assert block.rank() == 7
    leading_blocks.append(block)
assert sp.Matrix.vstack(*leading_blocks).rank() == 12


def canonical_subspace(basis: sp.Matrix) -> tuple[tuple[sp.Rational, ...], ...]:
    """Canonical key for the column space of an exact rational matrix."""

    reduced, _pivots = basis.T.rref()
    return tuple(tuple(row) for row in reduced.tolist())


def newton_face_pruning(
    nullspace: fmpz_mat,
    vector_monomials: list[tuple[int, int, int]],
    quotient_monomials: list[tuple[int, int, int]],
    trace: list[
        tuple[
            int,
            tuple[int | None, int | None, int | None],
            tuple[tuple[int, int, int], ...],
        ]
    ]
    | None = None,
) -> tuple[int, int]:
    """Exhaust exact-boundary candidates by recursive exposed-face tests.

    A node is a rational linear subspace of logarithmic vector fields.  Its
    componentwise maximal degrees define weights for ``x+H*V``.  If the
    exposed face of ``H`` is a single monomial, exact preservation forces at
    least one participating leading component to vanish.  Branching into
    those exact kernels and repeating is exhaustive.
    """

    block_size = len(vector_monomials)
    column_count = (
        3 * block_size
        + len(quotient_monomials)
        + 1
    )
    # python-flint stores a nullspace basis in the first ``rank`` columns of
    # a square container; the remaining columns are zero padding.
    nullity = nullspace.rank()
    assert all(
        int(nullspace[column_count - 1, index]) == 0
        for index in range(nullity)
    )
    vector_coefficient_matrix = sp.Matrix(
        [
            [
                int(nullspace[row, column])
                for column in range(nullity)
            ]
            for row in range(3 * block_size)
        ]
    )
    # Work in logarithmic parameter coordinates.  Multiplication by
    # ``vector_coefficient_matrix`` recovers monomial coefficients, while
    # canonicalization and kernel intersections stay in dimension at most
    # the logarithmic nullity.
    initial_basis = sp.eye(nullity)

    degree_rows = [
        [
            index
            for index, monomial in enumerate(vector_monomials)
            if sum(monomial) == degree
        ]
        for degree in range(
            max(sum(monomial) for monomial in vector_monomials) + 1
        )
    ]
    h_support = [monomial for monomial, _coefficient in primitive_H.terms()]

    queue = [initial_basis]
    seen: set[tuple[tuple[sp.Rational, ...], ...]] = set()
    terminal_nodes = 0
    while queue:
        basis = queue.pop()
        key = canonical_subspace(basis)
        if key in seen:
            continue
        seen.add(key)

        component_degrees = []
        leading_blocks_at_node: list[sp.Matrix | None] = []
        for component in range(3):
            component_offset = component * block_size
            degree = -1
            leading_block = None
            for candidate_degree in range(len(degree_rows) - 1, -1, -1):
                block = vector_coefficient_matrix[
                    [
                        component_offset + row
                        for row in degree_rows[candidate_degree]
                    ],
                    :,
                ] * basis
                if not block.is_zero_matrix:
                    degree = candidate_degree
                    leading_block = block
                    break
            component_degrees.append(degree)
            leading_blocks_at_node.append(leading_block)

        weights = tuple(
            1 if degree < 0 else primitive_H.total_degree() + degree
            for degree in component_degrees
        )
        face_weight = max(
            sum(
                exponent * weight
                for exponent, weight in zip(monomial, weights, strict=True)
            )
            for monomial in h_support
        )
        exposed_face = [
            monomial
            for monomial in h_support
            if sum(
                exponent * weight
                for exponent, weight in zip(monomial, weights, strict=True)
            )
            == face_weight
        ]
        if trace is not None:
            finite_degrees = [
                degree
                for degree in component_degrees
                if degree >= 0
            ]
            baseline = max(finite_degrees)
            normalized_degrees = tuple(
                None if degree < 0 else degree - baseline
                for degree in component_degrees
            )
            trace.append(
                (
                    basis.cols,
                    normalized_degrees,
                    tuple(exposed_face),
                )
            )
        if len(exposed_face) != 1:
            terminal_nodes += 1
            continue

        exposed_monomial = exposed_face[0]
        children = []
        active_components = 0
        for component, exponent in enumerate(exposed_monomial):
            leading_block = leading_blocks_at_node[component]
            if exponent == 0 or leading_block is None:
                continue
            active_components += 1
            kernel = leading_block.nullspace()
            if not kernel:
                continue
            child = basis * sp.Matrix.hstack(*kernel)
            assert 0 < child.cols < basis.cols
            children.append(child)

        if active_components == 0:
            terminal_nodes += 1
            continue
        if not children:
            # Every nonzero vector has all leading components required by the
            # unique exposed monomial, so its top term cannot cancel.
            continue
        queue.extend(children)

    return len(seen), terminal_nodes


newton_pruning_results = {}
newton_pruning_traces = {}
for degree_bound in (8, 9, 10, 11, 12):
    trace = []
    newton_pruning_results[degree_bound] = newton_face_pruning(
        *newton_data[degree_bound],
        trace=trace,
    )
    newton_pruning_traces[degree_bound] = trace
assert newton_pruning_results == {
    8: (10, 0),
    9: (20, 0),
    10: (33, 0),
    11: (56, 0),
    12: (81, 0),
}, newton_pruning_results
assert {
    face[0]
    for trace in newton_pruning_traces.values()
    for _dimension, _normalized_degrees, face in trace
} == {
    (12, 0, 4),
    (2, 5, 1),
}


# In the unstabilized target, T=id+H*(lambda*V_1+mu*V_2).  The degree-seven
# logarithmic space has kappa=0, so a polynomial automorphism must have
# determinant one.  Three exact evaluations already force lambda=mu=0.
generic_vector = LAMBDA * vector_basis[0] + MU * vector_basis[1]
gradient_H = sp.Matrix(
    [
        sp.diff(primitive_H_expression, variable)
        for variable in TARGET_VARIABLES
    ]
)
jacobian_vector = generic_vector.jacobian(TARGET_VARIABLES)
evaluation_points = (
    (1, 1, 0),
    (1, 0, 1),
    (1, 1, 1),
)
jacobian_equations: list[sp.Expr] = []
for point in evaluation_points:
    substitution = dict(zip(TARGET_VARIABLES, point, strict=True))
    derivative = (
        sp.eye(3)
        + generic_vector.subs(substitution)
        * gradient_H.subs(substitution).T
        + primitive_H_expression.subs(substitution)
        * jacobian_vector.subs(substitution)
    )
    equation = sp.primitive(
        sp.Poly(
            sp.expand(derivative.det() - 1),
            LAMBDA,
            MU,
            domain=sp.QQ,
        )
    )[1].as_expr()
    jacobian_equations.append(equation)

jacobian_groebner = sp.groebner(
    jacobian_equations,
    LAMBDA,
    MU,
    order="grevlex",
    domain=sp.QQ,
)
assert set(jacobian_groebner) == {LAMBDA, MU}


print("PASS: the fixed quintic coefficient-torus stabilizer is trivial")
print("PASS: the prime ramified discriminant has degree 16 and 59 terms")
print("PASS: an exact triple-root singular point forces kappa=0 in every degree")
print("PASS: logarithmic corrections vanish through vector degree six")
print("PASS: the degree-seven logarithmic frontier has dimension two")
print("PASS: its evaluation moves the selected target in rank two")
print("PASS: exact boundary preservation kills the stable degree-23 frontier")
print("PASS: three Jacobian evaluations kill the unstabilized frontier")
print("PASS: the degree-24 logarithmic frontier has dimension 21 and evaluation rank 3")
print("PASS: recursive Newton-face pruning kills target degrees 24 through 28")
print("PASS: every branch uses one of two exposed Newton vertices")
print("PASS: the positive upper Newton hull is the edge joining those vertices")
print("PASS: P^2*B^5*C is an exposed Newton intruder")
print("PASS: the intruder satisfies Kuroda's stable-invariant hypothesis")
print("PASS: the unstabilized target self-equivalence is trivial in every degree")
print("PASS: an exact logarithmic example makes the binomial wall cancellable")
print("PASS: the P-zero Koszul ladder first ties at target degree 50")
print("PASS: its degree-50 tie fails the UFD power condition")
print("PASS: its target-degree-55 rung has exact leading cancellation")
print("PASS: the stable target orbit is a point through target degree 28")
print("PASS: the unstabilized target self-equivalence is trivial through degree 28")
print(
    "SCOPE: full vertical wall classification remains open; "
    "Kuroda descent closes the standard marked orbit"
)
