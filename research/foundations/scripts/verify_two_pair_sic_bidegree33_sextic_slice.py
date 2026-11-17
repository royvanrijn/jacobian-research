#!/usr/bin/env python3
"""Exact first frontier for two-pair SIC in bidegree (3,3).

The full sixteen-coefficient one-sided nullcone is parametrized and
eliminated over Q.  The first thirteen full moments are proved
algebraically independent by an exact nonzero Jacobian minor, but an exact
Hilbert-series coefficient proves that their degrees cannot be those of a
homogeneous system of parameters.  Replacing moment 13 by moment 14 gives
another exact rank-thirteen Jacobian certificate and the first
Hilbert-compatible minimal candidate.  The script then closes each of the
three pure irreducible SL_2 summands.  On the binary-sextic summand, moments
2, 4, 6, and 10 cut out exactly the nullcone L^4 Q; on the binary-quartic
and binary-quadratic summands, moments (2,3) and moment 2 respectively
generate their nullcone ideals.
"""

from __future__ import annotations

from collections import defaultdict
import json
from math import factorial, gcd
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_frontier.json"
)
DEGREE = 3
MONOMIALS = tuple(
    (i, j) for i in range(DEGREE + 1) for j in range(DEGREE + 1)
)
POSITIVE_POSITIONS = tuple(
    (i, j)
    for i in range(DEGREE + 1)
    for j in range(DEGREE + 1)
    if i > j
)
SEXTIC_MOMENT_ORDERS = (2, 4, 6, 10)
SEXTIC_NULLCONE_GENERATORS = (
    "10*s3^2-15*s2*s4+6*s1*s5-s0*s6",
    (
        "5*s2*s3*s5-9*s1*s4*s5+5*s0*s5^2-10*s2^2*s6"
        "+15*s1*s3*s6-6*s0*s4*s6"
    ),
    (
        "25*s2*s3*s4-45*s1*s4^2-45*s2^2*s5+64*s1*s3*s5"
        "+5*s0*s4*s5+5*s1*s2*s6-9*s0*s3*s6"
    ),
    (
        "5*s1*s3*s4-10*s0*s4^2-9*s1*s2*s5+15*s0*s3*s5"
        "+5*s1^2*s6-6*s0*s2*s6"
    ),
    (
        "45*s2^2*s5^2-64*s1*s3*s5^2+20*s0*s4*s5^2"
        "-50*s2^2*s4*s6+120*s0*s4^2*s6+130*s1*s2*s5*s6"
        "-216*s0*s3*s5*s6-75*s1^2*s6^2+90*s0*s2*s6^2"
    ),
    (
        "75*s2^2*s4*s5-180*s0*s4^2*s5-192*s1*s2*s5^2"
        "+320*s0*s3*s5^2-100*s2^2*s3*s6+225*s1*s2*s4*s6"
        "-60*s0*s3*s4*s6-103*s0*s2*s5*s6+15*s0*s1*s6^2"
    ),
    (
        "75*s1*s2*s4*s5-100*s0*s3*s4*s5-192*s1^2*s5^2"
        "+225*s0*s2*s5^2-100*s1*s2*s3*s6+225*s1^2*s4*s6"
        "-50*s0*s2*s4*s6-101*s0*s1*s5*s6+18*s0^2*s6^2"
    ),
    (
        "1875*s2^2*s4^2-4500*s0*s4^3+1600*s0*s3*s4*s5"
        "-12288*s1^2*s5^2+14400*s0*s2*s5^2-4500*s2^3*s6"
        "+1600*s1*s2*s3*s6+14400*s1^2*s4*s6"
        "-8650*s0*s2*s4*s6-4864*s0*s1*s5*s6+927*s0^2*s6^2"
    ),
    (
        "75*s1*s2*s4^2-100*s0*s3*s4^2-192*s1^2*s4*s5"
        "+225*s0*s2*s4*s5-180*s1*s2^2*s6+320*s1^2*s3*s6"
        "-60*s0*s2*s3*s6-103*s0*s1*s4*s6+15*s0^2*s5*s6"
    ),
    (
        "45*s1^2*s4^2-50*s0*s2*s4^2-64*s1^2*s3*s5"
        "+130*s0*s1*s4*s5-75*s0^2*s5^2+20*s1^2*s2*s6"
        "+120*s0*s2^2*s6-216*s0*s1*s3*s6+90*s0^2*s4*s6"
    ),
)
SEXTIC_RADICAL_POWERS = (1,) + (5,) * 9
FULL_MOMENT_ORDERS = tuple(range(1, 14))
CORRECTED_MOMENT_ORDERS = tuple(range(1, 13)) + (14,)
FULL_JACOBIAN_POINT = (
    (2, 2, -4, 0),
    (4, 3, 2, 0),
    (3, 1, -1, 4),
    (-2, 0, -2, -3),
)
FULL_JACOBIAN_COLUMNS = tuple(range(3, 16))
FULL_JACOBIAN_DETERMINANT = -(
    2**256
    * 3**107
    * 5**48
    * 7**29
    * 11**13
    * 13**7
    * 17**3
    * 19
    * 139
    * 4493
    * 886069
    * 651443921434147
    * 108355984865758174686774716693198468303195999878781902931
)
CORRECTED_JACOBIAN_COLUMNS = tuple(range(13))
CORRECTED_JACOBIAN_DETERMINANT = int(
    "-177773463337872042258599314936100079730243854705203912887777623860"
    "541231537411339305382169268682266803601791521012996366373976876910"
    "055068237469775158074527906699099566496272363910115715195225686054"
    "006682543315888455588033482573001213880570668925999015526400000000"
    "000000000000000000000000000000000000000"
)
INVARIANT_WEIGHTS = (
    0,
    0,
    0,
    0,
    2,
    2,
    2,
    -2,
    -2,
    -2,
    4,
    4,
    -4,
    -4,
    6,
    -6,
)
HILBERT_SERIES_CUTOFF = 100


def singular_expression(expression: sp.Expr) -> str:
    """Serialize a SymPy polynomial for Singular."""

    return str(sp.expand(expression)).replace("**", "^").replace(" ", "")


def multiply_dense_bivariate(
    left: list[list[int]], right: tuple[tuple[int, ...], ...]
) -> list[list[int]]:
    """Multiply dense coefficient grids in two auxiliary variables."""

    answer = [
        [0] * (len(left[0]) + len(right[0]) - 1)
        for _ in range(len(left) + len(right) - 1)
    ]
    for left_i, row in enumerate(left):
        for left_j, left_coefficient in enumerate(row):
            if left_coefficient == 0:
                continue
            for right_i, right_row in enumerate(right):
                for right_j, right_coefficient in enumerate(right_row):
                    answer[left_i + right_i][left_j + right_j] += (
                        left_coefficient * right_coefficient
                    )
    return answer


def evaluated_moment_jacobian(moment_orders: tuple[int, ...]) -> sp.Matrix:
    """Evaluate a full moment Jacobian at the fixed integral point.

    If F(x,y)=sum c_(i,j)x^i y^j, then

      mu_m=sum_I (3m-I)! I! [x^I y^I] F(x,y)^m.

    Differentiation with respect to c_(a,b) replaces F^m by
    m*x^a*y^b*F^(m-1).  This evaluates the full Jacobian without expanding
    any sixteen-variable moment polynomial.
    """

    power: list[list[int]] = [[1]]
    rows: list[list[int]] = []
    for order in moment_orders:
        gradient: list[int] = []
        for coefficient_i, coefficient_j in MONOMIALS:
            value = 0
            for diagonal in range(DEGREE * order + 1):
                power_i = diagonal - coefficient_i
                power_j = diagonal - coefficient_j
                if (
                    0 <= power_i < len(power)
                    and 0 <= power_j < len(power[0])
                ):
                    value += (
                        factorial(DEGREE * order - diagonal)
                        * factorial(diagonal)
                        * power[power_i][power_j]
                    )
            gradient.append(order * value)
        rows.append(gradient)
        power = multiply_dense_bivariate(power, FULL_JACOBIAN_POINT)

    return sp.Matrix(rows)


def check_full_moment_jacobian() -> dict[str, object]:
    """Certify algebraic independence of the first thirteen full moments."""

    jacobian = evaluated_moment_jacobian(FULL_MOMENT_ORDERS)
    minor = jacobian[:, FULL_JACOBIAN_COLUMNS]
    determinant = int(minor.det(method="domain-ge"))
    assert determinant == FULL_JACOBIAN_DETERMINANT, (
        determinant,
        FULL_JACOBIAN_DETERMINANT,
    )
    assert jacobian.rank() == 13
    return {
        "moment_orders": list(FULL_MOMENT_ORDERS),
        "coefficient_point": [list(row) for row in FULL_JACOBIAN_POINT],
        "minor_columns_zero_based": list(FULL_JACOBIAN_COLUMNS),
        "minor_variables": [
            f"c_{MONOMIALS[index][0]}{MONOMIALS[index][1]}"
            for index in FULL_JACOBIAN_COLUMNS
        ],
        "determinant": str(determinant),
        "rank": 13,
        "conclusion": (
            "the first thirteen full moments are algebraically independent"
        ),
    }


def check_corrected_moment_jacobian() -> dict[str, object]:
    """Certify independence after replacing moment 13 by moment 14."""

    jacobian = evaluated_moment_jacobian(CORRECTED_MOMENT_ORDERS)
    minor = jacobian[:, CORRECTED_JACOBIAN_COLUMNS]
    determinant = int(minor.det(method="domain-ge"))
    assert determinant == CORRECTED_JACOBIAN_DETERMINANT, (
        determinant,
        CORRECTED_JACOBIAN_DETERMINANT,
    )
    return {
        "moment_orders": list(CORRECTED_MOMENT_ORDERS),
        "coefficient_point": [list(row) for row in FULL_JACOBIAN_POINT],
        "minor_columns_zero_based": list(CORRECTED_JACOBIAN_COLUMNS),
        "minor_variables": [
            f"c_{MONOMIALS[index][0]}{MONOMIALS[index][1]}"
            for index in CORRECTED_JACOBIAN_COLUMNS
        ],
        "determinant": str(determinant),
        "rank": 13,
        "conclusion": (
            "moments 1 through 12 and moment 14 are algebraically "
            "independent"
        ),
    }


def invariant_hilbert_coefficients(cutoff: int) -> list[int]:
    """Return Hilbert coefficients from the exact SL_2 weight formula.

    For an SL_2 module, the multiplicity of the trivial representation in
    a finite-dimensional module is its weight-zero dimension minus its
    weight-two dimension.  Expanding the symmetric-algebra character

      product_w (1-t*q^w)^(-1)

    through ``cutoff`` therefore gives the invariant-ring Hilbert function.
    """

    character = [defaultdict(int) for _ in range(cutoff + 1)]
    character[0][0] = 1
    for weight in INVARIANT_WEIGHTS:
        for degree in range(1, cutoff + 1):
            for previous_weight, multiplicity in tuple(
                character[degree - 1].items()
            ):
                character[degree][previous_weight + weight] += multiplicity
    return [
        weights[0] - weights.get(2, 0)
        for weights in character
    ]


def hilbert_numerator_prefix(
    hilbert: list[int], parameter_degrees: tuple[int, ...]
) -> list[int]:
    """Multiply a Hilbert-series prefix by the proposed denominator."""

    numerator = hilbert.copy()
    for degree in parameter_degrees:
        for index in range(len(numerator) - 1, degree - 1, -1):
            numerator[index] -= numerator[index - degree]
    return numerator


def check_moment_hsop_degree_obstruction() -> dict[str, object]:
    """Rule out degrees 1,...,13 and test the corrected degree candidate."""

    hilbert = invariant_hilbert_coefficients(HILBERT_SERIES_CUTOFF)
    assert hilbert[:14] == [
        1,
        1,
        4,
        8,
        26,
        53,
        146,
        305,
        704,
        1417,
        2920,
        5533,
        10500,
        18825,
    ]

    consecutive = hilbert_numerator_prefix(hilbert, FULL_MOMENT_ORDERS)
    first_negative = next(
        (index, value)
        for index, value in enumerate(consecutive)
        if value < 0
    )
    assert first_negative == (63, -2186)

    corrected = hilbert_numerator_prefix(
        hilbert, CORRECTED_MOMENT_ORDERS
    )
    assert min(corrected) == 0
    assert max(
        index for index, value in enumerate(corrected) if value
    ) == 76
    assert corrected[76] == 1
    assert all(value == 0 for value in corrected[77:])
    assert sum(corrected) == 9226602

    return {
        "representation_weights": list(INVARIANT_WEIGHTS),
        "hilbert_coefficients_0_through_13": hilbert[:14],
        "checked_through_degree": HILBERT_SERIES_CUTOFF,
        "consecutive_degrees_1_through_13": {
            "first_negative_numerator_coefficient": {
                "degree": first_negative[0],
                "coefficient": first_negative[1],
            },
            "conclusion": (
                "no homogeneous system of parameters can have degrees "
                "1 through 13"
            ),
        },
        "corrected_degrees_1_through_12_and_14": {
            "all_checked_coefficients_nonnegative": True,
            "last_checked_nonzero_degree": 76,
            "coefficient_at_degree_76": corrected[76],
            "coefficients_77_through_100_zero": True,
            "checked_numerator_coefficient_sum": sum(corrected),
            "status": (
                "passes this exact necessary Hilbert-series test; this "
                "does not prove the moments form a homogeneous system "
                "of parameters"
            ),
        },
    }


def check_global_quadratic_projection() -> dict[str, object]:
    """Recover the Sym^2 coordinates from all sixteen coefficients."""

    W, V, Z, Y = sp.symbols("W V Z Y")
    pairing = W * Z + V * Y

    def lowering(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            Y * sp.diff(polynomial, Z) - W * sp.diff(polynomial, V)
        )

    columns: list[list[sp.Expr]] = []
    names: list[str] = []
    for summand, prefix in ((3, "s"), (2, "t"), (1, "r"), (0, "u")):
        polynomial = sp.expand(
            V**summand
            * Z**summand
            * pairing ** (DEGREE - summand)
        )
        for order in range(2 * summand + 1):
            expanded = sp.Poly(polynomial, W, V, Z, Y)
            columns.append(
                [
                    expanded.coeff_monomial(
                        W ** (DEGREE - i)
                        * V**i
                        * Z ** (DEGREE - j)
                        * Y**j
                    )
                    for i, j in MONOMIALS
                ]
            )
            names.append(f"{prefix}{order}")
            polynomial = lowering(polynomial) / (order + 1)

    change_of_basis = sp.Matrix(
        len(MONOMIALS),
        len(columns),
        lambda row, column: columns[column][row],
    )
    assert change_of_basis.det() == -31104000
    inverse = change_of_basis.inv()
    expected = {
        "r0": {
            MONOMIALS.index((1, 0)): sp.Rational(3, 10),
            MONOMIALS.index((2, 1)): sp.Rational(1, 5),
            MONOMIALS.index((3, 2)): sp.Rational(3, 10),
        },
        "r1": {
            MONOMIALS.index((0, 0)): sp.Rational(-9, 20),
            MONOMIALS.index((1, 1)): sp.Rational(-1, 20),
            MONOMIALS.index((2, 2)): sp.Rational(1, 20),
            MONOMIALS.index((3, 3)): sp.Rational(9, 20),
        },
        "r2": {
            MONOMIALS.index((0, 1)): sp.Rational(-3, 10),
            MONOMIALS.index((1, 2)): sp.Rational(-1, 5),
            MONOMIALS.index((2, 3)): sp.Rational(-3, 10),
        },
    }
    for name, nonzero_entries in expected.items():
        row = inverse.row(names.index(name))
        assert {
            index: coefficient
            for index, coefficient in enumerate(row)
            if coefficient
        } == nonzero_entries

    return {
        "binary_quadratic": "r0*X^2+2*r1*X*T+r2*T^2",
        "r0": "(3*c10+2*c21+3*c32)/10",
        "r1": "(-9*c00-c11+c22+9*c33)/20",
        "r2": "-(3*c01+2*c12+3*c23)/10",
        "discriminant": "r1^2-r0*r2",
    }


def check_diagonal_moment_slice() -> dict[str, object]:
    """Close the maximal-torus fixed diagonal slice exactly."""

    auxiliary = sp.symbols("q")
    coefficients = sp.symbols("c0:4")
    polynomial = sum(
        coefficient * auxiliary**index
        for index, coefficient in enumerate(coefficients)
    )
    power = sp.Integer(1)
    moments: list[sp.Expr] = []
    for order in range(1, 5):
        power = sp.expand(power * polynomial)
        univariate = sp.Poly(power, auxiliary)
        moment = sp.expand(
            sum(
                factorial(DEGREE * order - diagonal)
                * factorial(diagonal)
                * univariate.coeff_monomial(auxiliary**diagonal)
                for diagonal in range(DEGREE * order + 1)
            )
        )
        _, primitive = sp.Poly(moment, *coefficients).primitive()
        moments.append(primitive.as_expr())

    basis = sp.groebner(moments, *coefficients, order="grevlex")
    assert basis.is_zero_dimensional
    powers: list[int] = []
    for coefficient in coefficients:
        for exponent in range(1, 8):
            if basis.reduce(coefficient**exponent)[1] == 0:
                powers.append(exponent)
                break
        else:
            raise AssertionError(f"no power certificate for {coefficient}")
    assert powers == [7, 7, 7, 7]

    return {
        "slice": "c_ij=0 for i!=j",
        "dimension": 4,
        "moment_orders": [1, 2, 3, 4],
        "groebner_basis_size": len(basis.polys),
        "radical_power_certificate": powers,
        "radical": "(c00,c11,c22,c33)",
        "conclusion": (
            "the torus-fixed diagonal slice has no nonzero moment-zero point"
        ),
    }


def normalized_branch_invariant_counts() -> dict[str, object]:
    """Count residual-torus invariant monomials after setting c=1."""

    weights = (6, 4, 2, 0, -2, -4, -6, 4, 2, 0, -2, -4)
    cutoff = 13
    counts: dict[tuple[int, int], int] = {(0, 0): 1}
    for weight in weights:
        updated: dict[tuple[int, int], int] = defaultdict(int)
        for (degree, total_weight), multiplicity in counts.items():
            for exponent in range(cutoff - degree + 1):
                updated[
                    (degree + exponent, total_weight + exponent * weight)
                ] += multiplicity
        counts = dict(updated)
    by_order = [
        sum(
            multiplicity
            for (degree, total_weight), multiplicity in counts.items()
            if degree <= order and total_weight == 0
        )
        for order in range(2, cutoff + 1)
    ]
    assert by_order == [
        15,
        57,
        192,
        564,
        1508,
        3692,
        8438,
        18146,
        37076,
        72396,
        135918,
        246354,
    ]
    return {
        "normal_form": "F2=2*X*T",
        "variables": 12,
        "equations": "mu_2 through mu_13",
        "residual_torus_weights": list(weights),
        "invariant_monomial_counts_by_order": {
            str(order): count
            for order, count in zip(range(2, cutoff + 1), by_order)
        },
        "interpretation": (
            "raw order-13 expansion has up to 246354 invariant monomials; "
            "use residual-weight sparse elimination"
        ),
    }


def check_sextic_embedding() -> list[tuple[int, int, int, int]]:
    """Construct the V_6 lowering chain and its coefficient map."""

    W, V, Z, Y = sp.symbols("W V Z Y")
    parameters = sp.symbols("s0:7")

    def lowering(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            Y * sp.diff(polynomial, Z) - W * sp.diff(polynomial, V)
        )

    basis: list[sp.Expr] = []
    polynomial = V**3 * Z**3
    for order in range(7):
        basis.append(sp.expand(polynomial))
        polynomial = lowering(polynomial) / (order + 1)

    general = sp.Poly(
        sp.expand(
            sum(
                parameter * vector
                for parameter, vector in zip(parameters, basis)
            )
        ),
        W,
        V,
        Z,
        Y,
    )
    coefficient_map: list[tuple[int, int, int, int]] = []
    for i, j in MONOMIALS:
        coefficient = sp.Poly(
            general.coeff_monomial(
                W ** (DEGREE - i)
                * V**i
                * Z ** (DEGREE - j)
                * Y**j
            ),
            *parameters,
        )
        terms = coefficient.terms()
        assert len(terms) == 1
        parameter_exponents, scalar = terms[0]
        assert sum(parameter_exponents) == 1
        parameter_index = parameter_exponents.index(1)
        coefficient_map.append((i, j, parameter_index, int(scalar)))

    assert coefficient_map == [
        (0, 0, 3, -1),
        (0, 1, 4, -3),
        (0, 2, 5, -3),
        (0, 3, 6, -1),
        (1, 0, 2, 3),
        (1, 1, 3, 9),
        (1, 2, 4, 9),
        (1, 3, 5, 3),
        (2, 0, 1, -3),
        (2, 1, 2, -9),
        (2, 2, 3, -9),
        (2, 3, 4, -3),
        (3, 0, 0, 1),
        (3, 1, 1, 3),
        (3, 2, 2, 3),
        (3, 3, 3, 1),
    ]
    return coefficient_map


def restricted_moments(
    coefficient_map: list[tuple[int, int, int, int]],
    parameter_count: int = 7,
    orders: tuple[int, ...] = SEXTIC_MOMENT_ORDERS,
    variable_prefix: str = "s",
) -> dict[int, str]:
    """Generate exact scalar moments on an irreducible summand."""

    polynomial = {(0, 0, (0,) * parameter_count): 1}
    answer: dict[int, str] = {}
    for order in range(1, max(orders) + 1):
        product: dict[tuple[int, int, tuple[int, ...]], int] = defaultdict(int)
        for (dual_v, coordinate_y, exponents), coefficient in polynomial.items():
            for i, j, parameter_index, scalar in coefficient_map:
                new_exponents = list(exponents)
                new_exponents[parameter_index] += 1
                product[
                    (dual_v + i, coordinate_y + j, tuple(new_exponents))
                ] += coefficient * scalar
        polynomial = dict(product)
        if order not in orders:
            continue

        terms: dict[tuple[int, ...], int] = defaultdict(int)
        for (dual_v, coordinate_y, exponents), coefficient in polynomial.items():
            if dual_v != coordinate_y:
                continue
            terms[exponents] += (
                coefficient
                * factorial(DEGREE * order - dual_v)
                * factorial(dual_v)
            )
        content = 0
        for coefficient in terms.values():
            content = gcd(content, abs(coefficient))
        assert content

        serialized = ""
        for exponents, coefficient in sorted(terms.items()):
            coefficient //= content
            factors = []
            for index, exponent in enumerate(exponents):
                if exponent == 1:
                    factors.append(f"{variable_prefix}{index}")
                elif exponent > 1:
                    factors.append(
                        f"{variable_prefix}{index}^{exponent}"
                    )
            monomial = "*".join(factors) or "1"
            if coefficient == 1:
                term = monomial
            elif coefficient == -1:
                term = f"-{monomial}"
            else:
                term = f"{coefficient}*{monomial}"
            serialized += (
                term
                if not serialized or term.startswith("-")
                else f"+{term}"
            )
        answer[order] = serialized
    return answer


def check_lower_summands(singular: str) -> None:
    """Close the pure binary-quartic and binary-quadratic summands."""

    W, V, Z, Y = sp.symbols("W V Z Y")

    def lowering(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            Y * sp.diff(polynomial, Z) - W * sp.diff(polynomial, V)
        )

    def coefficient_map(
        highest_weight: sp.Expr,
        parameter_count: int,
        expected: list[tuple[int, int, int, int]],
    ) -> list[tuple[int, int, int, int]]:
        parameters = sp.symbols(f"a0:{parameter_count}")
        basis = []
        polynomial = highest_weight
        for order in range(parameter_count):
            basis.append(sp.expand(polynomial))
            polynomial = lowering(polynomial) / (order + 1)
        general = sp.Poly(
            sp.expand(
                sum(
                    parameter * vector
                    for parameter, vector in zip(parameters, basis)
                )
            ),
            W,
            V,
            Z,
            Y,
        )
        answer = []
        for i, j in MONOMIALS:
            coefficient = sp.Poly(
                general.coeff_monomial(
                    W ** (DEGREE - i)
                    * V**i
                    * Z ** (DEGREE - j)
                    * Y**j
                ),
                *parameters,
            )
            if coefficient.is_zero:
                continue
            terms = coefficient.terms()
            assert len(terms) == 1
            exponents, scalar = terms[0]
            answer.append((i, j, exponents.index(1), int(scalar)))
        assert answer == expected
        return answer

    quartic_map = coefficient_map(
        W * V**2 * Z**3 + V**3 * Z**2 * Y,
        5,
        [
            (0, 0, 2, 1),
            (0, 1, 3, 2),
            (0, 2, 4, 1),
            (1, 0, 1, -2),
            (1, 1, 2, -3),
            (1, 3, 4, 1),
            (2, 0, 0, 1),
            (2, 2, 2, -3),
            (2, 3, 3, -2),
            (3, 1, 0, 1),
            (3, 2, 1, 2),
            (3, 3, 2, 1),
        ],
    )
    quadratic_map = coefficient_map(
        W**2 * V * Z**3
        + 2 * W * V**2 * Z**2 * Y
        + V**3 * Z * Y**2,
        3,
        [
            (0, 0, 1, -1),
            (0, 1, 2, -1),
            (1, 0, 0, 1),
            (1, 1, 1, -1),
            (1, 2, 2, -2),
            (2, 1, 0, 2),
            (2, 2, 1, 1),
            (2, 3, 2, -1),
            (3, 2, 0, 1),
            (3, 3, 1, 1),
        ],
    )
    quartic_moments = restricted_moments(
        quartic_map, 5, (2, 3), "t"
    )
    quadratic_moments = restricted_moments(
        quadratic_map, 3, (2,), "r"
    )

    X, T, q, u, v = sp.symbols("X T q u v")
    quartic = sp.Poly(
        sp.expand((X + q * T) ** 3 * (u * X + v * T)), X, T
    )
    quartic_coefficients = [
        sp.expand(
            quartic.coeff_monomial(X ** (4 - index) * T**index)
            / sp.binomial(4, index)
        )
        for index in range(5)
    ]
    quartic_relations = ",".join(
        f"t{index}-({singular_expression(coefficient)})"
        for index, coefficient in enumerate(quartic_coefficients)
    )
    quartic_generators = (
        "t0*t4-4*t1*t3+3*t2^2",
        (
            "-t0*t2*t4+t0*t3^2+t1^2*t4"
            "-2*t1*t2*t3+t2^3"
        ),
    )

    quadratic = sp.Poly(sp.expand(u * (X + q * T) ** 2), X, T)
    quadratic_coefficients = [
        sp.expand(
            quadratic.coeff_monomial(X ** (2 - index) * T**index)
            / sp.binomial(2, index)
        )
        for index in range(3)
    ]
    quadratic_relations = ",".join(
        f"r{index}-({singular_expression(coefficient)})"
        for index, coefficient in enumerate(quadratic_coefficients)
    )
    quadratic_generator = "-r0*r2+r1^2"

    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring quartic_elimination=0,(u,v,q,t0,t1,t2,t3,t4),(dp(3),dp(5));
option(redSB);
ideal QP={quartic_relations};
ideal QE=eliminate(QP,u*v*q);
ideal QG=std(QE);
ideal QW={",".join(quartic_generators)};
print(
  "QUARTIC_ELIM "
  +string(dim(QG))+" "+string(size(QG))+" "
  +string(size(reduce(QE,std(QW))))+" "
  +string(size(reduce(QW,QG)))
);
degree(QG);
ring quartic_moment=0,(t0,t1,t2,t3,t4),dp;
ideal QI={quartic_moments[2]},{quartic_moments[3]};
ideal QJ={",".join(quartic_generators)};
print(
  "QUARTIC_MOMENT "
  +string(size(reduce(QI,std(QJ))))+" "
  +string(size(reduce(QJ,std(QI))))
);
ring quadratic_elimination=0,(u,q,r0,r1,r2),(dp(2),dp(3));
ideal RP={quadratic_relations};
ideal RE=eliminate(RP,u*q);
ideal RG=std(RE);
ideal RW={quadratic_generator};
print(
  "QUADRATIC_ELIM "
  +string(dim(RG))+" "+string(size(RG))+" "
  +string(size(reduce(RE,std(RW))))+" "
  +string(size(reduce(RW,RG)))
);
degree(RG);
ring quadratic_moment=0,(r0,r1,r2),dp;
ideal RI={quadratic_moments[2]};
ideal RJ={quadratic_generator};
print(
  "QUADRATIC_MOMENT "
  +string(size(reduce(RI,std(RJ))))+" "
  +string(size(reduce(RJ,std(RI))))
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
    )
    quartic_marker = re.search(
        r"(?m)^QUARTIC_ELIM (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert quartic_marker is not None
    assert tuple(map(int, quartic_marker.groups())) == (6, 3, 0, 0), (
        quartic_marker.groups(),
        completed.stdout,
    )
    assert re.search(r"degree \(proj\.\)\s+= 6", completed.stdout)
    quartic_moment_marker = re.search(
        r"(?m)^QUARTIC_MOMENT (\d+) (\d+)$", completed.stdout
    )
    assert quartic_moment_marker is not None
    assert tuple(map(int, quartic_moment_marker.groups())) == (0, 0)
    quadratic_marker = re.search(
        r"(?m)^QUADRATIC_ELIM (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert quadratic_marker is not None
    assert tuple(map(int, quadratic_marker.groups())) == (4, 1, 0, 0), (
        quadratic_marker.groups(),
        completed.stdout,
    )
    degrees = re.findall(r"degree \(proj\.\)\s+= (\d+)", completed.stdout)
    assert degrees == ["6", "2"]
    quadratic_moment_marker = re.search(
        r"(?m)^QUADRATIC_MOMENT (\d+) (\d+)$", completed.stdout
    )
    assert quadratic_moment_marker is not None
    assert tuple(map(int, quadratic_moment_marker.groups())) == (0, 0)


def check_nonnull_quadratic_branches(singular: str) -> None:
    """Attack mixed branches after normalizing a non-null quadratic."""

    quartic_map = [
        (0, 0, 2, 1),
        (0, 1, 3, 2),
        (0, 2, 4, 1),
        (1, 0, 1, -2),
        (1, 1, 2, -3),
        (1, 3, 4, 1),
        (2, 0, 0, 1),
        (2, 2, 2, -3),
        (2, 3, 3, -2),
        (3, 1, 0, 1),
        (3, 2, 1, 2),
        (3, 3, 2, 1),
    ]
    sextic_map = [
        (0, 0, 3, -1),
        (0, 1, 4, -3),
        (0, 2, 5, -3),
        (0, 3, 6, -1),
        (1, 0, 2, 3),
        (1, 1, 3, 9),
        (1, 2, 4, 9),
        (1, 3, 5, 3),
        (2, 0, 1, -3),
        (2, 1, 2, -9),
        (2, 2, 3, -9),
        (2, 3, 4, -3),
        (3, 0, 0, 1),
        (3, 1, 1, 3),
        (3, 2, 2, 3),
        (3, 3, 3, 1),
    ]
    # In divided-power coordinates the non-null quadratic normal form is
    # 2*c*X*T.  Its biform embedding has these four diagonal terms.
    quartic_quadratic_map = quartic_map + [
        (0, 0, 5, -1),
        (1, 1, 5, -1),
        (2, 2, 5, 1),
        (3, 3, 5, 1),
    ]
    sextic_quadratic_map = sextic_map + [
        (0, 0, 7, -1),
        (1, 1, 7, -1),
        (2, 2, 7, 1),
        (3, 3, 7, 1),
    ]
    quartic_moments = restricted_moments(
        quartic_quadratic_map,
        6,
        (2, 3, 4, 5, 6),
        "t",
    )
    # Rename the last parameter t5 to the quadratic scale c.
    quartic_moments = {
        order: expression.replace("t5", "c")
        for order, expression in quartic_moments.items()
    }
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring quartic_quadratic=0,(t0,t1,t2,t3,t4,c),dp;
option(redSB);
ideal I={",".join(quartic_moments.values())};
ideal G=std(I);
print(
  "QUARTIC_QUADRATIC "
  +string(dim(G))+" "+string(size(G))+" "
  +string(size(reduce(c^5,G)))+" "
  +string(size(reduce(c^6,G)))
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
    )
    marker = re.search(
        r"(?m)^QUARTIC_QUADRATIC (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert marker is not None
    dimension, basis_size, power_five_size, power_six_size = map(
        int, marker.groups()
    )
    assert (dimension, basis_size, power_six_size) == (3, 10, 0)
    assert power_five_size > 0

    sextic_moments = restricted_moments(
        sextic_quadratic_map,
        8,
        (2, 4, 6, 8, 10, 12, 14),
        "s",
    )
    sextic_moments = {
        order: expression.replace("s7", "c")
        for order, expression in sextic_moments.items()
    }
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring sextic_quadratic=32003,(s0,s1,s2,s3,s4,s5,s6,c),dp;
option(redSB);
ideal I={",".join(sextic_moments.values())};
ideal G=std(I);
print(
  "SEXTIC_QUADRATIC_MODULAR "
  +string(dim(G))+" "+string(size(G))+" "
  +string(size(reduce(c^24,G)))+" "
  +string(size(reduce(c^25,G)))
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    marker = re.search(
        r"(?m)^SEXTIC_QUADRATIC_MODULAR "
        r"(\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert marker is not None
    dimension, basis_size, power_24_size, power_25_size = map(
        int, marker.groups()
    )
    assert (dimension, basis_size, power_25_size) == (4, 7576, 0)
    assert power_24_size > 0


def check_global_nullcone(singular: str) -> None:
    """Eliminate the full bidegree-(3,3) one-sided parametrization."""

    W, V, Z, Y, q = sp.symbols("W V Z Y q")
    parameters = sp.symbols("a10 a20 a21 a30 a31 a32")
    normal_form = sp.expand(
        sum(
            parameter
            * W ** (DEGREE - i)
            * (V - q * W) ** i
            * (Z + q * Y) ** (DEGREE - j)
            * Y**j
            for parameter, (i, j) in zip(
                parameters, POSITIVE_POSITIONS
            )
        )
    )
    polynomial = sp.Poly(normal_form, W, V, Z, Y)
    coefficients = [
        polynomial.coeff_monomial(
            W ** (DEGREE - i)
            * V**i
            * Z ** (DEGREE - j)
            * Y**j
        )
        for i, j in MONOMIALS
    ]
    variables = ",".join(
        [str(parameter) for parameter in parameters]
        + ["q"]
        + [f"x{index}" for index in range(16)]
    )
    relations = ",".join(
        f"x{index}-({singular_expression(coefficient)})"
        for index, coefficient in enumerate(coefficients)
    )
    elimination_monomial = "*".join(
        [str(parameter) for parameter in parameters] + ["q"]
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring global_ring=0,({variables}),(dp(7),dp(16));
option(redSB);
ideal parametrization={relations};
ideal eliminated=eliminate(parametrization,{elimination_monomial});
ideal G=std(eliminated);
print("GLOBAL33 "+string(dim(G))+" "+string(size(G)));
degree(G);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    marker = re.search(r"(?m)^GLOBAL33 (\d+) (\d+)$", completed.stdout)
    assert marker is not None
    # Seven eliminated parameters remain free in the 23-variable ring.
    assert tuple(map(int, marker.groups())) == (14, 148)
    assert re.search(r"degree \(proj\.\)\s+= 20", completed.stdout)


def check_sextic_nullcone(
    singular: str, moments: dict[int, str]
) -> None:
    """Prove radical equality on the binary-sextic summand."""

    X, T, q, u, v, w = sp.symbols("X T q u v w")
    sextic = sp.Poly(
        sp.expand((X + q * T) ** 4 * (u * X**2 + v * X * T + w * T**2)),
        X,
        T,
    )
    coefficients = [
        sp.expand(
            sextic.coeff_monomial(X ** (6 - index) * T**index)
            / sp.binomial(6, index)
        )
        for index in range(7)
    ]
    relations = ",".join(
        f"s{index}-({singular_expression(coefficient)})"
        for index, coefficient in enumerate(coefficients)
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring elimination_ring=0,(u,v,w,q,s0,s1,s2,s3,s4,s5,s6),(dp(4),dp(7));
option(redSB);
ideal parametrization={relations};
ideal eliminated=eliminate(parametrization,u*v*w*q);
ideal GE=std(eliminated);
ideal written={",".join(SEXTIC_NULLCONE_GENERATORS)};
ideal GW=std(written);
print(
  "SEXTIC_ELIM "
  +string(dim(GE))+" "+string(size(GE))+" "
  +string(dim(GW))+" "+string(size(GW))+" "
  +string(size(reduce(eliminated,GW)))+" "
  +string(size(reduce(written,GE)))
);
degree(GE);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
    )
    marker = re.search(
        r"(?m)^SEXTIC_ELIM (\d+) (\d+) (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert marker is not None
    assert tuple(map(int, marker.groups())) == (8, 10, 8, 10, 0, 0)
    assert re.search(r"degree \(proj\.\)\s+= 12", completed.stdout)

    ordered_moments = [moments[order] for order in SEXTIC_MOMENT_ORDERS]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring moment_ring=0,(s0,s1,s2,s3,s4,s5,s6),dp;
option(redSB);
ideal I={",".join(ordered_moments)};
ideal GI=std(I);
ideal J={",".join(SEXTIC_NULLCONE_GENERATORS)};
ideal GJ=std(J);
print(
  "SEXTIC_BASE "
  +string(dim(GI))+" "+string(size(GI))+" "
  +string(dim(GJ))+" "+string(size(GJ))+" "
  +string(size(reduce(I,GJ)))
);
int index;
int exponent;
for (index=1;index<=size(GJ);index++)
{{
  poly generator=GJ[index];
  poly remainder=reduce(generator,GI);
  exponent=1;
  while (exponent<5 && remainder!=0)
  {{
    remainder=reduce(remainder*generator,GI);
    exponent++;
  }}
  print(
    "SEXTIC_POWER "
    +string(index)+" "+string(exponent)+" "+string(size(remainder))
  );
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=30,
    )
    marker = re.search(
        r"(?m)^SEXTIC_BASE (\d+) (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    assert marker is not None
    assert tuple(map(int, marker.groups())) == (4, 65, 4, 10, 0)
    powers = tuple(
        int(exponent)
        for _, exponent, remainder_size in re.findall(
            r"(?m)^SEXTIC_POWER (\d+) (\d+) (\d+)$",
            completed.stdout,
        )
        if int(remainder_size) == 0
    )
    assert powers == SEXTIC_RADICAL_POWERS


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"

    coefficient_map = check_sextic_embedding()
    moments = restricted_moments(coefficient_map)
    full_moment_jacobian = check_full_moment_jacobian()
    corrected_moment_jacobian = check_corrected_moment_jacobian()
    moment_hsop_degree_obstruction = check_moment_hsop_degree_obstruction()
    global_quadratic_projection = check_global_quadratic_projection()
    diagonal_moment_slice = check_diagonal_moment_slice()
    normalized_branch_counts = normalized_branch_invariant_counts()
    check_global_nullcone(singular)
    check_sextic_nullcone(singular, moments)
    check_lower_summands(singular)
    check_nonnull_quadratic_branches(singular)

    artifact = {
        "format": "two-pair-sic-bidegree33-frontier-v4",
        "field": "characteristic zero",
        "full_bidegree_33": {
            "coefficient_count": 16,
            "one_sided_parameters": 7,
            "nullcone_dimension": 7,
            "nullcone_projective_degree": 20,
            "nullcone_groebner_basis_size": 148,
            "full_moment_jacobian": full_moment_jacobian,
            "corrected_moment_jacobian": corrected_moment_jacobian,
            "moment_hsop_degree_obstruction": (
                moment_hsop_degree_obstruction
            ),
            "global_quadratic_projection": global_quadratic_projection,
            "diagonal_moment_slice": diagonal_moment_slice,
            "normalized_nonnull_quadratic_branch": normalized_branch_counts,
            "status": (
                "exact nullcone elimination proved; equality with the full "
                "moment zero set remains open"
            ),
        },
        "binary_sextic_slice": {
            "dimension": 7,
            "decomposition_label": "Sym^6",
            "moment_orders": list(SEXTIC_MOMENT_ORDERS),
            "nullcone": "binary sextics of the form L^4*Q",
            "nullcone_dimension": 4,
            "nullcone_projective_degree": 12,
            "nullcone_generators": list(SEXTIC_NULLCONE_GENERATORS),
            "radical_power_certificate": list(SEXTIC_RADICAL_POWERS),
            "conclusion": (
                "the binary-sextic slice contains no SIC(2) counterexample"
            ),
        },
        "binary_quartic_slice": {
            "dimension": 5,
            "decomposition_label": "Sym^4",
            "moment_orders": [2, 3],
            "nullcone": "binary quartics of the form L^3*R",
            "nullcone_dimension": 3,
            "nullcone_projective_degree": 6,
            "conclusion": (
                "moments 2 and 3 generate the nullcone ideal"
            ),
        },
        "binary_quadratic_slice": {
            "dimension": 3,
            "decomposition_label": "Sym^2",
            "moment_orders": [2],
            "nullcone": "binary quadratics of the form L^2",
            "nullcone_dimension": 2,
            "nullcone_projective_degree": 2,
            "conclusion": "moment 2 generates the nullcone ideal",
        },
        "mixed_quartic_quadratic_branch": {
            "normal_form": "non-null quadratic 2*c*X*T",
            "moment_orders": [2, 3, 4, 5, 6],
            "field": "characteristic zero",
            "certificate": "c^6 lies in the homogeneous moment ideal",
            "minimal_power_in_checked_basis": 6,
            "conclusion": (
                "the Sym^4+Sym^2 branch with non-null quadratic "
                "contains no SIC(2) counterexample"
            ),
        },
        "mixed_sextic_quadratic_experiment": {
            "normal_form": "non-null quadratic 2*c*X*T",
            "moment_orders": [2, 4, 6, 8, 10, 12, 14],
            "field": "GF(32003)",
            "groebner_basis_size": 7576,
            "observed_power": 25,
            "certificate": (
                "c^25 reduces to zero while c^24 does not"
            ),
            "status": (
                "finite-field computation only; exact characteristic-zero "
                "membership remains open"
            ),
        },
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_two_pair_sic_bidegree33_sextic_slice.py"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS SIC(2) bidegree (3,3): exact seven-dimensional "
        "one-sided nullcone elimination"
    )
    print(
        "PASS SIC(2) bidegree (3,3): moments 1 through 13 have "
        "exact full Jacobian rank 13"
    )
    print(
        "PASS SIC(2) bidegree (3,3): Hilbert coefficient -2186 "
        "rules out hsop degrees 1 through 13"
    )
    print(
        "PASS SIC(2) bidegree (3,3): moments 1 through 12 and 14 "
        "have exact full Jacobian rank 13 and pass the Hilbert test"
    )
    print(
        "PASS SIC(2) bidegree (3,3): exact global Sym^2 projection "
        "and quadratic discriminant"
    )
    print(
        "PASS SIC(2) bidegree (3,3): moments 1 through 4 have only "
        "the origin on the torus-fixed diagonal slice"
    )
    print(
        "PASS SIC(2) binary-sextic slice: moments 2,4,6,10 have "
        "the L^4*Q nullcone radical"
    )
    print(
        "PASS SIC(2) lower pure slices: moments 2,3 cut out L^3*R "
        "and moment 2 cuts out L^2"
    )
    print(
        "PASS SIC(2) quartic+quadratic branch: c^6 lies in the "
        "characteristic-zero moment ideal"
    )
    print(
        "EVIDENCE SIC(2) sextic+quadratic branch: over GF(32003), "
        "c^25 lies in the moment ideal through order 14"
    )
    print(
        "OPEN SIC(2) bidegree (3,3): mixed Sym^6+Sym^4+Sym^2 "
        "components remain"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
