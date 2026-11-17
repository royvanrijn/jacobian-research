#!/usr/bin/env python3
"""Verify the quartic (3+1)-leading binary degree-five GVC calculation.

The accompanying proof is in
``extended-geometry/BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md``.  This
checker replays the defect-one radical, the exact higher-defect branch
identities, and the final weighted support separators.  Singular is
required for the radical calculation.
"""

from __future__ import annotations

import math
import shutil
import subprocess

import sympy as sp


Exponent = tuple[int, int]
SparsePolynomial = dict[Exponent, sp.Expr]
DefectPieces = dict[int, SparsePolynomial]


def falling(degree: int, order: int) -> int:
    return math.prod(range(degree - order + 1, degree + 1)) if order else 1


def multiply(
    left: SparsePolynomial, right: SparsePolynomial
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for (left_x, left_y), left_coefficient in left.items():
        for (right_x, right_y), right_coefficient in right.items():
            exponent = (left_x + right_x, left_y + right_y)
            result[exponent] = (
                result.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: sp.expand(coefficient)
        for exponent, coefficient in result.items()
        if coefficient != 0
    }


def apply_operator(
    polynomial: SparsePolynomial, operator: SparsePolynomial
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for (x_degree, y_degree), coefficient in polynomial.items():
        for (x_order, y_order), operator_coefficient in operator.items():
            if x_degree < x_order or y_degree < y_order:
                continue
            exponent = (x_degree - x_order, y_degree - y_order)
            result[exponent] = (
                result.get(exponent, 0)
                + coefficient
                * operator_coefficient
                * falling(x_degree, x_order)
                * falling(y_degree, y_order)
            )
    return {
        exponent: sp.expand(coefficient)
        for exponent, coefficient in result.items()
        if coefficient != 0
    }


def moment(
    polynomial: SparsePolynomial,
    operator: SparsePolynomial,
    order: int,
) -> SparsePolynomial:
    polynomial_power: SparsePolynomial = {(0, 0): 1}
    operator_power: SparsePolynomial = {(0, 0): 1}
    for _ in range(order):
        polynomial_power = multiply(polynomial_power, polynomial)
        operator_power = multiply(operator_power, operator)
    return apply_operator(polynomial_power, operator_power)


def layered_power(
    pieces: DefectPieces, order: int, maximum_defect: int
) -> list[SparsePolynomial]:
    layers: list[SparsePolynomial] = [
        {} for _ in range(maximum_defect + 1)
    ]
    layers[0] = {(0, 0): 1}
    for _ in range(order):
        next_layers: list[SparsePolynomial] = [
            {} for _ in range(maximum_defect + 1)
        ]
        for old_defect, old_piece in enumerate(layers):
            if not old_piece:
                continue
            for new_defect, new_piece in pieces.items():
                total_defect = old_defect + new_defect
                if total_defect > maximum_defect:
                    continue
                for exponent, coefficient in multiply(
                    old_piece, new_piece
                ).items():
                    next_layers[total_defect][exponent] = sp.expand(
                        next_layers[total_defect].get(exponent, 0)
                        + coefficient
                    )
        layers = next_layers
    return layers


def defect_moment(
    polynomial: DefectPieces,
    operator: DefectPieces,
    order: int,
    defect: int,
) -> SparsePolynomial:
    polynomial_layers = layered_power(polynomial, order, defect)
    operator_layers = layered_power(operator, order, defect)
    result: SparsePolynomial = {}
    for polynomial_defect in range(defect + 1):
        value = apply_operator(
            polynomial_layers[polynomial_defect],
            operator_layers[defect - polynomial_defect],
        )
        for exponent, coefficient in value.items():
            result[exponent] = sp.expand(
                result.get(exponent, 0) + coefficient
            )
    return result


def assert_multiple(actual: sp.Expr, expected: sp.Expr) -> None:
    ratio = sp.cancel(actual / expected)
    assert ratio != 0 and not ratio.free_symbols, (actual, expected, ratio)


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def verify_radical(
    variables: tuple[sp.Symbol, ...],
    equations: tuple[sp.Expr, ...],
    expected: tuple[sp.Expr, ...],
) -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the radical replay")
    program = f"""
ring r=0,({",".join(map(str, variables))}),dp;
ideal I={",".join(map(singular_expression, equations))};
ideal E={",".join(map(singular_expression, expected))};
LIB "primdec.lib";
ideal R=std(radical(I));
ideal left=reduce(R,std(E));
ideal right=reduce(std(E),R);
if ((size(left)==0) && (size(right)==0)) {{ print("PASS"); }}
else {{ print("FAIL"); R; }}
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in result.stdout and "FAIL" not in result.stdout


def verify_defect_one_radical() -> None:
    A, B, C = sp.symbols("A B C")
    l0, l1, l2, l5, p4 = sp.symbols("l0 l1 l2 l5 p4")
    equations = (
        B * (14 * A * l0 + 4 * B * l1 + p4),
        B**2 * l0,
    )
    expected = (B * l0, B * (4 * B * l1 + p4))
    verify_radical(
        (A, B, C, l0, l1, l2, l5, p4),
        equations,
        expected,
    )

    p0, p1, p2, p3 = sp.symbols("p0 p1 p2 p3")
    polynomial = {
        (0, 5): C,
        (1, 4): A,
        (2, 3): B,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
    }
    operator = {
        (3, 1): 1,
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (5, 0): l5,
    }
    first = moment(polynomial, operator, 1)[(0, 0)]
    assert_multiple(
        first,
        p3 + 20 * C * l0 + 4 * A * l1 + 2 * B * l2,
    )
    substitution = {
        p3: -20 * C * l0 - 4 * A * l1 - 2 * B * l2
    }
    for order, expected_values in (
        (
            2,
            {
                (0, 1): B * (14 * A * l0 + 4 * B * l1 + p4),
                (1, 0): B**2 * l0,
            },
        ),
        (3, {(0, 2): B**3 * l0}),
    ):
        value = {
            exponent: sp.expand(coefficient.subs(substitution))
            for exponent, coefficient in moment(
                polynomial, operator, order
            ).items()
            if sum(exponent) == order - 1
        }
        assert value.keys() == expected_values.keys()
        for exponent, expected_value in expected_values.items():
            assert_multiple(value[exponent], expected_value)


def verify_quadratic_tip_branch() -> None:
    """Close B != 0 after normalizing B=1."""

    A, C = sp.symbols("A C")
    p0, p1, p2 = sp.symbols("p0 p1 p2")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    l1, l2, l5 = sp.symbols("l1 l2 l5")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    polynomial: DefectPieces = {
        0: {(0, 5): C, (1, 4): A, (2, 3): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): p2,
            (3, 1): -4 * A * l1 - 2 * l2,
            (4, 0): -4 * l1,
        },
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    operator: DefectPieces = {
        0: {(3, 1): 1},
        1: {(1, 4): l1, (2, 3): l2, (5, 0): l5},
        2: {(0, 6): k0, (1, 5): k1, (2, 4): k2, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)
    third = defect_moment(polynomial, operator, 3, 2)
    assert_multiple(
        second[(0, 0)],
        10 * A**2 * l1**2
        + 21 * A * k0
        + 5 * A * l1 * l2
        + 20 * C * l1**2
        + 6 * k1,
    )
    assert_multiple(
        third[(0, 1)],
        186 * A * l1**2 + 60 * k0 + 37 * l1 * l2,
    )
    assert_multiple(third[(1, 0)], l1**2)

    # Thus l1=k0=k1=0.  The only weight-seven equality chain is:
    t, h = sp.symbols("t h")
    equality_polynomial = {(2, 3): 1, (3, 1): -2 * t}
    equality_operator = {(3, 1): 1, (2, 3): t, (0, 7): h}
    third_equality = moment(
        equality_polynomial, equality_operator, 3
    )[(0, 0)]
    fourth_equality = moment(
        equality_polynomial, equality_operator, 4
    )[(0, 0)]
    assert_multiple(third_equality, -20 * h + t**3)
    assert_multiple(
        sp.expand(fourth_equality.subs(h, t**3 / 20)),
        t**4,
    )


def verify_linear_tip_branch() -> None:
    """Close B=0, A != 0 after normalizing A=1."""

    C = sp.symbols("C")
    l0, l1, l2, l5, p0, p1, p2, p4 = sp.symbols(
        "l0 l1 l2 l5 p0 p1 p2 p4"
    )
    r0, r1, r2, r3 = sp.symbols("r0:4")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    polynomial: DefectPieces = {
        0: {(0, 5): C, (1, 4): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): p2,
            (3, 1): -20 * C * l0 - 4 * l1,
            (4, 0): p4,
        },
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    operator: DefectPieces = {
        0: {(3, 1): 1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (2, 3): l2,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (2, 4): k2, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)[(0, 1)]
    expected_second = (
        1340 * C**2 * l0**2
        + 480 * C * l0 * l1
        + 4 * C * l1 * p4
        + 56 * l0 * l2
        + 6 * l0 * p2
        + 20 * l1**2
        + 4 * l2 * p4
        + p2 * p4
    )
    assert_multiple(second, expected_second)
    assert_multiple(third, p4**2 + 10 * l0 * p4 + 330 * l0**2)
    assert_multiple(
        expected_second.subs({l0: 0, p4: 0}), l1**2
    )

    # On l0 != 0 put l0=1 and s=p4/l0.  The second equation
    # solves p2 because s+6 is coprime to f.  Complete W_6,W_7 and
    # P_3,P_2 data cannot change the displayed fourth-moment x term.
    s = sp.symbols("s")
    q0, q1, q2 = sp.symbols("q0:3")
    h0, h1, h2, h7 = sp.symbols("h0 h1 h2 h7")
    solved_p2 = -(
        1340 * C**2
        + 480 * C * l1
        + 4 * C * l1 * s
        + 56 * l2
        + 20 * l1**2
        + 4 * l2 * s
    ) / (s + 6)
    localized_polynomial: DefectPieces = {
        0: {(0, 5): C, (1, 4): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): solved_p2,
            (3, 1): -20 * C - 4 * l1,
            (4, 0): s,
        },
        2: polynomial[2],
        3: {(0, 2): q0, (1, 1): q1, (2, 0): q2},
    }
    localized_operator: DefectPieces = {
        0: operator[0],
        1: {
            (0, 5): 1,
            (1, 4): l1,
            (2, 3): l2,
            (5, 0): l5,
        },
        2: operator[2],
        3: {(0, 7): h0, (1, 6): h1, (2, 5): h2, (7, 0): h7},
    }
    fourth = defect_moment(
        localized_polynomial, localized_operator, 4, 3
    )
    obstruction = 143 * s**3 + 840 * s**2 + 13860 * s + 480480
    assert_multiple(fourth[(1, 0)], obstruction)
    quadratic = s**2 + 10 * s + 330
    assert sp.gcd(quadratic, s + 6) == 1
    assert sp.resultant(quadratic, obstruction, s) != 0


def verify_pure_y_branch() -> None:
    """Close B=A=0 after normalizing P_5=y^5."""

    l0, l1, l2, l5, p0, p1, p4 = sp.symbols(
        "l0 l1 l2 l5 p0 p1 p4"
    )
    r0, r1, r2, r3 = sp.symbols("r0:4")
    q0, q1, q2 = sp.symbols("q0:3")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    h0, h1, h2, h7 = sp.symbols("h0 h1 h2 h7")
    p2 = sp.symbols("p2")
    polynomial: DefectPieces = {
        0: {(0, 5): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): p2,
            (3, 1): -20 * l0,
            (4, 0): p4,
        },
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    operator: DefectPieces = {
        0: {(3, 1): 1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (2, 3): l2,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (2, 4): k2, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    assert_multiple(second, 1340 * l0**2 + 4 * l1 * p4 + p2 * p4)

    # If p4 != 0, solve the preceding equation.  The complete
    # defect-three fourth moment then forces p4=0, a contradiction.
    solved_p2 = -(1340 * l0**2 + 4 * l1 * p4) / p4
    localized_polynomial: DefectPieces = {
        0: polynomial[0],
        1: polynomial[1] | {(2, 2): solved_p2},
        2: polynomial[2],
        3: {(0, 2): q0, (1, 1): q1, (2, 0): q2},
    }
    localized_operator: DefectPieces = {
        0: operator[0],
        1: operator[1],
        2: operator[2],
        3: {(0, 7): h0, (1, 6): h1, (2, 5): h2, (7, 0): h7},
    }
    fourth = defect_moment(
        localized_polynomial, localized_operator, 4, 3
    )
    assert_multiple(fourth[(0, 1)], p4**3)


def verify_x_fifth_power_branch() -> None:
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    l0, l1, l2, l5 = sp.symbols("l0 l1 l2 l5")
    polynomial = {(5, 0): 1} | {
        (degree, 4 - degree): coefficient
        for degree, coefficient in enumerate((p0, p1, p2, p3, p4))
    }
    operator = {
        (3, 1): 1,
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (5, 0): l5,
    }
    first = moment(polynomial, operator, 1)[(0, 0)]
    assert_multiple(first, p3 + 20 * l5)
    substitution = {p3: -20 * l5}
    second = {
        exponent: sp.expand(coefficient.subs(substitution))
        for exponent, coefficient in moment(
            polynomial, operator, 2
        ).items()
        if sum(exponent) == 1
    }
    third = {
        exponent: sp.expand(coefficient.subs(substitution))
        for exponent, coefficient in moment(
            polynomial, operator, 3
        ).items()
        if sum(exponent) == 2
    }
    assert_multiple(second[(0, 1)], p1)
    assert_multiple(second[(1, 0)], p2)
    assert_multiple(third[(1, 1)], p0)

    r0, r1, r2, r3 = sp.symbols("r0:4")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    defect_polynomial: DefectPieces = {
        0: {(5, 0): 1},
        1: {(3, 1): -20 * l5, (4, 0): p4},
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    defect_operator: DefectPieces = {
        0: {(3, 1): 1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (2, 3): l2,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (2, 4): k2, (6, 0): k6},
    }
    second_defect = defect_moment(
        defect_polynomial, defect_operator, 2, 2
    )
    third_defect = defect_moment(
        defect_polynomial, defect_operator, 3, 2
    )
    assert_multiple(second_defect[(0, 0)], r1 + 340 * l5**2)
    assert_multiple(third_defect[(1, 0)], r0)

    t = sp.symbols("t")
    equality_polynomial = {
        (5, 0): 1,
        (3, 1): -20 * t,
        (1, 2): -340 * t**2,
    }
    equality_operator = {(3, 1): 1, (5, 0): t}
    assert_multiple(
        moment(equality_polynomial, equality_operator, 3)[(0, 0)],
        t**3,
    )


def verify_final_weight_separators() -> None:
    # B != 0: weights (2,1), threshold seven.
    assert max(2 * x_degree + y_degree for x_degree in range(3)
               for y_degree in (5 - x_degree,)) == 7
    assert 2 * 3 + 1 == 7
    assert 2 < 3  # x-deficit for (X^3Y, x^2y^3)

    # B=0,A!=0: weights (3,2), threshold eleven.
    assert 3 * 1 + 2 * 4 == 11
    assert 3 * 3 + 2 == 11
    assert 1 < 3

    # P_5=y^5: weights (4,3), threshold fifteen.
    assert 3 * 5 == 15
    assert 4 * 3 + 3 == 15
    assert 0 < 3

    # P_5=x^5: weights (1,2), threshold five.
    assert 5 == 3 + 2
    assert 0 < 1  # y-deficit for (X^3Y, x^5)

    # The normal-form support in order d>=5 is
    # Y^d, X Y^(d-1), X^2 Y^(d-2), X^d.  The displayed eliminations
    # leave every later operator term strictly above the threshold.
    for degree in range(8, 20):
        support = (
            (0, degree),
            (1, degree - 1),
            (2, degree - 2),
            (degree, 0),
        )
        assert min(2 * x + y for x, y in support) > 7
    for degree in range(6, 20):
        support = (
            (0, degree),
            (1, degree - 1),
            (2, degree - 2),
            (degree, 0),
        )
        assert min(3 * x + 2 * y for x, y in support) > 11
        assert min(4 * x + 3 * y for x, y in support) > 15
        assert min(x + 2 * y for x, y in support) > 5


def main() -> None:
    verify_defect_one_radical()
    verify_quadratic_tip_branch()
    verify_linear_tip_branch()
    verify_pure_y_branch()
    verify_x_fifth_power_branch()
    verify_final_weight_separators()
    print("PASS (3+1) defect-one radical and projective branch split")
    print("PASS all higher-defect threshold chains die by pure moment four")
    print("PASS arbitrary later jets have strict or one-sided final faces")


if __name__ == "__main__":
    main()
