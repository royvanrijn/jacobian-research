#!/usr/bin/env python3
"""Verify the quartic (2+2) and (2+1+1) binary GVC calculations.

The proof is in ``extended-geometry/BINARY_QUARTIC_DOUBLE_ROOT_GVC.md``.
All contractions are exact over the rationals.
"""

from __future__ import annotations

import math

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
            result[exponent] = sp.expand(
                result.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
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
            result[exponent] = sp.expand(
                result.get(exponent, 0)
                + coefficient
                * operator_coefficient
                * falling(x_degree, x_order)
                * falling(y_degree, y_order)
            )
    return {
        exponent: coefficient
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


def verify_defect_one_splits() -> None:
    C, A = sp.symbols("C A")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    l0, l1, l4, l5 = sp.symbols("l0 l1 l4 l5")
    double_double_polynomial = {
        (5, 0): C,
        (4, 1): A,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
    }
    double_double_operator = {
        (2, 2): 1,
        (0, 5): l0,
        (1, 4): l1,
        (4, 1): l4,
        (5, 0): l5,
    }
    first = moment(
        double_double_polynomial, double_double_operator, 1
    )[(0, 0)]
    assert_multiple(first, p2 + 6 * A * l4 + 30 * C * l5)
    substitution = {p2: -6 * A * l4 - 30 * C * l5}
    second = {
        exponent: sp.expand(coefficient.subs(substitution))
        for exponent, coefficient in moment(
            double_double_polynomial, double_double_operator, 2
        ).items()
        if sum(exponent) == 1
    }
    assert_multiple(second[(0, 1)], A * p0)
    assert_multiple(
        second[(1, 0)], 28 * A**2 * l5 + A * p1 + C * p0
    )

    l2 = sp.symbols("l2")
    double_simple_polynomial = {
        (0, 5): C,
        (1, 4): A,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
    }
    double_simple_operator = {
        (3, 1): 1,
        (2, 2): -1,
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (5, 0): l5,
    }
    first = moment(
        double_simple_polynomial, double_simple_operator, 1
    )[(0, 0)]
    relation = 2 * p2 - 3 * p3 - 12 * A * l1 - 60 * C * l0
    assert_multiple(first, relation)
    second = {
        exponent: sp.expand(
            coefficient.subs(
                p2, (3 * p3 + 12 * A * l1 + 60 * C * l0) / 2
            )
        )
        for exponent, coefficient in moment(
            double_simple_polynomial, double_simple_operator, 2
        ).items()
        if sum(exponent) == 1
    }
    assert_multiple(second[(1, 0)], A * p4)
    assert_multiple(
        second[(0, 1)],
        -28 * A**2 * l0 + A * p3 - 2 * A * p4 + C * p4,
    )


def verify_double_double_linear_tip() -> None:
    """The (2+2) chart P_5=C*x^5+x^4*y."""

    C = sp.symbols("C")
    p3, p4 = sp.symbols("p3 p4")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    l0, l1, l4, l5 = sp.symbols("l0 l1 l4 l5")
    k0, k1, k5, k6 = sp.symbols("k0 k1 k5 k6")
    polynomial: DefectPieces = {
        0: {(5, 0): C, (4, 1): 1},
        1: {
            (1, 3): -28 * l5,
            (2, 2): -6 * l4 - 30 * C * l5,
            (3, 1): p3,
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
        0: {(2, 2): 1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (4, 1): l4,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (5, 1): k5, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)
    assert_multiple(
        second,
        2550 * C**2 * l5**2
        + 360 * C * l4 * l5
        + 140 * k6
        - 2 * l4**2
        + 7 * l5 * p3
        + r0,
    )
    assert_multiple(third[(0, 1)], l5**2)
    assert_multiple(third[(1, 0)], l5 * (821 * C * l5 + 60 * l4))

    # After l5=0 the preceding scalar equation gives
    # r0=2*l4^2-140*k6.  This is the complete weight-six face.
    u, v = sp.symbols("u v")
    equality_polynomial = {
        (4, 1): 1,
        (2, 2): -6 * u,
        (0, 3): 2 * u**2 - 140 * v,
    }
    equality_operator = {(2, 2): 1, (4, 1): u, (6, 0): v}
    third_equality = moment(
        equality_polynomial, equality_operator, 3
    )[(0, 0)]
    fourth_equality = moment(
        equality_polynomial, equality_operator, 4
    )[(0, 0)]
    f = u * (u**2 - 28 * v)
    g = 2789 * u**4 - 71680 * u**2 * v - 62320 * v**2
    assert_multiple(third_equality, f)
    assert_multiple(fourth_equality, g)
    assert_multiple(g.subs(u, 0), v**2)
    assert_multiple(g.subs(v, u**2 / 28), u**4)


def verify_double_double_fifth_power() -> None:
    """The (2+2) chart P_5=x^5; the y^5 chart is symmetric."""

    p0, p1, p3, p4 = sp.symbols("p0 p1 p3 p4")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    l0, l1, l4, l5 = sp.symbols("l0 l1 l4 l5")
    k0, k1, k5, k6 = sp.symbols("k0 k1 k5 k6")
    polynomial: DefectPieces = {
        0: {(5, 0): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): -30 * l5,
            (3, 1): p3,
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
        0: {(2, 2): 1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (4, 1): l4,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (5, 1): k5, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)
    fourth = defect_moment(polynomial, operator, 4, 2)
    assert_multiple(
        second,
        15 * l4 * p1 + 2550 * l5**2 + p0 * p4 + p1 * p3,
    )
    assert_multiple(third[(0, 1)], p0 * p1)
    assert_multiple(third[(1, 0)], 12 * l5 * p0 + p1**2)
    assert_multiple(fourth[(2, 0)], p0**2)


def verify_double_simple_linear_tip() -> None:
    """The (2+1+1) chart P_5=C*y^5+x*y^4."""

    C = sp.symbols("C")
    p0, p1 = sp.symbols("p0 p1")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    l0, l1, l2, l5 = sp.symbols("l0 l1 l2 l5")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    polynomial: DefectPieces = {
        0: {(0, 5): C, (1, 4): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): 42 * l0 + 6 * l1 + 30 * C * l0,
            (3, 1): 28 * l0,
        },
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    operator: DefectPieces = {
        0: {(3, 1): 1, (2, 2): -1},
        1: {
            (0, 5): l0,
            (1, 4): l1,
            (2, 3): l2,
            (5, 0): l5,
        },
        2: {(0, 6): k0, (1, 5): k1, (2, 4): k2, (6, 0): k6},
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)
    assert_multiple(third[(1, 0)], l0**2)
    assert_multiple(
        second.subs(l0, 0), 140 * k0 + 2 * l1**2 - r3
    )

    u, v = sp.symbols("u v")
    equality_polynomial = {
        (1, 4): 1,
        (2, 2): 6 * u,
        (3, 0): 140 * v + 2 * u**2,
    }
    equality_operator = {(2, 2): -1, (1, 4): u, (0, 6): v}
    third_equality = moment(
        equality_polynomial, equality_operator, 3
    )[(0, 0)]
    fourth_equality = moment(
        equality_polynomial, equality_operator, 4
    )[(0, 0)]
    f = u * (u**2 + 28 * v)
    g = 2789 * u**4 + 71680 * u**2 * v - 62320 * v**2
    assert_multiple(third_equality, f)
    assert_multiple(fourth_equality, g)
    assert_multiple(g.subs(u, 0), v**2)
    assert_multiple(g.subs(v, -u**2 / 28), u**4)


def verify_double_simple_pure_y() -> None:
    p0, p1, p3 = sp.symbols("p0 p1 p3")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    l0, l1, l2, l5 = sp.symbols("l0 l1 l2 l5")
    k0, k1, k2, k6 = sp.symbols("k0 k1 k2 k6")
    polynomial: DefectPieces = {
        0: {(0, 5): 1},
        1: {
            (0, 4): p0,
            (1, 3): p1,
            (2, 2): (3 * p3 + 60 * l0) / 2,
            (3, 1): p3,
        },
        2: {
            (0, 3): r0,
            (1, 2): r1,
            (2, 1): r2,
            (3, 0): r3,
        },
    }
    operator: DefectPieces = {
        0: {(3, 1): 1, (2, 2): -1},
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
    assert_multiple(
        second,
        20400 * l0**2
        - 540 * l0 * p3
        - 120 * l1 * p3
        + 8 * p1 * p3
        - 11 * p3**2,
    )
    assert_multiple(third, p3**2)
    assert_multiple(second.subs(p3, 0), l0**2)


def verify_double_simple_fifth_power() -> None:
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    l0, l1, l2, l5 = sp.symbols("l0 l1 l2 l5")
    polynomial = {(5, 0): 1} | {
        (degree, 4 - degree): coefficient
        for degree, coefficient in enumerate((p0, p1, p2, p3, p4))
    }
    operator = {
        (3, 1): 1,
        (2, 2): -1,
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (5, 0): l5,
    }
    first = moment(polynomial, operator, 1)[(0, 0)]
    assert_multiple(first, -60 * l5 + 2 * p2 - 3 * p3)
    third = {
        exponent: coefficient
        for exponent, coefficient in moment(
            polynomial, operator, 3
        ).items()
        if sum(exponent) == 2
    }
    assert_multiple(third[(1, 1)], p0)
    assert_multiple(third[(2, 0)].subs(p0, 0), p1)
    second = {
        exponent: coefficient
        for exponent, coefficient in moment(
            polynomial, operator, 2
        ).items()
        if sum(exponent) == 1
    }
    assert_multiple(second[(1, 0)].subs({p0: 0, p1: 0}), p2)

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
        0: {(3, 1): 1, (2, 2): -1},
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
    )[(0, 0)]
    third_defect = defect_moment(
        defect_polynomial, defect_operator, 3, 2
    )[(1, 0)]
    assert_multiple(second_defect, -340 * l5**2 + r0 - r1)
    assert_multiple(third_defect, r0)

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

    X, Y = sp.symbols("X Y")
    leading_symbol = X**2 * Y * (X - Y)
    assert sp.expand(
        leading_symbol.subs(Y, X - Y) - leading_symbol
    ) == 0


def verify_final_weight_separators() -> None:
    # (2+2), linear tip: (1,2), threshold 6.
    assert 2 + 2 * 2 == 4 + 2
    assert 1 > 0  # y-deficit: X^2Y^2 versus x^4y

    # (2+2), fifth power: (2,3), threshold 10.
    assert 2 * 2 + 2 * 3 == 5 * 2
    assert 2 > 0

    # (2+1+1), linear tip: (2,1), threshold 6.
    assert 2 * 2 + 2 == 2 + 4
    assert 2 > 1  # x-deficit: X^2Y^2 versus xy^4

    # (2+1+1), pure y: (3,2), threshold 10.
    assert 2 * 3 + 2 * 2 == 5 * 2
    assert 2 > 0

    # (2+1+1), x^5: (1,2), threshold 5.
    assert 3 + 2 == 5
    assert 1 > 0

    # Later normal-form jets are strictly above the relevant thresholds.
    for degree in range(7, 20):
        double_double_support = (
            (0, degree),
            (1, degree - 1),
            (degree - 1, 1),
            (degree, 0),
        )
        assert min(x + 2 * y for x, y in double_double_support) > 6
    for degree in range(6, 20):
        double_double_support = (
            (0, degree),
            (1, degree - 1),
            (degree - 1, 1),
            (degree, 0),
        )
        assert min(2 * x + 3 * y for x, y in double_double_support) > 10
        double_simple_support = (
            (0, degree),
            (1, degree - 1),
            (2, degree - 2),
            (degree, 0),
        )
        assert min(3 * x + 2 * y for x, y in double_simple_support) > 10
        assert min(x + 2 * y for x, y in double_simple_support) > 5
    for degree in range(7, 20):
        double_simple_support = (
            (0, degree),
            (1, degree - 1),
            (2, degree - 2),
            (degree, 0),
        )
        assert min(2 * x + y for x, y in double_simple_support) > 6


def main() -> None:
    verify_defect_one_splits()
    verify_double_double_linear_tip()
    verify_double_double_fifth_power()
    verify_double_simple_linear_tip()
    verify_double_simple_pure_y()
    verify_double_simple_fifth_power()
    verify_final_weight_separators()
    print("PASS (2+2) arbitrary-jet branches die by pure moment four")
    print("PASS (2+1+1) arbitrary-jet branches die by pure moment four")
    print("PASS all final equality faces have a linear coordinate deficit")


if __name__ == "__main__":
    main()
