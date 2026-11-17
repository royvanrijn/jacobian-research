#!/usr/bin/env python3
"""Verify the r=5, deg(P)=6 binary GVC calculation for Lambda_5=X^5.

The proof is in
``extended-geometry/BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md``.  Singular is
required for the defect-one radical; all contractions are exact over Q.
"""

from __future__ import annotations

import sympy as sp

from verify_binary_quartic_triple_simple_root_gvc import (
    apply_operator,
    assert_multiple,
    defect_moment,
    moment,
    multiply,
    verify_radical,
)


def verify_defect_one_radical() -> None:
    A, B, C, D = sp.symbols("A B C D")
    l0, l1, l2, l3, l4 = sp.symbols("l0:5")
    variables = (A, B, C, D, l0, l1, l2, l3, l4)
    polynomial = {
        (0, 6): 1,
        (1, 5): A,
        (2, 4): B,
        (3, 3): C,
        (4, 2): D,
    }
    sextic_jet = {
        (0, 6): l0,
        (1, 5): l1,
        (2, 4): l2,
        (3, 3): l3,
        (4, 2): l4,
    }
    equations: list[sp.Expr] = []
    power = multiply(polynomial, polynomial)
    for order in range(2, 4):
        if order > 2:
            power = multiply(power, polynomial)
        value = apply_operator(
            apply_operator(power, sextic_jet),
            {(5 * (order - 1), 0): 1},
        )
        equations.extend(value.values())
    expected = (
        D * l0,
        D * l1,
        D * l2,
        D * l3,
        C * l0,
        C * l1,
    )
    verify_radical(variables, tuple(equations), expected)

    p5 = sp.symbols("p5")
    first = apply_operator(
        polynomial | {(5, 0): p5},
        {(5, 0): 1} | sextic_jet,
    )[(0, 0)]
    assert_multiple(
        first,
        p5
        + 6 * l0
        + A * l1
        + sp.Rational(2, 5) * B * l2
        + sp.Rational(3, 10) * C * l3
        + sp.Rational(2, 5) * D * l4,
    )


def verify_x4y2_terminal_branch() -> None:
    """Close D != 0 after D=1 and l0=l1=l2=l3=0."""

    A, B, C = sp.symbols("A B C")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    t = sp.symbols("t")
    h0, h1, h2, h3, h4 = sp.symbols("h0:5")
    polynomial = {
        0: {
            (0, 6): 1,
            (1, 5): A,
            (2, 4): B,
            (3, 3): C,
            (4, 2): 1,
        },
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): p4,
            (5, 0): -sp.Rational(2, 5) * t,
        },
        2: {
            (0, 4): q0,
            (1, 3): q1,
            (2, 2): q2,
            (3, 1): q3,
            (4, 0): q4,
        },
    }
    operator = {
        0: {(5, 0): 1},
        1: {(4, 2): t},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (3, 4): h3,
            (4, 3): h4,
        },
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)
    fourth = defect_moment(polynomial, operator, 4, 2)
    assert_multiple(fourth[(1, 1)], h0)
    assert_multiple(third[(1, 0)].subs(h0, 0), h1)
    assert_multiple(
        third[(0, 1)].subs({h0: 0, h1: 0}), h2
    )
    assert_multiple(
        second.subs({h0: 0, h1: 0, h2: 0}), 5 * h3 + t**2
    )

    r0, r1, r2, r3 = sp.symbols("r0:4")
    k0, k1, k2, k3, k4 = sp.symbols("k0:5")
    polynomial[3] = {
        (0, 3): r0,
        (1, 2): r1,
        (2, 1): r2,
        (3, 0): r3,
    }
    operator[2] = {
        (3, 4): -sp.Rational(1, 5) * t**2,
        (4, 3): h4,
    }
    operator[3] = {
        (0, 8): k0,
        (1, 7): k1,
        (2, 6): k2,
        (3, 5): k3,
        (4, 4): k4,
    }
    third_defect = defect_moment(polynomial, operator, 3, 3)[(0, 0)]
    fourth_defect = defect_moment(polynomial, operator, 4, 3)
    assert_multiple(fourth_defect[(1, 0)], k0)
    assert_multiple(fourth_defect[(0, 1)].subs(k0, 0), k1)
    assert_multiple(
        third_defect.subs({k0: 0, k1: 0}),
        12375 * k2 - 3971 * t**3,
    )

    s0, s1, s2 = sp.symbols("s0:3")
    z0, z1, z2, z3, z4 = sp.symbols("z0:5")
    polynomial[4] = {(0, 2): s0, (1, 1): s1, (2, 0): s2}
    operator[3] = {
        (2, 6): sp.Rational(3971, 12375) * t**3,
        (3, 5): k3,
        (4, 4): k4,
    }
    operator[4] = {
        (0, 9): z0,
        (1, 8): z1,
        (2, 7): z2,
        (3, 6): z3,
        (4, 5): z4,
    }
    fourth_defect = defect_moment(polynomial, operator, 4, 4)[(0, 0)]
    fifth_defect = defect_moment(polynomial, operator, 5, 4)[(0, 1)]
    assert_multiple(fifth_defect, z0)
    assert_multiple(
        fourth_defect.subs(z0, 0), 78750 * z1 + 38939 * t**4
    )

    u = sp.symbols("u")
    equality_polynomial = {
        (4, 2): 1,
        (5, 0): -sp.Rational(2, 5) * t,
    }
    equality_operator = {
        (5, 0): 1,
        (4, 2): t,
        (3, 4): -sp.Rational(1, 5) * t**2,
        (2, 6): sp.Rational(3971, 12375) * t**3,
        (1, 8): -sp.Rational(38939, 78750) * t**4,
        (0, 10): u,
    }
    fifth = moment(equality_polynomial, equality_operator, 5)[(0, 0)]
    sixth = moment(equality_polynomial, equality_operator, 6)[(0, 0)]
    assert_multiple(fifth, 4872527 * t**5 - 5906250 * u)
    solved_u = sp.Rational(4872527, 5906250) * t**5
    assert_multiple(sp.expand(sixth.subs(u, solved_u)), t**6)


def verify_x3y3_branch() -> None:
    """Close D=0, C != 0 after C=1 and l0=l1=0."""

    A, B = sp.symbols("A B")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    l2, l3, l4 = sp.symbols("l2 l3 l4")
    h0, h1, h2, h3, h4 = sp.symbols("h0:5")
    polynomial = {
        0: {(0, 6): 1, (1, 5): A, (2, 4): B, (3, 3): 1},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): p4,
            (5, 0): -sp.Rational(2, 5) * B * l2
            - sp.Rational(3, 10) * l3,
        },
        2: {
            (0, 4): q0,
            (1, 3): q1,
            (2, 2): q2,
            (3, 1): q3,
            (4, 0): q4,
        },
    }
    operator = {
        0: {(5, 0): 1},
        1: {(2, 4): l2, (3, 3): l3, (4, 2): l4},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (3, 4): h3,
            (4, 3): h4,
        },
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)[(0, 1)]
    assert_multiple(third, l2**2)
    assert_multiple(
        second.subs(l2, 0), 1400 * B * h0 + 600 * h1 + 321 * l3**2
    )

    r0, r1, r2, r3 = sp.symbols("r0:4")
    k0, k1, k2, k3, k4 = sp.symbols("k0:5")
    solved_h1 = (
        -sp.Rational(7, 3) * B * h0
        - sp.Rational(107, 200) * l3**2
    )
    polynomial[1][(5, 0)] = -sp.Rational(3, 10) * l3
    polynomial[3] = {
        (0, 3): r0,
        (1, 2): r1,
        (2, 1): r2,
        (3, 0): r3,
    }
    operator[1] = {(3, 3): l3, (4, 2): l4}
    operator[2][(1, 6)] = solved_h1
    operator[3] = {
        (0, 8): k0,
        (1, 7): k1,
        (2, 6): k2,
        (3, 5): k3,
        (4, 4): k4,
    }
    third_defect = defect_moment(polynomial, operator, 3, 3)[(0, 0)]
    assert_multiple(
        third_defect.subs(h0, 0), l3**3
    )

    # On h0 != 0, moment three solves p4.  Complete defect-four
    # data then give a positive-degree fifth-moment obstruction h0^2.
    solved_p4 = (
        565 * l3**3
        - 1484 * B * h0 * l3
        - 1344 * h0 * l4
    ) / (280 * h0)
    polynomial[1][(4, 1)] = solved_p4
    s0, s1, s2 = sp.symbols("s0:3")
    z0, z1, z2, z3, z4 = sp.symbols("z0:5")
    polynomial[4] = {(0, 2): s0, (1, 1): s1, (2, 0): s2}
    operator[4] = {
        (0, 9): z0,
        (1, 8): z1,
        (2, 7): z2,
        (3, 6): z3,
        (4, 5): z4,
    }
    fifth_defect = defect_moment(polynomial, operator, 5, 4)[(0, 1)]
    assert_multiple(fifth_defect, h0**2)


def verify_x2y4_branch() -> None:
    """Close C=D=0, B != 0 after B=1."""

    A = sp.symbols("A")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    l0, l1, l2, l3, l4 = sp.symbols("l0:5")
    h0, h1, h2, h3, h4 = sp.symbols("h0:5")
    polynomial = {
        0: {(0, 6): 1, (1, 5): A, (2, 4): 1},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): p4,
            (5, 0): -6 * l0 - A * l1 - sp.Rational(2, 5) * l2,
        },
        2: {
            (0, 4): q0,
            (1, 3): q1,
            (2, 2): q2,
            (3, 1): q3,
            (4, 0): q4,
        },
    }
    operator = {
        0: {(5, 0): 1},
        1: {
            (0, 6): l0,
            (1, 5): l1,
            (2, 4): l2,
            (3, 3): l3,
            (4, 2): l4,
        },
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (3, 4): h3,
            (4, 3): h4,
        },
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    third = defect_moment(polynomial, operator, 3, 2)
    assert_multiple(third[(1, 0)], l0**2)
    reduced_second = sp.expand(second.subs(l0, 0))
    assert_multiple(reduced_second.subs(l1, 0), l2**2)

    # On l1 != 0, solve the reduced second moment for p4.
    solved_p4 = -(
        1525 * A**2 * l1**2
        + 1590 * A * l1 * l2
        + 2100 * l1**2
        + 280 * l1 * l3
        + 196 * l2**2
    ) / (50 * l1)
    polynomial[1][(4, 1)] = solved_p4
    polynomial[1][(5, 0)] = -A * l1 - sp.Rational(2, 5) * l2
    polynomial[3] = {
        (0, 3): sp.symbols("r0"),
        (1, 2): sp.symbols("r1"),
        (2, 1): sp.symbols("r2"),
        (3, 0): sp.symbols("r3"),
    }
    operator[1].pop((0, 6))
    operator[3] = {
        (0, 8): sp.symbols("k0"),
        (1, 7): sp.symbols("k1"),
        (2, 6): sp.symbols("k2"),
        (3, 5): sp.symbols("k3"),
        (4, 4): sp.symbols("k4"),
    }
    fourth = defect_moment(polynomial, operator, 4, 3)[(0, 1)]
    assert_multiple(fourth, l1**3)


def verify_xy5_and_y6_branches() -> None:
    """Close C=D=B=0, including A != 0 and A=0."""

    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    l0, l1, l2, l3, l4 = sp.symbols("l0:5")
    h0, h1, h2, h3, h4 = sp.symbols("h0:5")
    polynomial = {
        0: {(0, 6): 1, (1, 5): 1},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): p4,
            (5, 0): -6 * l0 - l1,
        },
        2: {
            (0, 4): q0,
            (1, 3): q1,
            (2, 2): q2,
            (3, 1): q3,
            (4, 0): q4,
        },
    }
    operator = {
        0: {(5, 0): 1},
        1: {
            (0, 6): l0,
            (1, 5): l1,
            (2, 4): l2,
            (3, 3): l3,
            (4, 2): l4,
        },
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (3, 4): h3,
            (4, 3): h4,
        },
    }
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    expected = (
        3516 * l0**2
        + 1162 * l0 * l1
        + 84 * l0 * l2
        + 2 * l0 * p4
        + 61 * l1**2
    )
    assert_multiple(second, expected)
    assert_multiple(expected.subs(l0, 0), l1**2)

    # Localize at l0 and normalize l0=1.  Defects two and three solve
    # p4 and p3; a defect-four fifth-moment coefficient is then already
    # a nonzero constant.
    polynomial[1][(5, 0)] = -6 - l1
    solved_p4 = -(3516 + 1162 * l1 + 84 * l2 + 61 * l1**2) / 2
    polynomial[1][(4, 1)] = solved_p4
    polynomial[3] = {
        (0, 3): sp.symbols("r0"),
        (1, 2): sp.symbols("r1"),
        (2, 1): sp.symbols("r2"),
        (3, 0): sp.symbols("r3"),
    }
    operator[1][(0, 6)] = 1
    operator[3] = {
        (0, 8): sp.symbols("k0"),
        (1, 7): sp.symbols("k1"),
        (2, 6): sp.symbols("k2"),
        (3, 5): sp.symbols("k3"),
        (4, 4): sp.symbols("k4"),
    }
    third = defect_moment(polynomial, operator, 3, 3)[(0, 0)]
    solved_p3 = (
        3034 * l1**3
        - 21204 * l1**2
        - 10332 * l1 * l2
        - 685530 * l1
        - 131292 * l2
        - 9009 * l3
        - 1835172
    ) / 198
    assert sp.simplify(third.subs(p3, solved_p3)) == 0
    polynomial[1][(3, 2)] = solved_p3
    polynomial[4] = {
        (0, 2): sp.symbols("s0"),
        (1, 1): sp.symbols("s1"),
        (2, 0): sp.symbols("s2"),
    }
    operator[4] = {
        (0, 9): sp.symbols("z0"),
        (1, 8): sp.symbols("z1"),
        (2, 7): sp.symbols("z2"),
        (3, 6): sp.symbols("z3"),
        (4, 5): sp.symbols("z4"),
    }
    fifth = defect_moment(polynomial, operator, 5, 4)[(0, 1)]
    assert fifth == 9306726025998591590400000000

    # The pure-y chart has a complete two-term weight-thirty face.
    t = sp.symbols("t")
    pure_y_polynomial = {(0, 6): 1, (5, 0): -6 * t}
    pure_y_operator = {(5, 0): 1, (0, 6): t}
    assert_multiple(
        moment(pure_y_polynomial, pure_y_operator, 2)[(0, 0)],
        t**2,
    )


def verify_final_weight_separators() -> None:
    # The five projective tips use weights (2,1), (3,2), (4,3),
    # (5,4), and (6,5), all with threshold 5*w(x).
    for x_weight, y_weight, x_tip in (
        (2, 1, 4),
        (3, 2, 3),
        (4, 3, 2),
        (5, 4, 1),
        (6, 5, 0),
    ):
        threshold = 5 * x_weight
        assert x_tip * x_weight + (6 - x_tip) * y_weight == threshold
        assert x_tip < 5

    # After the order-ten equality chain on the x^4*y^2 chart dies,
    # every later Weierstrass remainder monomial has weight > 10.
    for degree in range(11, 25):
        support = tuple((x_degree, degree - x_degree) for x_degree in range(5))
        assert min(2 * x + y for x, y in support) > 10

    # On the other charts, all later jets are already strict at order 7.
    for x_weight, y_weight, threshold, first_strict_order in (
        (3, 2, 15, 8),
        (4, 3, 20, 7),
        (5, 4, 25, 7),
        (6, 5, 30, 7),
    ):
        for degree in range(first_strict_order, 20):
            support = tuple(
                (x_degree, degree - x_degree)
                for x_degree in range(5)
            )
            assert min(
                x_weight * x + y_weight * y for x, y in support
            ) > threshold


def main() -> None:
    verify_defect_one_radical()
    verify_x4y2_terminal_branch()
    verify_x3y3_branch()
    verify_x2y4_branch()
    verify_xy5_and_y6_branches()
    verify_final_weight_separators()
    print("PASS quintuple-root sextic defect-one radical and five branches")
    print("PASS terminal weight-ten equality chain dies at pure moment six")
    print("PASS arbitrary later jets have strict or one-sided final faces")


if __name__ == "__main__":
    main()
