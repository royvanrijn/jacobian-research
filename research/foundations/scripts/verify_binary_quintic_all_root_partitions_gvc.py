#!/usr/bin/env python3
"""Verify the complete r=5, deg(P)=6 binary GVC calculation.

The proof is in
``extended-geometry/BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md``.
This checker verifies the Hall leading-locus classification, the new
multiplicity-four, multiplicity-three, multiplicity-two, and simple-root
local calculations, and then replays the existing multiplicity-five
checker.  Singular is required for the three defect-one radicals.
"""

from __future__ import annotations

import itertools
import sys

import sympy as sp

from verify_binary_quartic_triple_simple_root_gvc import (
    assert_multiple,
    defect_moment,
    moment,
    verify_radical,
)
from verify_binary_quintic_quintuple_root_gvc import (
    verify_defect_one_radical as verify_multiplicity_five_radical,
    verify_final_weight_separators as verify_multiplicity_five_weights,
    verify_x2y4_branch as verify_multiplicity_five_x2y4,
    verify_x3y3_branch as verify_multiplicity_five_x3y3,
    verify_x4y2_terminal_branch as verify_multiplicity_five_x4y2,
    verify_xy5_and_y6_branches as verify_multiplicity_five_xy5_y6,
)


ROOT_PARTITIONS = (
    (5,),
    (4, 1),
    (3, 2),
    (3, 1, 1),
    (2, 2, 1),
    (2, 1, 1, 1),
    (1, 1, 1, 1, 1),
)


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def has_matching(
    derivative_directions: tuple[int, ...],
    polynomial_annihilators: tuple[int, ...],
) -> bool:
    """Can the five derivative copies hit five distinct P factors?"""

    def extend(position: int, used: frozenset[int]) -> bool:
        if position == len(derivative_directions):
            return True
        direction = derivative_directions[position]
        for factor, annihilator in enumerate(polynomial_annihilators):
            if factor in used or annihilator == direction:
                continue
            if extend(position + 1, used | {factor}):
                return True
        return False

    return extend(0, frozenset())


def verify_hall_leading_locus() -> None:
    """Brute-force the Hall classification for every quintic partition."""

    for partition in ROOT_PARTITIONS:
        directions = tuple(
            direction
            for direction, multiplicity in enumerate(partition)
            for _ in range(multiplicity)
        )
        number_of_roots = len(partition)
        # c_j P-factors are annihilated by operator direction j; the last
        # part counts P-factors at all other directions.
        for counts in compositions(6, number_of_roots + 1):
            annihilators = tuple(
                direction
                for direction, count in enumerate(counts[:-1])
                for _ in range(count)
            ) + (-1,) * counts[-1]
            expected_failure = any(
                counts[direction] >= 7 - multiplicity
                for direction, multiplicity in enumerate(partition)
            )
            assert has_matching(directions, annihilators) != expected_failure


def verify_multiplicity_four_family() -> None:
    """Local model X^4*Y with P_6=y^3*C_3."""

    C0, A, B, D = sp.symbols("C0 A B D")
    l0, l1, l2, l3, l6, p5 = sp.symbols("l0 l1 l2 l3 l6 p5")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    polynomial = {
        0: {(0, 6): C0, (1, 5): A, (2, 4): B, (3, 3): D},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): p4,
            (5, 0): p5,
        },
    }
    operator = {
        0: {(4, 1): 1},
        1: {
            (0, 6): l0,
            (1, 5): l1,
            (2, 4): l2,
            (3, 3): l3,
            (6, 0): l6,
        },
    }
    first = defect_moment(polynomial, operator, 1, 1)[(0, 0)]
    solved_p4 = -5 * A * l1 - 2 * B * l2 - 30 * C0 * l0 - sp.Rational(3, 2) * D * l3
    assert sp.simplify(first.subs(p4, solved_p4)) == 0
    polynomial[1][(4, 1)] = solved_p4
    equations = []
    for order in range(2, 5):
        equations.extend(defect_moment(polynomial, operator, order, 1).values())
    expected = (
        D * l1,
        D * l0,
        B * l0,
        D * (15 * D * l2 + 7 * p5),
    )
    verify_radical(
        (C0, A, B, D, l0, l1, l2, l3, l6, p5),
        tuple(equations),
        expected,
    )

    # D != 0: the sole migrating l2/p5 pair dies at defect two.
    h0, h1, h2, h3, h7 = sp.symbols("h0 h1 h2 h3 h7")
    q0, q1, q2, q3, q4 = sp.symbols("q0:5")
    d_polynomial = {
        0: {(0, 6): C0, (1, 5): A, (2, 4): B, (3, 3): 1},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): p3,
            (4, 1): -2 * B * l2 - sp.Rational(3, 2) * l3,
            (5, 0): -sp.Rational(15, 7) * l2,
        },
        2: {(i, 4 - i): q for i, q in enumerate((q0, q1, q2, q3, q4))},
    }
    d_operator = {
        0: {(4, 1): 1},
        1: {(2, 4): l2, (3, 3): l3, (6, 0): l6},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (3, 4): h3,
            (7, 0): h7,
        },
    }
    d4 = defect_moment(d_polynomial, d_operator, 4, 2)[(0, 2)]
    assert_multiple(d4, 196 * h0 + 151 * l2**2)
    solved_h0 = -sp.Rational(151, 196) * l2**2
    d3 = defect_moment(d_polynomial, d_operator, 3, 2)
    assert_multiple(d3[(1, 0)].subs(h0, solved_h0), l2**2)
    assert_multiple(d3[(0, 1)].subs({l2: 0, h0: 0}), h1)
    d2 = defect_moment(d_polynomial, d_operator, 2, 2)[(0, 0)]
    assert_multiple(
        d2.subs({l2: 0, h0: 0, h1: 0}),
        40 * h2 - l3**2,
    )

    # B != 0: either l1=p5=0, or a coprime ratio obstruction.
    s = sp.symbols("s")
    ratio_f = 11 * s**2 + 42 * s + 330
    ratio_g = 221 * s**3 + 546 * s**2 + 2970 * s + 21840
    assert sp.resultant(ratio_f, ratio_g, s) == 1909272615840

    b_l1, b_l2, b_l3, b_l6 = sp.symbols("b_l1 b_l2 b_l3 b_l6")
    b_p0, b_p1, b_p2, b_p3 = sp.symbols("b_p0:4")
    b_q = sp.symbols("b_q0:5")
    b_r = sp.symbols("b_r0:4")
    b_h = sp.symbols("b_h0:5")
    b_k = sp.symbols("b_k0:5")
    b_polynomial = {
        0: {(0, 6): C0, (1, 5): A, (2, 4): 1},
        1: {
            (0, 5): b_p0,
            (1, 4): b_p1,
            (2, 3): b_p2,
            (3, 2): b_p3,
            (4, 1): -5 * A - 2 * b_l2,
            (5, 0): s,
        },
        2: {(i, 4 - i): b_q[i] for i in range(5)},
        3: {(i, 3 - i): b_r[i] for i in range(4)},
    }
    b_operator = {
        0: {(4, 1): 1},
        1: {
            (1, 5): 1,
            (2, 4): b_l2,
            (3, 3): b_l3,
            (6, 0): b_l6,
        },
        2: {(i, 7 - i): b_h[i] for i in range(4)} | {(7, 0): b_h[4]},
        3: {(i, 8 - i): b_k[i] for i in range(4)} | {(8, 0): b_k[4]},
    }
    b_second = defect_moment(b_polynomial, b_operator, 2, 2)[(0, 0)]
    b_polynomial[1][(3, 2)] = sp.solve(b_second, b_p3)[0]
    b_fourth = defect_moment(b_polynomial, b_operator, 4, 3)[(1, 0)]
    assert_multiple(b_fourth, ratio_g)
    b_third = defect_moment(
        {
            0: {(0, 6): C0, (1, 5): A, (2, 4): 1},
            1: {(4, 1): -5 * A, (5, 0): s},
        },
        {0: {(4, 1): 1}, 1: {(1, 5): 1}},
        3,
        2,
    )[(0, 1)]
    assert_multiple(b_third, ratio_f)

    # A != 0: two extremal coefficients eliminate (l0,p5).
    u, v = sp.symbols("u v")
    extreme_polynomial = {0: {(1, 5): 1}, 1: {(5, 0): v}}
    extreme_operator = {0: {(4, 1): 1}, 1: {(0, 6): u}}
    cubic = defect_moment(extreme_polynomial, extreme_operator, 4, 3)[(0, 1)]
    quartic = defect_moment(extreme_polynomial, extreme_operator, 5, 4)[(1, 0)]
    expected_cubic = (
        13 * v**3 + 54 * u * v**2 + 1638 * u**2 * v + 302328 * u**3
    )
    expected_quartic = (
        323 * v**4
        + 680 * u * v**3
        + 8580 * u**2 * v**2
        + 465120 * u**3 * v
        + 98062800 * u**4
    )
    assert_multiple(cubic, expected_cubic)
    assert_multiple(quartic, expected_quartic)
    assert sp.resultant(
        expected_cubic.subs({u: 1, v: s}),
        expected_quartic.subs({u: 1, v: s}),
        s,
    ) == 340491886133329409922608332032

    # Pure y: p5^4 is extremal, then the complete defect-two scalar is l0^2.
    pure_y_polynomial = {
        0: {(0, 6): 1},
        1: {(4, 1): -30 * u, (5, 0): v},
    }
    pure_y_operator = {0: {(4, 1): 1}, 1: {(0, 6): u}}
    assert_multiple(
        defect_moment(pure_y_polynomial, pure_y_operator, 5, 4)[(0, 1)],
        v**4,
    )
    assert_multiple(
        defect_moment(pure_y_polynomial, pure_y_operator, 2, 2)[(0, 0)].subs(v, 0),
        u**2,
    )


def verify_multiplicity_three_family() -> None:
    """Local model X^3*Y^2 with P_6=y^4*C_2."""

    C0, A, B = sp.symbols("C0 A B")
    alpha, beta = sp.symbols("alpha beta")
    l0, l1, l2, l5, l6 = sp.symbols("l0 l1 l2 l5 l6")
    p0, p1, p2, p3, p4, p5 = sp.symbols("p0:6")
    polynomial = {
        0: {(0, 6): C0, (1, 5): A, (2, 4): B},
        1: {(i, 5 - i): p for i, p in enumerate((p0, p1, p2, p3, p4, p5))},
    }
    operator = {
        0: {(3, 2): 1, (4, 1): alpha, (5, 0): beta},
        1: {
            (0, 6): l0,
            (1, 5): l1,
            (2, 4): l2,
            (5, 1): l5,
            (6, 0): l6,
        },
    }
    first = defect_moment(polynomial, operator, 1, 1)[(0, 0)]
    solved_p3 = (
        -10 * A * l1
        - 4 * B * l2
        - 60 * C0 * l0
        - 2 * alpha * p4
        - 10 * beta * p5
    )
    assert sp.simplify(first.subs(p3, solved_p3)) == 0
    polynomial[1][(3, 2)] = solved_p3
    equations = []
    for order in range(2, 4):
        equations.extend(defect_moment(polynomial, operator, order, 1).values())
    expected = (
        B * p5,
        A * p5,
        B * l0,
        B * (56 * B * l1 + 5 * p4),
    )
    verify_radical(
        (C0, A, B, alpha, beta, l0, l1, l2, l5, l6, p4, p5),
        tuple(equations),
        expected,
    )

    # B != 0: l1 is an exposed square; the remaining jets are equality terms.
    h0, h1, h2, h6, h7 = sp.symbols("h0 h1 h2 h6 h7")
    q = sp.symbols("q0:5")
    b_polynomial = {
        0: {(0, 6): C0, (1, 5): A, (2, 4): 1},
        1: {
            (0, 5): p0,
            (1, 4): p1,
            (2, 3): p2,
            (3, 2): (
                -10 * A * l1
                - 4 * l2
                + sp.Rational(112, 5) * alpha * l1
            ),
            (4, 1): -sp.Rational(56, 5) * l1,
        },
        2: {(i, 4 - i): q[i] for i in range(5)},
    }
    b_operator = {
        0: {(3, 2): 1, (4, 1): alpha, (5, 0): beta},
        1: {(1, 5): l1, (2, 4): l2, (5, 1): l5, (6, 0): l6},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (6, 1): h6,
            (7, 0): h7,
        },
    }
    b_third = defect_moment(b_polynomial, b_operator, 3, 2)
    assert_multiple(b_third[(1, 0)], l1**2)
    assert_multiple(b_third[(0, 1)].subs(l1, 0), h0)
    b_second = defect_moment(b_polynomial, b_operator, 2, 2)[(0, 0)]
    assert_multiple(
        b_second.subs({l1: 0, h0: 0}),
        56 * h1 - 4 * l2**2 + q[4],
    )

    # A != 0: the l0/p4 ratio equations are coprime.
    s = sp.symbols("s")
    ratio_f = 7 * s**2 + 330 * s + 30030
    ratio_g = 13 * s**3 + 540 * s**2 + 32760 * s + 3023280
    assert sp.resultant(ratio_f, ratio_g, s) == 5436596606290200

    a_l1, a_l2, a_l5, a_l6 = sp.symbols("a_l1 a_l2 a_l5 a_l6")
    a_p0, a_p1, a_p2 = sp.symbols("a_p0:3")
    a_q = sp.symbols("a_q0:5")
    a_r = sp.symbols("a_r0:4")
    a_h = sp.symbols("a_h0:5")
    a_k = sp.symbols("a_k0:5")
    a_polynomial = {
        0: {(0, 6): C0, (1, 5): 1},
        1: {
            (0, 5): a_p0,
            (1, 4): a_p1,
            (2, 3): a_p2,
            (3, 2): -10 * a_l1 - 60 * C0 - 2 * alpha * s,
            (4, 1): s,
        },
        2: {(i, 4 - i): a_q[i] for i in range(5)},
        3: {(i, 3 - i): a_r[i] for i in range(4)},
    }
    a_operator = {
        0: {(3, 2): 1, (4, 1): alpha, (5, 0): beta},
        1: {
            (0, 6): 1,
            (1, 5): a_l1,
            (2, 4): a_l2,
            (5, 1): a_l5,
            (6, 0): a_l6,
        },
        2: {
            (0, 7): a_h[0],
            (1, 6): a_h[1],
            (2, 5): a_h[2],
            (6, 1): a_h[3],
            (7, 0): a_h[4],
        },
        3: {
            (0, 8): a_k[0],
            (1, 7): a_k[1],
            (2, 6): a_k[2],
            (7, 1): a_k[3],
            (8, 0): a_k[4],
        },
    }
    a_second = defect_moment(a_polynomial, a_operator, 2, 2)[(0, 0)]
    a_polynomial[1][(2, 3)] = sp.solve(a_second, a_p2)[0]
    a_fourth = defect_moment(a_polynomial, a_operator, 4, 3)[(1, 0)]
    assert_multiple(a_fourth, ratio_g)
    assert_multiple(
        defect_moment(
            {0: {(1, 5): 1}, 1: {(4, 1): s}},
            {0: {(3, 2): 1}, 1: {(0, 6): 1}},
            3,
            2,
        )[(0, 1)],
        ratio_f,
    )

    # Pure y: extremal coefficients kill p5, then p4, then l0.
    u, v, w = sp.symbols("u v w")
    pure_polynomial = {
        0: {(0, 6): 1},
        1: {(3, 2): -60 * u, (4, 1): v, (5, 0): w},
    }
    pure_operator = {0: {(3, 2): 1}, 1: {(0, 6): u}}
    assert_multiple(
        defect_moment(pure_polynomial, pure_operator, 3, 2)[(1, 0)],
        w**2,
    )
    assert_multiple(
        defect_moment(pure_polynomial, pure_operator, 4, 3)[(0, 1)].subs(w, 0),
        v**3,
    )
    assert_multiple(
        defect_moment(pure_polynomial, pure_operator, 2, 2)[(0, 0)].subs({v: 0, w: 0}),
        u**2,
    )


def verify_multiplicity_two_family() -> None:
    """Swapped local model X^2*Y^3 with P_6=y^5*C_1."""

    # It is computationally cheaper to use the swapped X^3*Y^2 chart,
    # where the family is P_6=x^5(A*y+B*x).
    A, B = sp.symbols("A B")
    alpha, beta, gamma = sp.symbols("alpha beta gamma")
    l0, l1, l2, l5, l6 = sp.symbols("l0 l1 l2 l5 l6")
    p0, p1, p2, p3, p4, p5 = sp.symbols("p0:6")
    polynomial = {
        0: {(5, 1): A, (6, 0): B},
        1: {(i, 5 - i): p for i, p in enumerate((p0, p1, p2, p3, p4, p5))},
    }
    operator = {
        0: {
            (3, 2): 1,
            (2, 3): alpha,
            (1, 4): beta,
            (0, 5): gamma,
        },
        1: {
            (0, 6): l0,
            (1, 5): l1,
            (2, 4): l2,
            (5, 1): l5,
            (6, 0): l6,
        },
    }
    first = defect_moment(polynomial, operator, 1, 1)[(0, 0)]
    solved_p3 = (
        -10 * A * l5
        - 60 * B * l6
        - alpha * p2
        - 2 * beta * p1
        - 10 * gamma * p0
    )
    assert sp.simplify(first.subs(p3, solved_p3)) == 0
    polynomial[1][(3, 2)] = solved_p3
    equations = []
    for order in range(2, 5):
        equations.extend(defect_moment(polynomial, operator, order, 1).values())
    expected = (
        B * p1,
        A * p1,
        B * p0,
        A * p0,
        A * (60 * A * l6 + p2),
    )
    verify_radical(
        (
            A,
            B,
            alpha,
            beta,
            gamma,
            l0,
            l1,
            l2,
            l5,
            l6,
            p0,
            p1,
            p2,
            p4,
            p5,
        ),
        tuple(equations),
        expected,
    )

    # A != 0: the migrating l6/p2 pair is an exposed square.
    q = sp.symbols("q0:5")
    h0, h1, h2, h6, h7 = sp.symbols("h0 h1 h2 h6 h7")
    a_polynomial = {
        0: {(5, 1): 1, (6, 0): B},
        1: {
            (2, 3): -60 * l6,
            (3, 2): -10 * l5 - 60 * B * l6 + 60 * alpha * l6,
            (4, 1): p4,
            (5, 0): p5,
        },
        2: {(i, 4 - i): q[i] for i in range(5)},
    }
    a_operator = {
        0: {
            (3, 2): 1,
            (2, 3): alpha,
            (1, 4): beta,
            (0, 5): gamma,
        },
        1: {(0, 6): l0, (1, 5): l1, (2, 4): l2, (5, 1): l5, (6, 0): l6},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (6, 1): h6,
            (7, 0): h7,
        },
    }
    a_third = defect_moment(a_polynomial, a_operator, 3, 2)
    assert_multiple(a_third[(0, 1)], l6**2)
    assert_multiple(a_third[(1, 0)].subs(l6, 0), q[0])
    a_second = defect_moment(a_polynomial, a_operator, 2, 2)[(0, 0)]
    assert_multiple(
        a_second.subs({l6: 0, q[0]: 0}),
        420 * h7 - 20 * l5**2 + q[1],
    )

    # Pure x: p2 is an exposed square; the equality chain dies at moment 3.
    t = sp.symbols("t")
    pure_polynomial = {
        0: {(6, 0): 1},
        1: {
            (2, 3): p2,
            (3, 2): -60 * t - alpha * p2,
            (4, 1): p4,
            (5, 0): p5,
        },
        2: {(i, 4 - i): q[i] for i in range(5)},
    }
    pure_operator = {
        0: {
            (3, 2): 1,
            (2, 3): alpha,
            (1, 4): beta,
            (0, 5): gamma,
        },
        1: {(0, 6): l0, (1, 5): l1, (2, 4): l2, (5, 1): l5, (6, 0): t},
        2: {
            (0, 7): h0,
            (1, 6): h1,
            (2, 5): h2,
            (6, 1): h6,
            (7, 0): h7,
        },
    }
    assert_multiple(
        defect_moment(pure_polynomial, pure_operator, 3, 2)[(1, 0)],
        p2**2,
    )
    pure_second = defect_moment(pure_polynomial, pure_operator, 2, 2)[(0, 0)]
    assert_multiple(pure_second.subs(p2, 0), 10620 * t**2 + q[0])
    equality_polynomial = {
        (6, 0): 1,
        (3, 2): -60 * t,
        (0, 4): -10620 * t**2,
    }
    equality_operator = {
        (3, 2): 1,
        (2, 3): alpha,
        (1, 4): beta,
        (0, 5): gamma,
        (6, 0): t,
    }
    assert_multiple(moment(equality_polynomial, equality_operator, 3)[(0, 0)], t**3)


def verify_simple_root_tip() -> None:
    """Universal simple-root sixth-power correction."""

    t, p5 = sp.symbols("t p5")
    alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")
    p = sp.symbols("p0:5")
    polynomial = {
        0: {(6, 0): 1},
        1: {(i, 5 - i): p[i] for i in range(5)} | {(5, 0): p5},
    }
    leading_operator = {
        (4, 1): 1,
        (3, 2): alpha,
        (2, 3): beta,
        (1, 4): gamma,
        (0, 5): delta,
    }
    operator = {0: leading_operator, 1: {(6, 0): t}}
    first = defect_moment(polynomial, operator, 1, 1)[(0, 0)]
    solved_p4 = (
        -sp.Rational(1, 2) * alpha * p[3]
        - sp.Rational(1, 2) * beta * p[2]
        - gamma * p[1]
        - 5 * delta * p[0]
        - 30 * t
    )
    assert sp.simplify(first.subs(p[4], solved_p4)) == 0
    polynomial[1][(4, 1)] = solved_p4
    third_first_defect = defect_moment(polynomial, operator, 3, 1)
    assert_multiple(third_first_defect[(0, 2)], p[0])
    assert_multiple(third_first_defect[(1, 1)].subs(p[0], 0), p[1])
    assert_multiple(
        third_first_defect[(2, 0)].subs({p[0]: 0, p[1]: 0}),
        p[2],
    )
    second_first_defect = defect_moment(polynomial, operator, 2, 1)
    assert_multiple(
        second_first_defect[(1, 0)].subs(
            {p[0]: 0, p[1]: 0, p[2]: 0}
        ),
        p[3],
    )
    for exponent in ((0, 5), (1, 4), (2, 3), (3, 2)):
        polynomial[1].pop(exponent)
    polynomial[1][(4, 1)] = -30 * t

    q = sp.symbols("q0:5")
    polynomial[2] = {(i, 4 - i): q[i] for i in range(5)}
    second = defect_moment(polynomial, operator, 2, 2)[(0, 0)]
    assert_multiple(
        second.subs({q[0]: 0, q[1]: 0}),
        720 * t**2 + q[2],
    )
    polynomial[2][(2, 2)] = -720 * t**2
    third_second_defect = defect_moment(polynomial, operator, 3, 2)
    assert_multiple(third_second_defect[(0, 1)], q[0])
    assert_multiple(third_second_defect[(1, 0)].subs(q[0], 0), q[1])
    polynomial[2].pop((0, 4))
    polynomial[2].pop((1, 3))
    r = sp.symbols("r0:4")
    polynomial[3] = {(i, 3 - i): r[i] for i in range(4)}
    third = defect_moment(polynomial, operator, 3, 3)[(0, 0)]
    assert_multiple(third, 154320 * t**3 + r[0])
    equality_polynomial = {
        (6, 0): 1,
        (4, 1): -30 * t,
        (2, 2): -720 * t**2,
        (0, 3): -154320 * t**3,
    }
    assert_multiple(
        moment(equality_polynomial, leading_operator | {(6, 0): t}, 4)[(0, 0)],
        t**4,
    )


def verify_weight_supports() -> None:
    """Audit every local projective chart and all later normalized jets."""

    for multiplicity in range(1, 6):
        for x_tip in range(multiplicity):
            k = multiplicity - x_tip
            x_weight, y_weight = k + 1, k
            threshold = 5 * k + multiplicity
            assert (
                multiplicity * x_weight
                + (5 - multiplicity) * y_weight
                == threshold
            )
            assert (
                x_tip * x_weight + (6 - x_tip) * y_weight
                == threshold
            )
            # The nonterminal leading-symbol monomials are strict.
            for shift in range(1, 6 - multiplicity):
                assert (
                    (multiplicity + shift) * x_weight
                    + (5 - multiplicity - shift) * y_weight
                    > threshold
                )
            # Normal remainders can be non-strict only at finitely many
            # orders; by threshold+1 every one is strict.
            for degree in range(threshold + 1, threshold + 8):
                support = (
                    (x_degree, degree - x_degree)
                    for x_degree in range(degree + 1)
                    if x_degree < multiplicity
                    or degree - x_degree < 5 - multiplicity
                )
                assert min(
                    x_weight * x + y_weight * y for x, y in support
                ) > threshold


def main() -> None:
    verify_hall_leading_locus()
    verify_multiplicity_four_family()
    verify_multiplicity_three_family()
    verify_multiplicity_two_family()
    verify_simple_root_tip()
    verify_weight_supports()
    verify_multiplicity_five_radical()
    verify_multiplicity_five_x4y2()
    verify_multiplicity_five_x3y3()
    verify_multiplicity_five_x2y4()
    verify_multiplicity_five_xy5_y6()
    verify_multiplicity_five_weights()
    print("PASS Hall leading-locus classification for all quintic partitions")
    print("PASS local correction systems for root multiplicities 1 through 5")
    print("PASS every r=5, deg(P)=6 binary root partition is GVC-safe")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
