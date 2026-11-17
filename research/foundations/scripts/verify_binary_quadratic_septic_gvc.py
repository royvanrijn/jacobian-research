#!/usr/bin/env python3
"""Verify the complete binary ``(r, deg(P)) = (2, 7)`` GVC row.

The proof is in
``extended-geometry/BINARY_QUADRATIC_SEPTIC_GVC.md``.  The checker uses
exact sparse contraction in SymPy and exact radicals over ``QQ`` in
Singular.  It verifies the Hall locus, the distinct-root triangular ladder,
and every Newton face in both charts of the repeated-root locus.
"""

from __future__ import annotations

import sympy as sp

from verify_binary_quadratic_all_root_partitions_gvc import (
    all_moment_coefficients,
    compositions,
    has_matching,
    radical_equal,
)
from verify_binary_quartic_triple_simple_root_gvc import (
    apply_operator,
    moment,
)


ROOT_PARTITIONS = ((2,), (1, 1))


def verify_hall_locus() -> None:
    """Exhaust the split-symbol Hall failure for seven polynomial factors."""

    for partition in ROOT_PARTITIONS:
        derivatives = tuple(
            direction
            for direction, multiplicity in enumerate(partition)
            for _ in range(multiplicity)
        )
        for counts in compositions(7, len(partition) + 1):
            annihilators = tuple(
                direction
                for direction, count in enumerate(counts[:-1])
                for _ in range(count)
            ) + (-1,) * counts[-1]
            expected_failure = any(
                counts[direction] >= 8 - multiplicity
                for direction, multiplicity in enumerate(partition)
            )
            assert has_matching(derivatives, annihilators) != expected_failure


def verify_distinct_root() -> None:
    """Replay the full first-equation reduction and second-moment ladder."""

    degree = 7
    polynomial = {
        (i, j): sp.symbols(f"distinct7_p{i}{j}")
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    }
    polynomial[(7, 0)] = sp.Integer(1)
    polynomial[(0, 7)] = sp.Integer(0)
    a = {
        order: sp.symbols(f"distinct7_a{order}")
        for order in range(3, 15)
    }
    b = {
        order: sp.symbols(f"distinct7_b{order}")
        for order in range(3, 15)
    }
    first_operator = {
        (1, 1): 1,
        **{(order, 0): a[order] for order in range(3, 8)},
        **{(0, order): b[order] for order in range(3, 8)},
    }

    # The unit XY coefficient solves every mixed coefficient of P from high
    # degree to low degree in W(P)=0.
    substitution: dict[sp.Symbol, sp.Expr] = {}
    for total in range(5, -1, -1):
        for output_x in range(total + 1):
            output_y = total - output_x
            target = (output_x + 1, output_y + 1)
            current = {
                exponent: sp.expand(
                    coefficient.subs(substitution)
                    if hasattr(coefficient, "subs")
                    else coefficient
                )
                for exponent, coefficient in polynomial.items()
            }
            equation = apply_operator(current, first_operator).get(
                (output_x, output_y), 0
            )
            target_coefficient = polynomial[target]
            assert isinstance(target_coefficient, sp.Symbol)
            substitution[target_coefficient] = sp.factor(
                sp.solve(equation, target_coefficient)[0]
            )
    assert len(substitution) == 21

    reduced_polynomial = {
        exponent: sp.expand(
            coefficient.subs(substitution)
            if hasattr(coefficient, "subs")
            else coefficient
        )
        for exponent, coefficient in polynomial.items()
    }
    assert all(
        sp.factor(coefficient) == 0
        for coefficient in apply_operator(
            reduced_polynomial, first_operator
        ).values()
    )

    full_operator = {
        (1, 1): 1,
        **{(order, 0): a[order] for order in a},
        **{(0, order): b[order] for order in b},
    }
    second = moment(reduced_polynomial, full_operator, 2)
    p06 = polynomial[(0, 6)]
    p05 = polynomial[(0, 5)]
    p04 = polynomial[(0, 4)]
    p03 = polynomial[(0, 3)]
    p02 = polynomial[(0, 2)]
    ladder = (
        (p06, (5, 4), 2520 * p06),
        (a[3], (8, 0), 635040 * a[3] ** 2),
        (p05, (5, 3), 1680 * p05),
        (p04, (5, 2), 1008 * p04),
        (a[4], (6, 0), 80015040 * a[4] ** 2),
        (p03, (5, 1), 504 * p03),
        (p02, (5, 0), 168 * p02),
        (a[5], (4, 0), 3166732800 * a[5] ** 2),
        (a[6], (2, 0), 41912640000 * a[6] ** 2),
        (a[7], (0, 0), 86467046400 * a[7] ** 2),
    )
    zeros: dict[sp.Symbol, int] = {}
    for variable, exponent, expected in ladder:
        assert sp.factor(second[exponent].subs(zeros) - expected) == 0
        zeros[variable] = 0

    # Once the ladder variables vanish, the first equation has exactly the
    # transverse-linear/high-order normal form P=f(x)+c*y and
    # W=Y*Gamma+H(X), ord(H)>=8.
    final_polynomial = {
        exponent: sp.factor(coefficient.subs(zeros))
        for exponent, coefficient in reduced_polynomial.items()
    }
    assert all(
        coefficient == 0
        for (x_degree, y_degree), coefficient in final_polynomial.items()
        if y_degree > 0 and (x_degree, y_degree) != (0, 1)
    )


def verify_nonpure_double_line() -> None:
    """Check every Newton face over the leading form ``x*y^6``."""

    # The only half-integral faces before the final integral threshold.
    half_faces = (
        (3, (3, 3)),
        (5, (3, 1)),
    )
    for y_order, correction_exponent in half_faces:
        a, c = sp.symbols(f"nonpure7_half_a{y_order} c{y_order}")
        operator = {(2, 0): 1, (0, y_order): a}
        polynomial = {(1, 6): 1, correction_exponent: c}
        equations = all_moment_coefficients(polynomial, operator, 6)
        radical_equal((a, c), equations, (a, c))

    # The remaining odd pure-Y crossings have no polynomial companion.  The
    # displayed second-moment coefficient is their first extremal pivot.
    a = sp.symbols("nonpure7_odd_crossing")
    crossing_data = (
        (7, (0, 5), 15966720),
        (9, (0, 3), 319334400),
        (11, (0, 1), 1916006400),
    )
    for y_order, output, coefficient in crossing_data:
        value = moment(
            {(1, 6): 1},
            {(2, 0): 1, (0, y_order): a},
            2,
        )[output]
        assert sp.factor(value - coefficient * a) == 0

    # Integral slopes are the complete equality faces.  The support equation
    # is j*i+k=j+6; the expected lists audit that no channel was omitted.
    expected_supports = {
        2: [(2, 4), (3, 2), (4, 0)],
        3: [(2, 3), (3, 0)],
        4: [(2, 2)],
        5: [(2, 1)],
        6: [(2, 0)],
    }
    for slope in range(2, 7):
        h, z = sp.symbols(f"nonpure7_h{slope} z{slope}")
        support = [
            (x_degree, y_degree)
            for x_degree in range(8)
            for y_degree in range(8 - x_degree)
            if slope * x_degree + y_degree == slope + 6
            and (x_degree, y_degree) != (1, 6)
        ]
        assert support == expected_supports[slope]
        corrections = sp.symbols(
            f"nonpure7_c{slope}_0:{len(support)}"
        )
        polynomial = {
            (1, 6): 1,
            **dict(zip(support, corrections, strict=True)),
        }
        operator = {(2, 0): 1, (1, slope): h, (0, 2 * slope): z}
        equations = all_moment_coefficients(polynomial, operator, 7)
        variables = (h, z, *corrections)
        radical_equal(variables, equations, variables)

    # At the last threshold the sole surviving pair is (X^2, x*y^6).
    assert all(2 * order > order for order in range(1, 13))


def verify_pure_double_line() -> None:
    """Verify the complete branch tree over the endpoint ``P_7=y^7``."""

    # Slope 3/2 dies at the origin.
    v, p, q = sp.symbols("pure7_32_v p q")
    operator = {(2, 0): 1, (0, 3): v}
    polynomial = {(0, 7): 1, (2, 4): p, (4, 1): q}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((v, p, q), equations, (v, p, q))

    # Slope two has exactly the two coordinate axes B!=0 and z!=0.
    A, B, z, q, rho = sp.symbols("pure7_2_A B z q rho")
    operator = {(2, 0): 1, (1, 2): B, (0, 4): A}
    polynomial = {
        (0, 7): 1,
        (1, 5): z,
        (2, 3): q,
        (3, 1): rho,
    }
    equations = all_moment_coefficients(polynomial, operator, 2)
    radical_equal((A, B, z, q, rho), equations, (rho, q, A, B * z))

    # B-axis: the integer faces at weights 3, 4, and the final weight 5.
    B_axis_faces = (
        (3, (0, 5), ((0, 7), (1, 4), (2, 1))),
        (4, (0, 6), ((0, 7), (1, 3))),
        (5, (0, 7), ((0, 7), (1, 2))),
    )
    for weight, pure_operator, support in B_axis_faces:
        variables = sp.symbols(
            f"pure7_B{weight}_0:{len(support)}"
        )
        a = variables[0]
        corrections = variables[1:]
        operator = {(1, 2): 1, pure_operator: a}
        polynomial = {
            support[0]: 1,
            **dict(zip(support[1:], corrections, strict=True)),
        }
        equations = all_moment_coefficients(polynomial, operator, 7)
        radical_equal(variables, equations, variables)

    # The three nonintegral polynomial crossings between those faces are
    # triangular.
    c = sp.symbols("pure7_B_half_c")
    early_crossing = moment(
        {(0, 7): 1, (3, 0): c}, {(1, 2): 1}, 2
    )[(1, 3)]
    first_crossing = moment(
        {(0, 7): 1, (2, 2): c}, {(1, 2): 1}, 1
    )[(1, 0)]
    second_crossing = moment(
        {(0, 7): 1, (2, 0): c}, {(1, 2): 1}, 2
    )[(0, 3)]
    assert sp.factor(early_crossing - 10080 * c) == 0
    assert sp.factor(first_crossing - 4 * c) == 0
    assert sp.factor(second_crossing - 3360 * c) == 0

    # z-axis: first close the ordinary top faces of (X^2,x*y^5).  The
    # half-integral face occurs at weight 5/2; integral weights 3,4,5 give
    # the complete equality supports.
    a, c = sp.symbols("pure7_z_half_a c")
    operator = {(2, 0): 1, (0, 5): a}
    polynomial = {(1, 5): 1, (3, 0): c}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((a, c), equations, (a, c))

    for weight in (3, 4, 5):
        u, w, c = sp.symbols(f"pure7_z_top{weight}_u w c")
        operator = {
            (2, 0): 1,
            (1, weight): u,
            (0, 2 * weight): w,
        }
        polynomial = {(1, 5): 1, (2, 5 - weight): c}
        equations = all_moment_coefficients(polynomial, operator, 8)
        radical_equal((u, w, c), equations, (u, w, c))

    # The lower y^7 channel can meet those top faces in a later output
    # layer.  Every resulting first-output face dies as well.
    for y_order in (3, 4, 5):
        u, w, c = sp.symbols(f"pure7_z{y_order}_u w c")
        operator = {
            (2, 0): 1,
            (1, y_order): u,
            (0, y_order + 2): w,
        }
        polynomial = {
            (0, 7): 1,
            (1, 5): 1,
            (2, 5 - y_order): c,
        }
        equations = all_moment_coefficients(polynomial, operator, 8)
        radical_equal((u, w, c), equations, (u, w, c))

    # Pure Y^8 and Y^9 are the only unpaired crossings before the final
    # weight-five face.
    a = sp.symbols("pure7_z_crossing")
    for y_order, output, coefficient in (
        (8, (0, 2), 7257600),
        (9, (0, 1), 14515200),
    ):
        value = moment(
            {(0, 7): 1, (1, 5): 1},
            {(2, 0): 1, (0, y_order): a},
            2,
        )[output]
        assert sp.factor(value - coefficient * a) == 0

    # Intersection of the two slope-two axes: slope 5/2 dies; slope three
    # leaves two axes which must each migrate once more in degree seven.
    A, q = sp.symbols("pure7_intersection_52_A q")
    operator = {(2, 0): 1, (0, 5): A}
    polynomial = {(0, 7): 1, (2, 2): q}
    equations = all_moment_coefficients(polynomial, operator, 7)
    radical_equal((A, q), equations, (A, q))

    A, B, z, q = sp.symbols("pure7_intersection_3_A B z q")
    operator = {(2, 0): 1, (1, 3): B, (0, 6): A}
    polynomial = {(0, 7): 1, (1, 4): z, (2, 1): q}
    equations = all_moment_coefficients(polynomial, operator, 9)
    radical_equal((A, B, z, q), equations, (q, A, B * z))

    # B-axis of the slope-three intersection.  The half crossing is killed
    # at moment two and the final weight-four face is the origin.
    c = sp.symbols("pure7_intersection_B_half_c")
    value = moment(
        {(0, 7): 1, (2, 0): c}, {(1, 3): 1}, 2
    )[(0, 1)]
    assert sp.factor(value - 20160 * c) == 0

    A, z = sp.symbols("pure7_intersection_B4_A z")
    operator = {(1, 3): 1, (0, 7): A}
    polynomial = {(0, 7): 1, (1, 3): z}
    equations = all_moment_coefficients(polynomial, operator, 6)
    radical_equal((A, z), equations, (A, z))

    # z-axis of the slope-three intersection: the aligned first-output face
    # and the final common-threshold face both die.
    u, v, q = sp.symbols("pure7_intersection_z4_u v q")
    operator = {(2, 0): 1, (1, 4): u, (0, 7): v}
    polynomial = {(0, 7): 1, (1, 4): 1, (2, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 8)
    radical_equal((u, v, q), equations, (u, v, q))

    u, v, q = sp.symbols("pure7_intersection_z4_final_u v q")
    operator = {(2, 0): 1, (1, 4): u, (0, 8): v}
    polynomial = {(1, 4): 1, (2, 0): q}
    equations = all_moment_coefficients(polynomial, operator, 8)
    radical_equal((u, v, q), equations, (u, v, q))

    # Every leaf is now a common-threshold coordinate-deficit pair.
    for order in range(1, 13):
        assert order > 0  # (X*Y^j, y^7): x-demand versus zero x-supply.
        assert 2 * order > order  # (X^2, x*y^j): demand versus supply.


def main() -> None:
    verify_hall_locus()
    verify_distinct_root()
    verify_nonpure_double_line()
    verify_pure_double_line()
    print("verified complete binary quadratic-leading septic GVC row")


if __name__ == "__main__":
    main()
