#!/usr/bin/env python3
"""Verify the complete binary (r,deg(P))=(4,6) GVC row.

The proof is in
``extended-geometry/BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md``.
Singular is required for the exact radical computations.
"""

from __future__ import annotations

import sympy as sp

from verify_binary_quartic_triple_simple_root_gvc import (
    assert_multiple,
    defect_moment,
    moment,
    verify_radical,
)


ROOT_PARTITIONS = (
    (4,),
    (3, 1),
    (2, 2),
    (2, 1, 1),
    (1, 1, 1, 1),
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
    def extend(position: int, used: frozenset[int]) -> bool:
        if position == len(derivative_directions):
            return True
        direction = derivative_directions[position]
        return any(
            factor not in used
            and annihilator != direction
            and extend(position + 1, used | {factor})
            for factor, annihilator in enumerate(polynomial_annihilators)
        )

    return extend(0, frozenset())


def verify_hall_locus() -> None:
    for partition in ROOT_PARTITIONS:
        directions = tuple(
            direction
            for direction, multiplicity in enumerate(partition)
            for _ in range(multiplicity)
        )
        for counts in compositions(6, len(partition) + 1):
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


def scalar_moments(polynomial, operator, first: int, last: int):
    equations = []
    for order in range(first, last + 1):
        value = moment(polynomial, operator, order)
        assert set(value).issubset({(0, 0)})
        equations.extend(value.values())
    return tuple(equations)


def verify_multiplicity_four_faces() -> None:
    # D != 0 terminal face.  Moments four and five kill (t,u).
    t, u = sp.symbols("t u")
    d_operator = {
        (4, 0): 1,
        (3, 3): t,
        (2, 6): -sp.Rational(33, 80) * t**2,
        (1, 9): sp.Rational(1029, 1920) * t**3,
        (0, 12): u,
    }
    d_polynomial = {(3, 3): 1, (4, 0): -sp.Rational(3, 2) * t}
    d4 = moment(d_polynomial, d_operator, 4)[(0, 0)]
    d5 = moment(d_polynomial, d_operator, 5)[(0, 0)]
    assert_multiple(d4, 44113 * t**4 + 51200 * u)
    assert sp.factor(
        d5.subs(u, -sp.Rational(44113, 51200) * t**4)
    ).as_poly(t).degree() == 5
    verify_radical((t, u), (d4, d5), (t, u))

    # B != 0 terminal face.  Moment one has already imposed q.
    ell, h, k, z, p = sp.symbols("ell h k z p")
    q = -2 * h - sp.Rational(1, 2) * ell * p
    b_operator = {
        (4, 0): 1,
        (3, 2): ell,
        (2, 4): h,
        (1, 6): k,
        (0, 8): z,
    }
    b_polynomial = {(2, 4): 1, (3, 2): p, (4, 0): q}
    equations = scalar_moments(b_polynomial, b_operator, 2, 6)

    # Fast exact cover of weighted projective space: ell != 0, p != 0,
    # and the boundary ell=p=0.  The three radicals glue to
    # (z,k,h,ell*p).
    second = equations[0]
    z_ell = sp.solve(second.subs(ell, 1), z)[0]
    ell_chart = tuple(
        sp.together(equation.subs({ell: 1, z: z_ell}))
        for equation in equations[1:]
    )
    verify_radical((h, k, p), ell_chart, (h, k, p))

    z_p = sp.solve(second.subs(p, 1), z)[0]
    p_chart = tuple(
        sp.together(equation.subs({p: 1, z: z_p}))
        for equation in equations[1:]
    )
    verify_radical((h, k, ell), p_chart, (h, k, ell))

    boundary = tuple(
        equation.subs({ell: 0, p: 0}) for equation in equations
    )
    verify_radical((h, k, z), boundary, (h, k, z))

    # A != 0 equality face is killed already by moment two.
    a = sp.symbols("a")
    assert_multiple(
        moment(
            {(1, 5): 1, (5, 0): -a},
            {(4, 0): 1, (0, 5): a},
            2,
        )[(2, 0)],
        a**2,
    )

    # Pure-y boundary: l0 first dies, then the l1/p3 ratio is impossible.
    l0 = sp.symbols("l0")
    assert_multiple(
        moment(
            {(0, 6): 1, (4, 1): -30 * l0},
            {(4, 0): 1, (0, 5): l0},
            2,
        )[(0, 2)],
        l0**2,
    )
    a1, b3, ratio = sp.symbols("a1 b3 ratio")
    pair_operator = {(4, 0): 1, (1, 4): a1}
    pair_polynomial = {(0, 6): 1, (3, 2): b3}
    pair2 = moment(pair_polynomial, pair_operator, 2)[(1, 0)]
    pair3 = moment(pair_polynomial, pair_operator, 3)[(0, 2)]
    assert_multiple(pair2, a1 * b3 * (14 * a1 + b3))
    assert_multiple(
        pair3,
        a1 * b3 * (2002 * a1**2 + 30 * a1 * b3 + b3**2),
    )
    assert sp.resultant(
        14 + ratio, 2002 + 30 * ratio + ratio**2, ratio
    ) != 0


def verify_multiplicity_three_faces() -> None:
    # A != 0: the two sub-threshold pairs die successively.
    u = sp.symbols("u")
    assert_multiple(
        moment(
            {(1, 5): 1, (4, 1): -5 * u},
            {(3, 1): 1, (0, 5): u},
            2,
        )[(2, 0)],
        u**2,
    )
    assert_multiple(
        moment(
            {(1, 5): 1, (3, 2): -10 * u},
            {(3, 1): 1, (1, 4): u},
            2,
        )[(0, 2)],
        u**2,
    )

    # The h0/q4 extremal ratio equations are coprime.
    h0, q4, ratio = sp.symbols("h0 q4 ratio")
    ratio_operator = {(3, 1): 1, (0, 6): h0}
    ratio_polynomial = {(1, 5): 1, (4, 0): q4}
    f = sp.Poly(
        moment(ratio_polynomial, ratio_operator, 3)[(0, 2)], h0, q4
    ).primitive()[1].as_expr()
    g = sp.Poly(
        moment(ratio_polynomial, ratio_operator, 4)[(1, 1)], h0, q4
    ).primitive()[1].as_expr()
    assert sp.resultant(
        f.subs({h0: 1, q4: ratio}),
        g.subs({h0: 1, q4: ratio}),
        ratio,
    ) == 257444582155895534400

    # Terminal cubic face.
    ell, h, k, p = sp.symbols("ell h k p")
    q = -20 * h - 2 * ell * p
    operator = {
        (3, 1): 1,
        (2, 3): ell,
        (1, 5): h,
        (0, 7): k,
    }
    polynomial = {(1, 5): 1, (2, 3): p, (3, 1): q}
    equations = scalar_moments(polynomial, operator, 2, 5)
    verify_radical((ell, h, k, p), equations, (k, h, ell * p))

    # Pure-sixth-power top face: this is a bounded face calculation only.
    # It records exactly what remains to be coupled to adjacent weight levels.
    a, b, q = sp.symbols("a b q")
    top_operator = {(3, 1): 1, (1, 4): a}
    top_polynomial = {(0, 6): 1, (2, 3): b, (4, 0): q}
    top_equations = tuple(
        coefficient
        for order in range(2, 5)
        for coefficient in moment(top_polynomial, top_operator, order).values()
    )
    verify_radical((a, b, q), top_equations, (q, a * b))

    # Complete adjacent-weight face at P_6=y^6.  The strict cofactor A and
    # all four matched spectator pairs are retained.
    A, ell, ell2, h0, h1, p1, p2, q3 = sp.symbols(
        "A ell ell2 h0 h1 p1 p2 q3"
    )
    endpoint_operator = {
        (3, 1): 1,
        (4, 0): A,
        (1, 4): ell,
        (2, 3): ell2,
        (0, 6): h0,
        (1, 5): h1,
    }
    endpoint_polynomial = {
        (0, 6): 1,
        (1, 4): p1,
        (2, 3): p2,
        (3, 1): q3,
    }
    endpoint_equations = tuple(
        coefficient
        for order in range(1, 7)
        for coefficient in moment(
            endpoint_polynomial, endpoint_operator, order
        ).values()
    )
    verify_radical(
        (A, ell, ell2, h0, h1, p1, p2, q3),
        endpoint_equations,
        (
            q3,
            h0,
            h1 * p2,
            ell2 * p2,
            ell * p2,
            ell * p1,
        ),
    )

    # On the only balanced tail component p2 != 0, moment three kills
    # the first surviving pure jet.
    tail = sp.symbols("tail")
    assert_multiple(
        moment(
            {(2, 3): 1},
            {(3, 1): 1, (0, 7): tail},
            3,
        )[(0, 0)],
        tail,
    )

    # On p1=p2=0, the only two successive cost-two balances are
    # ell1/q2 and then ell1/r3.
    tail_ell, tail_q2, tail_r3 = sp.symbols(
        "tail_ell tail_q2 tail_r3"
    )
    assert_multiple(
        moment(
            {(0, 6): 1, (2, 2): tail_q2},
            {(3, 1): 1, (1, 4): tail_ell},
            2,
        )[(0, 0)],
        tail_ell**2 * tail_q2,
    )
    assert_multiple(
        moment(
            {(0, 6): 1, (3, 0): tail_r3},
            {(3, 1): 1, (1, 4): tail_ell},
            3,
        )[(0, 0)],
        tail_ell**3 * tail_r3,
    )


def verify_multiplicity_two_faces() -> None:
    # Successive migrating pairs.
    u = sp.symbols("u")
    assert_multiple(
        moment(
            {(1, 5): 1, (2, 3): -10 * u},
            {(2, 2): 1, (1, 4): u},
            2,
        )[(0, 2)],
        u**2,
    )
    assert_multiple(
        moment(
            {(1, 5): 1, (3, 1): -420 * u},
            {(2, 2): 1, (0, 6): u},
            3,
        )[(1, 1)],
        u**2,
    )
    assert_multiple(
        moment(
            {(1, 5): 1, (3, 0): -2520 * u},
            {(2, 2): 1, (0, 7): u},
            4,
        )[(0, 2)],
        u**2,
    )

    # Terminal face.
    h, z = sp.symbols("h z")
    operator = {(2, 2): 1, (1, 5): h, (0, 8): z}
    polynomial = {(1, 5): 1, (2, 2): -30 * h}
    equations = scalar_moments(polynomial, operator, 2, 3)
    verify_radical((h, z), equations, (h, z))

    # Pure-sixth-power endpoint: defect one kills p4,p5; its remaining
    # l0/p2 tilted pair dies at moment two.
    u = sp.symbols("pure_u")
    assert_multiple(
        moment(
            {(0, 6): 1, (2, 3): -60 * u},
            {(2, 2): 1, (0, 5): u},
            2,
        )[(0, 2)],
        u**2,
    )

    # Complete weight-relevant face, including the generic quadratic
    # cofactor and the weight-14 and weight-16 spectator blocks.
    A, B, ell, h0, h1, k0, k1, z0, p, q2, q3, rho = sp.symbols(
        "A B ell h0 h1 k0 k1 z0 p q2 q3 rho"
    )
    endpoint_operator = {
        (2, 2): 1,
        (1, 4): ell,
        (0, 6): h0,
        (3, 1): A,
        (1, 5): h1,
        (0, 7): k0,
        (4, 0): B,
        (1, 6): k1,
        (0, 8): z0,
    }
    endpoint_polynomial = {
        (0, 6): 1,
        (1, 4): p,
        (2, 2): q2,
        (3, 1): q3,
        (3, 0): rho,
    }
    endpoint_equations = tuple(
        coefficient
        for order in range(1, 7)
        for coefficient in moment(
            endpoint_polynomial, endpoint_operator, order
        ).values()
    )
    verify_radical(
        (A, B, ell, h0, h1, k0, k1, z0, p, q2, q3, rho),
        endpoint_equations,
        (rho, q3, q2, h0, ell * p),
    )


def verify_simple_root_complete() -> None:
    # Defect one, with the complete quartic cofactor and normalized W5.
    A, B, G = sp.symbols("simple_A simple_B simple_G")
    t, l3, l4, l5 = sp.symbols("simple_t simple_l3 simple_l4 simple_l5")
    p = sp.symbols("simple_p0:6")
    polynomial = {
        0: {(0, 6): 1},
        1: {(i, 5 - i): p[i] for i in range(6)},
    }
    operator = {
        0: {(1, 3): 1, (2, 2): A, (3, 1): B, (4, 0): G},
        1: {(0, 5): t, (3, 2): l3, (4, 1): l4, (5, 0): l5},
    }
    first = defect_moment(polynomial, operator, 1, 1)
    solved = sp.solve(tuple(first.values()), (p[1], p[2]), dict=True)[0]
    polynomial[1][(1, 4)] = solved[p[1]]
    polynomial[1][(2, 3)] = solved[p[2]]
    equations = tuple(
        coefficient
        for order in range(2, 6)
        for coefficient in defect_moment(
            polynomial, operator, order, 1
        ).values()
    )
    verify_radical(
        (A, B, G, t, l3, l4, l5, p[3], p[4], p[5]),
        equations,
        (p[5], p[4], p[3]),
    )

    # The remaining Y5/xy4 pair.
    t = sp.symbols("t")
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 4): -30 * t},
            {(1, 3): 1, (0, 5): t},
            2,
        )[(0, 2)],
        t**2,
    )

    # Defect two, retaining all strict W5 and W6 coefficients.
    p0 = sp.symbols("simple_p0_surv")
    q = sp.symbols("simple_q0:5")
    s, h4, h5, h6 = sp.symbols("simple_s simple_h4 simple_h5 simple_h6")
    polynomial2 = {
        0: {(0, 6): 1},
        1: {(0, 5): p0},
        2: {(i, 4 - i): q[i] for i in range(5)},
    }
    operator2 = {
        0: {(1, 3): 1, (2, 2): A, (3, 1): B, (4, 0): G},
        1: {(3, 2): l3, (4, 1): l4, (5, 0): l5},
        2: {(0, 6): s, (4, 2): h4, (5, 1): h5, (6, 0): h6},
    }
    equations2 = tuple(
        coefficient
        for order in range(1, 5)
        for coefficient in defect_moment(
            polynomial2, operator2, order, 2
        ).values()
    )
    verify_radical(
        (
            A,
            B,
            G,
            l3,
            l4,
            l5,
            s,
            h4,
            h5,
            h6,
            q[1],
            q[2],
            q[3],
            q[4],
        ),
        equations2,
        (q[4], q[3], q[2], 120 * s + q[1]),
    )

    # The remaining Y6/xy3 pair.
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 3): -120 * s},
            {(1, 3): 1, (0, 6): s},
            2,
        )[(0, 0)],
        s**2,
    )

    # Defect three kills the x^2 and x^3 terms of P3.
    q0 = sp.symbols("simple_q0_surv")
    r = sp.symbols("simple_r0:4")
    k0, k5, k6, k7 = sp.symbols(
        "simple_k0 simple_k5 simple_k6 simple_k7"
    )
    polynomial3 = {
        0: {(0, 6): 1},
        1: {(0, 5): p0},
        2: {(0, 4): q0},
        3: {(i, 3 - i): r[i] for i in range(4)},
    }
    operator3 = {
        0: {(1, 3): 1, (2, 2): A, (3, 1): B, (4, 0): G},
        1: {(3, 2): l3, (4, 1): l4, (5, 0): l5},
        2: {(4, 2): h4, (5, 1): h5, (6, 0): h6},
        3: {(0, 7): k0, (5, 2): k5, (6, 1): k6, (7, 0): k7},
    }
    equations3 = tuple(
        coefficient
        for order in range(1, 5)
        for coefficient in defect_moment(
            polynomial3, operator3, order, 3
        ).values()
    )
    verify_radical(
        (
            A,
            B,
            G,
            l3,
            l4,
            l5,
            h4,
            h5,
            h6,
            k0,
            k5,
            k6,
            k7,
            r[1],
            r[2],
            r[3],
        ),
        equations3,
        (r[3], r[2]),
    )

    # Defect four kills the last x^2 term.
    r0, r1 = sp.symbols("simple_r0_surv simple_r1_surv")
    a = sp.symbols("simple_a0:3")
    z0, z6, z7, z8 = sp.symbols(
        "simple_z0 simple_z6 simple_z7 simple_z8"
    )
    polynomial4 = {
        0: {(0, 6): 1},
        1: {(0, 5): p0},
        2: {(0, 4): q0},
        3: {(0, 3): r0, (1, 2): r1},
        4: {(i, 2 - i): a[i] for i in range(3)},
    }
    operator4 = operator3 | {
        4: {(0, 8): z0, (6, 2): z6, (7, 1): z7, (8, 0): z8}
    }
    equations4 = tuple(
        coefficient
        for order in range(1, 4)
        for coefficient in defect_moment(
            polynomial4, operator4, order, 4
        ).values()
    )
    verify_radical(
        (
            A,
            B,
            G,
            l3,
            l4,
            l5,
            h4,
            h5,
            h6,
            k0,
            k5,
            k6,
            k7,
            z0,
            z6,
            z7,
            z8,
            r1,
            a[1],
            a[2],
        ),
        equations4,
        (a[2],),
    )


def main() -> None:
    verify_hall_locus()
    verify_multiplicity_four_faces()
    verify_multiplicity_three_faces()
    verify_multiplicity_two_faces()
    verify_simple_root_complete()
    print("verified complete binary quartic-leading sextic GVC row")


if __name__ == "__main__":
    main()
