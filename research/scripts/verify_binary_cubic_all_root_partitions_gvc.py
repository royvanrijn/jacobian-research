#!/usr/bin/env python3
"""Verify the complete binary (r,deg(P))=(3,6) GVC row.

The proof is in
``extended-geometry/BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md``.
Singular is used by the imported exact-radical helper; msolve is used for
the two weighted projective chart covers whose global characteristic-zero
primary decompositions have severe coefficient swell.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

import sympy as sp

from verify_binary_quartic_triple_simple_root_gvc import (
    assert_multiple,
    defect_moment,
    moment,
    singular_expression,
)


ROOT_PARTITIONS = ((3,), (2, 1), (1, 1, 1))


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def has_matching(derivatives, annihilators) -> bool:
    def extend(position: int, used: frozenset[int]) -> bool:
        if position == len(derivatives):
            return True
        direction = derivatives[position]
        return any(
            factor not in used
            and annihilator != direction
            and extend(position + 1, used | {factor})
            for factor, annihilator in enumerate(annihilators)
        )

    return extend(0, frozenset())


def verify_hall_locus() -> None:
    for partition in ROOT_PARTITIONS:
        derivatives = tuple(
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
            assert has_matching(derivatives, annihilators) != expected_failure


def msolve_empty(variables, equations) -> None:
    executable = shutil.which("msolve")
    if executable is None:
        raise RuntimeError("msolve is required for the affine-chart replay")
    source = (
        ",".join(map(str, variables))
        + "\n0\n"
        + ",\n".join(singular_expression(e) for e in equations)
        + "\n"
    )
    input_name = output_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msolve", delete=False
        ) as handle:
            handle.write(source)
            input_name = handle.name
        output_name = input_name + ".out"
        subprocess.run(
            [executable, "-f", input_name, "-o", output_name, "-t", "4"],
            check=True,
            capture_output=True,
            text=True,
        )
        with open(output_name, encoding="utf-8") as handle:
            assert handle.read().strip() == "[-1]:"
    finally:
        for name in (input_name, output_name):
            if name and os.path.exists(name):
                os.unlink(name)


def primitive_equations(polynomial, operator, first: int, last: int, variables):
    return tuple(
        sp.Poly(coefficient, *variables).primitive()[1].as_expr()
        for order in range(first, last + 1)
        for coefficient in moment(polynomial, operator, order).values()
    )


def verify_triple_root() -> None:
    # x^2 y^4 chart: the defect-one equations remove the two lower
    # slopes.  The weight-two face then collapses completely.
    a, b, c = sp.symbols("triple_a triple_b triple_c")
    p, q = sp.symbols("triple_p triple_q")
    operator = {(3, 0): 1, (2, 2): a, (1, 4): b, (0, 6): c}
    polynomial = {
        (2, 4): 1,
        (3, 2): -4 * a,
        (4, 0): 2 * a**2 - 2 * b,
    }
    second = moment(polynomial, operator, 2)
    assert_multiple(second[(0, 4)], a**2 - 27 * b)
    assert_multiple(second[(1, 2)], 7 * a**3 - 25 * a * b - 84 * c)
    assert_multiple(
        second[(2, 0)],
        8 * a**4 - 33 * a**2 * b + 24 * a * c + 43 * b**2,
    )
    assert (
        sp.factor(
            (8 * a**4 - 33 * a**2 * b + 24 * a * c + 43 * b**2).subs(
                {b: a**2 / 27, c: sp.Rational(41, 567) * a**3}
            )
        )
        == sp.Rational(43744, 5103) * a**4
    )

    # The two last x^2y^4 slopes, at weights (3,1) and (4,1).
    t = sp.symbols("triple_t")
    for weighted_operator, weighted_polynomial, order, expected in (
        (
            {(3, 0): 1, (2, 3): a, (1, 6): b, (0, 9): c},
            {(2, 4): 1, (3, 1): t},
            3,
            a**4,
        ),
        (
            {(3, 0): 1, (2, 4): a, (1, 8): b, (0, 12): c},
            {(2, 4): 1, (3, 0): t},
            4,
            a**4,
        ),
    ):
        equations = tuple(
            moment(weighted_polynomial, weighted_operator, m).values()
            for m in range(1, order + 1)
        )
        flattened = tuple(value for group in equations for value in group)
        solved_t = sp.solve(flattened[0], t)[0]
        solved_b = sp.solve(flattened[1].subs(t, solved_t), b)[0]
        solved_c = sp.solve(
            flattened[2].subs({t: solved_t, b: solved_b}), c
        )[0]
        final = flattened[-1].subs({t: solved_t, b: solved_b, c: solved_c})
        assert_multiple(final, expected)

    # xy^5 weight-two face: radical (h,k,q,l*p), checked by affine
    # saturations over Q.  q is eliminated by moment one.
    ell, h, k, p, z = sp.symbols("ell h k p z")
    q = -20 * h - 2 * ell * p
    op = {(3, 0): 1, (2, 2): ell, (1, 4): h, (0, 6): k}
    pol = {(1, 5): 1, (2, 3): p, (3, 1): q}
    equations = primitive_equations(pol, op, 2, 7, (ell, h, k, p))
    for target in (h, k, p):
        msolve_empty(
            (z, h, k, p),
            tuple(e.subs(ell, 1) for e in equations) + (z * target - 1,),
        )
    for target in (ell, h, k):
        msolve_empty(
            (z, ell, h, k),
            tuple(e.subs(p, 1) for e in equations) + (z * target - 1,),
        )
    msolve_empty((ell, k, p), tuple(e.subs(h, 1) for e in equations))
    msolve_empty((ell, h, p), tuple(e.subs(k, 1) for e in equations))

    # Pure-y endpoint: first two slopes.
    u, v = sp.symbols("pure_u pure_v")
    assert_multiple(
        moment(
            {(0, 6): 1, (3, 2): -60 * u},
            {(3, 0): 1, (0, 4): u},
            2,
        )[(0, 4)],
        u**2,
    )
    face = moment(
        {(0, 6): 1, (2, 3): v, (4, 0): -u * v / 2},
        {(3, 0): 1, (1, 3): u},
        2,
    )
    assert_multiple(face[(0, 3)], u * v * (328 * u + 7 * v))
    assert_multiple(face[(2, 0)], u**2 * v * (12 * u - 13 * v))

    # Complete weight-(2,1) pure-y face.  Moment one eliminates r.
    ell, h, k, p, q, z = sp.symbols("E H K P Q Z")
    r = -(12 * h * p + 360 * k + 2 * ell * q) / 3
    op = {(3, 0): 1, (2, 2): ell, (1, 4): h, (0, 6): k}
    pol = {(0, 6): 1, (1, 4): p, (2, 2): q, (3, 0): r}
    equations = primitive_equations(pol, op, 2, 8, (ell, h, k, p, q))
    for target in (k, p, q):
        msolve_empty(
            (z, k, p, q, ell),
            tuple(e.subs(h, 1) for e in equations) + (z * target - 1,),
        )
    for target in (k, q, h * p):
        msolve_empty(
            (z, h, k, p, q),
            tuple(e.subs(ell, 1) for e in equations)
            + (z * target - 1,),
        )
    msolve_empty(
        (ell, h, p, q), tuple(e.subs(k, 1) for e in equations)
    )

    # XY^3-anchored branch: Y^5, then Y^6, then x^2y.
    b = sp.symbols("anchor_b")
    anchored_operator = {(1, 3): 1, (0, 5): b}
    anchored_polynomial = {
        (0, 6): 1,
        (1, 4): -30 * b,
        (2, 2): -720 * b**2,
        (3, 0): 13680 * b**3,
    }
    assert_multiple(
        moment(anchored_polynomial, anchored_operator, 3)[(0, 3)], b**3
    )
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 3): -120 * b, (2, 0): -113040 * b**2},
            {(1, 3): 1, (0, 6): b},
            3,
        )[(0, 0)],
        b**3,
    )
    assert_multiple(
        defect_moment(
            {0: {(0, 6): 1}, 3: {(2, 1): b}},
            {0: {(3, 0): 1}, 1: {(1, 3): 1}},
            2,
            5,
        )[(0, 1)],
        b,
    )


def verify_double_root() -> None:
    # Non-pure xy^5 chart: the six successive faces.
    a, b, c, p, q = sp.symbols("double_a double_b double_c double_p double_q")
    assert_multiple(
        moment(
            {(1, 5): 1, (3, 2): -10 * a},
            {(2, 1): 1, (0, 4): a},
            2,
        )[(0, 5)],
        a,
    )
    op = {(2, 1): 1, (1, 3): b, (0, 5): c}
    pol = {
        (1, 5): 1,
        (2, 3): -10 * b,
        (3, 1): 20 * b**2 - 20 * c,
    }
    assert_multiple(moment(pol, op, 2)[(0, 4)], c)
    assert_multiple(moment(pol, op, 2)[(1, 2)].subs(c, 0), b**3)
    assert_multiple(
        moment(
            {(1, 5): 1, (3, 0): -2520 * a},
            {(2, 1): 1, (0, 6): a},
            3,
        )[(1, 2)],
        a**2,
    )
    for y_order, x_order in ((4, 7), (5, 9)):
        h = sp.symbols(f"double_h_{y_order}")
        k = sp.symbols(f"double_k_{y_order}")
        polynomial = {(1, 5): 1, (2, 6 - y_order): -30 * h}
        operator = {(2, 1): 1, (1, y_order): h, (0, x_order): k}
        first = moment(polynomial, operator, 2)[(0, 2 if y_order == 4 else 0)]
        solved_k = sp.solve(first, k)[0]
        final = moment(polynomial, operator, 3)[
            (0, 3 if y_order == 4 else 0)
        ].subs(k, solved_k)
        assert_multiple(final, h**3)

    # Pure-y weight-two face: radical (c,r,b*p), after eliminating q.
    b, c, p, r, z = sp.symbols("db dc dp dr dz")
    q = -6 * b * p - 180 * c
    op = {(2, 1): 1, (1, 3): b, (0, 5): c}
    pol = {(0, 6): 1, (1, 4): p, (2, 2): q, (3, 0): r}
    equations = primitive_equations(pol, op, 2, 8, (b, c, p, r))
    for target in (c, p, r):
        msolve_empty(
            (z, c, p, r),
            tuple(e.subs(b, 1) for e in equations) + (z * target - 1,),
        )
    for target in (b, c, r):
        msolve_empty(
            (z, b, c, r),
            tuple(e.subs(p, 1) for e in equations) + (z * target - 1,),
        )
    msolve_empty((b, p, r), tuple(e.subs(c, 1) for e in equations))
    assert_multiple(
        moment(
            {(0, 6): 1, (2, 1): -360 * a},
            {(2, 1): 1, (0, 6): a},
            2,
        )[(0, 0)],
        a**2,
    )


def verify_simple_root() -> None:
    A, B = sp.symbols("simple_A simple_B")
    t, l3, l4 = sp.symbols("simple_t simple_l3 simple_l4")
    p = sp.symbols("simple_p0:6")
    polynomial = {
        0: {(0, 6): 1},
        1: {(i, 5 - i): p[i] for i in range(6)},
    }
    operator = {
        0: {(1, 2): 1, (2, 1): A, (3, 0): B},
        1: {(0, 4): t, (3, 1): l3, (4, 0): l4},
    }
    first = defect_moment(polynomial, operator, 1, 1)
    solved_p1 = sp.solve(first[(0, 2)], p[1])[0]
    second = defect_moment(polynomial, operator, 2, 1)
    assert_multiple(second[(3, 2)], p[5])
    assert_multiple(second[(2, 3)].subs(p[5], 0), p[4])
    assert_multiple(second[(1, 4)].subs({p[5]: 0, p[4]: 0}), p[3])
    assert_multiple(
        second[(0, 5)].subs({p[5]: 0, p[4]: 0, p[3]: 0}), p[2]
    )
    assert sp.factor(
        solved_p1.subs({p[2]: 0, p[3]: 0})
    ) == -30 * t
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 4): -30 * t},
            {(1, 2): 1, (0, 4): t},
            2,
        )[(0, 4)],
        t**2,
    )

    # Defects two and three have the same triangular shape.
    s = sp.symbols("simple_s")
    q = sp.symbols("simple_q0:5")
    polynomial2 = {
        0: {(0, 6): 1},
        1: {(0, 5): sp.Symbol("simple_p0_surv")},
        2: {(i, 4 - i): q[i] for i in range(5)},
    }
    operator2 = {
        0: {(1, 2): 1, (2, 1): A, (3, 0): B},
        1: {(3, 1): l3, (4, 0): l4},
        2: {(0, 5): s},
    }
    second_defect = defect_moment(polynomial2, operator2, 2, 2)
    assert_multiple(second_defect[(2, 2)], q[4])
    assert_multiple(second_defect[(1, 3)].subs(q[4], 0), q[3])
    assert_multiple(
        second_defect[(0, 4)].subs({q[4]: 0, q[3]: 0}), q[2]
    )
    first_defect = defect_moment(polynomial2, operator2, 1, 2)
    assert sp.factor(
        sp.solve(first_defect[(0, 1)], q[1])[0].subs(
            {q[2]: 0, q[3]: 0}
        )
    ) == -120 * s
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 3): -120 * s},
            {(1, 2): 1, (0, 5): s},
            3,
        )[(0, 3)],
        s**3,
    )
    u = sp.symbols("simple_u")
    r = sp.symbols("simple_r0:4")
    polynomial3 = {
        0: {(0, 6): 1},
        3: {(i, 3 - i): r[i] for i in range(4)},
    }
    operator3 = {
        0: {(1, 2): 1, (2, 1): A, (3, 0): B},
        3: {(0, 6): u},
    }
    third_defect = defect_moment(polynomial3, operator3, 2, 3)
    assert_multiple(third_defect[(1, 2)], r[3])
    assert_multiple(third_defect[(0, 3)].subs(r[3], 0), r[2])
    first_third = defect_moment(polynomial3, operator3, 1, 3)
    assert sp.factor(
        sp.solve(first_third[(0, 0)], r[1])[0].subs(
            {r[2]: 0, r[3]: 0}
        )
    ) == -360 * u
    assert_multiple(
        moment(
            {(0, 6): 1, (1, 2): -360 * u},
            {(1, 2): 1, (0, 6): u},
            2,
        )[(0, 0)],
        u**2,
    )
    a2 = sp.symbols("simple_a2")
    assert_multiple(
        defect_moment(
            {0: {(0, 6): 1}, 4: {(2, 0): a2}},
            {0: {(1, 2): 1, (2, 1): A, (3, 0): B}},
            2,
            4,
        )[(0, 2)],
        a2,
    )


def main() -> None:
    verify_hall_locus()
    verify_triple_root()
    verify_double_root()
    verify_simple_root()
    print("verified complete binary cubic-leading sextic GVC row")


if __name__ == "__main__":
    main()
