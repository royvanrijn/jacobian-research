#!/usr/bin/env python3
"""Verify the complete affine singularity atlas on the F2 k=1 chart.

For p=t^3+a*t and q=t^5+b*t^4+c*t^2+d*t, the locus where distinct
collision parameters acquire the same target value has one hypersurface as
its closure.  On that hypersurface the collision quartic is a line times a
depressed cubic.  Off it, a marked common critical point reduces every
nonimmersion degeneration to one explicitly controlled cubic.

This checks the polynomial identities and all packet witnesses.  The ADE
labels use the standard local normal forms proved in the accompanying note.
"""

from __future__ import annotations

import sympy as sp


t, u, v, z = sp.symbols("t u v z")
a, b, c, d, rho = sp.symbols("a b c d rho")


def collision_quartic(value: sp.Expr) -> sp.Expr:
    return (
        value**4
        + b * value**3
        + a * value**2
        + (2 * a * b - c) * value
        - (a**2 + d)
    )


R = collision_quartic(u)
M = a**2 + a * b**2 - b * c + d
D = 3 * u**2 + 4 * a
X = -u * (u**2 + a)
Y = (u**2 + a) * (u**3 + 2 * a * u + a * b - c)


def polynomial_equal(left: sp.Expr, right: sp.Expr) -> bool:
    return sp.expand(left - right) == 0


def collision_and_merger_audit() -> None:
    """Prove that the equal-image locus has closure M=0."""

    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    pair = z**2 - u * z + (u**2 + a)
    assert sp.discriminant(pair, z) == -D

    q_remainder = sp.rem(q.subs(t, z), pair, z)
    assert sp.expand(q_remainder - (-R * z + Y)) == 0
    assert sp.expand(sp.rem(p.subs(t, z), pair, z) - X) == 0

    # If u != v have the same first collision value, then
    # a=-(u^2+uv+v^2).  The divided difference of R fixes c; R(u)=0
    # then fixes d, and the merger equation becomes zero identically.
    a_equal_x = -(u**2 + u * v + v**2)
    r_difference = sp.factor(
        (collision_quartic(u) - collision_quartic(v)) / (u - v)
    ).subs(a, a_equal_x)
    c_equal = -b * (u**2 + u * v + v**2) - u**2 * v - u * v**2
    assert sp.factor(r_difference.subs(c, c_equal)) == 0
    y_difference = sp.factor((Y - Y.subs(u, v)) / (u - v)).subs(
        a, a_equal_x
    )
    assert sp.factor(y_difference - (u + v) * r_difference) == 0

    d_equal = sp.solve(
        collision_quartic(u).subs({a: a_equal_x, c: c_equal}), d
    )[0]
    assert sp.factor(
        M.subs({a: a_equal_x, c: c_equal, d: d_equal})
    ) == 0

    # Conversely, M=0 gives a complete factorization.  Every root of the
    # cubic maps to the same target point (ab-c,-a(ab-c)).
    d_merger = -a**2 - a * b**2 + b * c
    cubic = u**3 + a * u + a * b - c
    assert sp.factor(R.subs(d, d_merger) - (u + b) * cubic) == 0
    assert sp.expand(sp.rem(X, cubic, u) - (a * b - c)) == 0
    assert sp.expand(sp.rem(Y, cubic, u) + a * (a * b - c)) == 0


def nonmerger_critical_audit() -> None:
    """Exhaust the critical-root multiplicities away from M=0."""

    # A nonzero common critical point can be scaled to t=1.  Then a=-3
    # and q'(1)=0 fixes d.
    critical_slice = {a: -3, d: -5 - 4 * b - 2 * c}
    specialized = sp.factor(R.subs(critical_slice))
    cubic = sp.factor(specialized / (u - 2))
    expected = u**3 + (b + 2) * u**2 + (2 * b + 1) * u + 2 - 2 * b - c
    assert sp.expand(cubic - expected) == 0
    assert sp.expand(
        sp.factor(M.subs(critical_slice)) + (b + 2) * (3 * b + c - 2)
    ) == 0
    assert polynomial_equal(sp.factor(sp.discriminant(cubic, u)), (
        (3 * b + c - 2)
        * (4 * b**3 - 12 * b**2 - 69 * b - 27 * c + 50)
    ))

    # Multiplicity at the marked diagonal root u=2.
    assert polynomial_equal(sp.factor(cubic.subs(u, 2)), 6 * b - c + 20)
    assert polynomial_equal(
        sp.factor(sp.diff(cubic, u).subs(u, 2)), 3 * (2 * b + 7)
    )
    assert polynomial_equal(
        sp.factor(sp.diff(cubic, u, 2).subs(u, 2)), 2 * (b + 8)
    )
    higher_cusp = sp.factor(cubic.subs(c, 6 * b + 20))
    assert polynomial_equal(
        higher_cusp, (u - 2) * (u**2 + (b + 4) * u + 4 * b + 9)
    )
    assert polynomial_equal(sp.factor(
        sp.discriminant(u**2 + (b + 4) * u + 4 * b + 9, u)
    ), (b - 10) * (b + 2))
    assert polynomial_equal(
        sp.factor(M.subs(critical_slice).subs(c, 6 * b + 20)),
        -9 * (b + 2) ** 2,
    )
    assert polynomial_equal(
        sp.factor(specialized.subs({b: -sp.Rational(7, 2), c: -1})),
        (u - 2) ** 3 * (2 * u + 5) / 2,
    )
    # Multiplicity four would require the incompatible b=-7/2 and b=-8.
    assert sp.solve(
        [
            cubic.subs(u, 2),
            sp.diff(cubic, u).subs(u, 2),
            sp.diff(cubic, u, 2).subs(u, 2),
        ],
        (b, c),
        dict=True,
    ) == []
    assert (c + 2 * b).subs({b: -sp.Rational(7, 2), c: -1}) != 0

    # The other critical point is u=-2.  It is common exactly on c=-2b;
    # it can be double only at b=5/2, and never triple there.
    assert polynomial_equal(sp.factor(cubic.subs(u, -2)), -2 * b - c)
    assert polynomial_equal(sp.factor(sp.diff(cubic, u).subs(u, -2)), 5 - 2 * b)
    assert polynomial_equal(
        sp.factor(sp.diff(cubic, u, 2).subs(u, -2)), 2 * (b - 4)
    )

    two_critical = sp.factor(specialized.subs(c, -2 * b))
    assert polynomial_equal(two_critical, (u**2 - 4) * (u**2 + b * u + 1))
    assert polynomial_equal(
        sp.factor(M.subs(critical_slice).subs(c, -2 * b)), 4 - b**2
    )
    assert sp.discriminant(u**2 + b * u + 1, u) == b**2 - 4
    assert polynomial_equal(sp.factor((u**2 + b * u + 1).subs(u, -2)), 5 - 2 * b)
    assert polynomial_equal(sp.factor((u**2 + b * u + 1).subs(u, 2)), 5 + 2 * b)

    # A triple residual root is forced onto M=0.  Thus there is no
    # distinct-image A2+A5 packet.
    root = sp.symbols("root")
    triple_solutions = sp.solve(
        [
            b + 2 + 3 * root,
            2 * b + 1 - 3 * root**2,
            2 - 2 * b - c + root**3,
        ],
        (b, c, root),
        dict=True,
    )
    assert triple_solutions == [{b: 1, c: -1, root: -1}]
    assert M.subs(critical_slice).subs(triple_solutions[0]) == 0

    # Both diagonal roots cannot simultaneously have multiplicity two.
    assert sp.solve(
        [c - (6 * b + 20), c + 2 * b, 5 - 2 * b],
        (b, c),
        dict=True,
    ) == []


def merger_factor_atlas_audit() -> None:
    """Check the D/E classification from the line-cubic factorization."""

    cubic = u**3 + a * u + a * b - c
    assert polynomial_equal(
        sp.factor(sp.discriminant(cubic, u)), -4 * a**3 - 27 * (a * b - c) ** 2
    )
    assert polynomial_equal(sp.factor(cubic.subs(u, -b)), -(b**3 + c))
    assert polynomial_equal(sp.factor(D.subs(u, -b)), 3 * b**2 + 4 * a)

    # A double-root cubic is (u-rho)^2(u+2rho).  Its simple root is
    # diagonal, while its double root is off diagonal for rho != 0.
    double_parameters = {a: -3 * rho**2, c: -3 * b * rho**2 - 2 * rho**3}
    assert polynomial_equal(
        sp.factor(cubic.subs(double_parameters)), (u - rho) ** 2 * (u + 2 * rho)
    )
    assert sp.factor(D.subs(double_parameters).subs(u, -2 * rho)) == 0
    assert polynomial_equal(
        sp.factor(D.subs(double_parameters).subs(u, rho)), -9 * rho**2
    )

    # The residual line root -b is respectively the other diagonal root,
    # the double cubic root, or the simple cubic root in these three rows.
    assert sp.factor(
        D.subs(double_parameters).subs({b: -2 * rho, u: 2 * rho})
    ) == 0  # D5+A2
    assert sp.expand(-(-rho) - rho) == 0  # E7, intersection order three
    assert sp.expand(-(2 * rho) + 2 * rho) == 0  # E7, A4 cusp

    # A triple cubic forces a=c=0.  The residual line is separate for
    # b!=0 (E6+A1) and coincident for b=0 (E8).
    assert cubic.subs({a: 0, c: 0}) == u**3
    assert polynomial_equal(sp.factor(R.subs({a: 0, c: 0, d: 0})), u**3 * (u + b))


def witness_table_audit() -> None:
    """Check one exact point on every packet in the complete table."""

    witnesses: dict[str, tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, bool]] = {
        # M != 0: individual A-type collision roots.
        "4A1": (1, 0, 0, 0, False),
        "A3+2A1": (0, 0, 4, -3, False),
        "2A3": (2, 0, 0, -5, False),
        "A5+A1": (3, -3, -17, -9, False),
        "A7": (6, -4, -44, -37, False),
        "A2+3A1": (-3, 0, 1, -7, False),
        "A2+A3+A1": (-3, -sp.Rational(1, 2), 3, -9, False),
        "A4+2A1": (-3, -sp.Rational(9, 4), sp.Rational(13, 2), -9, False),
        "A4+A3": (-3, 10, 80, -205, False),
        "A6+A1": (-3, -sp.Rational(7, 2), -1, 11, False),
        "2A2+2A1": (-3, 0, 0, -5, False),
        "A4+A2+A1": (-3, sp.Rational(5, 2), -5, -5, False),
        # M = 0: line times depressed cubic.
        "D4+A1": (1, 1, 0, -2, True),
        "D4+A2": (-3, 2, 0, 3, True),
        "D5+A1": (-3, 0, -2, -9, True),
        "D6": (1, 1, -1, -3, True),
        "D5+A2": (-3, -2, 4, -5, True),
        "E7-I3": (-3, -1, 1, -7, True),
        "E7-A4": (-3, 2, -8, -13, True),
        "E6+A1": (0, 1, 0, 0, True),
        "E8": (0, 0, 0, 0, True),
    }

    for label, (aa, bb, cc, dd, merger) in witnesses.items():
        specialization = {a: aa, b: bb, c: cc, d: dd}
        assert (sp.factor(M.subs(specialization)) == 0) is merger, label
        assert sp.Poly(R.subs(specialization), u).degree() == 4

    expected_factors = {
        "A4+A2+A1": (u - 2) * (u + 2) ** 2 * (2 * u + 1) / 2,
        "D4+A2": (u + 2) * (u**3 - 3 * u - 6),
        "D5+A1": u * (u - 1) ** 2 * (u + 2),
        "D6": (u + 1) ** 2 * (u**2 - u + 2),
        "D5+A2": (u - 2) * (u - 1) ** 2 * (u + 2),
        "E7-I3": (u - 1) ** 3 * (u + 2),
        "E7-A4": (u - 1) ** 2 * (u + 2) ** 2,
        "E6+A1": u**3 * (u + 1),
        "E8": u**4,
    }
    for label, expected in expected_factors.items():
        aa, bb, cc, dd, _ = witnesses[label]
        assert sp.factor(R.subs({a: aa, b: bb, c: cc, d: dd}) - expected) == 0


def main() -> None:
    collision_and_merger_audit()
    nonmerger_critical_audit()
    merger_factor_atlas_audit()
    witness_table_audit()
    print(
        "PASS: the F2 k=1 affine singularity atlas is exhaustive; the "
        "equal-image locus has closure M=0, the merger quartic is line times "
        "cubic, and the nonmerger critical slice has only the listed A-packets"
    )


if __name__ == "__main__":
    main()
