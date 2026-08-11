#!/usr/bin/env python3
"""Verify the algebraic k=1 loci used by the Chau tangency obstruction.

The external topological input is Nguyen Van Chau's theorem that the
exceptional-value set of a nonsingular polynomial map of C^2 cannot, in one
common affine target coordinate, be a union of curves parametrized by
``(t^n, q(t))``.  This script checks the exact normal-form reductions to that
theorem and isolates the first target-image merger outside its scope.
"""

from __future__ import annotations

import sympy as sp


t, s = sp.symbols("t s")
a, b, c, d, lam, h = sp.symbols("a b c d lam h")
P, Q = sp.symbols("P Q")


def pure_projection_loci() -> None:
    """Classify pure-power linear projections of the k=1 normal form."""

    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t

    # The degree-three target directions are precisely the multiples of p.
    # Its two finite critical points merge exactly on a=0, when p=t^3.
    assert sp.discriminant(sp.diff(p, t), t) == -12 * a
    assert sp.expand(p.subs(a, 0) - t**3) == 0

    # Every degree-five linear target direction is q+lam*p after scaling.
    # Matching it with (t-h)^5 plus a constant gives the complete locus.
    pure_fifth = (t - h) ** 5 + h**5
    difference = sp.Poly(sp.expand(q + lam * p - pure_fifth), t)
    equations = [difference.coeff_monomial(t**power) for power in range(1, 5)]
    solution = {
        h: -b / 5,
        lam: 2 * b**2 / 5,
        c: 2 * b**3 / 25,
        d: b**4 / 125 - 2 * a * b**2 / 5,
    }
    assert all(sp.factor(equation.subs(solution)) == 0 for equation in equations)

    # Conversely the coefficient comparison is triangular: t^4 fixes h,
    # t^3 fixes lam, t^2 fixes c, and t fixes d.
    solved = sp.solve(equations, (h, lam, c, d), dict=True)
    assert len(solved) == 1
    assert all(sp.factor(solved[0][key] - solution[key]) == 0 for key in solution)

    pure_fifth_conditions = {
        c: 2 * b**3 / 25,
        d: b**4 / 125 - 2 * a * b**2 / 5,
        lam: 2 * b**2 / 5,
    }
    shifted = sp.expand(
        (q + lam * p).subs(pure_fifth_conditions).subs(t, s - b / 5)
    )
    assert sp.diff(shifted, s) == 5 * s**4


def named_closed_faces() -> None:
    """Check the E6, E8, and cusp/triple witnesses on the cubic locus."""

    witnesses = {
        "E6+A1": (t**3, t**5 + t**4),
        "E8": (t**3, t**5),
        "A2+ordinary-triple": (t**3, t**5 + t**2),
    }
    for p, _q in witnesses.values():
        assert p == t**3

    # The E8 curve is also explicitly the image P^5-Q^3=0.
    assert sp.expand((P**5 - Q**3).subs({P: t**3, Q: t**5})) == 0


def triple_image_merger_frontier() -> None:
    """Derive the merger hypersurface and one witness beyond the Chau locus."""

    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    value = sp.symbols("value")

    # Reduce q modulo p-value.  Three distinct points of a p-fibre have the
    # same q-value exactly when the t^2 and t coefficients vanish.
    remainder = sp.rem(q, p - value, t)
    expected = (
        (value - a * b + c) * t**2
        + (a**2 + b * value + d) * t
        - a * value
    )
    assert sp.expand(remainder - expected) == 0
    triple_value = a * b - c
    merger_equation = sp.factor(
        (a**2 + b * value + d).subs(value, triple_value)
    )
    assert merger_equation == a**2 + a * b**2 - b * c + d

    # The first a!=0 witness is p=t^3+t, q=t^5-t.  Its fibre over P=0
    # consists of 0,+i,-i and all three map to (0,0).  It lies outside both
    # pure linear-projection loci, so Chau's theorem does not decide it.
    witness = {a: 1, b: 0, c: 0, d: -1}
    assert merger_equation.subs(witness) == 0
    assert witness[a] != 0
    assert (25 * c - 2 * b**3).subs(witness) == 0
    assert (125 * d + 50 * a * b**2 - b**4).subs(witness) != 0

    witness_p = sp.expand(p.subs(witness))
    witness_q = sp.expand(q.subs(witness))
    assert sp.factor(witness_p) == t * (t**2 + 1)
    assert all(
        sp.expand(poly.subs(t, root)) == 0
        for root in (0, sp.I, -sp.I)
        for poly in (witness_p, witness_q)
    )

    implicit = sp.factor(
        sp.resultant(P - witness_p, Q - witness_q, t)
    )
    expected_implicit = -P**5 + 4 * P**3 + 8 * P**2 * Q + 5 * P * Q**2 + Q**3
    assert implicit == expected_implicit
    assert sp.factor(P**5 - (Q + P) * (Q + 2 * P) ** 2) == -expected_implicit
    assert sp.factor(sp.discriminant(expected_implicit, Q)) == -P**8 * (27 * P**2 + 4)


def main() -> None:
    pure_projection_loci()
    named_closed_faces()
    triple_image_merger_frontier()
    print(
        "PASS: the k=1 pure cubic/quintic projection loci are exact; "
        "E6, E8, and the a=0 triple merger fall under Chau's theorem; "
        "the first a!=0 triple-image witness is isolated explicitly"
    )


if __name__ == "__main__":
    main()
