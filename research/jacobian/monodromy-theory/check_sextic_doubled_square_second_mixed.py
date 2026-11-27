"""Check the finite weight and determinant calculations for the k branch."""

from fractions import Fraction as Q
from itertools import product

from sympy import expand, symbols

from check_sextic_doubled_square_third_mixed import (
    competing_monomials, cut, jacobian, weight,
)


P_SQUARED = {(4, 0, 2), (3, 2, 1), (2, 4, 0)}
X_P_SQUARED = {(3, 0, 2), (2, 2, 1), (1, 4, 0)}
X_SIXTH = {(6, 0, 0)}
XY_REFERENCE = ((2, 4, 0), (1, 4, 0), (6, 0, 0))


def lower_x_ties(supports, y_weight, z_weight, x_start, reference):
    """Enumerate ties strictly inside the indicated interval."""
    return tuple(
        {s: monomials for s, monomials in rows.items() if s < x_start}
        for rows in competing_monomials(
            supports, (0, y_weight, z_weight), 0, 0, x_start, reference,
        )
    )


def check_weights():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    original = cut(
        (monomials,) * 3, (2, 3, 4),
        (P_SQUARED, X_P_SQUARED | {(5, 0, 1), (4, 2, 0)}, X_SIXTH),
    )
    assert all(
        weight(m, (2, 3, 4)) <= weight(ref, (2, 3, 4))
        for allowed, ref in zip(original, XY_REFERENCE) for m in allowed
    )
    assert lower_x_ties(original, 3, 4, 2, XY_REFERENCE) == (
        {Q(3, 2): {(0, 1, 3), (0, 5, 0)}, 1: {(0, 2, 2), (1, 3, 1)},
         Q(1, 2): {(0, 3, 1)}},
        {1: {(0, 3, 1)}},
        {
            Q(11, 6): {(0, 1, 2)}, Q(9, 5): {(1, 3, 0)},
            Q(7, 4): {(2, 1, 1)}, Q(5, 3): {(0, 2, 1)},
            Q(8, 5): {(1, 0, 2)}, Q(3, 2): {(0, 3, 0), (2, 2, 0), (4, 1, 0)},
            Q(7, 5): {(1, 1, 1)}, Q(4, 3): {(0, 0, 2), (3, 0, 1)},
            Q(6, 5): {(1, 2, 0)}, Q(7, 6): {(0, 1, 1)},
            1: {(0, 2, 0), (2, 0, 1), (3, 1, 0)},
            Q(4, 5): {(1, 0, 1)}, Q(3, 4): {(2, 1, 0)},
            Q(2, 3): {(0, 0, 1)}, Q(3, 5): {(1, 1, 0)}, Q(1, 2): {(0, 1, 0)},
        },
    )

    quadratic_two = cut(
        original, (1, 2, 3),
        (P_SQUARED, X_P_SQUARED, X_SIXTH | {(3, 0, 1), (2, 2, 0)}),
    )
    assert lower_x_ties(quadratic_two, 2, 3, 1, (*XY_REFERENCE[:2], (2, 2, 0))) == (
        {Q(1, 2): {(0, 0, 3), (0, 3, 1)}}, {}, {Q(1, 2): {(0, 1, 1)}},
    )

    quadratic_one = cut(
        original, (6, 15, 24),
        (P_SQUARED, X_P_SQUARED, X_SIXTH | {(2, 0, 1), (1, 2, 0)}),
    )
    assert lower_x_ties(
        quadratic_one, 3, Q(24, 5), Q(6, 5), (*XY_REFERENCE[:2], (1, 2, 0)),
    ) == (
        {Q(9, 10): {(0, 3, 1)}, Q(3, 5): {(1, 1, 2)}, Q(3, 10): {(0, 1, 2)}},
        {Q(3, 5): {(0, 1, 2)}},
        {},
    )

    shear = cut(
        original, (1, 3, 4),
        ({(2, 4, 0), (1, 3, 1)}, {(1, 4, 0), (0, 3, 1)}, X_SIXTH),
    )
    excluded_third = {
        (0, 1, 2), (1, 3, 0), (2, 1, 1), (0, 2, 1), (1, 0, 2),
        (0, 3, 0), (2, 2, 0), (4, 1, 0), (1, 1, 1), (0, 0, 2),
        (3, 0, 1), (1, 2, 0), (0, 1, 1),
    }
    shear = (*shear[:2], shear[2] - excluded_third)
    assert lower_x_ties(shear, 3, 4, 1, ((1, 3, 1), (0, 3, 1), (6, 0, 0))) == (
        {}, {},
        {Q(4, 5): {(1, 0, 1)}, Q(3, 4): {(2, 1, 0)}, Q(2, 3): {(0, 0, 1)},
         Q(3, 5): {(1, 1, 0)}, Q(1, 2): {(0, 1, 0)}},
    )


def check_determinants():
    x, y, z = symbols("x y z")
    a, b, c, d, e, f = symbols("a b c d e f")
    at_one = jacobian(
        (x**2*y**4+a*y**2*z**2+b*x*y**3*z,
         x*y**4+c*y**3*z,
         x**6+d*y**2+e*x**2*z+f*x**3*y), (x, y, z),
    )
    equations = [
        at_one.coeff_monomial((7, 6, 0)) + 24*(b-c),
        at_one.coeff_monomial((6, 5, 1)) + 48*a,
        at_one.coeff_monomial((4, 7, 0)).subs(b, c) - (4*e-c*f),
        at_one.coeff_monomial((2, 5, 2)).subs({b: c, a: 0}) - 3*c**2*e,
        at_one.coeff_monomial((1, 8, 0)).subs(b, c) + 2*c*d,
    ]
    assert all(expand(equation) == 0 for equation in equations)
    at_half = jacobian(
        (x**2*y**4+a*z**3+b*y**3*z, x*y**4, x**2*y**2+e*y*z), (x, y, z),
    )
    assert at_half.coeff_monomial((2, 5, 2)) == -18*a
    assert expand(at_half.coeff_monomial((2, 8, 0)) + 6*b - 4*e) == 0
    assert at_half.coeff_monomial((0, 7, 1)) == -2*b*e


if __name__ == "__main__":
    check_weights()
    check_determinants()
    print("Doubled-square second mixed correction: support and determinants verified.")
