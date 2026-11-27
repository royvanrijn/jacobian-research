"""Check the axial support and determinant coefficients of the pure frame."""

from fractions import Fraction as Q
from itertools import product

from sympy import expand, symbols

from check_sextic_doubled_square_third_mixed import (
    competing_monomials, cut, jacobian, weight,
)


def check_weights():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    initial = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    pure = cut((monomials,) * 3, (2, 3, 4), initial)
    reference = ((4, 0, 2), (5, 0, 1), (6, 0, 0))
    assert all(
        weight(m, (1, 1, 2)) <= weight(ref, (1, 1, 2))
        for allowed, ref in zip(pure, reference) for m in allowed
    )
    assert competing_monomials(pure, (1, 1, 0), 2, 2, 3, reference) == (
        {3: {(1, 0, 3), (0, 1, 3)}},
        {Q(5, 2): {(0, 0, 3)}, 3: {(2, 0, 2), (1, 1, 2)}},
        {Q(5, 2): {(1, 0, 2), (0, 1, 2)}, 3: {(0, 0, 2), (3, 0, 1), (2, 1, 1)}},
    )


def check_coefficients():
    x, y, z = symbols("x y z")
    a, b, c, d, e, f, g = symbols("a b c d e f g")
    determinant = jacobian(
        (x**4*z**2+(a*x+b*y)*z**3,
         x**5*z+(c*x**2+d*x*y)*z**2,
         x**6+e*z**2+(f*x**3+g*x**2*y)*z),
        (x, y, z),
    )
    coefficient = determinant.coeff_monomial
    equations = [
        coefficient((10, 0, 3)) - 6*(b-2*d+g),
        coefficient((6, 1, 4)) + 6*(b*d-2*b*g+d*g),
        coefficient((0, 1, 6)).subs({d: b, g: b}) + 2*b**2*e,
        coefficient((7, 0, 4)).subs({d: b, g: b}) + 4*b*(a-2*c+f),
        coefficient((4, 0, 5)).subs({d: b, g: b, e: 0}) - 4*b*(a*c-2*a*f+c*f),
    ]
    assert all(expand(equation) == 0 for equation in equations)


if __name__ == "__main__":
    check_weights()
    check_coefficients()
    print("Uncorrected doubled-square frame: axial support and determinants verified.")
