"""Exact support and binomial checks for the mixed square-cusp reduction."""

from itertools import product

from sympy import Function, Matrix, Poly, Rational, expand, series, symbols


def check_first_axial_weight():
    # On this edge the nonreference terms have strictly lower (2,3,4)-weight.
    degrees = (18, 14, 12)
    bounds = (5, 5, 5)
    for index, (degree, bound) in enumerate(zip(degrees, bounds)):
        allowed = []
        for i, j, l in product(range(6), repeat=3):
            if i + j + l > bound or 2*i + 3*j + 4*l >= degree:
                continue
            allowed.append((2*i + 3*j + 5*l, (i, j, l)))
        maximum = max(weight for weight, _ in allowed)
        top = {monomial for weight, monomial in allowed if weight == maximum}
        expected = (
            (20, {(0, 0, 4), (1, 1, 3)}),
            (15, {(0, 0, 3), (1, 1, 2)}),
            (13, {(0, 1, 2)}),
        )
        assert (maximum, top) == expected[index]


def check_fractional_coefficient():
    x, y, eps = symbols("x y eps")
    b, c, d, f, g = symbols("b c d f g")
    ratio = 1 + eps*b*y/x**2 + eps**2*c/x + eps**3*d*y/x**3
    ratio += eps**4*f/x**2 + eps**6*g/x**3
    power = expand(x**2 * series(ratio**Rational(2, 3), eps, 0, 4).removeO())
    assert power.coeff(eps, 1) == 2*b*y/3
    assert expand(power.coeff(eps, 2) - (2*c*x/3 - b**2*y**2/(9*x**2))) == 0
    assert power.coeff(eps, 3).subs(b, 0) == 2*d*y/(3*x)

    # Centered cubic and its matched monic quadratic at the last axial stage.
    first, second = symbols("f1 g1")
    cubic = x**3 + first*x + second
    quadratic = x**2 + 2*first/3
    difference = Poly(cubic**2 - quadratic**3, x)
    assert difference.nth(3) == 2*second
    assert difference.nth(2) == -first**2/3


def check_parallel_plane():
    y, z, a, b = symbols("y z a b")
    p2, p1, p0, q = (Function(name)(y) for name in ("p2", "p1", "p0", "q"))
    pair = Matrix([a*z**3 + p2*z**2 + p1*z + p0, b*z**2 + q])
    determinant = Poly(pair.jacobian((y, z)).det(), z)
    assert determinant.nth(0) == -p1*q.diff(y)
    assert determinant.nth(2) == 2*b*p1.diff(y) - 3*a*q.diff(y)


if __name__ == "__main__":
    check_first_axial_weight()
    check_fractional_coefficient()
    check_parallel_plane()
    print("Checked the first axial supports, fractional coefficients and parallel-plane bracket.")
