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


def check_zero_e_binary_ratios():
    from sympy import cancel
    x, t, theta, sigma, rho = symbols("x t theta sigma rho")
    z = t*x**3
    first = z**3*(x**3 + theta*z)
    second = x**2*z*(x**3 + sigma*z)
    third = x*z*(x**3 + rho*z)
    ratio_one = t**3*(1 + theta*t)**2/(1 + sigma*t)**3
    ratio_two = (1 + sigma*t)**7/(t*(1 + rho*t)**8)
    assert cancel(first**2/second**3 - ratio_one) == 0
    assert cancel(second**7/third**8 - ratio_two) == 0


def check_generic_fibre():
    x, y, a, h2, h1, t, sigma, delta, theta = symbols(
        "x y a h2 h1 t sigma delta theta"
    )
    z_level = (t - h1*y - h2*y**2)/a
    p_level = x*z_level + y**2
    leading_a = Poly(expand(p_level**3 + delta*y**4*z_level), y).nth(6)
    leading_b = Poly(expand(sigma*(p_level**2 - y**4)), y).nth(4)
    assert expand(leading_a - ((1-x*h2/a)**3-delta*h2/a)) == 0
    assert expand(leading_b - sigma*(x*h2/a)*(x*h2/a-2)) == 0

    # The final signed face is computed after the residue forces a=k*x**4.
    pair = Matrix([theta**3/x**9, sigma*theta**2/x**6 - sigma*y**4])
    assert expand(pair.jacobian((x, y)).det()) == 36*sigma*theta**3*y**3/x**10


def check_zero_k_axial_binomial():
    from sympy import Rational
    x, y, z, e, b, c, s, xi, gamma = symbols("x y z e b c s xi gamma")
    toric = x*z
    phi = b*x*y*z**3 + c*x**2*z**3 + s*y**3*z**2
    norm = (
        x**6 + xi*x*z**2 + gamma*x**2*y*z + Rational(2, 3)*e*b*y*z**2
        - Rational(2, 3)*e*phi*(1/toric - y**2/toric**2)
        + e*phi**2/(9*toric**4)
    )
    terms = expand(norm).as_ordered_terms()
    face = sum(
        term for term in terms
        if sum(w*term.as_powers_dict().get(v, 0)
               for v, w in ((x, 2), (y, 3), (z, 5))) == 12
    )
    expected = (
        x**6 + (xi-Rational(2, 3)*e*c)*x*z**2 + gamma*x**2*y*z
        - Rational(2, 3)*e*(s-b)*y**3*z/x + e*b**2*y**2*z**2/(9*x**2)
    )
    assert expand(face - expected) == 0


def check_both_zero_projective_index():
    from sympy import Rational, factor, rem
    t, theta, sigma, r = symbols("t theta sigma r")
    f = 1 + (3*sigma-theta)*t/2 + Rational(3, 8)*(sigma-theta)**2*t**2
    u = series(t**3*(1+theta*t)/f**2, t, 0, 6).removeO().expand()
    ratio = series((1+sigma*t)**3/((1+theta*t)*f**2), t, 0, 6).removeO().expand()
    leftover = (ratio - ratio.coeff(t, 3)*u).expand()
    quadratic = 9*sigma**2+38*theta*sigma-95*theta**2
    cubic = 27*sigma**3+105*theta*sigma**2-259*theta**2*sigma+31*theta**3
    assert factor(leftover.coeff(t, 4) + (sigma-theta)**2*quadratic/64) == 0
    assert factor(leftover.coeff(t, 5) - (sigma-theta)**2*cubic/64) == 0
    assert rem(cubic.subs({sigma:r, theta:1}), quadratic.subs({sigma:r, theta:1}), r) == 64*r-64


def check_nonzero_axial_slope_fibres():
    from sympy import diff
    x, y, u, theta, rho, sigma, beta, d, a20, a13 = symbols(
        "x y u theta rho sigma beta d a20 a13"
    )

    # Substitution uses the lowest x-order of each scalar coefficient.
    # For a double root and beta != 0 the top of A is its z^2 term.
    z2 = theta/x**2
    p2 = x*z2 + y**2
    b2 = sigma*(p2**2-y**4) + beta*z2
    a2 = a20*z2**2
    det = expand(Matrix([a2, b2]).jacobian((x, y)).det())
    assert det == -16*a20*sigma*theta**3*y/x**6
    a2 = p2**3 + d*x*z2**2
    assert expand(x**3*a2.subs(y, u/x**Rational(1, 2))
                  - ((theta+u**2)**3 + d*theta**2)) == 0
    assert expand(x**2*b2.subs({y:u/x**Rational(1, 2), beta:0})
                  - sigma*(theta**2+2*theta*u**2)) == 0

    # For a triple root the pure y^3 term of A_1 can also tie.
    z3 = (theta-rho*(x*y)**2)/x**3
    p3 = x*z3 + y**2
    a3 = p3**3 + a20*z3**2 + a13*y**3*z3
    b3 = sigma*(p3**2-y**4)
    expected_a = (
        (theta+(1-rho)*u**2)**3
        + a20*(theta-rho*u**2)**2 + a13*u**3*(theta-rho*u**2)
    )
    expected_b = sigma*(theta-rho*u**2)*(theta+(2-rho)*u**2)
    assert expand(x**6*a3.subs(y, u/x) - expected_a) == 0
    assert expand(x**4*b3.subs(y, u/x) - expected_b) == 0

    first, second = Function("U")(u), Function("V")(u)
    for exponent, denominator, slope in ((3, 2, Rational(1, 2)), (6, 4, 1)):
        pair = Matrix([first/x**exponent, second/x**denominator])
        pair = pair.subs(u, x**slope*y)
        bracket = pair.jacobian((x, y)).det()
        expected = (
            denominator*diff(first, u)*second - exponent*first*diff(second, u)
        ).subs(u, x**slope*y)
        expected *= x**(slope-1-exponent-denominator)
        assert expand(bracket - expected) == 0


if __name__ == "__main__":
    check_first_axial_weight()
    check_fractional_coefficient()
    check_parallel_plane()
    check_zero_e_binary_ratios()
    check_generic_fibre()
    check_zero_k_axial_binomial()
    check_both_zero_projective_index()
    check_nonzero_axial_slope_fibres()
    print("Checked axial supports, binomial expansion, fibre brackets, and projective field ratios.")
