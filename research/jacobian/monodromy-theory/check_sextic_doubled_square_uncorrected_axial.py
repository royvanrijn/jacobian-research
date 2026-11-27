"""Check the axial support and determinant coefficients of the pure frame."""

from fractions import Fraction as Q
from itertools import product

from sympy import Function, Poly, diff, expand, symbols

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


def check_y_lifts():
    expected = {
        1: {},
        2: {2: {(0, 1, 0)}},
        3: {Q(3, 2): {(0, 2, 0)}, 2: {(1, 1, 0)}},
        4: {Q(4, 3): {(0, 3, 0)}, Q(3, 2): {(1, 2, 0)},
            2: {(0, 2, 0), (2, 1, 0)}},
        5: {Q(5, 3): {(0, 3, 0)}, 2: {(1, 2, 0), (0, 1, 1)}},
    }
    for target, desired in expected.items():
        possible = {}
        for i, j, k in product(range(4), repeat=3):
            if j and i+j+k <= 3:
                y_weight = Q(target-i-3*k, j)
                if 1 < y_weight <= 2:
                    possible.setdefault(y_weight, set()).add((i, j, k))
        assert possible == desired


def check_veronese_fibre():
    x, y, z, u, v = symbols("x y z u v")
    c, mu, nu = symbols("c mu nu")
    U = x**2*z + x*y**2 + mu*y*z
    V = x**3 + c*z + nu*x*y
    D = c-mu*nu
    K = mu*v-(mu+nu)*x**3
    rhs = x**2*(v-x**3)+K*y+D*x*y**2
    assert expand((c*U).subs(z, (v-x**3-nu*x*y)/c) - rhs) == 0
    Q = expand(K**2 + 4*D*(x**6-v*x**3+c*u*x))
    assert expand(Q.coeff(x, 6) - ((mu-nu)**2+4*c)) == 0
    assert expand(Q.coeff(x, 3) + 2*(mu**2-mu*nu+2*c)*v) == 0
    exceptional = {nu: -mu, c: -mu**2}
    assert expand((c*U-(x**2+mu*y)*V+x**5).subs(exceptional)) == 0
    assert expand((jacobian((U, V, x), (x, y, z)).as_expr()-mu*V).subs(exceptional)) == 0


def check_isolated_third_correction():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    start = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    pure = cut((monomials,) * 3, (2, 3, 4), start)
    top = (*start[:2], start[2] | {(2, 2, 0), (3, 0, 1)})
    corrected = cut(pure, (1, 2, 3), top)
    ties = competing_monomials(
        corrected, (0, 2, 3), 0, 0, 1, ((2, 4, 0), (4, 2, 0), (2, 2, 0)),
    )
    assert tuple({k: v for k, v in rows.items() if k < 1} for rows in ties) == (
        {Q(1, 2): {(0, 0, 3), (0, 3, 1)}},
        {Q(3, 4): {(0, 2, 1)}, Q(2, 3): {(1, 0, 2), (1, 3, 0)},
         Q(1, 2): {(0, 0, 2), (0, 3, 0), (2, 1, 1)},
         Q(1, 3): {(1, 1, 1)}, Q(1, 4): {(0, 1, 1)}},
        {Q(1, 2): {(0, 1, 1)}},
    )
    x, y, z = symbols("x y z")
    a, b, p, d, e, ell = symbols("a b p d e ell")
    determinant = jacobian(
        (x**2*y**4 + a*z**3 + b*y**3*z,
         x**4*y**2 + p*z**2 + d*y**3 + e*x**2*y*z,
         x**2*y**2 + ell*y*z), (x, y, z),
    )
    coef = determinant.coeff_monomial
    e_value = 3*ell-b
    p_value = 3*ell**2-3*b*ell+b**2
    equations = [
        coef((5, 3, 2))-12*a,
        coef((5, 6, 0))-4*(b+e-3*ell),
        coef((3, 5, 1)).subs(e, e_value)-8*(p-p_value),
        coef((1, 4, 2)).subs({a: 0, e: e_value, p: p_value})-12*(b-ell)**3,
    ]
    assert all(expand(equation) == 0 for equation in equations)


def check_isolated_second_correction():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    start = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    square = {(2, 0, 2), (1, 2, 1), (0, 4, 0)}
    pure = cut((monomials,) * 3, (2, 3, 4), start)
    axial = cut(pure, (1, 2, 3), (start[0], start[1] | square, start[2]))
    ties = competing_monomials(
        axial, (0, 2, 3), 0, 0, 1, ((2, 4, 0), (0, 4, 0), (6, 0, 0)),
    )
    assert tuple({s: terms for s, terms in rows.items() if s < 1} for rows in ties) == (
        {Q(1, 2): {(0, 0, 3), (0, 3, 1)}}, {},
        {Q(5, 6): {(0, 1, 1)}, Q(4, 5): {(1, 2, 0)}, Q(3, 4): {(2, 0, 1)},
         Q(2, 3): {(0, 2, 0), (3, 1, 0)}, Q(3, 5): {(1, 0, 1)},
         Q(1, 2): {(0, 0, 1), (2, 1, 0)}, Q(2, 5): {(1, 1, 0)},
         Q(1, 3): {(0, 1, 0)}},
    )
    refined = cut(
        axial, (4, 10, 16), (start[0], square, start[2] | {(1, 2, 0), (2, 0, 1)}),
    )
    ties = competing_monomials(
        refined, (0, 2, Q(16, 5)), 0, 0, Q(4, 5),
        ((2, 4, 0), (0, 4, 0), (1, 2, 0)),
    )
    assert tuple({s: terms for s, terms in rows.items() if s < Q(4, 5)}
                 for rows in ties) == (
        {Q(3, 5): {(0, 3, 1)}, Q(2, 5): {(1, 1, 2)}, Q(1, 5): {(0, 1, 2)}},
        {}, {},
    )


def check_simultaneous_corrections():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    start = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    pure = cut((monomials,) * 3, (2, 3, 4), start)
    axial = cut(
        pure, (1, 2, 3),
        (start[0], start[1] | {(2, 0, 2), (1, 2, 1), (0, 4, 0)},
         start[2] | {(3, 0, 1), (2, 2, 0)}),
    )
    ties = competing_monomials(
        axial, (0, 2, 3), 0, 0, 1, ((2, 4, 0), (0, 4, 0), (2, 2, 0)),
    )
    assert tuple({s: terms for s, terms in rows.items() if s < 1} for rows in ties) == (
        {Q(1, 2): {(0, 0, 3), (0, 3, 1)}}, {},
        {Q(1, 2): {(0, 1, 1)}},
    )
    x, y, z = symbols("x y z")
    a, b, ell, u, v, w, q = symbols("a b ell u v w q")
    A = x**2*y**4+a*z**3+b*y**3*z
    C = x**2*y**2+ell*y*z
    assert expand(
        diff(A,x)*diff(C,z)-diff(A,z)*diff(C,x)
        -2*x*y**2*((ell-b)*y**3-3*a*z**2)
    ) == 0
    Y = y+v/ell
    Z = ell*z+q*Y+u-2*q*v/ell
    assert expand(Y*Z+w-u*v/ell+q*v**2/ell**2-(q*y**2+u*y+w+(ell*y+v)*z)) == 0


def check_axial_four_reduction():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    first = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    ref = ((4, 0, 2), (5, 0, 1), (6, 0, 0))
    pure = cut((monomials,) * 3, (2, 3, 4), first)
    through_three = cut(pure, (1, 1, 3), tuple({m} for m in ref))
    assert competing_monomials(through_three, (1, 1, 0), 2, 3, 4, ref) == (
        {4: {(0, 0, 3)}},
        {4: {(1, 0, 2), (0, 1, 2)}},
        {4: {(2, 0, 1), (1, 1, 1), (0, 2, 1)}},
    )

    x, y, z = symbols("x y z")
    a, b, beta, c, d, e = symbols("a b beta c d e")
    determinant = jacobian(
        (x**4*z**2+a*z**3, x**5*z+(b*x+beta*y)*z**2,
         x**6+(c*x**2+d*x*y+e*y**2)*z), (x, y, z),
    )
    expected = {
        (9, 0, 3): 6*(d-2*beta),
        (8, 1, 3): 12*e,
        (4, 1, 4): 6*(5*a*e-2*b*e-beta*d),
    }
    for monomial, value in expected.items():
        assert expand(determinant.coeff_monomial(monomial)-value) == 0

    lifts = {
        1: {},
        3: {Q(3, 2): {(0, 2, 0)}, 2: {(1, 1, 0)}},
        4: {Q(4, 3): {(0, 3, 0)}, Q(3, 2): {(1, 2, 0)},
            2: {(0, 2, 0), (2, 1, 0)}},
        6: {2: {(0, 3, 0), (0, 1, 1)}, Q(5, 2): {(1, 2, 0)}},
    }
    for target, desired in lifts.items():
        actual = {}
        for i, j, k in product(range(4), repeat=3):
            if j and i+j+k <= 3:
                y_weight = Q(target-i-4*k, j)
                if 1 < y_weight <= Q(5, 2):
                    actual.setdefault(y_weight, set()).add((i, j, k))
        assert actual == desired


def check_four_three_bracket():
    u, y = symbols("u y")
    b, c, d, e, f = (Function(name)(u) for name in ("b", "c", "d", "e", "f"))
    d0, e0 = symbols("d0 e0")
    P, Q = y**4+d*y**2+e*y+f, y**3+b*y+c
    coefficients = Poly(expand(diff(P,u)*diff(Q,y)-diff(P,y)*diff(Q,u)),y)
    expected = {
        4: 3*diff(d,u)-4*diff(b,u),
        3: 3*diff(e,u)-4*diff(c,u),
        2: b*diff(d,u)-2*d*diff(b,u)+3*diff(f,u),
        1: b*diff(e,u)-2*d*diff(c,u)-e*diff(b,u),
        0: b*diff(f,u)-e*diff(c,u),
    }
    for degree, value in expected.items():
        assert expand(coefficients.nth(degree)-value) == 0
    normalized = {
        d: 4*b/3+d0,
        e: 4*c/3+e0,
        f: 2*b**2/9+2*d0*b/3,
    }
    substitute = lambda value: expand(value.subs(normalized).doit())
    for degree in (4,3,2):
        assert substitute(expected[degree]) == 0
    B, C = b+3*d0/2, c+3*e0/4
    assert expand(substitute(expected[1])+4*diff(B*C,u)/3) == 0
    integral = 4*B**3/27-d0*B**2/3-2*C**2/3
    assert expand(substitute(expected[0])-diff(integral,u)) == 0


def check_axial_five_reduction():
    monomials = {m for m in product(range(7), repeat=3) if sum(m) <= 6}
    initial = (
        {(4, 0, 2), (3, 2, 1), (2, 4, 0)},
        {(5, 0, 1), (4, 2, 0)},
        {(6, 0, 0)},
    )
    ref = ((4, 0, 2), (5, 0, 1), (6, 0, 0))
    pure = cut((monomials,) * 3, (2, 3, 4), initial)
    through_four = cut(pure, (1, 1, 4), tuple({m} for m in ref))
    assert competing_monomials(through_four, (1, 1, 0), 2, 4, 5, ref) == (
        {},
        {5: {(0, 0, 2)}},
        {5: {(1, 0, 1), (0, 1, 1)}},
    )
    through_five = cut(through_four, (1, 1, 5), tuple({m} for m in ref))
    assert tuple(max(m[2] for m in allowed) for allowed in through_five) == (2, 1, 1)
    assert {m for m in through_five[2] if m[2]} == {(0, 0, 1)}
    for target, expected in ((2, {2: {(0, 1, 0)}}), (6, {3: {(0, 2, 0)}})):
        lifts = {}
        for i, j, k in product(range(3), repeat=3):
            if j and i+j+k <= 2:
                y_weight = Q(target-i-5*k, j)
                if 1 < y_weight <= 3:
                    lifts.setdefault(y_weight, set()).add((i, j, k))
        assert lifts == expected


if __name__ == "__main__":
    check_weights()
    check_coefficients()
    check_y_lifts()
    check_veronese_fibre()
    check_isolated_third_correction()
    check_isolated_second_correction()
    check_simultaneous_corrections()
    check_axial_four_reduction()
    check_four_three_bracket()
    check_axial_five_reduction()
    print("Uncorrected doubled-square frame: axial support, determinants and fibre identities verified.")
