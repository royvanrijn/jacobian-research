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


if __name__ == "__main__":
    check_weights()
    check_coefficients()
    check_y_lifts()
    check_veronese_fibre()
    check_isolated_third_correction()
    print("Uncorrected doubled-square frame: axial support, determinants and fibre identities verified.")
