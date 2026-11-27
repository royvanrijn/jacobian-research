"""Replay the cubic–octic conductor equations and all highest-weight gap tests."""

from sympy import Poly, Rational, diff, expand, rem, resultant, symbols


a, b, c, d, e, f, z = symbols("a b c d e f z")
D, E, F = symbols("D E F")
da, dc, dd, de, df = symbols("da dc dd de df")

U = z**3 + a*z
V = z**8 + b*z**7 + c*z**5 + d*z**4 + e*z**2 + f*z

d0 = D - 28*a**2/9 - 49*a*b**2/24 + 7*b*c/8
e0 = E + 4*D*a/3 - 176*a**3/81 - 49*a**2*b**2/48 - 7*a*b*c/24 + 5*c**2/16
f0 = (
    F - 7*D*a*b/6 + D*c/2 + 98*a**3*b/81 + 343*a**2*b**3/288
    - 10*a**2*c/9 - 49*a*b**2*c/48 + 7*b*c**2/32
)
R = (
    7*D*a**2*b/6 - D*a*c/2 + 7*E*a*b/4 - 3*E*c/4 - F*a
    - 28*a**4*b/81 - 343*a**3*b**3/1728 + 5*a**3*c/27
    - 49*a**2*b**2*c/192 + 21*a*b*c**2/64 - 5*c**3/64
)
G = (
    4*D*a**3/27 - 49*D*a**2*b**2/96 + 7*D*a*b*c/16 - 3*D*c**2/32
    + E*a**2/3 + 7*F*a*b/8 - 3*F*c/8 - 8*a**5/243
    + 343*a**4*b**2/1296 + 2401*a**3*b**4/6912 - 77*a**3*b*c/216
    - 343*a**2*b**3*c/768 + 5*a**2*c**2/48 + 49*a*b**2*c**2/256
    - 7*b*c**3/256
)


def derivative(expression):
    return sum(diff(expression, coefficient)*rate for coefficient, rate in
               zip((a, c, d, e, f), (da, dc, dd, de, df)))


def check_parameter_equations():
    db = symbols("db")
    left = diff(U, z)*(derivative(V) + db*z**7) - derivative(U)*diff(V, z)
    assert Poly(left, z).nth(9) == 3*db
    H = (
        -8*da*U**3/9 + (3*dc - 7*b*da)*V/8 + (a*dc - 5*c*da)*U**2/6
        + (8*a**3*da/9 + (a*de - 2*e*da)/3)*U
    )
    reduced = Poly(expand(left.subs(db, 0) - diff(H, z)), z)
    assert expand(reduced.nth(6) - 3*derivative(d-d0)) == 0
    known_d = {d: d0, dd: derivative(d0)}
    assert expand(reduced.nth(4).subs(known_d, simultaneous=True)
                  - 3*derivative(e-e0)) == 0
    assert expand(reduced.nth(3).subs(known_d, simultaneous=True)
                  - 3*derivative(f-f0)) == 0
    substitutions = {**known_d, e: e0, de: derivative(e0), f: f0, df: derivative(f0)}
    assert expand(reduced.as_expr().subs(substitutions, simultaneous=True)
                  - z*derivative(R) - derivative(G)) == 0

    weights = {a: 2, b: 1, c: 3, D: 4, E: 6, F: 7}
    s = symbols("s")
    for polynomial, weight in ((d0, 4), (e0, 6), (f0, 7), (R, 9), (G, 10)):
        rescale = {variable: variable*s**power for variable, power in weights.items()}
        assert expand(polynomial.subs(rescale, simultaneous=True) - polynomial*s**weight) == 0


def check_collisions():
    S, p, A, C, t, tau = symbols("S p A C t tau")

    # The divided difference of z^i is the complete symmetric polynomial
    # of degree i-1 in the two endpoints.
    divided_monomials = [0, 1]
    for _ in range(2, 9):
        divided_monomials.append(expand(S*divided_monomials[-1] - p*divided_monomials[-2]))

    def divided_difference(polynomial):
        polynomial = Poly(polynomial, z)
        return sum(polynomial.nth(i) * divided_monomials[i] for i in range(9))

    assert expand(divided_difference(U) - (S**2-p+a)) == 0
    f_S = (
        S**7 + b*S**6 + 2*a*S**5 + (4*a*b-c)*S**4 + (-2*a**2-d)*S**3
        + (3*a**2*b-a*c)*S**2 + (-4*a**3-2*a*d+e)*S - a**3*b + a**2*c + f
    )
    assert expand(divided_difference(V).subs(p, S**2+a) - f_S) == 0
    assert expand((S**2-4*p).subs(p, S**2+a) + 3*S**2 + 4*a) == 0
    collision = f_S.subs({d: d0, e: e0, f: f0})

    poles = {a: A/t**2, c: C/t**3, S: S/t}
    initial_R = Poly(expand(t**9*R.subs(poles, simultaneous=True)), t).nth(0)
    initial_G = Poly(expand(t**10*G.subs(poles, simultaneous=True)), t).nth(0)
    assert expand(initial_R - 5*A**3*C/27 + 5*C**3/64) == 0
    assert expand(initial_G + 8*A**5/243 - 5*A**2*C**2/48) == 0
    initial_collision = Poly(expand(t**7*collision.subs(poles, simultaneous=True)), t).nth(0)
    leading = (
        S**7 + 2*A*S**5 - C*S**4 + 10*A**2*S**3/9 - A*C*S**2
        + (4*A**3/81 + 5*C**2/16)*S - A**2*C/9
    )
    assert expand(initial_collision - leading) == 0
    gap = 3*S**2 + 4*A
    balanced = C**2 - 64*A**3/27
    assert rem(resultant(leading, gap, S), balanced, C) == 20800*A**7/9
    assert rem(leading.subs(C, 0), gap, S) == -20*A**3*S/81
    assert initial_G.subs(C, 0) == -8*A**5/243
    assert rem(initial_G, balanced, C) == 52*A**5/243

    homogeneous = {b: 0, D: 0, E: 0, F: 0}
    assert expand(R.subs(homogeneous) - c*(5*a**3/27 - 5*c**2/64)) == 0
    for shape in ({A: 1, C: 0}, {A: Rational(27, 64), C: Rational(27, 64)}):
        scale = {a: shape[A]*tau**2, c: shape[C]*tau**3, S: S*tau}
        assert expand(collision.subs(homogeneous).subs(scale, simultaneous=True)
                      - tau**7*leading.subs(shape)) == 0
        g0 = initial_G.subs(shape)
        assert g0 != 0
        assert expand(G.subs(homogeneous).subs(scale, simultaneous=True)
                      - tau**10*g0) == 0


if __name__ == "__main__":
    check_parameter_equations()
    check_collisions()
    print("Verified the cubic–octic coefficient, collision, pole and rescaling identities.")
