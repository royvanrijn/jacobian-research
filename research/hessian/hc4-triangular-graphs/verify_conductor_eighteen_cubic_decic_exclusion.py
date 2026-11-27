"""Replay the cubic–decic parameter integrations and finite-pair pole tests."""

from sympy import Poly, Rational, diff, expand, rem, resultant, symbols


z, a, b, c, d, e, f, g = symbols("z a b c d e f g")
B, D, E, F, M = symbols("B D E F M")
da, db, dc, dd, de, df, dg = symbols("da db dc dd de df dg")


def derivative(expression):
    return sum(diff(expression, coefficient) * rate for coefficient, rate in
               zip((a, b, c, d, e, f, g), (da, db, dc, dd, de, df, dg)))


lower = z**3 + a*z
upper = z**10 + b*z**8 + c*z**7 + d*z**5 + e*z**4 + f*z**2 + g*z
b_integral = B + 10*a/3
d_integral = D + 4*B*c/5 + 7*a*c/3
e_integral = E - 28*B*a**2/9 - 490*a**3/81 + 7*c**2/20
f_integral = (
    F - 176*B*a**3/81 + B*c**2/5 + D*c/2 + 4*E*a/3
    - 910*a**4/243 + 7*a*c**2/15
)
g_integral = (
    M - 8*B*a**2*c/9 - 10*D*a**2/9 + 2*E*c/5
    - 112*a**3*c/81 + 7*c**3/150
)
first_integral = (
    1200*B*a**3*c - 324*B*c**3 + 1500*D*a**3 - 1215*D*c**2
    - 3240*E*a*c - 4860*F*c - 8100*M*a + 700*a**4*c - 378*a*c**3
)/8100
scale_primitive = -(
    144000*B*a**5 - 291600*B*a**2*c**2 - 729000*D*a**2*c - 648000*E*a**3
    + 262440*E*c**2 - 1458000*F*a**2 + 1312200*M*c
    + 70000*a**6 - 226800*a**3*c**2 + 15309*c**4
)/4374000


def check_parameter_equations():
    left = diff(lower, z)*derivative(upper) - derivative(lower)*diff(upper, z)
    assert Poly(left, z).nth(10) == 3*db - 10*da
    substitutions = {b: b_integral, db: derivative(b_integral)}
    hamiltonian = (
        3*dc*upper/10 - 2*da*(12*B + 35*a)*lower**3/27
        + (a*dd - 5*d*da)*lower**2/6
        + (24*B*a**3*da + 70*a**4*da + 9*a*df - 18*da*f)*lower/27
    )
    residual = Poly((left - diff(hamiltonian, z)).subs(substitutions, simultaneous=True), z)
    assert residual.degree() == 7
    integrations = (
        (7, d, dd, d_integral),
        (6, e, de, e_integral),
        (4, f, df, f_integral),
        (3, g, dg, g_integral),
    )
    for degree, variable, rate, integral in integrations:
        assert expand(residual.nth(degree) - 3*derivative(variable - integral)) == 0
        residual = Poly(residual.as_expr().subs(
            {variable: integral, rate: derivative(integral)}, simultaneous=True), z)
    assert expand(residual.as_expr() - z*derivative(first_integral)
                  - derivative(scale_primitive)) == 0
    weights = {a: 2, c: 3, B: 2, D: 5, E: 6, F: 8, M: 9}
    t = symbols("t")
    scale = {variable: variable*t**weight for variable,weight in weights.items()}
    for expression, weight in ((b_integral, 2), (d_integral, 5), (e_integral, 6),
                               (f_integral, 8), (g_integral, 9), (first_integral, 11),
                               (scale_primitive, 12)):
        assert expand(expression.subs(scale, simultaneous=True) - expression*t**weight) == 0


def check_collision_and_pole_tests():
    S, p, A, C, t = symbols("S p A C t")
    divided_monomials = [0, 1]
    for degree in range(2, 11):
        divided_monomials.append(expand(S*divided_monomials[-1] - p*divided_monomials[-2]))

    def divided_difference(expression):
        poly = Poly(expression, z)
        return sum(poly.nth(i)*divided_monomials[i] for i in range(11))

    product = S**2 + a
    assert expand(divided_difference(lower) - product + p) == 0
    assert expand(S**2 - 4*product + 3*S**2 + 4*a) == 0
    collision = -divided_difference(upper).subs(p, product).subs(
        {b: b_integral, d: d_integral, e: e_integral, f: f_integral, g: g_integral},
        simultaneous=True,
    )
    assert Poly(collision, S).LC() == 1 and Poly(collision, S).degree() == 9
    parametric_pole = {a: A/t**2, c: C/t**3, S: S/t}

    def initial(expression, weight):
        return Poly(expand(t**weight*expression.subs(parametric_pole, simultaneous=True)), t).nth(0)

    initial_r = initial(first_integral, 11)
    initial_g = initial(scale_primitive, 12)
    assert expand(initial_r - 7*A*C*(50*A**3 - 27*C**2)/4050) == 0
    assert expand(initial_g + 7*(10000*A**6 - 32400*A**3*C**2 + 2187*C**4)/4374000) == 0
    leading = (
        S**9 + 8*A*S**7/3 - C*S**6 + 7*A**2*S**5/3 - 5*A*C*S**4/3
        + (50*A**3/81 + 7*C**2/20)*S**3 - 2*A**2*C*S**2/3
        + (-5*A**4/243 + 7*A*C**2/30)*S + 4*A**3*C/81 - 7*C**3/150
    )
    assert expand(initial(collision, 9) - leading) == 0
    degenerate_gap = 3*S**2 + 4*A
    assert resultant(leading.subs(C, 0), degenerate_gap, S) == 4900*A**9/9
    assert resultant(leading.subs(A, 0), degenerate_gap.subs(A, 0), S) == 107163*C**6/2500
    balanced = C**2 - 50*A**3/27
    assert rem(resultant(leading, degenerate_gap, S), balanced, C) == 41650*A**9/9
    assert rem(initial_g, balanced, C) == 595*A**6/8748
    assert initial_g.subs(A, 0) == -7*C**4/2000
    assert initial_g.subs(C, 0) == -35*A**6/2187

    homogeneous = {B: 0, D: 0, E: 0, F: 0, M: 0}
    for A0, C0 in ((1, 0), (0, 1), (Rational(27, 50), Rational(27, 50))):
        scale = {a: A0*t**2, c: C0*t**3, S: S*t}
        assert expand(collision.subs(homogeneous).subs(scale, simultaneous=True)
                      - t**9*leading.subs({A: A0, C: C0})) == 0
        g0 = initial_g.subs({A: A0, C: C0})
        assert g0 != 0
        assert expand(scale_primitive.subs(homogeneous).subs(scale, simultaneous=True)
                      - t**12*g0) == 0


if __name__ == "__main__":
    check_parameter_equations()
    check_collision_and_pole_tests()
    print("Verified the cubic–decic conductor parameter equations and all three pole regimes.")
