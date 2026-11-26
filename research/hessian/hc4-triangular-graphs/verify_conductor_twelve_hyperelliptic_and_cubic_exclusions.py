"""Replay the exact coefficient identities for two conductor-twelve branches."""

from sympy import Poly, Rational, diff, expand, integrate, rem, resultant, symbols


def check_parameter_equations():
    z, a, b, c, d, e = symbols("z a b c d e")
    da, db, dc, dd, de = symbols("da db dc dd de")
    B, D, E = symbols("B D E")

    def derivative(expression):
        return sum(diff(expression, value) * rate for value, rate in
                   zip((a, b, c, d, e), (da, db, dc, dd, de)))

    U = z**3 + a*z
    V = z**7 + b*z**5 + c*z**4 + d*z**2 + e*z
    jacobian = Poly(expand(diff(U, z)*derivative(V) - derivative(U)*diff(V, z)), z)
    H = 3*dc*V/7 + (a*db - 5*b*da)*U**2/6 + (a*dd - 2*d*da)*U/3
    reduced = Poly(expand(jacobian.as_expr() - diff(H, z)), z)
    assert jacobian.nth(8) == 0
    assert reduced.nth(7) == 3*derivative(b - 7*a/3)
    assert all(reduced.nth(j) == 0 for j in (2, 5, 6))

    b0 = B + 7*a/3
    d0 = D + (4*a/3 + 5*B/7)*c
    e0 = E - 112*a**3/81 - 10*B*a**2/9 + 2*c**2/7
    with_b = {b: b0, db: derivative(b0)}
    assert expand(reduced.nth(4).subs(with_b, simultaneous=True)
                  - 3*derivative(d - d0)) == 0
    assert expand(reduced.nth(3).subs(with_b, simultaneous=True)
                  - 3*derivative(e - e0)) == 0

    R = 7*a**4/81 + 5*B*a**3/27 - (2*a/7 + 15*B/49)*c**2 - E*a - 6*D*c/7
    G = (4*a**3/27 + 5*B*a**2/21 - 3*E/7)*c + D*a**2/3 - 2*c**3/49
    substitutions = dict(with_b)
    substitutions.update({d: d0, e: e0, dd: derivative(d0), de: derivative(e0)})
    assert expand(reduced.as_expr().subs(substitutions, simultaneous=True)
                  - z*derivative(R) - derivative(G)) == 0


def check_collision_equations():
    z, S, P, a, b, c, d, e = symbols("z S P a b c d e")
    A, C, B, D, E, t, tau = symbols("A C B D E t tau")
    U = z**3 + a*z
    V = z**7 + b*z**5 + c*z**4 + d*z**2 + e*z

    # The divided difference of z^j is the complete symmetric polynomial of
    # degree j-1 in the endpoints, also on their diagonal by polynomial continuation.
    differences = [0, 1]
    for _ in range(2, 8):
        differences.append(expand(S*differences[-1] - P*differences[-2]))

    def divided_difference(polynomial):
        polynomial = Poly(polynomial, z)
        return sum(polynomial.nth(j)*differences[j] for j in range(8))

    assert expand(divided_difference(U) - (S**2 - P + a)) == 0
    collision = (
        S**6 + (4*a-b)*S**4 - c*S**3 + (3*a**2-a*b)*S**2
        + (d-2*a*c)*S - a**3 + a**2*b + e
    )
    assert expand(divided_difference(V).subs(P, S**2+a) - collision) == 0
    collision = expand(collision.subs({
        b: B + 7*a/3,
        d: D + (4*a/3 + 5*B/7)*c,
        e: E - 112*a**3/81 - 10*B*a**2/9 + 2*c**2/7,
    }))
    G = (4*a**3/27 + 5*B*a**2/21 - 3*E/7)*c + D*a**2/3 - 2*c**3/49

    first = {a: A/t**2, c: C/t**3, S: S/t}
    balanced = Poly(expand(t**6*collision.subs(first, simultaneous=True)), t).nth(0)
    ratio = C**2 - 49*A**3/162
    leading = S**6 + 5*A*S**4/3 - C*S**3 + 2*A**2*S**2/3 - 2*A*C*S/3 + A**3/27
    assert rem(balanced-leading, ratio, C) == 0
    assert rem(resultant(leading, 3*S**2+4*A, S), ratio, C) == 539*A**6/3
    initial_G = Poly(expand(t**9*G.subs(first, simultaneous=True)), t).nth(0)
    assert rem(initial_G - 11*A**3*C/81, ratio, C) == 0

    second = {c: C/t**3, S: S/t}
    pure_c = Poly(expand(t**6*collision.subs(second, simultaneous=True)), t).nth(0)
    assert pure_c == S**6 - C*S**3 + 2*C**2/7
    assert resultant(pure_c, 3*S**2, S) == 2916*C**4/49
    assert Poly(expand(t**9*G.subs(second, simultaneous=True)), t).nth(0) == -2*C**3/49

    # Inertia of order two forces c and D to vanish at the limiting equation.
    d2_residual = E - 28*a**3/81 - 5*B*a**2/9
    assert expand(rem(collision.subs({c: 0, D: 0}), 3*S**2+4*a, S) - d2_residual) == 0
    assert expand(Poly(G, c).nth(1) + 3*d2_residual/7) == 0
    d5_equation = collision.subs({a: 0, c: 0, B: 0, E: 0})
    assert d5_equation == S**6 + D*S
    assert rem(diff(d5_equation, S), S**5+D, S) == -5*D

    # The two homogeneous branches give the pole-scale equations exactly.
    homogeneous = {B: 0, D: 0, E: 0}
    pure = {a: 0, c: tau**3, S: S*tau}
    assert expand(collision.subs(homogeneous).subs(pure, simultaneous=True)
                  - tau**6*pure_c.subs(C, 1)) == 0
    balanced_a = Rational(162, 49)
    balanced_scale = {a: balanced_a*tau**2, c: balanced_a*tau**3, S: S*tau}
    assert expand(collision.subs(homogeneous).subs(balanced_scale, simultaneous=True)
                  - tau**6*leading.subs({A: balanced_a, C: balanced_a})) == 0


def check_hyperelliptic_resonance_and_node():
    r, u = symbols("r u")
    resonances = [Rational(12-2*k, 2*k+3) for k in range(1, 6)]
    assert resonances == [2, Rational(8, 7), Rational(2, 3),
                          Rational(4, 11), Rational(2, 13)]
    T = u**2
    V = u*(u**2-r**2)
    assert integrate(T*diff(V, u), (u, -r, r)) == 8*r**5/15


if __name__ == "__main__":
    check_parameter_equations()
    check_collision_equations()
    check_hyperelliptic_resonance_and_node()
    print("Verified the conductor-twelve parameter equations, collision limits and nodal action.")
