"""Exact Laurent and local profile identities for the (9,6,8) exclusion."""

from sympy import Matrix, Poly, Rational, expand, series, simplify, symbols


def check_parameter_equations():
    z, u, k, l, r, s, t, X, v = symbols("z u k l r s t X v")
    sigma, tau, K, B = symbols("sigma tau K B")
    L = z**3 + k*z + l
    M = r*z**2 + s*z + t
    quotient, remainder = divmod(Poly(M**2, z), Poly(L, z))
    assert quotient.as_expr() == r**2*z + 2*r*s
    tail = 3*remainder.as_expr()/(8*L) - M**3/(16*L**3)
    expansion = expand(series(tail.subs(z, 1/u), u, 0, 6).removeO())
    N1, N2, N3, N4, N5 = (expansion.coeff(u, j) for j in range(1, 6))

    substitution = {
        s: r*X, t: r*(k-X**2)/2-2*sigma/(3*r),
        l: 4*(tau-sigma*X)/(3*r**2)-X*(k+X**2),
    }
    assert simplify(N1.subs(substitution, simultaneous=True) + sigma/2) == 0
    assert simplify(N2.subs(substitution, simultaneous=True) + tau/2) == 0

    E3 = 9*r**4*v**2 + 8*sigma*r**2*v + 96*sigma*r**2*X**2 - 96*tau*r**2*X
    E3 += -6*r**5 + 16*sigma**2 - 16*K*r**2
    E4 = (48*sigma*X-24*tau)*v - 96*sigma*X**3 + 72*tau*X**2 + 27*r**3*X
    E4 += 64*sigma*(sigma*X-tau)/r**2 + 24*B
    reduced = lambda expression: expression.subs(substitution, simultaneous=True).subs(k, v-3*X**2)
    assert simplify(reduced(6*N3-k*sigma-K) - E3/(16*r**2)) == 0
    assert simplify(reduced(6*N4-2*k*tau-l*sigma-B) + E4/24) == 0

    G = -6*N5-K*k-Rational(7, 6)*k**2*sigma+2*l*tau+(k**2+r)*sigma/2
    recorded_G = (
        3*sigma*X**4 - 2*tau*X**3 + 9*r**3*X**2/8 - 2*sigma*X**2*v
        + tau*X*v + 4*sigma*tau*X/(3*r**2) - 3*r**3*v/16
        - r*sigma/4 + sigma*v**2/3 - 4*tau**2/(3*r**2)
    )
    assert simplify(reduced(G)-recorded_G-(v-3*X**2)*E3/(16*r**2)) == 0

    # If r vanishes identically, sigma, tau, K and B fix s, t, k and l in order.
    zero = {r: 0}
    assert expand((-2*N1).subs(zero) + 3*s**2/4) == 0
    assert expand((-2*N2).subs(zero) + 3*s*t/2) == 0
    assert expand((6*N3+2*k*N1).subs(zero) - (9*t**2/4-3*k*s**2/2)) == 0
    assert expand((6*N4+4*k*N2+2*l*N1).subs(zero) - (-3*l*s**2/2-3*k*s*t/2)) == 0


def check_removable_zero():
    r, u, D, C, tau, K, B = symbols("r u D C tau K B", nonzero=True)
    sigma = -3*C**2/4
    equations = (
        9*D**2+8*sigma*D+96*sigma*u**2-96*tau*r*u-6*r**5+16*sigma**2-16*K*r**2,
        (48*sigma*u-24*tau*r)*D - 96*sigma*u**3 + 72*tau*r*u**2
        + 27*r**5*u + 64*sigma*(sigma*u-tau*r) + 24*B*r**3,
    )
    jacobian = Matrix(equations).jacobian((u, D))
    assert jacobian.det().subs({r: 0, u: C, D: 3*C**2}) == -1728*C**6
    u0 = C + 2*tau*r/(3*C**2) + (K*C**2-3*tau**2)*r**2/(3*C**5)
    D0 = 3*C**2 + 4*tau*r/C + 4*(K*C**2-3*tau**2)*r**2/(3*C**4)
    for equation in equations:
        substituted = expand(equation.subs({u: u0, D: D0}, simultaneous=True))
        assert all(substituted.coeff(r, j) == 0 for j in range(3))
    k_scaled = expand(D0-3*u0**2)
    l_scaled = expand(4*(tau*r-sigma*u0)/3-u0*(D0-2*u0**2))
    assert all(k_scaled.coeff(r, j) == 0 for j in range(2))
    assert all(l_scaled.coeff(r, j) == 0 for j in range(3))


def check_big_pole_midpoint():
    z, epsilon = symbols("z epsilon")
    eta = -3*epsilon**2/7
    L = (z-1)**2*(z+2)+epsilon*(z-1)+9*epsilon**2/56
    f = (z-1)*(z+2)+epsilon/2+9*epsilon**2/112
    near_minus_two = -2+epsilon/3
    assert expand(L.subs(z, near_minus_two)).coeff(epsilon, 1) == 0
    assert expand(f.subs(z, near_minus_two)).coeff(epsilon, 1) == -Rational(1, 2)
    near_minus_one = -1+epsilon/6+23*epsilon**2/1008
    first = expand((2*L*L.diff(z)+eta*f.diff(z)).subs(z, near_minus_one))
    assert all(first.coeff(epsilon, j) == 0 for j in range(3))
    numerator = expand((L/2-f*f.diff(z)).subs(z, near_minus_one))
    assert all(numerator.coeff(epsilon, j) == 0 for j in range(2))
    assert numerator.coeff(epsilon, 2) == Rational(3, 16)
    assert simplify(series(numerator/(2*L.subs(z, near_minus_one)), epsilon, 0, 3)
                    .removeO()-3*epsilon**2/128) == 0


if __name__ == "__main__":
    check_parameter_equations()
    check_removable_zero()
    check_big_pole_midpoint()
    print("Verified the conductor-twenty nonic equations, removable chart and pole midpoint.")
