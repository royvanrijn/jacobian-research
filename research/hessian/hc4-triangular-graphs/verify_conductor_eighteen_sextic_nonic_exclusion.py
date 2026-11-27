"""Exact Laurent, divided-difference and limiting-profile checks for (9,6,7)."""

from sympy import Poly, Rational, expand, factor, rem, resultant, series, simplify, symbols


def check_laurent_equations():
    z, u, k, l, r, s, t, x = symbols("z u k l r s t x")
    sigma, K = symbols("sigma K")
    L = z**3 + k*z + l
    M = r*z**2 + s*z + t
    polynomial_part, remainder = divmod(Poly(M**2, z), Poly(L, z))
    assert polynomial_part.as_expr() == r**2*z + 2*r*s
    tail = 3*remainder.as_expr()/(8*L) - M**3/(16*L**3)
    expansion = expand(series(tail.subs(z, 1/u), u, 0, 6).removeO())
    N1, N2, N3, N4, N5 = (expansion.coeff(u, j) for j in range(1, 6))
    assert expand(N1 + 3*(k*r**2 - 2*r*t - s**2)/8) == 0
    assert expand(N2 + 3*(2*k*r*s + l*r**2 - 2*s*t)/8) == 0
    substitute = {s: r*x, t: r*(k-x**2)/2, l: 4*sigma/(3*r**2)-x*(k+x**2)}
    assert simplify(N1.subs(substitute, simultaneous=True)) == 0
    assert simplify(N2.subs(substitute, simultaneous=True) + sigma/2) == 0
    residual = 3*r**2*(k+3*x**2)**2 - 2*r**3 - 32*sigma*x + 16*K
    assert simplify(N3.subs(substitute, simultaneous=True) + K/2 - residual/32) == 0
    B = -6*N4 + 2*k*sigma
    assert simplify(B.subs(substitute, simultaneous=True) - (9*r**3*x/8-k*sigma)) == 0
    G = -6*N5 + 9*k*K + 2*l*sigma
    reduced_G = (
        6*K*k - 3*k*r**3/16 + sigma*k*x + 9*r**3*x**2/16 + sigma*x**3
        - 4*sigma**2/(3*r**2)
    )
    assert simplify(G.subs(substitute, simultaneous=True) - reduced_G - 3*k*residual/16) == 0


def check_divided_differences_and_rescaling():
    z, w, g, eta, k, l, x = symbols("z w g eta k l x")
    L = z**3 + k*z + l
    f = z**2 + x*z + (k-x**2)/2
    A = z + 2*x
    W = L**2 + eta*f
    V = L**3 + 3*eta*L*f/2 + 3*eta**2*A/8

    def dd(polynomial):
        return expand((polynomial.subs(z, w+g/2)-polynomial.subs(z, w-g/2))/g)

    Q = dd(L)
    S = 2*L.subs(z, w) + Rational(3, 2)*w*g**2
    f0 = f.subs(z, w)
    assert expand(Q - L.diff(z).subs(z, w) - g**2/4) == 0
    assert expand(dd(W) - Q*S - eta*f.diff(z).subs(z, w)) == 0
    assert expand(4*(dd(V)-3*S*dd(W)/4)
                  - (Q**3*g**2+3*Q*eta*(2*f0+g**2/2)+3*eta**2/2)) == 0

    difference = Poly(expand(V**2-W**3), eta)
    assert difference.nth(0) == difference.nth(1) == 0
    assert expand(difference.nth(2) - 3*L**2*(A*L-f**2)/4) == 0
    assert expand(difference.nth(3) - 9*L*f*A/8 + f**3) == 0
    assert expand(A*L-f**2 - ((l+x*(k+x**2))*z + 2*x*l - (k-x**2)**2/4)) == 0


def check_limiting_profiles():
    z, epsilon = symbols("z epsilon")
    # The nonsquare profile, scaled by sqrt(k).
    L = z**3 + z
    M = Rational(3, 2)*(z**2 + Rational(1, 2))
    W = L**2 + M
    V = L**3 + 3*L*M/2 + Rational(27, 32)*z
    difference = expand(V**2-W**3)
    assert Poly(difference, z).degree() == 4
    assert factor(difference.diff(z)) == -27*z*(16*z**2+13)/512
    assert V.diff(z).subs(z, 0) == Rational(63, 32)
    assert rem(W.diff(z), z**2+Rational(13, 16), z) == Rational(315, 128)*z

    # The r=0, unbounded-x profile, in units of its nonzero x residue.
    assert resultant(z**3+Rational(1, 8), z**2+z-Rational(1, 2), z) == -Rational(27, 64)
    assert Rational(1, 4) - (-Rational(1, 2))/(2*Rational(1, 8)) == Rational(9, 4)

    # At the square-leading r-pole profile, eta ~ -3 epsilon²/10 and
    # the constant correction to L is 3 epsilon²/20.
    L = (z-1)**2*(z+2) + epsilon*(z-1) + 3*epsilon**2/20
    f = (z-1)*(z+2) + epsilon/2
    eta = -3*epsilon**2/10
    near_minus_two = -2 + epsilon/3
    assert expand(L.subs(z, near_minus_two)).coeff(epsilon, 1) == 0
    assert expand(f.subs(z, near_minus_two)).coeff(epsilon, 1) == -Rational(1, 2)

    near_minus_one = -1 + epsilon/6 + 29*epsilon**2/1440
    first = (2*L*L.diff(z)+eta*f.diff(z)).subs(z, near_minus_one)
    assert all(expand(first).coeff(epsilon, j) == 0 for j in range(3))
    numerator = (L/2-f*f.diff(z)).subs(z, near_minus_one)
    assert all(expand(numerator).coeff(epsilon, j) == 0 for j in range(2))
    assert expand(numerator).coeff(epsilon, 2) == Rational(3, 32)
    assert simplify(series((numerator/(2*L.subs(z, near_minus_one))), epsilon, 0, 3)
                    .removeO() - 3*epsilon**2/256) == 0


if __name__ == "__main__":
    check_laurent_equations()
    check_divided_differences_and_rescaling()
    check_limiting_profiles()
    print("Verified the sextic–nonic parameter equations, midpoint identities and limiting profiles.")
