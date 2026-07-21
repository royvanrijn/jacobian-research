"""Exact regressions for the universal global affine-rigidity identity."""

import sympy as sp

W = sp.symbols("W")


def monic(seed, degree):
    """A deterministic generic-looking monic polynomial over the integers."""
    return W**degree + sum(
        (seed + 2 * index + 1) * W**index for index in range(degree)
    )


# Check every 2/3 degree allocation pair through degree twelve on exact integer
# samples. The prose proof is the universal polynomial identity; keeping this
# regression univariate prevents symbolic coefficient swell.
for n in range(3, 13):
    allocations = [
        (a, b)
        for b in range(n // 3 + 1)
        for a in range(n // 2 + 1)
        if 2 * a + 3 * b == n
    ]
    for left_index, (a, b) in enumerate(allocations):
        for right_index, (c, d) in enumerate(allocations):
            Q = monic(11 + n + left_index, a)
            R = monic(23 + n + left_index, b)
            S = monic(37 + n + right_index, c)
            T = monic(53 + n + right_index, d)
            M = sp.expand(Q**2 * R**3)
            N = sp.expand(S**2 * T**3)
            wronskian = sp.expand(N * sp.diff(M, W) - M * sp.diff(N, W))
            divisor = sp.expand(Q * R**2 * S * T**2)
            quotient = sp.expand(
                S * T * (2 * sp.diff(Q, W) * R + 3 * Q * sp.diff(R, W))
                - Q * R * (2 * sp.diff(S, W) * T + 3 * S * sp.diff(T, W))
            )
            assert sp.expand(wronskian - divisor * quotient) == 0
            assert (b - d) % 2 == 0
            assert sp.degree(divisor, W) == n + (b + d) // 2


# Once the Wronskian vanishes and M-N=lambda*W+mu, the two leading
# coefficients are units times lambda and mu over every Q-algebra.
for n in range(3, 30):
    lam, mu = sp.symbols(f"lambda_{n} mu_{n}")
    lower = sp.symbols(f"v{n}_0:{n}")
    N = W**n + sum(lower[index] * W**index for index in range(n))
    affine_wronskian = sp.Poly(
        sp.expand(N * lam - sp.diff(N, W) * (lam * W + mu)), W
    )
    assert affine_wronskian.coeff_monomial(W**n) == (1 - n) * lam
    assert sp.expand(
        affine_wronskian.coeff_monomial(W ** (n - 1)).subs(lam, 0)
    ) == -n * mu


print("PASS: the monic Wronskian divisor is audited for 2/3 pairs through degree 12")
print("PASS: the shared affine coefficients vanish over every Q-algebra")
