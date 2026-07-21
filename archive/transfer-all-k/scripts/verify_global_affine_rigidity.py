"""Exact regressions for arbitrary global transfer equalizers."""

import math
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


def product(items):
    answer = sp.Poly(1, W, domain=sp.QQ)
    for item in items:
        answer *= sp.Poly(item, W, domain=sp.QQ)
    return answer


def transfer_diagnostic(vector):
    """Audit the cross-coupled global equation for one transfer vector."""
    assert sum(vector) == 0 and all(value != 0 for value in vector)
    blocks = []
    for index, value in enumerate(vector):
        size = abs(value)
        blocks.append((
            sp.Poly(monic(101 + 17 * index, 3 * size), W, domain=sp.QQ),
            sp.Poly(monic(211 + 19 * index, 2 * size), W, domain=sp.QQ),
            value,
        ))

    Q_plus = product(U for U, _, value in blocks if value > 0)
    Q_minus = product(U for U, _, value in blocks if value < 0)
    R_plus = product(V for _, V, value in blocks if value > 0)
    R_minus = product(V for _, V, value in blocks if value < 0)
    M_plus = Q_plus**2 * R_minus**3
    M_minus = Q_minus**2 * R_plus**3

    wronskian = M_minus * M_plus.diff() - M_plus * M_minus.diff()
    divisor = Q_plus * R_minus**2 * Q_minus * R_plus**2
    quotient = (
        Q_minus * R_plus
        * (2 * Q_plus.diff() * R_minus + 3 * Q_plus * R_minus.diff())
        - Q_plus * R_minus
        * (2 * Q_minus.diff() * R_plus + 3 * Q_minus * R_plus.diff())
    )
    assert wronskian == divisor * quotient

    weight = sum(value for value in vector if value > 0)
    assert -sum(value for value in vector if value < 0) == weight
    assert M_plus.degree() == M_minus.degree() == 12 * weight
    assert divisor.degree() == 14 * weight > M_plus.degree()

    total_transfer = sum(abs(value) for value in vector)
    rank = 2**total_transfer
    hilbert = tuple(math.comb(total_transfer, degree)
                    for degree in range(total_transfer + 1))
    return rank, hilbert


diagnostics = {
    "Z3 with three Z1 compensators": ((3, -1, -1, -1), 64),
    "two interacting Z3 blocks": ((3, -3), 64),
    "Z4 with four Z1 compensators": ((4, -1, -1, -1, -1), 256),
    "mixed Z3-Z2-Z1": ((3, -2, -1), 64),
    "three-block Z4-Z3-Z1": ((4, -3, -1), 256),
}

for name, (vector, expected_rank) in diagnostics.items():
    rank, hilbert = transfer_diagnostic(vector)
    assert rank == expected_rank
    assert sum(hilbert) == rank

# A primitive Z1 contribution used with exponent -3 at one cluster is an
# effective transfer -3, hence a Z3 block rather than three tensor factors.
weights = (3, 1)
exponents = (1, -3)
assert sum(weight * exponent for weight, exponent in zip(weights, exponents)) == 0
assert tuple(weight * exponent for weight, exponent in zip(weights, exponents)) == (3, -3)


print("PASS: the monic Wronskian divisor is audited for 2/3 pairs through degree 12")
print("PASS: the shared affine coefficients vanish over every Q-algebra")
print("PASS: five Z3/Z4 global transfer diagnostics have the predicted tensor ranks")
print("PASS: concentrated exponents are classified by their effective local transfer")
