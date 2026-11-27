"""Exact finite-field certificates for the remaining powered triangle windows.

Run from the repository root with SageMath. The optional argument selects one D.
"""

import argparse
from math import comb, factorial

from sage.all import GF, PolynomialRing, TermOrder

PRIME = 100003


def power_sums(s, h, r, limit):
    sums = [3 + 0 * s, s, s**2 - 2 * h]
    for _ in range(3, limit + 1):
        sums.append(s * sums[-1] - h * sums[-2] + r * sums[-3])
    return sums


def orbit_sum(i, j, k, r, sums, half):
    """Sum the distinct permutations of a^i b^j c^k, where i >= j >= k."""
    if i == k:
        return r**k
    if j == k:
        return r**k * sums[i - k]
    result = sums[i - k] * sums[j - k] - sums[i + j - 2 * k]
    return r**k * (half * result if i == j else result)


def verify_exponent(exponent):
    stop = 127 if exponent == 2 else 47
    assert PRIME > exponent * (stop + 2)
    field = GF(PRIME)
    symmetric = PolynomialRing(field, "h,r", order=TermOrder("wdegrevlex", (2, 3)))
    boundary = PolynomialRing(field, "t")
    h, r = symmetric.gens()
    t = boundary.gen()
    half = field(2)**(-1)
    ordinary_sums = power_sums(symmetric.one(), h, r, stop + 2)
    zero_sum_sums = power_sums(boundary.zero(), t, t + 1, stop + 2)

    fac = [1]
    for i in range(1, exponent * (stop + 2) + 1):
        fac.append(fac[-1] * i % PRIME)
    invfac = [pow(value, -1, PRIME) for value in fac[:stop + 3]]

    moments = {}
    for order in range(3, stop + 3):
        weight = [fac[exponent * (order - i)] * invfac[i] % PRIME
                  for i in range(order + 1)]
        affine = symmetric.zero()
        zero_sum = boundary.zero()
        for k in range(order // 3 + 1):
            for j in range(k, (order - k) // 2 + 1):
                i = order - j - k
                scalar = weight[i] * weight[j] * weight[k] % PRIME
                affine += scalar * orbit_sum(i, j, k, r, ordinary_sums, half)
                zero_sum += scalar * orbit_sum(i, j, k, t + 1, zero_sum_sums, half)
        binomial = boundary([weight[0] * weight[j] * weight[order - j] % PRIME
                             for j in range(order + 1)])
        moments[order] = (affine, zero_sum, binomial)
        if order < 5:
            continue
        start = order - 2
        (f0, g0, b0), (f1, g1, b1), (f2, g2, _) = (
            moments.pop(start), moments[start + 1], moments[order])
        assert b0.gcd(b1) == 1, (exponent, start, "two-term boundary")
        assert g0.gcd(g1).gcd(g2) == 1, (exponent, start, "zero-sum chart")
        assert symmetric.ideal([f0, f1, f2]).groebner_basis() == [symmetric.one()], (
            exponent, start, "unit sum chart")
    print(f"D={exponent}: exact certificates for all starts 3 through {stop}.", flush=True)


def verify_cutoff_constants():
    """Exact rational inequalities at the two starts of the analytic ranges."""
    from fractions import Fraction as Q

    assert sum(Q(3)**i / factorial(i) for i in range(9)) > 20
    assert sum(Q(12, 5)**i / factorial(i) for i in range(10)) > 11
    for m, theta in ((128, Q(9, 10)), (48, Q(4, 5))):
        B = [None] + [Q(comb(m, l), comb(2 * m, 2 * l)) for l in range(1, 4)]
        assert 64 * m * theta**m < 1
        if m == 128:
            tail = Q(17010, (2*m-7)**4) + Q(2*m, 20**6)
            low = 16*B[1] + 510*B[2] + 17*158*B[3]
        else:
            tail = Q(3, 2) * Q(10, 3)**4 / comb(m, 4)**2
            tail += Q(2*m*(m+1), 11**9)
            low = 20*B[1]/m + 840*B[2]/comb(m, 2) + 28*158*B[3]/comb(m, 3)
        assert 4 * (low + 8500 * Q(51, 50)**2 * tail) < Q(3, 4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exponent", type=int, nargs="?", choices=range(2, 36))
    args = parser.parse_args()
    verify_cutoff_constants()
    for d in [args.exponent] if args.exponent is not None else range(2, 36):
        verify_exponent(d)
