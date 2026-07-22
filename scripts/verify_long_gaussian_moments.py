#!/usr/bin/env python3
"""Exact, dependency-free regression for Long's three-Gaussian witness.

Source: Christopher D. Long, "Small Counterexamples to the Gaussian
Moments Conjecture", arXiv:2607.18186v1 (submitted 20 July 2026).

X1, X2, T are independent N(0,1) real Gaussians and
Z=(X1+iX2)/sqrt(2), W=(X1-iX2)/sqrt(2).  The source formula is

    P = (1+Z) * (W - (2+Z)T^2/2),   Q = Z.

The finite Wick computations below are exact regressions.  The separate
coefficient checks encode finite truncations of the all-m identity printed in
the paper; they are not represented as a new proof replacing Long's written
general argument.
"""

from fractions import Fraction
from math import comb, factorial


Poly = dict[tuple[int, int, int], Fraction]  # exponents of W, Z, T


def add(*polys: Poly) -> Poly:
    out: Poly = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def multiply(left: Poly, right: Poly) -> Poly:
    out: Poly = {}
    for (aw, az, at), ac in left.items():
        for (bw, bz, bt), bc in right.items():
            monomial = (aw + bw, az + bz, at + bt)
            out[monomial] = out.get(monomial, Fraction(0)) + ac * bc
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def power(poly: Poly, exponent: int) -> Poly:
    out: Poly = {(0, 0, 0): Fraction(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def odd_double_factorial(n: int) -> int:
    """Return (n-1)!! for even n, including (-1)!!=1 when n=0."""
    value = 1
    for factor in range(1, n, 2):
        value *= factor
    return value


def gaussian_expectation(poly: Poly) -> Fraction:
    """Use E[W^a Z^b]=delta_(a,b) a! and standard-real T moments."""
    total = Fraction(0)
    for (w_degree, z_degree, t_degree), coefficient in poly.items():
        if w_degree != z_degree or t_degree % 2:
            continue
        total += coefficient * factorial(w_degree) * odd_double_factorial(t_degree)
    return total


P_expanded: Poly = {
    (1, 0, 0): Fraction(1),
    (1, 1, 0): Fraction(1),
    (0, 0, 2): Fraction(-1),
    (0, 1, 2): Fraction(-3, 2),
    (0, 2, 2): Fraction(-1, 2),
}
ONE: Poly = {(0, 0, 0): Fraction(1)}
W: Poly = {(1, 0, 0): Fraction(1)}
Z: Poly = {(0, 1, 0): Fraction(1)}
P = multiply(
    add(ONE, Z),
    add(W, {(0, 0, 2): Fraction(-1)}, {(0, 1, 2): Fraction(-1, 2)}),
)
assert P == P_expanded  # compact source formula equals the five-term expansion
Q: Poly = Z


def binomial_coefficient(degree: int, index: int) -> int:
    return comb(degree, index) if 0 <= index <= degree else 0


def coefficient_master_lhs(m: int, j: int) -> Fraction:
    """Coefficient sum for A(z)=z^j in Long's all-m identity."""
    total = Fraction(0)
    for k in range(m + 1):
        wanted = m - k - j
        coefficient = 0
        # [z^wanted] (1+z)^m (2+z)^k
        for b in range(k + 1):
            a = wanted - b
            coefficient += binomial_coefficient(m, a) * comb(k, b) * 2 ** (k - b)
        total += Fraction((-1) ** k * comb(2 * k, k), 4**k) * coefficient
    return total


def coefficient_master_rhs(m: int, j: int) -> Fraction:
    return Fraction(binomial_coefficient(m - 1, m - j))


def substitute_central_binomial_series(bound: int) -> list[Fraction]:
    """Coefficients through z^bound of sum c_k(z(2+z))^k."""
    result = [Fraction(0) for _ in range(bound + 1)]
    base = [Fraction(0), Fraction(2), Fraction(1)]
    current = [Fraction(1)]
    for k in range(bound + 1):
        scalar = Fraction((-1) ** k * comb(2 * k, k), 4**k)
        for degree, coefficient in enumerate(current[: bound + 1]):
            result[degree] += scalar * coefficient
        next_poly = [Fraction(0) for _ in range(min(bound, len(current) + 1) + 1)]
        for a, ac in enumerate(current):
            for b, bc in enumerate(base):
                if a + b <= bound:
                    next_poly[a + b] += ac * bc
        current = next_poly
    return result


def main() -> None:
    # Direct exact Wick regression of the transcribed witness.
    p_power: Poly = {(0, 0, 0): Fraction(1)}
    for m in range(1, 11):
        p_power = multiply(p_power, P)
        assert gaussian_expectation(p_power) == 0
        assert gaussian_expectation(multiply(Q, p_power)) == factorial(m)

    # The master coefficient identity for a basis A(z)=z^j.
    for m in range(1, 15):
        for j in range(m + 3):
            assert coefficient_master_lhs(m, j) == coefficient_master_rhs(m, j)

    # The formal branch is (1+z(2+z))^(-1/2)=(1+z)^(-1).
    coefficients = substitute_central_binomial_series(20)
    assert coefficients == [Fraction((-1) ** degree) for degree in range(21)]

    print("PASS Long GMC: exact Wick moments m=1..10")
    print("PASS Long GMC: master coefficient identity m=1..14")
    print("PASS Long GMC: central-binomial formal identity through z^20")


if __name__ == "__main__":
    main()
