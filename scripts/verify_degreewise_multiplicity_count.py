#!/usr/bin/env python3
"""Exact arithmetic regression for the short M1 assembly.

This checker deliberately uses only the unconditional signature
  (e_Delta, mu) = (r+1, (N-1)((N-1)/r-1))
and the number N-1-r of parameter roots for each proper divisor r of N-1.
It has no dependency on the refined contact-resultant diagram.
"""


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


for degree in range(4, 201):
    n = degree - 1
    proper = divisors(n)[:-1]
    signatures: set[tuple[int, int]] = set()
    cancellation_count = 0

    for r in proper:
        m = n // r - 1
        assert m >= 1
        assert r * (m + 1) + 1 == degree

        parameter_roots = m * r
        assert parameter_roots == n - r
        cancellation_count += parameter_roots

        signature = (r + 1, n * (n // r - 1))
        assert signature == (r + 1, m * n)
        assert signature not in signatures
        signatures.add(signature)

    tau = len(divisors(n))
    sigma = sum(divisors(n))
    assert cancellation_count == n * tau - sigma
    assert 1 + cancellation_count == 1 + n * tau - sigma


print("PASS M1: unconditional signatures separate every divisor type")
print("PASS M1: parameter-root counting gives 1+(N-1)tau(N-1)-sigma(N-1)")
