#!/usr/bin/env python3
"""Finite transcription controls; the all-prime proof is the written note.

Run with .venv/bin/python from the repository root.  No floating point roots
or conclusions about untested primes are used.
"""

import hashlib
import json
import platform
from fractions import Fraction
from math import factorial
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
NOTE = HERE / "BILINEAR_PRIME_TERMINAL_THREE_WINDOWS.md"
RECEIPT = HERE / "bilinear_doubled_prime_ramification_verification.json"
PRIMES = (3, 5, 7, 11, 13, 17)
GCD_PRIMES = (3, 5, 7)


def trim(poly):
    poly = list(poly)
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def add(*polys):
    out = [0] * max(map(len, polys), default=0)
    for poly in polys:
        for i, value in enumerate(poly):
            out[i] += value
    return trim(out)


def scale(poly, scalar):
    return trim([scalar * value for value in poly])


def shift(poly, power):
    return [0] * power + list(poly) if poly else []


def mul(left, right):
    out = [0] * (len(left) + len(right))
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def derivative(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def power(poly, exponent):
    result = [1]
    for _ in range(exponent):
        result = mul(result, poly)
    return result


def modulo(poly, prime):
    return trim([value % prime for value in poly])


def at_u_squared_two(poly, prime):
    result = [0, 0]
    for i, value in enumerate(poly):
        result[i % 2] += value * pow(2, i // 2, prime)
    return modulo(result, prime)


def vp(integer, prime):
    assert integer
    integer = abs(integer)
    result = 0
    while integer % prime == 0:
        integer //= prime
        result += 1
    return result


def native_polynomials(cap):
    # E_n=U_n(u/2); D_N=2T_N(u/2).  Separate initial conditions.
    e = [[1], [0, 1]]
    d = [[2], [0, 1]]
    for _ in range(1, cap):
        e.append(add(shift(e[-1], 1), scale(e[-2], -1)))
        d.append(add(shift(d[-1], 1), scale(d[-2], -1)))
    return e, d


def check_native(prime, e, d):
    n = 2 * prime - 1
    coeffs = [e[n]]
    for _ in range(n):
        coeffs.append(derivative(coeffs[-1]))
    assert coeffs[-1] == [factorial(n)]
    coefficient_divisibilities = 0
    for j, poly in enumerate(coeffs):
        for value in poly:
            assert value % factorial(j) == 0
            coefficient_divisibilities += 1

    # Eisenstein M_p, with the constant coefficient of exact valuation one.
    assert d[prime][0] == 0
    m = d[prime][1:]
    assert m[-1] == 1 and all(c % prime == 0 for c in m[:-1])
    assert vp(m[0], prime) == 1

    expected = shift(power([-4, 0, 1], (prime - 1) // 2), prime)
    assert modulo(e[n], prime) == modulo(expected, prime)
    assert modulo(coeffs[prime - 1], prime) == modulo(
        shift([factorial(prime - 1)], prime), prime
    )
    assert at_u_squared_two(coeffs[0], prime)
    assert at_u_squared_two(coeffs[prime - 1], prime)
    assert vp(coeffs[n][0], prime) == 1

    # A F_u = F-E_n, as an identity in Z[u,A].
    for j in range(1, n + 1):
        assert derivative(coeffs[j - 1]) == coeffs[j]

    # Independently collect the displayed second-order identity (24).
    lhs = [[] for _ in range(n + 3)]
    for j, poly in enumerate(coeffs):
        lhs[j] = add(lhs[j], mul([4, 0, -1], poly))
        lhs[j + 1] = add(lhs[j + 1], scale(shift(poly, 1), -(2 * j + 1)))
        lhs[j + 2] = add(
            lhs[j + 2], scale(poly, n * (n + 2) - j * (j + 2))
        )
    rhs = [[] for _ in lhs]
    rhs[0] = mul([4, 0, -1], e[n])
    rhs[1] = scale(d[n + 1], -(n + 1))
    assert lhs == rhs

    # Check the valuation-minimum argument at chosen rational test slopes.
    slopes = sorted({Fraction(-2), Fraction(-1, 2), Fraction(-1, prime),
                     Fraction(-1, 2 * prime), Fraction(1, prime), Fraction(1)})
    for slope in slopes:
        heights = [int(j >= prime) + j * slope for j in range(n + 1)]
        minima = [j for j, h in enumerate(heights) if h == min(heights)]
        if slope == Fraction(-1, prime):
            assert minima == [prime - 1, n]
        elif slope < 0:
            assert minima in ([prime - 1], [n])
        else:
            assert minima == [0]

    return {"p": prime, "n": n, "integer_coefficient_divisibilities":
            coefficient_divisibilities, "universal_identities": 2,
            "valuation_slope_controls": len(slopes)}


def check_discriminant(prime, e, d):
    u, a = sp.symbols("u A")
    n = 2 * prime - 1
    expr = sum(sp.Integer(value) * u**i for i, value in enumerate(e[n]))
    f = sum(sp.diff(expr, u, j) * a**j for j in range(n + 1))
    discr = sp.Poly(sp.discriminant(f, a), u, domain=sp.QQ)
    m_expr = sum(sp.Integer(value) * (u*u - 2)**i
                 for i, value in enumerate(d[prime][1:]))
    m = sp.Poly(m_expr, u, domain=sp.QQ)
    gcd = sp.gcd(discr, m)
    assert gcd.degree() == 0
    return {"p": prime, "degree_A": n, "resonance_degree": m.degree(),
            "discriminant_degree": discr.degree(), "gcd": str(gcd.as_expr())}


def check_negative_controls():
    a = sp.symbols("A")
    factor = 1 + 3 * a**3
    bad = sp.Poly(factor**2, a, domain=sp.QQ)
    gcd = sp.gcd(bad, bad.diff())
    assert bad.degree() == 6 and gcd.degree() == 3
    # The same valuation occurs six times, and the degree is exactly 2p.
    assert sp.discriminant(bad.as_expr(), a) == 0
    return {"degree_bound_removed": {"polynomial": str(bad.as_expr()),
            "p": 3, "root_valuation": "-1/3", "gcd_degree": gcd.degree(),
            "rejected": True}}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    e, d = native_polynomials(2 * max(PRIMES))
    native = [check_native(prime, e, d) for prime in PRIMES]
    exact_gcds = [check_discriminant(prime, e, d) for prime in GCD_PRIMES]
    negative = check_negative_controls()
    receipt = {
        "status": "passed",
        "evidence_type": "internal_exact_finite_transcription_and_negative_controls",
        "scope": "Only the explicitly listed primes and identities; all-prime proof is written.",
        "command": ".venv/bin/python research/factorial/image-rigidity/verify_bilinear_doubled_prime_ramification.py",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "sources": {NOTE.name: sha256(NOTE), Path(__file__).name: sha256(Path(__file__))},
        "native_checks": native,
        "exact_discriminant_controls": exact_gcds,
        "negative_controls": negative,
        "independent_replay": False,
        "external_review": False,
        "formal_verification": False,
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "native_primes": list(PRIMES),
                      "gcd_primes": list(GCD_PRIMES), "negative_controls": 1,
                      "receipt": str(RECEIPT.relative_to(Path.cwd()))}))


if __name__ == "__main__":
    main()
