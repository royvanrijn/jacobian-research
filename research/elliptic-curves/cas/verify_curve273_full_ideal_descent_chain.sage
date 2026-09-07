#!/usr/bin/env sage
"""Verify the first full-ideal/CRT descent chain for ICARM curve 273.

The certificate uses arbitrary integral elements of products of cubic prime
ideals and the older one-dimensional elements ``m-theta``.  It checks
every declared degree-one ideal valuation and proves that the undeclared norm
cofactor is supported on the ordinary factor base and the curve's bad-prime
set.  The resulting support transitions are exact; the search bounds which
found them remain experiments.
"""

from __future__ import annotations

from pathlib import Path
import sys

from sage.all import NumberField, PolynomialRing, QQ, ZZ, prime_range


sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients
from curve273_full_ideal_chain import SUPPORTS, build_relations, prime_ideal as raw_prime_ideal


PROTOCOL = "R30IDEALCHAIN"
FACTOR_BASE_BOUND = 1_000_000
S_RATIONAL = {
    2,
    3,
    5,
    7,
    13,
    31,
    41,
    47,
    53,
    67,
    379,
    4349,
    25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
}


T0, T1, T2, T3, T4, T5, T6, T7, T8, T9 = SUPPORTS


def sage_q(value):
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return QQ(ZZ(numerator)) / QQ(ZZ(denominator))


def support_text(support):
    return ",".join(f"{q}:{residue}" for q, residue in support)


def xor_support(left, right):
    return tuple(sorted(set(left).symmetric_difference(right)))


coefficients = short_coefficients()
A = ZZ(sage_q(coefficients[3]))
B = ZZ(sage_q(coefficients[4]))
polynomial_ring = PolynomialRing(QQ, "x")
x = polynomial_ring.gen()
defining_polynomial = x**3 + A * x + B
field = NumberField(defining_polynomial, "theta")
theta = field.gen()
small_primes = tuple(prime_range(2, FACTOR_BASE_BOUND + 1))


def prime_ideal(label):
    q, residue = map(ZZ, label)
    assert q.is_prime(proof=True)
    assert defining_polynomial(residue) % q == 0
    ideal = raw_prime_ideal(field, theta, label)
    assert ideal.is_prime() and ideal.norm() == q
    return ideal

def strip_allowed_support(norm, declared):
    cofactor = abs(ZZ(norm))
    for q, _ in declared:
        q = ZZ(q)
        assert cofactor % q == 0
        cofactor //= q
    for prime in small_primes:
        while cofactor % prime == 0:
            cofactor //= prime
    for prime in S_RATIONAL:
        prime = ZZ(prime)
        if prime <= FACTOR_BASE_BOUND:
            continue
        while cofactor % prime == 0:
            cofactor //= prime
    return cofactor


def verify_relation(name, alpha, declared):
    assert alpha != 0 and alpha.is_integral()
    norm = abs(ZZ(alpha.norm()))
    valuations = []
    for label in declared:
        ideal = prime_ideal(label)
        valuation = ZZ(alpha.valuation(ideal))
        assert valuation == 1
        valuations.append(valuation)
    assert strip_allowed_support(norm, declared) == 1
    print(
        f"{PROTOCOL}|relation={name}|declared={support_text(declared)}"
        + f"|valuations={','.join(map(str, valuations))}|norm_bits={norm.nbits()}"
        + "|remaining=FB_PLUS_S",
        flush=True,
    )


for relation in build_relations(field, theta):
    verify_relation(
        relation["name"],
        relation["alpha"],
        relation["declared"],
    )
    assert xor_support(
        relation["before"],
        relation["declared"],
    ) == tuple(sorted(relation["after"]))

print(
    f"{PROTOCOL}|chain={','.join(str(len(support)) for support in SUPPORTS)}"
    + f"|final={support_text(SUPPORTS[-1])}|status=PASS",
    flush=True,
)
