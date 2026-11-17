#!/usr/bin/env python3
"""Certify the full-sieve weighted survival constant D(1)^2.

The accompanying SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md defines D(s).
For primes p=5,7 mod12, its squared local factor is
(p^2-1)^4/[p^4(p^2+1)^2]; for p=11 mod12 it is (p^2-1)^2/p^4.
Primes 1 mod12 contribute one; primes 2,3 were already removed in h_6.

All arithmetic below is integer or Fraction arithmetic.  At scale S, the
updates floor(lo*num/den), ceil(hi*num/den) enclose each finite product.
For q=p^-2, -log of either squared forbidden-prime factor is at most
6q/(1-q)=6/(p^2-1), using -log(1-q)<=q/(1-q) and log(1+q)<=q.
Thus the infinite tail after an integer T is at least
exp(-6*sum_{n>T}1/(n^2-1)) >= 1-3*(1/T+1/(T+1)).
The sum telescopes.  This controls every omitted prime, without a prime
distribution assumption.  The result concerns an Euler-product constant;
it proves neither the averaged formula nor pointwise binary coverage.
"""
import argparse
import hashlib
import json
import math
import platform
from fractions import Fraction as F
from pathlib import Path

PRIME_CUTOFF = 1_000_000
SCALE = 10 ** 30


def require_equal(actual, expected, label):
    if actual != expected:
        raise RuntimeError(f"{label}: {actual!r} != {expected!r}")


def primes_up_to(limit):
    prime = bytearray(b"\x01") * (limit + 1)
    prime[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if prime[p]:
            prime[p*p::p] = b"\x00" * ((limit - p*p) // p + 1)
    return [p for p in range(2, limit + 1) if prime[p]]


def factor_squared(p):
    """Numerator and denominator for a prime p>=5."""
    q = p * p
    if p % 12 in (5, 7):
        return (q - 1) ** 4, q * q * (q + 1) ** 2
    if p % 12 == 11:
        return (q - 1) ** 2, q * q
    if p >= 5 and p % 12 == 1:
        return 1, 1
    raise ValueError("the argument must be a prime at least five")


def finite_product(primes, scale):
    lo = hi = scale
    counts = {str(r): 0 for r in (1, 5, 7, 11)}
    for p in primes:
        if p < 5:
            continue
        counts[str(p % 12)] += 1
        num, den = factor_squared(p)
        lo = lo * num // den
        hi = (hi * num + den - 1) // den
    return lo, hi, counts


def require_enclosure(exact, lo, hi, scale, label):
    if not F(lo, scale) <= exact <= F(hi, scale):
        raise RuntimeError(f"{label}: directed interval misses the exact value")


def controls():
    primes = primes_up_to(50)
    require_equal(primes, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
                  "small prime enumeration")
    for p, expected in ((5, F(82944, 105625)), (7, F(1327104, 1500625)),
                        (11, F(14400, 14641)), (13, F(1))):
        require_equal(F(*factor_squared(p)), expected, f"local factor at {p}")
    exact = math.prod(F(*factor_squared(p)) for p in primes if p >= 5)
    lo, hi, _ = finite_product(primes, SCALE)
    require_enclosure(exact, lo, hi, SCALE, "small exact product")
    t, n = 7, 53
    tail = sum((F(1, j*j - 1) for j in range(t + 1, n + 1)), F(0))
    telescoped = (F(1, t) + F(1, t+1) - F(1, n) - F(1, n+1)) / 2
    require_equal(tail, telescoped, "finite telescoping identity")
    corrupt_lo = exact.numerator * SCALE // exact.denominator + 1
    rejected = []
    trials = (
        ("missing_p5_denominator", lambda: require_equal(
            F(*factor_squared(5)), F(24, 25) ** 4, "false p5 formula")),
        ("unsquared_p11_factor", lambda: require_equal(
            F(*factor_squared(11)), F(120, 121), "false p11 factor")),
        ("doubled_telescoping_tail", lambda: require_equal(
            tail, 2 * telescoped, "false tail formula")),
        ("corrupt_directed_interval", lambda: require_enclosure(
            exact, corrupt_lo, hi, SCALE, "false interval")),
    )
    for label, trial in trials:
        try:
            trial()
        except RuntimeError:
            rejected.append(label)
        else:
            raise RuntimeError(f"negative control accepted: {label}")
    return {"exact_small_product": str(exact), "local_factor_primes": [5, 7, 11, 13],
            "small_prime_cutoff": 50, "negative_controls_rejected": rejected}


def run():
    checked_controls = controls()
    lo, hi, counts = finite_product(primes_up_to(PRIME_CUTOFF), SCALE)
    tail_lower = 1 - 3 * (F(1, PRIME_CUTOFF) + F(1, PRIME_CUTOFF + 1))
    lower, upper = F(lo, SCALE) * tail_lower, F(hi, SCALE)
    claimed_lower, claimed_upper = F(633187, 1_000_000), F(633191, 1_000_000)
    if not 0 < tail_lower <= 1 or not claimed_lower < lower <= upper < claimed_upper:
        raise RuntimeError("the claimed infinite-product interval was not certified")
    return {
        "schema": "simultaneous-norm-survival-checks/v1", "status": "pass",
        "scope": "Certified Euler-product constant only; no pointwise binary theorem.",
        "arithmetic": "integers and Fraction; no floating-point arithmetic",
        "python": platform.python_version(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "prime_cutoff": PRIME_CUTOFF, "scale": SCALE,
        "prime_counts_modulo_12": counts,
        "forbidden_prime_factor_count": sum(counts[str(r)] for r in (5, 7, 11)),
        "finite_product_lower_numerator": lo, "finite_product_upper_numerator": hi,
        "tail_factor_lower_bound": str(tail_lower),
        "infinite_product_lower_bound": str(lower),
        "infinite_product_upper_bound": str(upper),
        "claimed_strict_lower_bound": str(claimed_lower),
        "claimed_strict_upper_bound": str(claimed_upper),
        "controls": checked_controls,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the exact deterministic receipt")
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "prime_cutoff": PRIME_CUTOFF,
                      "forbidden_prime_factors": result["forbidden_prime_factor_count"],
                      "certified_D_squared_interval": "0.633187 < D(1)^2 < 0.633191",
                      "negative_controls_rejected": len(result["controls"]["negative_controls_rejected"])}))
