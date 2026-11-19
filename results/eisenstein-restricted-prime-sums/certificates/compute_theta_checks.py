#!/usr/bin/env python3
"""Exact checks of the A2+A2 theta identity and candidate local bounds.

No floating point is used in the inequalities.  Square roots are removed
only after checking the sign.  These finite checks are not a proof of the
proposed estimates.
"""

import argparse
import hashlib
import json
import math
import pathlib
import time


def arithmetic_tables(limit):
    sigma = [0] * (limit + 1)
    tau = [0] * (limit + 1)
    eisenstein_r = [0] * (limit + 1)
    eisenstein_r[0] = 1
    for d in range(1, limit + 1):
        chi = 0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
        for n in range(d, limit + 1, d):
            sigma[n] += d
            tau[n] += 1
            eisenstein_r[n] += 6 * chi
    return sigma, tau, eisenstein_r


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=3000)
    parser.add_argument("--max-prime", type=int, default=53)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("theta_checks.json"))
    args = parser.parse_args()
    start = time.monotonic()
    sigma, tau, r = arithmetic_tables(args.limit)
    primes = [p for p in range(2, args.max_prime + 1) if p % 3 == 2 and is_prime(p)]
    theta = [sum(r[a] * r[n-a] for a in range(n + 1)) for n in range(args.limit + 1)]
    theta_errors = []
    local_errors = []
    tested_local = 0
    for n in range(1, args.limit + 1):
        expected = 12 * sigma[n] - (36 * sigma[n // 3] if n % 3 == 0 else 0)
        if theta[n] != expected:
            theta_errors.append({"n": n, "actual": theta[n], "expected": expected})
        if n % 3 != 2:
            continue
        for p in primes:
            single_pair = sum(r[u] * r[n-p*p*u] for u in range(n // (p*p) + 1))
            both_pairs = theta[n // (p*p)] if n % (p*p) == 0 else 0
            bad = 2 * single_pair - both_pairs
            tested_local += 1
            if n % p == 0:
                expected_bad = theta[n // (p*p)] if n % (p*p) == 0 else 0
                if bad != expected_bad:
                    local_errors.append({"n": n, "p": p, "bad": bad, "expected_bad": expected_bad})
            else:
                # bad/(12sigma) <= 2/[p(p-1)] + tau sqrt(n)/sigma.
                # Multiply by 12sigma p(p-1), then square positive sides.
                denominator = p * (p - 1)
                left = bad * denominator - 24 * sigma[n]
                right_factor = 12 * denominator * tau[n]
                if left > 0 and left * left > right_factor * right_factor * n:
                    local_errors.append({"n": n, "p": p, "bad": bad, "theta": theta[n], "left_squared": left * left, "right_squared": right_factor * right_factor * n})
    # Independent direct lattice enumeration checks the input binary theta
    # coefficients, so a mistaken divisor-sum convention is not hidden.
    small_limit = min(args.limit, 200)
    direct = [0] * (small_limit + 1)
    bound = math.isqrt(2 * small_limit) + 2
    for x in range(-bound, bound + 1):
        for y in range(-bound, bound + 1):
            value = x*x - x*y + y*y
            if value <= small_limit:
                direct[value] += 1
    binary_errors = [{"n": n, "direct": direct[n], "formula": r[n]} for n in range(small_limit + 1) if direct[n] != r[n]]
    receipt = {
        "scope": "Finite exact checks, not unrestricted proof",
        "theta_form": "x1^2-x1*x2+x2^2+x3^2-x3*x4+x4^2",
        "theta_formula": "r_Q(n)=12 sigma(n)-36 sigma(n/3), with sigma(noninteger)=0",
        "limit": args.limit,
        "theta_identity_checks": args.limit,
        "binary_theta_direct_checks": small_limit + 1,
        "candidate_bound": "b_p/r_Q <= 2/[p(p-1)] + tau(n)*sqrt(n)/sigma(n) when p does not divide n; b_p=r_Q(n/p^2) when p^2 divides n and zero when p divides n exactly once",
        "candidate_bound_target_scope": "n=2 mod3, p=2 mod3 prime",
        "tested_primes": primes,
        "local_bound_checks": tested_local,
        "theta_identity_errors": theta_errors,
        "binary_theta_errors": binary_errors,
        "local_bound_errors": local_errors,
        "script_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    if theta_errors or binary_errors or local_errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
