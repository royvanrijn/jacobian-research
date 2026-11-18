#!/usr/bin/env python3
"""Portable exact finite controls for well-factorable support and Rosser switching.

Uses only the standard library and exact integer/rational arithmetic.
Run from any working directory; write a receipt only with --output.
These controls do not verify infinite claims, source theorems,
asymptotic prime sums, the critical remainder, or binary coverage.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt, prod
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("Run the finite controls without Python -O.")

ROOT = Path(__file__).resolve().parents[3]
SOURCE_PATHS = (
    "research/arithmetic/ASYMMETRIC_WELL_FACTORABLE_SWITCH.md",
    "research/arithmetic/ASYMMETRIC_SEMILINEAR_SIEVE.md",
    "research/arithmetic/ASYMMETRIC_CRITICAL_SWITCHING.md",
    "research/arithmetic/ASYMMETRIC_AUXILIARY_AVERAGE.md",
    "research/arithmetic/ASYMMETRIC_SIGN_TRANSFER.md",
    "research/arithmetic/sources/iwaniec-half-dimensional-sieve-1976.pdf",
    "research/arithmetic/sources/fouvry-bombieri-vinogradov-II-1987.pdf",
    "research/arithmetic/sources/wu-chens-double-sieve-2004.pdf",
)
EXPECTED_PRIMARY = {
    SOURCE_PATHS[5]: "ca6e84aafb32ee62dcd5ce76679c65652bf677ccd4e36d200576024fac1eca67",
    SOURCE_PATHS[6]: "13dec04a3f215c0b5809dcefe77cd11fcad8a8982ddb9c471f1b6b08fad61c7f",
    SOURCE_PATHS[7]: "ebe75ef36478213bde937dce7d163bdd14b5b0352be54e085c9f696afc6e8a30",
}
COUNTS = Counter()


def require(ok, label):
    if not ok:
        raise RuntimeError(label)
    COUNTS[label] += 1


@lru_cache(maxsize=None)
def factors(n):
    if n < 1:
        raise ValueError("Only positive integers are factored.")
    result = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            result.append((p, e))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def primes_below(z):
    return [n for n in range(2, z) if factors(n) == ((n, 1),)]


def mu_from_primes(ps):
    return (-1) ** len(ps)


def lower_weight(ps, D):
    """Even-index squared conditions, written directly from source (4)."""
    ps = sorted(ps, reverse=True)
    for index in range(1, len(ps), 2):
        if prod(ps[:index]) * ps[index] ** 2 >= D:
            return 0
    return mu_from_primes(ps)


def upper_weight(ps, y):
    """Odd-index squared conditions, written independently from source (16)."""
    ps = sorted(ps, reverse=True)
    for index in range(0, len(ps), 2):
        if prod(ps[:index]) * ps[index] ** 2 >= y:
            return 0
    return mu_from_primes(ps)


def check_large_prime_switch():
    cases = 0
    for D, z in ((1000, 300), (10000, 700), (1000000, 2000)):
        pp = [p for p in primes_below(z) if p % 12 == 7]

        def walk(index, value, ps):
            nonlocal cases
            if ps and max(ps) ** 2 > D:
                p = max(ps)
                rest = [q for q in ps if q != p]
                low = lower_weight(ps, D)
                up = upper_weight(rest, F(D, p))
                require(low == -up, "large_prime_rosser_switch")
                if up:
                    require(all(q * q * p < D for q in rest),
                            "upper_cofactor_prime_bound")
                    require(prod(rest) * p < D, "upper_cofactor_product_bound")
                cases += 1
            for j in range(index, len(pp)):
                p = pp[j]
                if value * p >= D:
                    break
                walk(j + 1, value * p, ps + [p])

        walk(0, 1, [])
    require(cases == 1426, "switch_case_count")


def check_balanced_support():
    # The finite factors are arbitrary signed functions of balanced support.
    for Q in (36, 100, 256):
        h = isqrt(Q)
        left = {a: (-1) ** a for a in range(1, h + 1)}
        right = {b: 1 if b % 3 else -1 for b in range(1, h + 1)}
        conv = Counter()
        for a, x in left.items():
            for b, y in right.items():
                conv[a * b] += x * y
        for n in range(2, Q + 1):
            if max(p for p, _ in factors(n)) > h:
                require(conv[n] == 0, "balanced_large_prime_support")
    # Composite 7*43 is below level 1000 but has no balanced factorization.
    require(7 * 43 < 1000 and 43 * 43 > 1000,
            "large_prime_composite_inside_level")
    require(not any(301 % a == 0 and (301 // a) ** 2 <= 1000
                    for a in range(1, isqrt(1000) + 1)),
            "large_prime_composite_no_balanced_factors")
    # Source lower weight has a prime coefficient forbidden by WF support.
    require(lower_weight([43], 1000) == -1 and 43 * 43 > 1000,
            "control_lower_weight_not_balanced_factorable")


def check_actual_pairs():
    # Fixed alpha=1/10: for these primes the exact ceiling gives c=13.
    # Test the real parity-respecting pair with epsilon=1/100, B=0.
    for p in (7, 19, 31, 43, 67, 79, 103, 127):
        require(p % 12 == 7 and factors(p) == ((p, 1),),
                "actual_pair_forbidden_prime")
        c = 13
        require(1 < p and p <= 13 ** 10, "alpha_ceiling_gives_13")
        b = c * p * p
        N = b + 1
        require(N % 12 == 2, "actual_pair_admissible_target")
        require(all(q % 12 == 1 for q, _ in factors(c)),
                "actual_pair_good_cofactor")
        require(all(q % 3 == 1 for q, _ in factors(b)),
                "actual_pair_second_in_A3")
        require(sum(e for q, e in factors(b) if q % 12 == 7) == 2,
                "actual_pair_even_forbidden_parity")
        require([q for q, _ in factors(b) if q % 12 == 7] == [p],
                "actual_pair_only_sieve_divisors_1_p")
        require(p ** 4 > N and p ** 100 < N ** 49,
                "actual_pair_sqrtD_less_p_less_z")
        h = isqrt(isqrt(N))
        def balanced_coefficient(n):
            return sum((-1) ** a * (-1) ** (n // a)
                       for a in range(1, h + 1)
                       if n % a == 0 and n // a <= h)
        candidate = sum(balanced_coefficient(d) for d in (1, p))
        require(balanced_coefficient(1) == 1 and balanced_coefficient(p) == 0,
                "actual_pair_balanced_candidate_coefficients")
        require(candidate > 0 and p ** 100 < N ** 49,
                "control_positive_lambda1_fails_bad_pair")
    require((7 + 1) % 12 != 2, "control_single_prime_parity_test_not_actual")


def check_parameters():
    # Exact linear combination behind the full C1/C2 range ceiling.
    require(F(1, 9) * 2 + F(4, 9) * F(7, 4) == 1,
            "lp_r_coefficient")
    require(F(1, 9) + F(4, 9) * 2 == 1,
            "lp_s_coefficient")
    require(-F(1, 9) + F(4, 9) * F(3, 2) == F(5, 9),
            "lp_C1_nu_coefficient")
    require(F(4, 9) * F(5, 4) == F(5, 9),
            "lp_C2_nu_coefficient")
    for nu in (F(1, 100), F(1, 19), F(1, 12), F(1, 10)):
        r, s = (4 + 5 * nu) / 9, (1 - 10 * nu) / 9
        require(2 * r + s == 1 and F(7, 4) * r + 2 * s + F(5, 4) * nu == 1,
                "fouvry_C2_boundary_factorization")
        require(r + s == 5 * (1 - nu) / 9,
                "fouvry_large_residue_ceiling")
        require(F(3, 4) > r + s, "control_three_quarters_outside_source")
    delta = F(1, 1000)
    for nu in (F(1, 100), F(1, 19), F(1, 12),
               F(499, 5000), F(1, 10)):
        s = max(F(0), (1 - 10 * nu) / 9 - delta / 2)
        r = 5 * (1 - nu) / 9 - delta - s
        require(r >= 0 and s >= 0, "clamped_factors_at_least_one")
        require(2 * r + s <= 1 - delta
                and F(7, 4) * r + 2 * s + F(5, 4) * nu <= 1 - delta,
                "clamped_C2_source_conditions")
        require(r + s == 5 * (1 - nu) / 9 - delta,
                "clamped_total_level_unchanged")
    require((1 - 10 * F(1, 10)) / 9 - delta / 2 < 0,
            "control_unclamped_endpoint_gives_S_below_one")
    nu = F(1, 19)
    require(F(1, 2) + nu / 2 == F(4, 7) - 6 * nu / 7 == F(10, 19),
            "arbitrary_weight_optimal_exponent")
    require(F(4, 7) - 6 * F(1, 12) / 7 == F(1, 2),
            "arbitrary_weight_beyondhalf_endpoint")
    require(F(4, 7) - 6 * F(1, 10) / 7 < F(1, 2),
            "control_short_factor_too_long")
    eps = F(1, 100)
    for theta, alpha in ((F(1, 2), F(1, 10)),
                         (F(9, 10), F(1, 10)),
                         (F(97, 100), F(1, 20))):
        low, high = 4 * eps / (1 - 2 * eps), 2 / theta - 2
        require(low < alpha < high, "actual_pair_alpha_interval")
        require(theta / 2 < 1 / (2 + alpha) < F(1, 2) - eps,
                "actual_pair_power_separation")
        require((1 - 2 * eps) / theta > 1, "local_L1_ratio_positive")
    theta = F(99, 100)
    require(4 * eps / (1 - 2 * eps) >= 2 / theta - 2,
            "control_padding_at_or_above_z_squared_fails")
    require(F(1, 100) < F(1, 42), "fixed_alpha_epsilon_range")
    # Both the sign and D/p in the switched upper weight are indispensable.
    require(lower_weight([43], 1000) != upper_weight([], F(1000, 43)),
            "control_missing_switch_sign_fails")
    require(lower_weight([43, 7], 1000) != -upper_weight([7], 1000),
            "control_missing_division_by_p_fails")
    require(lower_weight([127, 7], 127 * 7 ** 2) == 0
            and upper_weight([7], F(127 * 7 ** 2, 127)) == 0,
            "strict_support_boundary")
    require((-1) ** 2 != lower_weight([127, 7], 127 * 7 ** 2),
            "control_non_strict_boundary_fails")


def provenance():
    hashes = {rel: hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
              for rel in SOURCE_PATHS}
    for rel, expected in EXPECTED_PRIMARY.items():
        require(hashes[rel] == expected, "immutable_primary_pin")
    return {
        "checker": Path(__file__).resolve().relative_to(ROOT).as_posix(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashes,
        "python_version": sys.version.split()[0],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write an optional JSON receipt.")
    args = parser.parse_args()
    if args.output is not None:
        protected = {(ROOT / rel).resolve() for rel in SOURCE_PATHS}
        protected.add(Path(__file__).resolve())
        if args.output.resolve() in protected:
            raise ValueError("Output would overwrite a proof, primary source, or checker.")
    pins = provenance()
    check_large_prime_switch()
    check_balanced_support()
    check_actual_pairs()
    check_parameters()
    receipt = {
        "result": "PASS",
        "scope": ("Finite exact support, target-pair, Rosser-switch and rational "
                  "parameter controls only; no infinite proposition, source theorem, "
                  "asymptotic prime sum, critical remainder or binary coverage is verified."),
        "provenance": pins,
        "checks": dict(sorted(COUNTS.items())),
        "switch_parameters": [[1000, 300], [10000, 700], [1000000, 2000]],
        "negative_controls": [
            "lower prime coefficient lies outside balanced WF support",
            "single forbidden prime does not satisfy the actual target parity",
            "positive lambda1 fails the actual bad even-parity pair",
            "three-quarters modulus exponent exceeds the primary C1/C2 range",
            "short exponent 1/10 loses the arbitrary-weight beyond-half range",
            "padding at or above z squared destroys the alpha interval",
            "omitting the switch sign or the level division by p fails",
            "replacing strict source inequalities by non-strict ones fails",
        ],
    }
    serialized = json.dumps(receipt, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")


if __name__ == "__main__":
    main()
