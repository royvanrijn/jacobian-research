#!/usr/bin/env python3
"""Portable finite controls for the asymmetric auxiliary-average identities.

Run from any directory without Python -O. This uses only the standard
library. It checks selected exact family decompositions, local factors,
CRT substitutions, and full/truncated finite-kernel identities. It does
not verify dispersion, error exponents, distribution, or global asymptotics.
No file is written unless --output is explicitly supplied.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import gcd, isqrt, prod
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("These finite controls require assertions; run without Python -O.")

ROOT = Path(__file__).resolve().parents[3]
SOURCE_PATHS = (
    "research/arithmetic/ASYMMETRIC_AUXILIARY_AVERAGE.md",
    "research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md",
    "research/arithmetic/ASYMMETRIC_SEMILINEAR_SIEVE.md",
    "research/arithmetic/FIXED_PERIODIC_DIVISOR_BOXES.md",
    "research/arithmetic/sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf",
)


@lru_cache(maxsize=None)
def factorization(n):
    if n < 1:
        raise ValueError("Factorization is used only for positive integers.")
    result = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            exponent += 1
            n //= p
        if exponent:
            result.append((p, exponent))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


@lru_cache(maxsize=None)
def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def tau(n):
    return prod(e + 1 for _, e in factorization(n))


def omega(n):
    return len(factorization(n))


def mobius(n):
    return 0 if any(e > 1 for _, e in factorization(n)) else (-1) ** omega(n)


def u_weight(n):
    return tau(n) ** 2 if all(p % 12 == 1 for p, _ in factorization(n)) else 0


def v_weight(n):
    return tau(n) if all(p % 3 == 1 for p, _ in factorization(n)) else 0


def chi3_mod6(n):
    return 1 if n % 6 == 1 else -1 if n % 6 == 5 else 0


def eisenstein_core(n):
    return sum(chi3_mod6(d) for d in divisors(n)) if gcd(n, 6) == 1 else 0


def strict_hyperbola(n):
    assert n % 6 == 1
    root = isqrt(n)
    return 2 * sum(chi3_mod6(q) for q in divisors(n) if q * q < n) + (
        chi3_mod6(root) if root * root == n else 0
    )


@lru_cache(maxsize=None)
def shifted_second_kernel(g, bound):
    assert u_weight(g) > 0 and bound >= 1
    kernel = [F(0) for _ in range(bound + 1)]
    for n in range(1, bound + 1):
        kernel[n] = F(v_weight(g * n), v_weight(g)) - sum(
            kernel[d] * eisenstein_core(n // d) for d in divisors(n) if d < n
        )
        if kernel[n]:
            assert n % 12 == 1
            assert all(p % 12 != 7 for p, _ in factorization(n))
    return tuple(kernel)


@lru_cache(maxsize=None)
def family_count(N, c):
    assert N % 12 == 2 and all(p % 12 == 7 for p, _ in factorization(c))
    return sum(
        u_weight(N - c * t) * v_weight(t) for t in range(1, (N - 1) // c + 1)
    )


def h_factor(c):
    assert all(p % 12 == 7 for p, _ in factorization(c))
    return prod(F(p * p, p * p - p + 1) for p, _ in factorization(c))


def family_model(N, c, formal_X):
    return formal_X * h_factor(c) / c if gcd(c, N) == 1 else F(0)


def squarefree_p7(n):
    return mobius(n) != 0 and all(p % 12 == 7 for p, _ in factorization(n))


def sieve_density(d, N):
    assert squarefree_p7(d)
    if gcd(d, N) > 1:
        return F(0)
    return prod(F(2 * p - 1, p * p - p + 1) for p, _ in factorization(d))


def check_family_decompositions():
    records = []
    for N in (14, 26, 38, 50, 62, 74, 86, 98, 170, 182, 338, 506):
        for c in (1, 7, 19, 31, 49):
            direct = family_count(N, c)
            if gcd(c, N) > 1:
                assert direct == 0, (N, c)
            else:
                expanded = F(0)
                for g in divisors(N):
                    if not u_weight(g):
                        continue
                    assert gcd(g, c) == 1
                    M = N // g
                    kernel = shifted_second_kernel(g, M)
                    for m in range(1, M):
                        if gcd(m, M) != 1 or (M - m) % c:
                            continue
                        first = F(u_weight(g * m), u_weight(g))
                        if not first:
                            continue
                        s = (M - m) // c
                        second = F(v_weight(g * s), v_weight(g))
                        hyperbola_expansion = F(0)
                        for v in divisors(s):
                            if not kernel[v]:
                                continue
                            assert gcd(c * v, M) == 1
                            t = s // v
                            assert t % 6 == 1
                            assert strict_hyperbola(t) == eisenstein_core(t)
                            hyperbola_expansion += kernel[v] * strict_hyperbola(t)
                        assert second == hyperbola_expansion, (N, c, g, m)
                        expanded += (
                            u_weight(g) * v_weight(g) * first * hyperbola_expansion
                        )
                assert direct == expanded, (N, c, direct, expanded)
            records.append(
                {"N": N, "c": c, "count": direct,
                 "zero_local_obstruction": gcd(c, N) > 1}
            )
    return records


def check_local_prime_powers():
    records = []
    for p in (7, 19, 31, 43, 67):
        x = F(1, p)
        q = 1 - x
        first_mean, second_mean = q, 1 / q
        unit_target_integral = first_mean + second_mean - 1
        expected_h = 1 / (1 - x + x * x)
        assert second_mean / unit_target_integral == expected_h
        for exponent in range(1, 9):
            completed_q_factor = x ** exponent / q ** 2
            assert first_mean * completed_q_factor == x ** exponent * second_mean
            ratio = first_mean * completed_q_factor / unit_target_integral
            assert ratio == x ** exponent * expected_h
            records.append({"p": p, "exponent": exponent, "ratio": str(ratio)})
    return records


def check_auxiliary_crt():
    crt_cases = 0
    injective_assignments = 0
    for v in (1, 13, 25):
        for auxiliary in (1, 7, 19, 49):
            images = {}
            for c in (1, 7, 19, 31, 49, 133, 343):
                support = {p for p, _ in factorization(c * v)}
                if any(p not in support for p, _ in factorization(auxiliary)):
                    continue
                d = c * v * auxiliary
                assert d not in images
                images[d] = c
                assert tau(c) <= tau(d)
                injective_assignments += 1
                for q_prime in (1, 5, 11, 13, 17, 23):
                    if gcd(q_prime, d) != 1:
                        continue
                    full_modulus = d * q_prime
                    for inverse_multiplier in (1, 13, 25, 121):
                        if gcd(inverse_multiplier, full_modulus) != 1:
                            continue
                        M = 2 * full_modulus + 1
                        residue = (
                            M * pow(inverse_multiplier, -1, full_modulus)
                        ) % full_modulus if full_modulus > 1 else 0
                        for n in (residue, residue + 1, residue + full_modulus):
                            full = (inverse_multiplier * n - M) % full_modulus == 0
                            at_d = (
                                (n - M * pow(inverse_multiplier, -1, d)) % d == 0
                            ) if d > 1 else True
                            at_q = (
                                (n - M * pow(inverse_multiplier, -1, q_prime))
                                % q_prime == 0
                            ) if q_prime > 1 else True
                            assert full == (at_d and at_q)
                            crt_cases += 1
    return {"crt_cases": crt_cases, "injective_assignments": injective_assignments}


def check_negative_controls():
    # Leaving q_aux out of d accepts a spurious hyperbola progression.
    N, c, v, auxiliary, q_prime, m = 506, 7, 1, 7, 1, 37
    assert u_weight(m) > 0 and m <= N - c * v * (auxiliary * q_prime) ** 2
    wrong_accepts = (m - N) % (c * v) == 0
    correct_accepts = (m - N) % (c * v * auxiliary * q_prime) == 0
    assert wrong_accepts and not correct_accepts

    # A v sharing M has zero actual progression and a nonunit source w.
    M, bad_v = 50, 25
    assert shifted_second_kernel(1, M)[bad_v] == -1
    actual = [
        a for a in range(1, M - bad_v + 1)
        if (a - M) % bad_v == 0 and gcd(a, M) == 1
    ]
    assert actual == []
    assert M - bad_v == 25 and M % bad_v == 0
    return [
        {"omission": "q_aux factor in auxiliary d", "N": N, "c": c,
         "v": v, "q_aux": auxiliary, "m": m,
         "incorrect_accepts": wrong_accepts, "correct_rejects": not correct_accepts},
        {"omission": "(v,M)=1 restriction", "M": M, "v": bad_v,
         "kernel_coefficient": "-1", "actual_progression_count": len(actual),
         "incorrect_principal_volume": M - bad_v, "source_residue": M % bad_v,
         "source_unit_condition_satisfied": gcd(M, bad_v) == 1},
    ]


def check_finite_kernel_compression():
    identity_checks = 0
    records = []
    for N in (14, 38, 98, 170, 338, 506):
        for D in sorted({min(50, N), min(200, N)}):
            ds = tuple(d for d in range(1, D) if squarefree_p7(d))
            # These are formal rational scalars, not estimates of the infinite
            # singular series. The exact identities are linear in this scalar.
            for formal_X in (F(0), F(7, 3), F(N)):
                remainders = {}
                for d in ds:
                    actual = sum(
                        u_weight(N - b) * v_weight(b) for b in range(d, N, d)
                    )
                    exact_kernel = sum(
                        F(2 ** omega(d) * mobius(k), 2 ** omega(k))
                        * family_count(N, d * k) for k in divisors(d)
                    )
                    exact_model = sum(
                        F(2 ** omega(d) * mobius(k), 2 ** omega(k))
                        * family_model(N, d * k, formal_X) for k in divisors(d)
                    )
                    assert actual == exact_kernel
                    assert exact_model == formal_X * sieve_density(d, N)
                    remainders[d] = actual - formal_X * sieve_density(d, N)
                    identity_checks += 1
                lambdas = {
                    d: 1 if r > 0 else -1 if r < 0 else 0
                    for d, r in remainders.items()
                }
                absolute_sum = sum(abs(r) for r in remainders.values())
                assert absolute_sum == sum(
                    lambdas[d] * remainders[d] for d in ds
                )
                full = F(0)
                for d in ds:
                    for k in divisors(d):
                        ell = d // k
                        full += (
                            mobius(k) * 2 ** omega(ell) * lambdas[d]
                            * (family_count(N, d * k)
                               - family_model(N, d * k, formal_X))
                        )
                assert full == absolute_sum

                for L in sorted({2, 7, D, 2 * D}):
                    truncated, tail = F(0), F(0)
                    preimages = {}
                    omitted_pairs = 0
                    for d in ds:
                        if gcd(d, N) > 1:
                            assert remainders[d] == 0
                            continue
                        for k in divisors(d):
                            ell = d // k
                            assert squarefree_p7(k) and squarefree_p7(ell)
                            assert gcd(k, ell) == 1
                            c = k * k * ell
                            coefficient = mobius(k) * 2 ** omega(ell) * lambdas[d]
                            term = coefficient * (
                                family_count(N, c) - family_model(N, c, formal_X)
                            )
                            if k <= L:
                                assert c < D * L
                                assert c not in preimages
                                preimages[c] = (k, ell)
                                fs = factorization(c)
                                assert all(e in (1, 2) for _, e in fs)
                                assert prod(p for p, e in fs if e == 2) == k
                                assert prod(p for p, e in fs if e == 1) == ell
                                assert abs(coefficient) <= tau(c)
                                truncated += term
                            else:
                                omitted_pairs += 1
                                tail += term
                    assert truncated + tail == absolute_sum
                    if L > D:
                        assert omitted_pairs == 0 and tail == 0
                        assert truncated == full
                    records.append(
                        {"N": N, "D": D, "L": L, "formal_X": str(formal_X),
                         "L_greater_than_D": L > D,
                         "full_identity": True, "unique_preimages": True,
                         "coefficient_tau_bound": True,
                         "omitted_pairs": omitted_pairs, "tail": str(tail)}
                    )
    return {
        "model_and_kernel_identity_checks": identity_checks,
        "cutoff_case_checks": len(records),
        "L_greater_than_D_cases": sum(r["L_greater_than_D"] for r in records),
        "records": records,
    }


def run_checks():
    family = check_family_decompositions()
    local = check_local_prime_powers()
    return {
        "family_exact_case_count": len(family),
        "family_exact_cases": family,
        "local_prime_power_case_count": len(local),
        "local_prime_power_cases": local,
        "auxiliary_crt": check_auxiliary_crt(),
        "negative_controls": check_negative_controls(),
        "finite_kernel_compression": check_finite_kernel_compression(),
    }


def provenance():
    checker = Path(__file__).resolve()
    source_hashes = {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in SOURCE_PATHS
    }
    return {
        "checker": checker.relative_to(ROOT).as_posix(),
        "checker_sha256": hashlib.sha256(checker.read_bytes()).hexdigest(),
        "source_sha256": source_hashes,
        "python_version": sys.version.split()[0],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Exact finite arithmetic controls; no dispersion or asymptotic verification."
    )
    parser.add_argument("--output", type=Path,
                        help="Optional JSON receipt path; otherwise print only.")
    args = parser.parse_args()
    receipt = {
        "scope": (
            "Finite exact family decompositions, local factors, auxiliary CRT, "
            "negative controls, and full/truncated kernel identities only. "
            "No dispersion, analytic error estimate, distribution theorem, "
            "or global asymptotic is verified."
        ),
        "provenance": provenance(),
        "checks": run_checks(),
        "result": "PASS",
    }
    serialized = json.dumps(receipt, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")


if __name__ == "__main__":
    main()
