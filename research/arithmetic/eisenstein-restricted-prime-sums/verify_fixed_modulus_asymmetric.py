#!/usr/bin/env python3
"""Finite exact controls for fixed-modulus asymmetric sums.

This checks selected character, local-density, kernel and fibre identities.
It does not verify the analytic estimates, imported theorems or asymptotics.
Run without -O. By default the receipt is printed; only --output writes it.
The optional --source override is for reviewing a source before promotion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb, gcd, isqrt, lcm, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
F = Fraction
DEFAULT_SOURCE = "research/arithmetic/ASYMMETRIC_FIXED_MODULUS_THEORY.md"
RETAINED_INPUTS = (
    "research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md",
    "research/arithmetic/ASYMMETRIC_UNIFORM_COEFFICIENTS.md",
    "research/arithmetic/ASYMMETRIC_SEMIGROUP_SCALE.md",
    "research/arithmetic/FIXED_PERIODIC_DIVISOR_BOXES.md",
    "research/arithmetic/sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf",
    "research/arithmetic/sources/heath-brown-d3-progressions-1986.pdf",
)
FIXED_QS = (3, 4, 5, 6, 8, 9, 12, 15, 18, 20, 21, 24, 27, 30, 36, 42, 60)
COUNTS: Counter[str] = Counter()
NEGATIVE_CONTROLS: list[dict[str, str]] = []


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def reject_wrong(label: str, expected: object, wrong: object) -> None:
    require(expected != wrong, f"negative control was not rejected: {label}")
    NEGATIVE_CONTROLS.append(
        {"label": label, "expected": str(expected), "wrong": str(wrong)}
    )


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * ((limit - p * p) // p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


PRIMES = primes_up_to(2500)


@lru_cache(maxsize=None)
def factorization(n: int) -> tuple[tuple[int, int], ...]:
    require(n >= 1, f"factorization requires a positive integer: {n}")
    factors: list[tuple[int, int]] = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            exponent += 1
            n //= p
        if exponent:
            factors.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def divisors(n: int) -> list[int]:
    answer = [1]
    for p, exponent in factorization(n):
        answer = [d * p**e for d in answer for e in range(exponent + 1)]
    return answer


def totient(n: int) -> int:
    return sum(gcd(a, n) == 1 for a in range(n))


def primitive_quadratic(discriminant: int, n: int) -> int:
    """Explicit small primary characters used in the finite case list."""
    if discriminant == -3:
        return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)
    if discriminant == -4:
        return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
    if discriminant == 5:
        return 0 if n % 5 == 0 else (1 if n % 5 in (1, 4) else -1)
    if discriminant == 8:
        return 0 if n % 2 == 0 else (1 if n % 8 in (1, 7) else -1)
    if discriminant == -8:
        return primitive_quadratic(-4, n) * primitive_quadratic(8, n)
    if discriminant == 12:
        return primitive_quadratic(-4, n) * primitive_quadratic(-3, n)
    raise CheckFailure(f"unsupported finite character: {discriminant}")


@dataclass(frozen=True)
class CharacterCase:
    q: int
    discriminant: int

    @property
    def period(self) -> int:
        return lcm(2, self.q)

    @property
    def kappa(self) -> int:
        return totient(self.period)

    def chi(self, n: int) -> int:
        if gcd(n, self.period) != 1:
            return 0
        return primitive_quadratic(self.discriminant, n)

    @property
    def label(self) -> str:
        return f"q={self.q}, D={self.discriminant}"


CASES = tuple(
    CharacterCase(q, d)
    for q, d in (
        (3, -3), (4, -4), (5, 5), (8, -4), (8, 8), (8, -8),
        (12, 12), (12, -3), (12, -4), (15, 5), (20, 5),
        (21, -3), (24, 8), (24, -8),
    )
)


def d_kappa(kappa: int, exponent: int) -> int:
    return comb(exponent + kappa - 1, kappa - 1)


@lru_cache(maxsize=None)
def first_weight(period: int, n: int) -> int:
    kappa = totient(period)
    factors = factorization(n)
    if any(p % period != 1 for p, _ in factors):
        return 0
    return prod(d_kappa(kappa, e) for _, e in factors)


@lru_cache(maxsize=None)
def second_weight(case: CharacterCase, n: int) -> int:
    factors = factorization(n)
    if any(case.chi(p) != 1 for p, _ in factors):
        return 0
    return prod(e + 1 for _, e in factors)


def second_base(case: CharacterCase, n: int) -> int:
    if gcd(n, case.period) != 1:
        return 0
    return prod(e + 1 if case.chi(p) == 1 else int(e % 2 == 0)
                for p, e in factorization(n))


def second_kernel(case: CharacterCase, shift: int, n: int) -> Fraction:
    shift_factors = dict(factorization(shift))
    answer = F(1)
    for p, exponent in factorization(n):
        if gcd(p, case.period) != 1:
            return F(0)
        if p in shift_factors:
            answer *= -F(shift_factors[p], shift_factors[p] + 1) if exponent == 1 else 0
        elif case.chi(p) == -1:
            answer *= -1 if exponent == 2 else 0
        else:
            return F(0)
    return answer


def shifted_second(case: CharacterCase, shift: int, n: int) -> Fraction:
    denominator = second_weight(case, shift)
    require(denominator > 0, f"inadmissible shift in control: {case.label}, {shift}")
    return F(second_weight(case, shift * n), denominator)


def multiplicative_order(p: int, period: int) -> int:
    value = p % period
    exponent = 1
    while value != 1:
        value = value * p % period
        exponent += 1
        require(exponent <= totient(period), f"order did not close: p={p}, L={period}")
    return exponent


def check_first_kernels() -> None:
    periods = sorted({lcm(2, q) for q in FIXED_QS})
    for q in FIXED_QS:
        period = lcm(2, q)
        require(totient(q) == totient(period), f"totient mismatch at q={q}")
        for p in PRIMES:
            require((p % q == 1) == (p % period == 1), f"identity class mismatch: q={q}, p={p}")
            COUNTS["identity_prime_classes"] += 1
    for period in periods:
        kappa = totient(period)
        require(kappa >= 2 and kappa % 2 == 0, f"unexpected kappa: L={period}")
        for p in PRIMES:
            if gcd(p, period) != 1:
                continue
            order = multiplicative_order(p, period)
            require(kappa % order == 0, f"order does not divide kappa: {period}, {p}")
            if order == 1:
                continue
            multiplicity = kappa // order
            for degree in range(16):
                coefficient = sum(
                    (-1)**a * comb(multiplicity, a)
                    * comb((degree - a * order) // order + multiplicity - 1, multiplicity - 1)
                    for a in range(multiplicity + 1)
                    if degree >= a * order and (degree - a * order) % order == 0
                )
                require(coefficient == int(degree == 0), f"bad Euler cancellation: L={period}, p={p}, degree={degree}")
            require(pow(p, order, period) == 1, "first kernel residue is not 1")
            COUNTS["bad_order_kernel_primes"] += 1
        for shift in (0, 1, 2, 10, 10000):
            for exponent in range(20):
                ratio = F(d_kappa(kappa, shift + exponent), d_kappa(kappa, shift))
                require(ratio <= d_kappa(kappa, exponent) <= kappa**exponent, "shifted divisor bound failed")
                correction = sum(
                    (-1)**a * comb(kappa, a)
                    * F(d_kappa(kappa, shift + exponent - a), d_kappa(kappa, shift))
                    for a in range(min(kappa, exponent) + 1)
                )
                if exponent >= kappa:
                    require(correction == 0, "shift correction did not terminate")
                COUNTS["shifted_first_polynomial_and_weight_cases"] += 1
        # Full target deletion replaces the shifted series, not its reciprocal.
        for degree in range(20):
            coefficient = sum(
                (-1)**a * comb(kappa, a) * d_kappa(kappa, degree - a)
                for a in range(min(kappa, degree) + 1)
            )
            require(coefficient == int(degree == 0), "deleted first local factor did not cancel")
            COUNTS["deleted_first_local_cases"] += 1


def check_characters_and_hyperbola() -> list[dict[str, object]]:
    details: list[dict[str, object]] = []
    for case in CASES:
        period, kappa = case.period, case.kappa
        allowed = [a for a in range(period) if case.chi(a) == 1]
        require(len(allowed) == kappa // 2, f"quadratic kernel size: {case.label}")
        for a in range(period):
            for b in range(period):
                require(case.chi(a * b) == case.chi(a) * case.chi(b), f"character multiplicativity: {case.label}")
                COUNTS["character_homomorphism_cases"] += 1
        for v in allowed:
            density = period * sum((F(2, kappa) for b in allowed if (1 + v * b - 2) % period == 0), F(0))
            require(density == F(2 * period, kappa), f"small-modulus density: {case.label}")
            COUNTS["fixed_residue_density_cases"] += 1
        for p in PRIMES:
            if p % period == 1:
                require(case.chi(p) == 1, f"Cq not contained in Achi: {case.label}")
        for t in range(1, 3001):
            if case.chi(t) != 1:
                continue
            strict = 2 * sum(case.chi(d) for d in range(1, isqrt(t) + 1) if d * d < t and t % d == 0)
            square = case.chi(isqrt(t)) if isqrt(t)**2 == t else 0
            require(second_base(case, t) == strict + square, f"strict pairing failed: {case.label}, t={t}")
            COUNTS["strict_hyperbola_cases"] += 1
        shifts = [1] + [p for p in PRIMES if p < 30 and case.chi(p) == 1][:2]
        if case == CharacterCase(8, -8):
            shifts += [3, 9]
        for shift in sorted(set(shifts)):
            for n in range(1, 501):
                reconstructed = sum((second_kernel(case, shift, d) * second_base(case, n // d) for d in divisors(n)), F(0))
                require(reconstructed == shifted_second(case, shift, n), f"second kernel failed: {case.label}, C={shift}, n={n}")
                coefficient = second_kernel(case, shift, n)
                if coefficient:
                    require(gcd(n, period) == 1 and case.chi(n) == 1, "second kernel has forbidden support")
                    COUNTS["nonzero_second_kernel_support_cases"] += 1
                COUNTS["shifted_second_kernel_cases"] += 1
        details.append({"q": case.q, "L": period, "kappa": kappa,
                        "discriminant": case.discriminant, "chi_minus_one": case.chi(-1),
                        "chi_three": case.chi(3), "small_factor": str(F(2 * period, kappa))})
    for period in sorted({case.period for case in CASES}):
        for v in range(period):
            if gcd(v, period) != 1:
                require(not any((1 + v * b - 2) % period == 0 for b in range(period)), "nonunit coefficient survived")
                COUNTS["nonunit_second_coefficient_zero_cases"] += 1
    return details


def good_local_integral(p: int, kappa: int, target_valuation: int) -> Fraction:
    x, q = F(1, p), 1 - F(1, p)
    r = target_valuation
    a = lambda e: d_kappa(kappa, e)
    b = lambda e: e + 1
    tail_a = q**(-kappa) - sum(a(i) * x**i for i in range(r + 1))
    tail_b = q**(-2) - sum(b(i) * x**i for i in range(r + 1))
    return (q * sum(a(i) * b(i) * x**i for i in range(r))
            + x**r * (1 - 2 * x) * a(r) * b(r)
            + q * (b(r) * tail_a + a(r) * tail_b))


def direct_good_integral(p: int, kappa: int, r: int) -> Fraction:
    """Enumerate residues and integrate the two zero cylinders exactly."""
    modulus, target = p**(r + 1), p**r
    x, q = F(1, p), 1 - F(1, p)
    def valuation(n: int) -> int:
        exponent = 0
        while n and n % p == 0:
            exponent += 1
            n //= p
        return exponent
    interior = sum(d_kappa(kappa, valuation(a)) * (valuation((target - a) % modulus) + 1)
                   for a in range(1, modulus) if (target - a) % modulus)
    tail_a = q**(-kappa) - sum(d_kappa(kappa, i) * x**i for i in range(r + 1))
    tail_b = q**(-2) - sum((i + 1) * x**i for i in range(r + 1))
    return F(interior, modulus) + q * ((r + 1) * tail_a + d_kappa(kappa, r) * tail_b)


def direct_three_integral(split: bool, r: int) -> Fraction:
    modulus, target = 3**(r + 1), 3**r
    interior = 0
    for a in range(1, modulus):
        if a % 3 == 0:
            continue
        b = (target - a) % modulus
        if b:
            exponent = 0
            while b % 3 == 0:
                exponent += 1
                b //= 3
            interior += exponent + 1 if split else int(exponent == 0)
    x, q = F(1, 3), F(2, 3)
    h_b = q**(-2) if split else F(1)
    tail_b = h_b - sum((i + 1 if split or i == 0 else 0) * x**i for i in range(r + 1))
    return F(interior, modulus) + (q * tail_b if r == 0 else 0)


def check_local_integrals() -> list[dict[str, object]]:
    for period in (6, 12, 18, 30, 42):
        kappa = totient(period)
        p = next(p for p in PRIMES if p % period == 1)
        for r in range(3):
            require(direct_good_integral(p, kappa, r) == good_local_integral(p, kappa, r), "independent good-prime integral failed")
            COUNTS["independent_good_prime_integrals"] += 1
    three_details: list[dict[str, object]] = []
    for case in CASES:
        sign = case.chi(3)
        if not sign:
            continue
        unit = F(7, 6) if sign == 1 else F(1, 3)
        for r in range(6):
            expected = unit if r == 0 else F(2, 3)
            require(direct_three_integral(sign == 1, r) == expected, f"p3 integral failed: {case.label}, r={r}")
            COUNTS["generic_prime_three_integrals"] += 1
        ratio = F(2, 3) / unit
        require(ratio == (F(4, 7) if sign == 1 else F(2)), "p3 target ratio failed")
        require(F(2) - (1 + F(2, 3)) == F(1, 3), "p3 allowed-prime compensation failed")
        three_details.append({"q": case.q, "discriminant": case.discriminant,
                              "I0": str(unit), "Ipositive": "2/3", "target_ratio": str(ratio)})
    for period in sorted({lcm(2, q) for q in FIXED_QS}):
        kappa = totient(period)
        require(F(kappa + 2, kappa) + 2 * (F(1, 2) - F(1, kappa)) == 2, "rough dimension failed")
        require(2 - (F(2, kappa) + F(1, 2) - F(1, kappa)) == F(3, 2) - F(1, kappa), "restoration exponent failed")
        require(2 - F(2, kappa) == 2 * (1 - F(1, kappa)), "restricted restoration exponent failed")
        for p in PRIMES:
            if p < 3 or gcd(p, period) != 1:
                continue
            e, a = F(p - 1, p - 2), F(p - 1, p * (p - 2))
            require(F(p - 1, p * p) / (1 - F(2, p)) == a, "directive normalization failed")
            if p % period == 1:
                for r in (1, 2, 5, 20):
                    previous = F(1) if r == 1 else e
                    require(e + F(1, p) * previous >= 1 + 2 * a, "common-divisor target compensation failed")
                    COUNTS["good_target_compensation_cases"] += 1
            else:
                require(e - (1 + a) == F(1, p * (p - 2)), "nongood allowed compensation failed")
                COUNTS["nongood_target_compensation_cases"] += 1
    return three_details


def check_fixed_prime_stripping() -> None:
    def old_eisenstein(n: int) -> int:
        factors = factorization(n)
        return 0 if any(p % 3 != 1 for p, _ in factors) else prod(e + 1 for _, e in factors)
    literal = CharacterCase(21, -3)
    for shift, stripped in ((7, 1), (49 * 13, 13)):
        for n in range(1, 501):
            old = F(old_eisenstein(shift * n), old_eisenstein(shift)) if gcd(n, 42) == 1 else F(0)
            require(old == shifted_second(literal, stripped, n), "fixed-prime stripping failed")
            COUNTS["fixed_prime_shift_stripping_cases"] += 1
    require(second_weight(literal, 7) == 0, "a literal fixed-prime shift must be inadmissible")
    reject_wrong("undeleted shifted-7 Euler factor at L42", F(1), F(91, 72))


def check_restricted_good_fibres() -> list[dict[str, object]]:
    """Independent finite tuple recovery for Section10, at fixed w,z."""
    summaries: list[dict[str, object]] = []
    target_max = 3000
    for case, w, z in ((CharacterCase(5, 5), 11, 13),
                       (CharacterCase(8, -4), 17, 19),
                       (CharacterCase(8, -8), 17, 19)):
        period, kappa = case.period, case.kappa
        good_prime = next(p for p in PRIMES if p % period == 1)
        require(good_prime == w, "fibre cases assume one small good prime")
        small = [p for p in PRIMES if p <= z]
        small_w = [p for p in small if p <= w]
        first_values = [n for n in range(1, target_max + 1) if first_weight(period, n)]
        seen: dict[tuple[int, int], tuple[int, int, int]] = {}
        count = high_target = with_excluded = 0
        for target in range(2, target_max + 1, period):
            for d in (1, good_prime):
                if target % d:
                    continue
                residual_target = target // d
                for u in (1, good_prime):
                    for v in (1, good_prime):
                        if gcd(d, u * v) != 1 or gcd(u, v) != 1 or gcd(u * v, target) != 1:
                            continue
                        for x in first_values:
                            if x >= residual_target:
                                break
                            y = residual_target - x
                            if not second_weight(case, y):
                                continue
                            fx, fy = dict(factorization(x)), dict(factorization(y))
                            if u > 1 and fx.get(good_prime, 0) != 1:
                                continue
                            if v > 1 and fy.get(good_prime, 0) != 1:
                                continue
                            if any((x % p == 0 or y % p == 0) for p in small if p not in factorization_primes(u * v)):
                                continue
                            a, b = d * x, d * y
                            key = (a, b)
                            require(key not in seen, f"duplicate restricted fibre: {case.label}, {key}")
                            seen[key] = (d, u, v)
                            fa, fb = dict(factorization(a)), dict(factorization(b))
                            rd = ru = rv = 1
                            for p in small_w:
                                ea, eb = fa.get(p, 0), fb.get(p, 0)
                                if ea and eb:
                                    require(ea == eb == 1 and p % period == 1, "invalid common valuation pattern")
                                    rd *= p
                                elif ea:
                                    require(ea == 1 and p % period == 1, "invalid first directive pattern")
                                    ru *= p
                                elif eb:
                                    require(eb == 1 and p % period == 1, "invalid second directive pattern")
                                    rv *= p
                            require((rd, ru, rv) == (d, u, v), "restricted tuple recovery failed")
                            denominator = kappa**len(factorization(u)) * 2**len(factorization(v))
                            actual_weight = F(first_weight(period, x) * second_weight(case, y), denominator)
                            residual_weight = first_weight(period, x // u) * second_weight(case, y // v)
                            require(actual_weight == residual_weight, "restricted weight removal failed")
                            require(all(p > z for p, _ in factorization(x // u)), "first residual is not rough")
                            require(all(p > z for p, _ in factorization(y // v)), "second residual is not rough")
                            wrong = [(p, e) for p, e in factorization(b) if p % case.q != 1]
                            require(all(p > z for p, _ in wrong), "excluded small prime survived")
                            omega = sum(e for _, e in wrong)
                            if omega:
                                require(z**omega < b, "multiplicity product bound failed")
                                with_excluded += 1
                            if target % (good_prime**2) == 0 and d == good_prime:
                                high_target += 1
                            count += 1
        require(count > 0 and with_excluded > 0 and high_target > 0, "vacuous restricted fibre case")
        COUNTS["restricted_good_fibres"] += count
        summaries.append({"q": case.q, "discriminant": case.discriminant,
                          "target_max": target_max, "w": w, "z": z,
                          "fibres": count, "with_excluded_primes": with_excluded,
                          "higher_good_target_common_d": high_target})
    return summaries


def factorization_primes(n: int) -> set[int]:
    return {p for p, _ in factorization(n)}


def check_negative_controls() -> None:
    even = CharacterCase(5, 5)
    require(even.chi(3) == -1 and second_base(even, 3) == 0, "even-character negative-control setup")
    wrong_pairing = 2 * sum(even.chi(d) for d in range(1, isqrt(3) + 1) if 3 % d == 0)
    reject_wrong("dropping chi(t)=1 in divisor pairing", second_base(even, 3), wrong_pairing)
    reject_wrong("using constant 6 for q5 fixed-place factor", F(5), F(6))
    # At kappa=4, shift valuation 2, retaining the old local correction
    # after target deletion leaves A_1=(2+4)/(2+1)=2 instead of zero.
    reject_wrong("retaining the old shifted correction after target deletion", F(0), F(2 + 4, 2 + 1))
    reject_wrong("omitting inert target prime3", F(2), F(1))
    reject_wrong("omitting split nongood target prime3 from Jq", F(2), F(1))
    reject_wrong("using old 3/5 unit mass at generic prime3", F(1, 3), F(3, 5))
    split_three = CharacterCase(8, -8)
    coefficient = second_kernel(split_three, 3, 3)
    require(coefficient == -F(1, 2) and 3 % 8 != 1 and split_three.chi(3) == 1,
            "strong second-kernel-support negative control failed")
    reject_wrong("claiming second kernels are 1 modulo L", "nonzero at 3 modulo 8", "only 1 modulo 8")
    reject_wrong("omitting inert-square kernel at prime3", second_weight(even, 9), second_base(even, 9))


def pin(path: Path) -> dict[str, str]:
    resolved = path.resolve()
    require(resolved.is_file(), f"missing retained input: {resolved}")
    try:
        relative = str(resolved.relative_to(ROOT)).replace("\\", "/")
    except ValueError as exc:
        raise CheckFailure(f"receipt input is outside ROOT: {resolved}") from exc
    return {"path": relative, "sha256": hashlib.sha256(resolved.read_bytes()).hexdigest()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=DEFAULT_SOURCE,
                        help="general source relative to ROOT (scratch override before promotion)")
    parser.add_argument("--output", type=Path, help="write the successful JSON receipt here")
    args = parser.parse_args()
    source = Path(args.source)
    if not source.is_absolute():
        source = ROOT / source
    try:
        source_pins = [pin(source)] + [pin(ROOT / item) for item in RETAINED_INPUTS]
        check_first_kernels()
        models = check_characters_and_hyperbola()
        prime_three = check_local_integrals()
        check_fixed_prime_stripping()
        fibres = check_restricted_good_fibres()
        check_negative_controls()
        receipt = {
            "status": "passed",
            "scope": "Finite exact character, local-density, kernel and fibre controls; not a verifier of analytic inputs, distribution estimates or unrestricted asymptotics.",
            "source_binding": "Hashes identify the inspected bytes; the arithmetic checks do not parse or certify the proof text.",
            "python": {"implementation": sys.implementation.name, "version": sys.version.split()[0]},
            "checker": pin(Path(__file__)),
            "sources": source_pins,
            "ranges": {"primes_through": 2500, "first_shift_values": [0, 1, 2, 10, 10000],
                       "first_degrees_below": 20, "hyperbola_t_through": 3000,
                       "second_kernel_n_through": 500, "prime3_target_valuations": list(range(6)),
                       "good_integral_target_valuations": [0, 1, 2]},
            "counts": dict(sorted(COUNTS.items())),
            "character_models": models,
            "generic_prime3": prime_three,
            "restricted_good_fibres": fibres,
            "negative_controls": NEGATIVE_CONTROLS,
        }
        rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
        if args.output is not None:
            protected = {source.resolve(), Path(__file__).resolve()}
            protected.update((ROOT / item).resolve() for item in RETAINED_INPUTS)
            require(args.output.resolve() not in protected, "output would overwrite a retained input")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(rendered, encoding="utf-8")
        print(rendered, end="")
    except (CheckFailure, OSError, ValueError) as exc:
        print(f"finite controls failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    if not __debug__:
        print("refusing optimized Python: run this checker without -O", file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(main())
