#!/usr/bin/env python3
"""Finite normalization checks for the pair-primitive quaternary deduction.

These checks do not prove its analytic input, limiting estimates, or
unrestricted local-global theorem. Standard library only; exact integers.
"""
from __future__ import annotations


if not __debug__:
    raise RuntimeError(
        "This certificate checker requires active Python assertions; "
        "run without -O, -OO, or PYTHONOPTIMIZE."
    )

import argparse
import hashlib
import json
from math import gcd, isqrt
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FORMS = {
    "-4": (1, 0, 1),
    "-8": (1, 0, 2),
    "-7": (1, 1, 2),
    "-11": (1, 1, 3),
    "-20": (1, 0, 5),
    "-27": (1, 1, 7),
    "nonprincipal_-20": (2, 2, 3),
}


def discriminant(q):
    a, b, c = q
    return b * b - 4 * a * c


def qvalue(q, x, y):
    a, b, c = q
    return a * x * x + b * x * y + c * y * y


def binary_residues(q, p, h):
    modulus = p**h
    counts = [0] * modulus
    for x in range(modulus):
        for y in range(modulus):
            if x % p or y % p:
                counts[qvalue(q, x, y) % modulus] += 1
    return counts


def convolution_mod(a, b):
    modulus = len(a)
    result = [0] * modulus
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    result[(i + j) % modulus] += ai * bj
    return result


def legendre(a, p):
    residue = pow(a % p, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def expected_good_numerator(p, e1, e2, target_zero):
    base = p**3 - (2 + e1 * e2) * p + e1 + e2
    return base + ((p * e1 - 1) * (p * e2 - 1) if target_zero else 0)


def mobius(n):
    answer = 1
    trial = 2
    while trial * trial <= n:
        if n % trial == 0:
            n //= trial
            answer = -answer
            if n % trial == 0:
                return 0
        trial += 1
    return -answer if n > 1 else answer


def binary_integer_counts(q, bound):
    a, b, c = q
    determinant = 4 * a * c - b * b
    # lambda_min(A) >= det(A)/trace(A), and Q = x^t A x/2.
    radius = isqrt((4 * (a + c) * bound) // determinant) + 1
    raw = [0] * (bound + 1)
    primitive = [0] * (bound + 1)
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            n = qvalue(q, x, y)
            if 0 <= n <= bound:
                raw[n] += 1
                if gcd(x, y) == 1:
                    primitive[n] += 1
    return raw, primitive



def eta32_coefficients(bound):
    product = [0] * bound
    product[0] = 1
    for step in (4, 8):
        for shift in range(step, bound, step):
            for _ in range(2):
                for n in range(bound - 1, shift - 1, -1):
                    product[n] -= product[n - shift]
    return [0] + product


def divisor_sums(bound):
    result = [0] * (bound + 1)
    for d in range(1, bound + 1):
        for n in range(d, bound + 1, d):
            result[n] += d
    return result


def dyadic_sturm_checks(bound):
    limit = max(8, bound)
    q = FORMS["-8"]
    raw, _ = binary_integer_counts(q, limit)
    dyadic = [0] * (limit + 1)
    radius = isqrt(limit) + 1
    for x in range(-radius, radius + 1):
        for y in range(-radius, radius + 1):
            value = qvalue(q, x, y)
            if value <= limit and (x % 2 or y % 2):
                dyadic[value] += 1
    eta = eta32_coefficients(limit)
    sigma = divisor_sums(limit)

    def sd(n, denominator=1):
        return sigma[n // denominator] if n and n % denominator == 0 else 0

    raw_shell = [sum(raw[a] * raw[n - a] for a in range(n + 1))
                 for n in range(limit + 1)]
    scaled_shell = [sum(raw[a] * raw[(n - a) // 4]
                        for a in range(n + 1) if (n - a) % 4 == 0)
                    for n in range(limit + 1)]
    dyadic_shell = [sum(dyadic[a] * dyadic[n - a] for a in range(n + 1))
                    for n in range(limit + 1)]
    r_formula = [1] + [4 * sd(n) - 4 * sd(n, 2) + 8 * sd(n, 4) - 32 * sd(n, 8)
                       for n in range(1, limit + 1)]
    s_formula = [1] + [sd(n) - sd(n, 2) + 8 * sd(n, 16) - 32 * sd(n, 32) + eta[n]
                       for n in range(1, limit + 1)]
    assert raw_shell[:3] == r_formula[:3] == [1, 4, 8]
    assert scaled_shell[:9] == s_formula[:9] == [1, 2, 2, 4, 4, 4, 8, 8, 8]
    assert raw_shell == r_formula
    assert scaled_shell == s_formula
    for n in range(1, limit + 1):
        e, m = 0, n
        while m % 2 == 0:
            e += 1
            m //= 2
        expected = (2 * sigma[m] - 2 * eta[m] if e == 0 else
                    4 * sigma[m] if e == 1 else
                    20 * sigma[m] if e == 2 else
                    16 * sigma[m] if e == 3 else 0)
        inclusion = raw_shell[n] - 2 * scaled_shell[n]
        if n % 4 == 0:
            inclusion += raw_shell[n // 4]
        assert dyadic_shell[n] == expected == inclusion, n
        if n % 2 == 0:
            assert eta[n] == 0
        assert eta[n] * eta[n] <= sum(1 for d in range(1, n + 1) if n % d == 0)**2 * n
    return {
        "completeness": "modularity plus sourced Sturm bounds proves unrestricted theta identities from the listed coefficients",
        "r_level": 8,
        "r_index": 12,
        "r_sturm_bound": 2,
        "r_coefficients_through_bound": raw_shell[:3],
        "s_level": 32,
        "s_index": 48,
        "s_sturm_bound": 8,
        "s_coefficients_through_bound": scaled_shell[:9],
        "dyadic_direct_comparison_through": limit,
        "eta_coefficients_through_17": eta[:min(18, len(eta))],
    }



def minus8_explicit_constants():
    from itertools import product

    def prime(n):
        return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))

    factor_certificates = {}
    products = {}
    for k, split, constant in ((6, True, 6), (6, False, 7), (4, True, 2)):
        assert Fraction(k + 2, k + 1)**k < 3
        factors = []
        value = Fraction(1)
        for p in range(3, 2**k, 2):
            if not prime(p) or (p % 8 in (1, 3)) != split:
                continue
            local, exponent = max((Fraction((a + 1)**k, p**a), a)
                                  for a in range(k + 1))
            factors.append({"prime": p, "max_exponent": exponent, "factor": str(local)})
            value *= local
        assert value < constant**k
        key = f"k{k}_{'split' if split else 'inert'}"
        products[key] = value
        factor_certificates[key] = {
            "power_product": str(value), "root_bound": constant,
            "prime_cutoff_exclusive": 2**k, "exponent_bound": k,
            "local_factors": factors,
        }
    assert products["k4_split"] == Fraction(4096, 297)
    assert products["k6_split"] * products["k6_inert"] < 42**6

    small_inert = [p for p in range(5, 65) if prime(p) and p % 8 in (5, 7)]
    wheel = sum((Fraction(2, p * (p - 1)) for p in small_inert), Fraction(0))
    wheel += Fraction(5, 2) * sum((Fraction(1, a*a) + Fraction(1, 8*a)
                                 for a in (69, 71)), Fraction(0))
    assert wheel == Fraction(49287299505759574293, 269836951018231737140)
    assert wheel < Fraction(3, 16)
    inequalities = {
        "all_odd_m_at_2^75": Fraction(3, 16) + Fraction(135, 448) + Fraction(159, 2**25),
        "inert_m_at_2^60": Fraction(3, 16) + Fraction(19, 64) + Fraction(69, 2**21),
        "odd_prime_at_2^22": Fraction(3, 16) + Fraction(3, 10) + Fraction(7, 2048),
    }
    assert 7**4 < 2 * 6**4
    assert 20**8 < 2**22 * 3**8
    assert all(value < Fraction(1, 2) for value in inequalities.values())
    explicit_survival_margin = Fraction(1,2)-inequalities["all_odd_m_at_2^75"]
    assert explicit_survival_margin == Fraction(2620327,234881024)
    assert explicit_survival_margin > Fraction(1,100)

    def multiply(v, w):
        a, b, c, d = v
        e, f, g, h = w
        return (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)

    def hurwitz_residue(v):
        a, b, c, d = v
        coordinates = (a-d, b-d, c-d, 2*d)
        assert all(x.denominator == 1 for x in coordinates)
        return tuple(int(x) % 2 for x in coordinates)

    def norm(v):
        return sum((x*x for x in v), Fraction(0))

    def in_order(v):
        return all(x.denominator == 1 for x in v) and int(v[1]-v[2]) % 2 == 0

    local_representatives = []
    for a, b, c, d in product((0, 1), repeat=4):
        local_representatives.append((Fraction(a)+Fraction(d, 2), Fraction(b)+Fraction(d, 2),
                                      Fraction(c)+Fraction(d, 2), Fraction(d, 2)))
    local_units = [v for v in local_representatives if norm(v).denominator == 1
                   and int(norm(v)) % 2]
    local_order_units = [v for v in local_units if in_order(v)]
    assert len(local_units) == 12 and len(local_order_units) == 2
    hurwitz_units = []
    for i in range(4):
        for sign in (-1, 1):
            v = [Fraction(0)] * 4
            v[i] = Fraction(sign)
            hurwitz_units.append(tuple(v))
    hurwitz_units += [tuple(Fraction(s, 2) for s in signs)
                     for signs in product((-1, 1), repeat=4)]
    assert len(hurwitz_units) == 24 and all(norm(v) == 1 for v in hurwitz_units)
    assert sum(in_order(v) for v in hurwitz_units) == 4
    cosets = {frozenset(hurwitz_residue(multiply(u, v)) for v in local_order_units)
              for u in hurwitz_units}
    assert len(cosets) == 6
    assert set().union(*cosets) == {hurwitz_residue(v) for v in local_units}
    return {
        "scope": "exact finite factor and local-unit certificates used by the explicit -8 Hecke sieve proof",
        "factor_certificates": factor_certificates,
        "wheel_primes_through_64": small_inert,
        "wheel_upper_bound": str(wheel),
        "cutoff_inequalities": {key: str(value) for key, value in inequalities.items()},
        "global_N_onset": 2**78,
        "global_surviving_shell_margin": str(explicit_survival_margin),
        "distinct_value_lower_bound": "N^(2/3)/4800 at the global onset",
        "reduced_inert_odd_m_onset": 2**60,
        "odd_prime_onset": 2**22,
        "hurwitz_units": 24,
        "order_units": 4,
        "local_hurwitz_unit_residues_mod_2H": len(local_units),
        "local_order_unit_residues_mod_2H": len(local_order_units),
        "unit_cosets": len(cosets),
    }



def minus8_prime_certificate():
    bound = 2**22
    short_bound = 5000
    primes = bytearray(b"\x01") * (bound + 1)
    primes[0:2] = b"\x00\x00"
    for p in range(2, isqrt(bound) + 1):
        if primes[p]:
            start = p*p
            primes[start::p] = b"\x00" * ((bound - start) // p + 1)

    split_core = bytearray(b"\x01") * (bound + 1)
    split_core[0] = 0
    for p in range(2, bound + 1):
        if primes[p] and (p == 2 or p % 8 in (5, 7)):
            split_core[p::p] = b"\x00" * (bound // p)
    primitive_norm = bytearray(split_core)
    for n in range(1, bound // 2 + 1):
        if split_core[n]:
            primitive_norm[2*n] = 1
    assert primitive_norm[0] == 0
    assert [bool(primitive_norm[n]) for n in (1, 2, 3, 5, 6, 9, 16)] == [
        True, True, True, False, True, True, False]
    _, enumerated = binary_integer_counts(FORMS["-8"], 10000)
    assert all(bool(primitive_norm[n]) == bool(enumerated[n]) for n in range(10001))

    packed = bytearray((bound + 8) // 8)
    for n in range(1, bound + 1):
        if primitive_norm[n]:
            packed[n >> 3] |= 1 << (n & 7)
    norm_bits = int.from_bytes(packed, "little")
    short_summands = [n for n in range(1, short_bound + 1) if primitive_norm[n]]
    covered = 0
    for a in short_summands:
        covered |= norm_bits << a
    covered &= (1 << (bound + 1)) - 1
    covered_bytes = covered.to_bytes((bound + 8) // 8, "little")

    count = 0
    fallback_targets = []
    fallback_witnesses = []
    exceptions = []
    for p in range(2, bound + 1):
        if not primes[p]:
            continue
        count += 1
        if covered_bytes[p >> 3] & (1 << (p & 7)):
            continue
        fallback_targets.append(p)
        witness = next((a for a in range(1, p // 2 + 1)
                        if primitive_norm[a] and primitive_norm[p-a]), None)
        if witness is None:
            exceptions.append(p)
        else:
            fallback_witnesses.append([p, witness, p-witness])
    assert not exceptions, exceptions
    return {
        "scope": "exhaustive finite proof certificate for every prime target up to 2^22, combined with the analytic prime threshold for the all-prime theorem",
        "bound_inclusive": bound,
        "prime_targets": count,
        "prime_exceptions": exceptions,
        "norm_membership": "E3: 2^epsilon times a product of primes 1 or 3 modulo 8, epsilon=0 or 1",
        "independent_norm_regression_through": 10000,
        "short_summand_bound": short_bound,
        "short_summand_count": len(short_summands),
        "fallback_targets": fallback_targets,
        "fallback_witnesses": fallback_witnesses,
        "norm_flags_sha256": hashlib.sha256(primitive_norm).hexdigest(),
        "sumset_coverage_sha256": hashlib.sha256(covered_bytes).hexdigest(),
    }


def run(bound):
    prime_formula_cases = 0
    lifted_cases = 0
    forms = list(FORMS.items())
    for name1, q1 in forms:
        for name2, q2 in forms:
            for p in (3, 5, 7, 11, 13, 17, 19):
                if discriminant(q1) * discriminant(q2) % p == 0:
                    continue
                e1 = legendre(discriminant(q1), p)
                e2 = legendre(discriminant(q2), p)
                roots = convolution_mod(binary_residues(q1, p, 1),
                                        binary_residues(q2, p, 1))
                for target, count in enumerate(roots):
                    expected = expected_good_numerator(p, e1, e2, target == 0)
                    assert count == expected, (name1, name2, p, target, count, expected)
                    assert count > 0
                    prime_formula_cases += 1
                if p <= 5 and name1 in ("-8", "-20") and name2 in ("-11", "-27"):
                    roots2 = convolution_mod(binary_residues(q1, p, 2),
                                             binary_residues(q2, p, 2))
                    for target, count in enumerate(roots2):
                        assert count == roots[target % p] * p**3
                        lifted_cases += 1

    integer_counts = {name: binary_integer_counts(q, bound) for name, q in forms}
    mobius_cases = 0
    for name, (raw, primitive) in integer_counts.items():
        for n in range(1, bound + 1):
            inverted = sum(mobius(d) * raw[n // (d * d)]
                           for d in range(1, isqrt(n) + 1) if n % (d * d) == 0)
            assert inverted == primitive[n], (name, n)
            mobius_cases += 1

    # An independent comparison of coordinate enumeration and two divisibility sums.
    double_mobius_cases = 0
    for name1, name2 in (("-8", "-8"), ("-11", "-11"), ("-20", "-7"),
                         ("nonprincipal_-20", "-27")):
        raw1, prim1 = integer_counts[name1]
        raw2, prim2 = integer_counts[name2]
        for n in range(2, bound + 1):
            direct = sum(prim1[a] * prim2[n - a] for a in range(1, n))
            expanded = 0
            for d in range(1, isqrt(n - 1) + 1):
                md = mobius(d)
                if not md:
                    continue
                for e in range(1, isqrt(n - 1) + 1):
                    me = mobius(e)
                    if not me:
                        continue
                    for a in range(d * d, n, d * d):
                        if (n - a) % (e * e) == 0:
                            expanded += md * me * raw1[a // (d * d)] * raw2[(n - a) // (e * e)]
            assert expanded == direct, (name1, name2, n)
            double_mobius_cases += 1

    # Minimal obstruction moduli for the two principal self-sum test cases.
    local_samples = {}
    for name, p, h, forbidden_modulus in (("-8", 2, 7, 16), ("-11", 11, 3, 121)):
        q = FORMS[name]
        roots = convolution_mod(binary_residues(q, p, h), binary_residues(q, p, h))
        support = [n for n, count in enumerate(roots) if count]
        assert all(bool(roots[n]) == (n % forbidden_modulus != 0)
                   for n in range(p**h)), (name, p, h)
        local_samples[name] = {
            "prime": p, "exponent": h, "target_residues": p**h,
            "soluble_residues": len(support),
            "observed_support": f"not divisible by {forbidden_modulus}",
        }

    # Unrestricted four-square normalization at n=1: det(2 I_4)=16.
    # The odd-prime density product is 8/pi^2; the real factor is
    # 4*pi^2/sqrt(16). Thus the exact product must equal r_4(1)=8.
    raw_gauss, _ = integer_counts["-4"]
    r4_one = raw_gauss[0] * raw_gauss[1] * 2
    counts_mod8 = [0] * 8
    for x in range(8):
        for y in range(8):
            counts_mod8[(x*x+y*y) % 8] += 1
    beta2_one = Fraction(convolution_mod(counts_mod8, counts_mod8)[1], 8**3)
    assert beta2_one == 1
    assert Fraction(4, isqrt(16)) * 8 * beta2_one == r4_one == 8

    source_pdf = ROOT / "sources" / "rouse-1802.03437v1.pdf"
    return {
        "scope": "finite exact checks; listed Sturm certificates prove the two theta identities under sourced modularity; no computational proof of the general analytic theorem or independent mathematical review",
        "python": "standard library; exact integers",
        "integer_bound": bound,
        "dyadic_sturm_certificate": dyadic_sturm_checks(bound),
        "minus8_explicit_constants": minus8_explicit_constants(),
        "minus8_all_prime_certificate": minus8_prime_certificate(),
        "four_square_normalization": {"gram_determinant": 16, "r4_one": r4_one,
                                      "beta2_one": str(beta2_one)},
        "good_prime_formula_residues": prime_formula_cases,
        "nonsingular_lift_residues": lifted_cases,
        "binary_mobius_cases": mobius_cases,
        "double_mobius_targets": double_mobius_cases,
        "local_sample_checks": local_samples,
        "rouse_pdf_sha256": hashlib.sha256(source_pdf.read_bytes()).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert args.bound >= 2
    result = run(args.bound)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
