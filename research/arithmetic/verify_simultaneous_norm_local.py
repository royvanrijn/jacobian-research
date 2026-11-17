#!/usr/bin/env python3
"""Exact checks for simultaneous Gaussian--Eisenstein local densities.

The proof is in SIMULTANEOUS_NORM_LOCAL_DENSITIES.md.  This script enumerates
binary norm value histograms, not eight-dimensional solution tuples.  Its finite
checks test the formulas; they do not establish a global integer-point count.
"""
import argparse
import hashlib
import itertools
import json
import math
import platform
from fractions import Fraction as F
from pathlib import Path


def require_equal(actual, expected, label):
    if actual != expected:
        raise RuntimeError(f"{label}: {actual!r} != {expected!r}")


def valuation(n, p, cutoff):
    if n == 0:
        return cutoff
    result = 0
    while result < cutoff and n % p == 0:
        n //= p
        result += 1
    return result


def signs(p):
    return (1 if p % 4 == 1 else -1, 1 if p % 3 == 1 else -1)


def binary_count(p, k, v, sign):
    """Number of solutions to a binary norm at valuation min(v,k)."""
    if v == k:
        if sign == 1:
            return (p + k * (p - 1)) * p ** (k - 1)
        return p ** (2 * (k // 2))
    if sign == 1:
        return (v + 1) * (p - 1) * p ** (k - 1)
    return (p + 1) * p ** (k - 1) if v % 2 == 0 else 0


def histograms(modulus):
    gaussian, eisenstein = [0] * modulus, [0] * modulus
    for x in range(modulus):
        xx = x * x
        for y in range(modulus):
            yy = y * y
            gaussian[(xx + yy) % modulus] += 1
            eisenstein[(xx - x * y + yy) % modulus] += 1
    return gaussian, eisenstein


def finite_from_valuations(p, k, n):
    """Independent valuation-stratum formula, including the zero residue."""
    eps, eta = signs(p)
    values = [binary_count(p, k, v, eps) * binary_count(p, k, v, eta)
              for v in range(k + 1)]
    r = valuation(n, p, k)
    strata = [(p - 1) * p ** (k - e - 1) for e in range(k)]
    if r == k:
        total = sum(strata[e] * values[e] ** 2 for e in range(k))
        total += values[k] ** 2
    else:
        total = sum(strata[e] * values[e] ** 2 for e in range(r))
        total += (p - 2) * p ** (k - r - 1) * values[r] ** 2
        total += 2 * sum(strata[e] * values[e] * values[r]
                         for e in range(r + 1, k))
        total += 2 * values[k] * values[r]
    good = (p - 2 if r == 0 else p - 1) * p ** (k - 1) * values[0] ** 2
    return total, total - good


def even_tail(p, d, first):
    """Sum_{j>=first} (2*j+1)^d / p^(2*j), for d=0 or1."""
    q = F(1, p * p)
    if d == 0:
        return q ** first / (1 - q)
    return q ** first * (F(2 * first + 1, 1) / (1 - q)
                         + 2 * q / (1 - q) ** 2)


def limiting_integral(p, d, r):
    """Integral after removing the constant block factor h_p^2."""
    unit_measure = F(p - 1, p)
    result = sum((unit_measure * F((e + 1) ** (2 * d), p ** e)
                  for e in range(0, r, 2)), F(0))
    if r % 2 == 0:
        result += F((p - 2) * (r + 1) ** (2 * d), p ** (r + 1))
        result += (2 * unit_measure * (r + 1) ** d
                   * even_tail(p, d, r // 2 + 1))
    return result


def limiting_bad(p, r):
    residue = p % 12
    if residue not in (5, 7, 11):
        raise ValueError("the prime must be forbidden")
    d = 0 if residue == 11 else 1
    integral = limiting_integral(p, d, r)
    good = F(p - 2 if r == 0 else p - 1, p)
    return (integral - good) / integral


def rank_mod_p(rows, p):
    rows = [[value % p for value in row] for row in rows]
    rank = 0
    for column in range(3):
        pivot = next((j for j in range(rank, len(rows))
                      if rows[j][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, p)
        rows[rank] = [(x * inverse) % p for x in rows[rank]]
        for j in range(rank + 1, len(rows)):
            multiplier = rows[j][column]
            rows[j] = [(x - multiplier * y) % p
                       for x, y in zip(rows[j], rows[rank])]
        rank += 1
    return rank


def singular_from_histograms(p, n, gaussian, eisenstein):
    columns = [(1, 0, 1), (-1, 0, 0), (0, 1, 1), (0, -1, 0)]
    singular = 0
    for a in range(p):
        b = (n - a) % p
        options = []
        for t, counts in ((a, gaussian), (a, eisenstein),
                          (b, gaussian), (b, eisenstein)):
            options.append([(False, 1), (True, counts[t] - 1)]
                           if t == 0 else [(True, counts[t])])
        for choice in itertools.product(*options):
            support = [columns[j] for j in range(4) if choice[j][0]]
            if rank_mod_p(support, p) < 3:
                singular += math.prod(item[1] for item in choice)
    return singular


def primitive_controls(limit=2000):
    gaussian, eisenstein = set(), set()
    bound = math.isqrt(4 * limit // 3)
    for x in range(-bound, bound + 1):
        for y in range(-bound, bound + 1):
            if math.gcd(x, y) != 1:
                continue
            g, e = x * x + y * y, x * x - x * y + y * y
            if 0 < g <= limit:
                gaussian.add(g)
            if 0 < e <= limit:
                eisenstein.add(e)
    intersection = gaussian & eisenstein
    expected = set()
    for n in range(1, limit + 1):
        remainder, allowed, p = n, True, 2
        while p * p <= remainder:
            if remainder % p == 0:
                allowed &= p % 12 == 1
                while remainder % p == 0:
                    remainder //= p
            p += 1
        if remainder > 1:
            allowed &= remainder % 12 == 1
        if allowed:
            expected.add(n)
    require_equal(intersection, expected, "primitive norm intersection")
    for n in (1, 13, 37, 169, 481):
        require_equal(n in intersection, True, f"positive value {n}")
    for n in (2, 3, 25, 49, 121):
        require_equal(n in intersection, False, f"negative value {n}")
    return {"limit": limit, "intersection_count": len(intersection),
            "positive_values": [1, 13, 37, 169, 481],
            "negative_values": [2, 3, 25, 49, 121]}


def ramified_controls():
    cases = []
    for p, exponents, congruence, expected_density in (
            (2, range(2, 8), 4, F(9, 4)), (3, range(1, 5), 3, F(64, 27))):
        for k in exponents:
            modulus = p ** k
            gaussian, eisenstein = histograms(modulus)
            for a in range(modulus):
                if a % p == 0:
                    continue
                if p == 2:
                    expected_g = 2 ** (k + 1) if a % 4 == 1 else 0
                    expected_e = 3 * 2 ** (k - 1)
                else:
                    expected_g = 4 * 3 ** (k - 1)
                    expected_e = 2 * 3 ** k if a % 3 == 1 else 0
                require_equal(gaussian[a], expected_g, "ramified Gaussian unit count")
                require_equal(eisenstein[a], expected_e, "ramified Eisenstein unit count")
            block = [g * e for g, e in zip(gaussian, eisenstein)]
            targets = list(range(2, modulus, congruence))
            for n in targets:
                total = sum(block[a] * block[(n - a) % modulus]
                            for a in range(modulus))
                bad = sum(block[a] * block[(n - a) % modulus]
                          for a in range(modulus) if a % p == 0 or (n - a) % p == 0)
                require_equal(F(total, p ** (5*k)), expected_density,
                              "ramified total density")
                require_equal(bad, 0, "automatic ramified-prime exclusion")
            cases.append({"p": p, "k": k, "target_residues_checked": len(targets),
                          "normalized_total_density": str(expected_density)})
    return cases


def run():
    plan = [(5, 4), (7, 3), (11, 3), (13, 2), (17, 2), (19, 2), (23, 2)]
    finite, local_limits = [], []
    for p, maximum_k in plan:
        eps, eta = signs(p)
        h_field = (p - eps) * (p - eta)
        z_field = (p + (p - 1) * eps) * (p + (p - 1) * eta)
        for k in range(1, maximum_k + 1):
            modulus = p ** k
            gaussian, eisenstein = histograms(modulus)
            for t in range(modulus):
                v = valuation(t, p, k)
                require_equal(gaussian[t], binary_count(p, k, v, eps),
                              f"Gaussian histogram p={p},k={k},t={t}")
                require_equal(eisenstein[t], binary_count(p, k, v, eta),
                              f"Eisenstein histogram p={p},k={k},t={t}")
            block = [g * e for g, e in zip(gaussian, eisenstein)]
            cases = []
            for n in sorted({j % modulus for j in (0, 1, 2, p - 1, p, p * p)}):
                total = sum(block[a] * block[(n - a) % modulus]
                            for a in range(modulus))
                bad = sum(block[a] * block[(n - a) % modulus]
                          for a in range(modulus) if a % p == 0 or (n - a) % p == 0)
                require_equal((total, bad), finite_from_valuations(p, k, n),
                              f"valuation strata p={p},k={k},N={n}")
                cases.append({"N_mod_modulus": n, "total": total, "bad": bad})
            if k == 1:
                require_equal(finite_from_valuations(p, 1, 1),
                              (h_field * ((p - 2) * h_field + 2 * z_field),
                               2 * h_field * z_field), "finite-field unit count")
                require_equal(singular_from_histograms(p, 1, gaussian, eisenstein),
                              2 * h_field, "unit-fiber singular count")
                l_g, l_e = (p - 1) * (1 + eps), (p - 1) * (1 + eta)
                singular_zero = 1 + 2 * l_g + 2 * l_e + l_g*l_g + 4*l_g*l_e + l_e*l_e
                require_equal(singular_from_histograms(p, 0, gaussian, eisenstein),
                              singular_zero, "zero-fiber singular count")
            digest = hashlib.sha256(json.dumps([gaussian, eisenstein],
                                               separators=(",", ":")).encode()).hexdigest()
            finite.append({"p": p, "k": k, "modulus": modulus,
                           "histograms_sha256": digest, "cases": cases})
        if p % 12 in (5, 7, 11):
            q = F(1, p * p)
            d = 0 if p % 12 == 11 else 1
            cauchy_bound = q if d == 0 else (9*q - 2*q*q + q**3)/(1 + 6*q + q*q)
            for r in range(13):
                beta = limiting_bad(p, r)
                if not 0 <= beta < 9*q:
                    raise RuntimeError("uniform local upper bound failed")
                if r and beta > cauchy_bound:
                    raise RuntimeError("Cauchy--Schwarz upper bound failed")
                if r == 0:
                    expected = (F(2, p * (p - 1)) if d == 0 else
                                F(2*(3*p*p - 1), p*(p**3 - p*p + 3*p + 1)))
                    require_equal(beta, expected, "unit p-adic fraction")
                if r == 1:
                    require_equal(beta, F(0), "valuation-one bad fraction")
            for r in (0, 1, 2, 3, 4):
                k = 20
                total, bad = finite_from_valuations(p, k, p ** r)
                beta = limiting_bad(p, r)
                error = abs(F(bad, total) - beta)
                if error >= F(1, p ** 12):
                    raise RuntimeError("prime-power fraction did not approach p-adic limit")
                local_limits.append({"p": p, "v_p_N": r, "bad_fraction": str(beta),
                                     "modulus_exponent": k, "finite_fraction_error": str(error)})
    # This rational integral-test bound already proves the all-prime sum <1.
    rational_upper = 9 * (F(1, 25) + F(1, 49) + F(1, 121) + F(1, 169)
                          + F(1, 66) + F(1, 78))
    if not rational_upper < 1:
        raise RuntimeError("strict all-prime margin failed")
    rejected = []
    for label, actual, false_expected in (
            ("raw_mod_p_is_not_p_adic", limiting_bad(11, 0), F(1, 649)),
            ("odd_inert_valuation", binary_count(5, 2, 1, -1), 30),
            ("corrupt_histogram", binary_count(7, 2, 0, 1), 43)):
        try:
            require_equal(actual, false_expected, label)
        except RuntimeError:
            rejected.append(label)
        else:
            raise RuntimeError(f"negative control was accepted: {label}")
    return {"schema": "simultaneous-norm-local-checks/v1", "status": "pass",
            "scope": "Exact finite histogram checks and formula checks; no global counting theorem.",
            "python": platform.python_version(),
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "finite_histogram_cases": finite, "p_adic_formula_checks": local_limits,
            "ramified_histogram_cases": ramified_controls(),
            "rational_all_prime_upper_bound": str(rational_upper),
            "primitive_intersection": primitive_controls(), "negative_controls_rejected": rejected}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("simultaneous_norm_local_checks.json"))
    args = parser.parse_args()
    receipt = run()
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"],
                      "histogram_prime_powers": len(receipt["finite_histogram_cases"]),
                      "ramified_prime_powers": len(receipt["ramified_histogram_cases"]),
                      "p_adic_formula_cases": len(receipt["p_adic_formula_checks"]),
                      "primitive_values_checked": receipt["primitive_intersection"]["limit"],
                      "negative_controls_rejected": len(receipt["negative_controls_rejected"])}))
