#!/usr/bin/env python3
"""Finite exact constants for the ordinary balanced Eisenstein bridge.

Standard library only; run from any working directory.
No new target scan, Hecke theorem, or complete exception list is checked.
Use --output explicitly to write a compact JSON receipt.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt, prod
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("Run without Python -O.")

ROOT = Path(__file__).resolve().parents[3]
NOTE = "research/arithmetic/EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md"
SOURCES = (
    NOTE,
    "research/arithmetic/EISENSTEIN_POWER_SHRINKING_WINDOWS.md",
    "unpublished/eisenstein-restricted-prime-sums/sections/elementary-estimates.tex",
    "unpublished/eisenstein-restricted-prime-sums/sections/mixed-local.tex",
    "research/arithmetic/eisenstein-restricted-prime-sums/research_notes/literature_balanced_polynomial_receipt.json",
    "research/arithmetic/eisenstein-restricted-prime-sums/receipts/elementary_constants_receipt.json",
    "research/arithmetic/eisenstein-restricted-prime-sums/receipts/python_finite_receipt.json",
)
COUNTS = {}


def need(ok, label):
    if not ok:
        raise RuntimeError(label)
    COUNTS[label] = COUNTS.get(label, 0) + 1


def primes_through(limit):
    flag = bytearray(b"\1") * (limit + 1)
    flag[:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if flag[p]:
            start = p * p
            flag[start:limit + 1:p] = b"\0" * ((limit - start) // p + 1)
    return [p for p in range(2, limit + 1) if flag[p]]


def root_ceiling(n, exponent):
    lo, hi = 0, n + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** exponent >= n:
            hi = mid
        else:
            lo = mid
    need((hi - 1) ** exponent < n <= hi ** exponent, "integer_root_bracket")
    return hi


def integrate(poly, x, extra=0):
    return sum(v * x ** (j + extra + 1) / (j + extra + 1)
               for j, v in enumerate(poly))


def trapezoid_certificate():
    h, a = F(1, 3), F(1, 4)
    delta = h - a
    coeffs = {}
    previous, current = [F(1)], [F(0), F(1)]
    for d in range(101):
        if d == 0:
            poly = previous
        elif d == 1:
            poly = current
        else:
            poly = [F(0)] * (d + 1)
            for j, v in enumerate(current):
                poly[j + 1] += F(2 * d - 1, d) * v
            for j, v in enumerate(previous):
                poly[j] -= F(d - 1, d) * v
            previous, current = current, poly
        if d % 2 == 0:
            coeffs[d] = (2 * d + 1) * (
                integrate(poly, a)
                + (h * (integrate(poly, h) - integrate(poly, a))
                   - (integrate(poly, h, 1) - integrate(poly, a, 1))) / delta
            )
    partial = sum(abs(v) for d, v in coeffs.items() if d)
    need(coeffs[0] == F(7, 24), "trapezoid_mean")
    need(partial < F(23, 4), "finite_legendre_l1_bound")
    variation = (2 * (1 - h * h) + 2 * (1 - a * a)
                 + 2 * (h * h - a * a)) / delta
    need(variation == 45, "derivative_measure_variation")
    need(variation * F(3, 2) / 10 == F(27, 4), "even_degree_tail_bound")
    need(partial + F(27, 4) < F(25, 2), "total_legendre_l1_bound")
    payload = json.dumps({str(d): str(c) for d, c in coeffs.items()},
                         sort_keys=True, separators=(",", ":")).encode()
    return {
        "mean": "7/24",
        "degrees": "even 0 through100",
        "coefficient_sha256": hashlib.sha256(payload).hexdigest(),
        "nonconstant_partial_l1": str(partial),
        "partial_l1_upper_bound": "23/4",
        "even_tail_upper_bound": "27/4",
        "total_nonconstant_l1_upper_bound": "25/2",
    }


def local_maxima(primes, corrected):
    value = F(1)
    rows = []
    for p in primes:
        end = 0
        while (end + 2) ** 5 > p * (end + 1) ** 5:
            end += 1
        choices = []
        for e in range(end + 1):
            ratio = (F((e + 1) * (p - 1) * p ** e, p ** (e + 1) - 1)
                     if corrected else F(e + 1))
            choices.append(ratio ** 5 / p ** e)
        e = max(range(len(choices)), key=choices.__getitem__)
        value *= choices[e]
        rows.append([p, end, e])
    return value, rows


def envelope_constants():
    tau5, tau_rows = local_maxima([2, 5, 11, 17, 23, 29], True)
    m5, m_rows = local_maxima([7, 13, 19, 31], False)
    need(tau5 < F(5, 2) ** 5, "corrected_inert_tau_constant")
    need(m5 < F(15, 8) ** 5, "split_divisor_constant")
    need(F(2) * F(5, 2) * F(4, 3) == F(20, 3),
         "candidate_20_over_3_normalization")
    return {"tau_product_fifth_power": str(tau5), "tau_maxima": tau_rows,
            "M_product_fifth_power": str(m5), "M_maxima": m_rows}


def prime_family_blocks():
    R = 10000
    inert = [p for p in primes_through(R) if p >= 5 and p % 3 == 2]
    tail = F(2, 15 * R) + F(8, 5 * R * R)
    result = []
    old_c0 = F(3824117519863, 17459442000000)
    old_B = F(69478795301277343, 20180620746000000)
    for j in range(4):
        L, U = 10 ** 8 * 2 ** j, 10 ** 8 * 2 ** (j + 1)
        a = F(44, 3 * isqrt(L))
        m = F(15, 8) * root_ceiling(U, 5)
        loss = sum(min(F(2, p * (p - 1)) + a, 4 * m / (p * p))
                   for p in inert) + 4 * m * tail
        margin = F(7, 24) - F(25, 2) * a - loss
        need(margin > F(137, 10000), "positive_uniform_prime_family_block")
        floor_micro = margin.numerator * 10 ** 6 // margin.denominator
        need(F(floor_micro, 10 ** 6) <= margin, "certified_margin_rounding")
        if j == 0:
            need(old_c0 - old_B * a - loss < 0,
                 "control_old_polynomial_does_not_certify_first_block")
        result.append({"L": L, "U": U, "a_upper": str(a), "M_upper": str(m),
                       "margin_lower_bound": str(F(floor_micro, 10 ** 6))})
    for p in [7, 13, 19, 31, 37]:
        for n, tau, S, bound in [
            (8 * p, 8, 12 * (p + 1), F(32, 3)),
            (44 * p, 12, 72 * (p + 1), F(44, 3)),
        ]:
            need(F(24 * tau * n, 12 * S) < bound,
                 "thin_family_spectral_normalization")
    need(22 ** 20 < 10 ** 27, "global_power_lower_bound")
    need(30000 ** 2 < 10 ** 9, "global_square_root_lower_bound")
    A, B = F(44, 3), F(25, 2)
    need(2 * A / 15 + 1 == F(133, 45), "global_three_twentieths_coefficient")
    need(17 * A / 15 + 12 + B * A == F(9538, 45),
         "global_half_power_coefficient")
    global_margin = F(7, 24) - F(1, 7) - F(133, 45 * 22) - F(9538, 45 * 30000)
    need(global_margin > 0, "positive_global_thin_family_tail")
    return {"prime_cutoff": R, "odd_inert_prime_count": len(inert),
            "wheel_tail_bound": str(tail), "blocks": result,
            "global_tail_onset": 10 ** 9, "global_tail_margin": str(global_margin)}


def inert_onsets():
    old_c0 = F(3824117519863, 17459442000000)
    old = old_c0 - F(1, 7) - F(17, 9 * 25) - F(386, 9 * 25 ** 6)
    conservative = F(7, 24) - F(1, 7) - F(17, 9 * 13) - F(926, 9 * 13 ** 6)
    refined = F(7, 24) - F(1, 7) - F(11, 80) - (F(111, 20) + F(250, 3)) / 8 ** 6
    need(old > 0, "fixed_polynomial_inert_onset_25_power20")
    need(conservative > 0, "trapezoid_inert_onset_13_power20")
    need(refined > 0, "refined_scalar_inert_onset_8_power20")
    return {"fixed_polynomial_onset": 25 ** 20, "fixed_polynomial_margin": str(old),
            "conservative_trapezoid_onset": 13 ** 20, "conservative_margin": str(conservative),
            "refined_scalar_onset": 8 ** 20, "refined_scalar_margin": str(refined)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.resolve() in (
            {Path(__file__).resolve()} | {(ROOT / r).resolve() for r in SOURCES}):
        raise ValueError("Output would overwrite a checker or proof dependency.")
    hashes = {r: hashlib.sha256((ROOT / r).read_bytes()).hexdigest() for r in SOURCES}
    receipt = {
        "result": "PASS",
        "scope": ("Exact finite coefficient, normalization, envelope and inequality checks. "
                  "No new target scan, no Hecke theorem verification, no inert-core "
                  "interval coverage and no complete exception theorem."),
        "provenance": {
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_sha256": hashes, "python_version": sys.version.split()[0],
        },
        "trapezoid": trapezoid_certificate(),
        "elementary_envelopes": envelope_constants(),
        "thin_families": prime_family_blocks(),
        "inert_onsets": inert_onsets(),
    }
    receipt["checks"] = COUNTS
    text = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()

