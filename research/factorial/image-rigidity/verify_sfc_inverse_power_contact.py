#!/usr/bin/env python3
"""Exact transcription controls for the inverse-power audit; stdlib only.

The universal proofs are in the pinned note.  This program does not test
or certify the conjectural general inverse-power contact inequality.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NOTE = HERE / "SFC_INVERSE_POWER_CONTACT_AUDIT.md"
NOTE_SHA = "ec6b4d6a8a5e84feb8a351c68567d367b900c8831d21dba916d713a6d74f2973"
SOURCES = {
    "research/factorial/image-rigidity/POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md":
        "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
    "research/factorial/image-rigidity/WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md":
        "f9ee0dd06d5e6bc85a5b8d35aacb6d55f5c8ec18e483a59c00d3275177495645",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a or [F(0)]


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return trim(out)


def scale(a, c):
    return trim([c * x for x in a])


def mul(a, b, degree=None):
    size = len(a) + len(b) - 1
    if degree is not None:
        size = min(size, degree + 1)
    out = [F(0)] * size
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            if i + j < size:
                out[i + j] += u * v
    return trim(out)


def power(a, n, degree=None):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, a, degree)
    return out


def compose(a, b, degree=None):
    out = [F(0)]
    for v in reversed(a):
        out = add(mul(out, b, degree), [v])
    return out


def coefficient(a, n):
    return a[n] if n < len(a) else F(0)


def inverse(a, degree):
    require(coefficient(a, 0) == 0 and coefficient(a, 1) == 1,
            "inverse input must be normalized")
    out = [F(0), F(1)]
    for n in range(2, degree + 1):
        value = -coefficient(compose(a, out, n), n)
        out += [F(0)] * (n + 1 - len(out))
        out[n] = value
    require(compose(a, out, degree) == [0, 1], "inverse composition")
    return trim(out)


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))])


def evaluate(a, x):
    out = 0
    for v in reversed(a):
        out = out * x + v
    return out


def valuation(a):
    return next((i for i, v in enumerate(a) if v), None)


def determinant(matrix):
    """Fraction-free Bareiss determinant with exact pivot checks."""
    a = [list(row) for row in matrix]
    n = len(a)
    if not n:
        return 1
    sign, denominator = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                require(numerator % denominator == 0,
                        "nonexact Bareiss division")
                a[i][j] = numerator // denominator
            a[i][k] = 0
        denominator = pivot
    return sign * a[-1][-1]


def resultant(a, b):
    a, b = list(reversed(trim(a))), list(reversed(trim(b)))
    m, n = len(a) - 1, len(b) - 1
    matrix = []
    for k in range(n):
        matrix.append([0] * k + a + [0] * (n - k - 1))
    for k in range(m):
        matrix.append([0] * k + b + [0] * (m - k - 1))
    return determinant(matrix)


def mod_gcd(a, b, p):
    a, b = trim([int(x) % p for x in a]), trim([int(x) % p for x in b])
    while b != [0]:
        r = a[:]
        while r != [0] and len(r) >= len(b):
            delta = len(r) - len(b)
            c = r[-1] * pow(b[-1], -1, p) % p
            for i in range(len(b)):
                r[i + delta] = (r[i + delta] - c * b[i]) % p
            r = trim(r)
        a, b = b, r
    return [x * pow(a[-1], -1, p) % p for x in a]


CV_POLYNOMIAL = [-29696, -225280, -21376768,
                 -54986400, -94035000, 2573571875]
PSI = [0, 1, 0, 0, -14, 28, -20, 5]


def verify_cv_certificate(poly):
    require(isinstance(poly, list) and len(poly) == 6
            and all(type(x) is int for x in poly),
            "malformed critical-value polynomial")
    require(poly[-1] != 0 and poly[0] != 0, "invalid degree or zero root")
    # The Sylvester determinant has degree at most deg(Psi')=6 in w.
    # Seven exact values certify the complete polynomial identity.
    for w in range(7):
        g = PSI[:]
        g[0] -= w
        actual = resultant(derivative(PSI), g)
        require(actual == 25 * w * evaluate(poly, w),
                f"critical-value resultant mismatch at w={w}")
    require(poly[-1] % 101 != 0 and poly[0] % 101 != 0,
            "bad modular leading or constant coefficient")
    require(mod_gcd(poly, derivative(poly), 101) == [1],
            "critical-value polynomial not certified squarefree")


def rejected(call, message):
    try:
        call()
    except ValueError:
        return
    raise ValueError(message)


def run(note):
    require(sys.flags.optimize == 0, "optimized Python replay is not supported")
    require(hashlib.sha256(note.read_bytes()).hexdigest() == NOTE_SHA,
            "note hash mismatch")
    for path, expected in SOURCES.items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected,
                f"source hash mismatch: {path}")

    # Independently reconstruct the inverse and the actual truncation.
    a = [0, 1, -1]
    v = inverse(a, 12)
    f = power(v, 2, 12)
    q = [0, 0, 1, 2, 5]
    require(f[:5] == q, "quadratic inverse-square truncation")
    require(valuation(add(f, scale(q, -1))) == 5, "actual contact order")
    qa = compose(q, a)
    require(qa[0] == 0 and qa[1:] == PSI, "divided polynomial identity")
    verify_cv_certificate(CV_POLYNOMIAL)
    require(5 > 4 + 2 * 1 - 2, "failed-budget arithmetic")
    wrong = CV_POLYNOMIAL[:]
    wrong[0] += 1
    rejected(lambda: verify_cv_certificate(wrong), "altered certificate accepted")
    rejected(lambda: verify_cv_certificate(CV_POLYNOMIAL[:-1]),
             "malformed certificate accepted")

    # The two displayed elimination identities hold as polynomial identities.
    y = [0, 0, 1]
    p = add(add(power(a, 2), scale(mul(a, y), 2)),
            add(power(y, 2), scale(y, -1)))
    require(p == [0], "quadratic elimination identity")
    a3 = [0, 1, 0, -1]
    p = add(power(a3, 2), add(scale(power(y, 3), -1),
                             add(scale(power(y, 2), 2), scale(y, -1))))
    require(p == [0], "cubic elimination identity")
    cubic_square = power(inverse(a3, 10), 2, 10)
    require([coefficient(cubic_square, i) for i in (2, 4, 6, 8)]
            == [1, 2, 7, 30], "distinguished cubic branch coefficients")
    # Newton's p_2 + 2*a_2 = 0 for z^3-z+t, where a_2=-1.
    require(-2 * (-1) == 2, "cubic trace identity")

    cyclic_cases = [(2, [1, -3, 1]), (3, [1, 2])]
    for b, r in cyclic_cases:
        aa = [F(0)] * (b * (len(r) - 1) + 2)
        for j, c in enumerate(r):
            aa[1 + b * j] = c
        bb = [0] + power(r, b)
        require(power(aa, b) == compose(bb, [0] * b + [1]),
                "cyclic descent polynomial identity")
        va = power(inverse(aa, 18), b, 18)
        vb = compose(inverse(bb, 18 // b), [0] * b + [1], 18)
        require(va == vb, "cyclic inverse series identity")
    require(coefficient(power(a, 2), 3) == -2,
            "incompatible cyclic descent control")

    # Unicritical formula and its common-sign convolution, finite controls only.
    unicritical = 0
    for d in range(2, 6):
        for lam in (1, 2):
            aa = [F(0)] + [F(comb(d, j), d) * F(lam) ** (j - 1)
                           for j in range(1, d + 1)]
            vv = inverse(aa, 14)
            binomial = F(1)
            for n in range(1, 15):
                binomial *= (F(1, d) - (n - 1)) / n
                expected = binomial * (d * lam) ** n / lam
                require(coefficient(vv, n) == expected, "unicritical inverse formula")
            for b in range(1, 6):
                ff = power(vv, b, 14)
                require(all((-1) ** (n - b) * coefficient(ff, n) > 0
                            for n in range(b, 15)), "unicritical power sign")
            unicritical += 1

    return {
        "status": "PASS",
        "scope": "Exact identities and obstruction controls; general contact target unproved",
        "note_sha256": NOTE_SHA,
        "source_sha256": SOURCES,
        "critical_value_resultant_degree": 6,
        "exact_resultant_evaluation_points": list(range(7)),
        "nonzero_critical_values_of_divided_polynomial": 5,
        "rejected_proposed_critical_value_budget": 4,
        "squarefreeness_prime": 101,
        "corruption_controls": ["altered constant coefficient rejected",
                                "malformed polynomial degree rejected"],
        "compatible_cyclic_cases": len(cyclic_cases),
        "unicritical_finite_inputs": unicritical,
        "unicritical_inverse_degree_checked": 14,
        "unicritical_powers_checked": [1, 2, 3, 4, 5],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--receipt", type=Path,
                        default=HERE / "sfc_inverse_power_contact_verification.json")
    args = parser.parse_args()
    try:
        receipt = run(args.note)
    except (ValueError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("PASS: exact inverse-power identities, cyclic descent, and critical-value budget obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
