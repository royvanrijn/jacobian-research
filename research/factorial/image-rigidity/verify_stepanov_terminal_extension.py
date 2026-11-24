#!/usr/bin/env python3
"""Exact terminal-11/12 factorial-moment certificates; stdlib replay.

Generation optionally uses python-flint. Verification uses only Python's
standard library, reconstructs the entire finite domain, and checks every
full polynomial Bezout identity. The written proof supplies the unbounded
exponent complements; this program verifies their exact numerical inputs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path


PARAMETERS = ((11, 14, 16, 57, 9, 14316), (12, 16, 17, 65, 10, 20745))
PRIMES = (65537, 1000003, 1000033)
LAST = 8436
TOTAL = 35061
PREFIX = Path("research/factorial/image-rigidity")
SOURCES = {
    "lower_exponent_initial_and_boundary.tex": "e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45",
    "GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md": "22700eab6cdb655dc676d75158f17ca5350fe0b7684fc6fe933c5a1195f2db9c",
    "GAMMA_PRODUCT_ROOT_MODULUS.tex": "b6cb212dd2b461a12c5f687375089da7839495031737ef04d0de2cc8151950b9",
    "GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md": "5ddbe130e436416808dda833a06b6512a0d2b62eda682ccd5029e8414f32f473",
    "ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex": "604e1687ea759785fe5d97a2cf214c84732be10ac0df0cae208d7bd01435dcd7",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def domain():
    for n, _, _, gap, constant_gap, count in PARAMETERS:
        seen = 0
        for d in range(1, gap):
            for a in range(1, d * (n - 2)):
                seen += 1
                yield n, a, d
        for d in range(1, constant_gap):
            seen += 1
            yield n, 0, d
        require(seen == count, "Domain count")


def is_prime(p):
    return type(p) is int and p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def complement_checks():
    rows = []
    for n, c0, c1, gap, constant_gap, count in PARAMETERS:
        require((c0 + 1) ** (n - 1) < 2 * c0 ** (n - 1), "Lower radius")
        require((c1 + 1) ** n < 2 * c1 ** n, "Upper radius")
        b0 = math.factorial((gap + 1) * (n - 1)) // math.factorial(n - 1)
        b1 = math.factorial((gap + 1) * n) // math.factorial(n)
        left = b1 ** (n - 1)
        right = (c0 * c1) ** (n * (n - 1)) * b0 ** n
        require(left > right, "Disjoint rational annuli")
        require(2 ** constant_gap >= 4 * n * n, "Constant-family large gap")
        require((n - 2) * (gap - 1) * gap // 2 - (gap - 1) + constant_gap - 1 == count,
                "Closed domain count")
        rows.append({"terminal": n, "C_previous": c0, "C_upper": c1,
                     "gap_cutoff": gap, "constant_gap_cutoff": constant_gap,
                     "cases": count, "annulus_inequality": True,
                     "annulus_difference_bits": (left - right).bit_length(),
                     "annulus_difference_sha256": hashlib.sha256(
                         (left - right).to_bytes(((left - right).bit_length() + 7) // 8, "big")
                     ).hexdigest()})
    require(sum(row[5] for row in PARAMETERS) == TOTAL, "Total domain count")
    require(max(n * (a + d) for n, a, d in domain()) == LAST, "Maximum factorial")
    for p in PRIMES:
        require(is_prime(p) and p > LAST, "Good prime")
    return rows


def factorials(p):
    table = [1] * (LAST + 1)
    for j in range(1, LAST + 1):
        table[j] = table[j - 1] * j % p
    return table


def coefficients(m, a, d, p, table):
    # Ascending powers of x; literal reversed factorial-moment formula.
    denominator = pow(table[a * m], -1, p)
    return [math.comb(m, j) * table[a * m + d * (m - j)] * denominator % p
            for j in range(m + 1)]


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def product(a, b, p):
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
    return trim(result)


def bezout_sum(u, upper, v, lower, p):
    one, two = product(u, upper, p), product(v, lower, p)
    result = [0] * max(len(one), len(two))
    for poly in (one, two):
        for j, c in enumerate(poly):
            result[j] = (result[j] + c) % p
    return trim(result)


def coefficient_controls():
    # Derive coefficients independently by multiplying the two-term phase,
    # then applying k! coefficient by coefficient over the rationals.
    controls = 0
    for n in (2, 3, 11, 12):
        for a, d in ((0, 3), (1, 2), (3, 2), (4, 3)):
            for m in (n - 1, n):
                phase = {a: 2, a + d: -3}
                current = {0: 1}
                for _ in range(m):
                    updated = {}
                    for k, v in current.items():
                        for l, w in phase.items():
                            updated[k + l] = updated.get(k + l, 0) + v * w
                    current = updated
                original = sum(c * math.factorial(k) for k, c in current.items())
                direct = sum(math.comb(m, k) * math.factorial(a * m + d * k)
                             * 2 ** (m - k) * (-3) ** k for k in range(m + 1))
                require(original == direct, "Original factorial functional")
                for p in PRIMES[:1]:
                    table = factorials(p)
                    c = coefficients(m, a, d, p, table)
                    value = sum(cj * pow(2 * pow(-3, -1, p), j, p) for j, cj in enumerate(c)) % p
                    require(original % p == pow(-3, m, p) * table[a * m] * value % p,
                            "Reversal and factorial normalization")
                controls += 1
    # An actual common factor must never be accepted by an alleged identity.
    p = PRIMES[0]
    require(bezout_sum([1], [1, 2, 1], [1], [1, 1], p) != [1], "Common-factor negative")
    # Ordinary xgcd witnesses for x and x+1; alteration must be rejected.
    require(bezout_sum([p - 1], [0, 1], [1], [1, 1], p) == [1], "Bezout positive control")
    require(bezout_sum([p - 2], [0, 1], [1], [1, 1], p) != [1], "Altered-witness negative")
    return {"literal_moment_controls": controls, "negative_controls": 2}


def generate(path):
    from flint import nmod_poly
    tables = {p: factorials(p) for p in PRIMES}
    exceptions = []
    count = 0
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            for n, a, d in domain():
                for p in PRIMES:
                    lower = nmod_poly(coefficients(n - 1, a, d, p, tables[p]), p)
                    upper = nmod_poly(coefficients(n, a, d, p, tables[p]), p)
                    g, u, v = upper.xgcd(lower)
                    require(u * upper + v * lower == g, "Generation identity")
                    if g == nmod_poly([1], p):
                        record = {"N": n, "a": a, "d": d, "p": p,
                                  "u": [int(c) for c in u], "v": [int(c) for c in v]}
                        compressed.write((json.dumps(record, separators=(",", ":")) + "\n").encode())
                        count += 1
                        break
                    exceptions.append({"N": n, "a": a, "d": d, "p": p, "gcd_degree": g.degree()})
                else:
                    raise ValueError(f"Unresolved support {(n, a, d)}")
    require(count == TOTAL, "Generation coverage")
    return exceptions


def verify(path):
    tables = {p: factorials(p) for p in PRIMES}
    per_prime = {str(p): 0 for p in PRIMES}
    per_terminal = {str(n): 0 for n, *_ in PARAMETERS}
    nondefault = []
    count = 0
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for n, a, d in domain():
            line = stream.readline()
            require(bool(line), "Truncated certificate")
            row = json.loads(line)
            require(set(row) == {"N", "a", "d", "p", "u", "v"}, "Record schema")
            require(all(type(row[k]) is int for k in ("N", "a", "d", "p")), "Integer parameters")
            require((row["N"], row["a"], row["d"]) == (n, a, d), "Ordered complete domain")
            p = row["p"]
            require(p in PRIMES, "Permitted good prime")
            for key, limit in (("u", n - 1), ("v", n)):
                vector = row[key]
                require(type(vector) is list and 0 < len(vector) <= limit, "Witness degree")
                require(all(type(c) is int and 0 <= c < p for c in vector), "Canonical residues")
            lower = coefficients(n - 1, a, d, p, tables[p])
            upper = coefficients(n, a, d, p, tables[p])
            require(len(lower) == n and len(upper) == n + 1, "Actual full degrees")
            require(lower[-1] == upper[-1] == 1 and lower[0] and upper[0], "Unit endpoints")
            require(bezout_sum(row["u"], upper, row["v"], lower, p) == [1], "Full Bezout identity")
            if a == 0:
                require([(j + 1) * upper[j + 1] % p for j in range(n)] ==
                        [n * c % p for c in lower], "Appell boundary identity")
            if p != PRIMES[0]:
                nondefault.append({"N": n, "a": a, "d": d, "p": p})
            per_prime[str(p)] += 1
            per_terminal[str(n)] += 1
            count += 1
        require(stream.read() == "", "Extra certificate records")
    require(count == TOTAL, "Verified full count")
    return {"full_Bezout_identities": count, "per_terminal": per_terminal,
            "per_prime": per_prime, "nondefault_prime_certificates": nondefault}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root", required=True, type=Path)
    ap.add_argument("--certificate", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--generate", action="store_true")
    ns = ap.parse_args()
    source_hashes = {name: file_hash(ns.repo_root / PREFIX / name) for name in SOURCES}
    require(source_hashes == SOURCES, "Changed pinned mathematical proof input")
    complements = complement_checks()
    controls = coefficient_controls()
    exceptions = None
    ns.certificate.parent.mkdir(parents=True, exist_ok=True)
    if ns.generate:
        exceptions = generate(ns.certificate)
    result = verify(ns.certificate)
    require(source_hashes == {name: file_hash(ns.repo_root / PREFIX / name) for name in SOURCES},
            "Proof sources changed during replay")
    receipt = {"schema": 1, "status": "PASS", "terminals": [11, 12],
               "conclusion": "all ordinary univariate binomials, all characteristic-zero fields, terminal 11 or 12",
               "all_orders": False, "multivariate_extension": False,
               "formal_verification": False, "uniform_eventual_theorem_premise": False,
               "source_sha256": source_hashes,
               "checker_sha256": file_hash(Path(__file__)),
               "certificate_sha256": file_hash(ns.certificate),
               "certificate_bytes": ns.certificate.stat().st_size,
               "maximum_factorial_argument": LAST,
               "complement_inequalities": complements, "finite": result, "controls": controls}
    if exceptions is not None:
        receipt["generation_exceptional_reductions"] = exceptions
    ns.output.parent.mkdir(parents=True, exist_ok=True)
    ns.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "cases": TOTAL, "finite": result,
                      "certificate_bytes": receipt["certificate_bytes"], "output": str(ns.output)}))


if __name__ == "__main__":
    main()
