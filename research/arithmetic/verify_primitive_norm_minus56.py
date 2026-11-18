#!/usr/bin/env python3
"""Finite arithmetic checks for principal primitive norms x^2+14*y^2.

Default scope: all integers through 10^6, direct gcd-restricted coordinates
versus the exact C4 class formula; finite local and absorber-profile checks.
This does not prove eventual additive representability. Standard library only.

python3 research/arithmetic/verify_primitive_norm_minus56.py \
  --output research/arithmetic/primitive_norm_minus56_checks.json
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
from hashlib import sha256
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import json


def prime_tables(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start:limit + 1:p] = b"\x00" * ((limit - start) // p + 1)
    primes = [p for p in range(2, limit + 1) if flags[p]]
    spf = array("I", range(limit + 1))
    for p in primes:
        if p * p > limit:
            break
        for k in range(p * p, limit + 1, p):
            if spf[k] == k:
                spf[k] = p
    return flags, primes, spf


def splitting(p):
    assert p not in (2, 7)
    a = pow((-14) % p, (p - 1) // 2, p)
    assert a in (1, p - 1)
    return 1 if a == 1 else -1


def psi8(n):
    assert n % 2
    return 1 if n % 8 in (1, 7) else -1


def reduced_prime_labels(limit, flags, primes):
    """1=principal, 2=nonzero class in 2G, 3=either odd class."""
    labels = bytearray(limit + 1)

    def record(n, code):
        if n <= 1 or n > limit or n in (2, 7) or not flags[n]:
            return
        assert labels[n] in (0, code), (n, labels[n], code)
        labels[n] = code

    for y in range(isqrt(limit // 14) + 1):
        for x in range(isqrt(limit - 14 * y * y) + 1):
            record(x * x + 14 * y * y, 1)
    for y in range(isqrt(limit // 7) + 1):
        for x in range(isqrt((limit - 7 * y * y) // 2) + 1):
            record(2 * x * x + 7 * y * y, 2)
    # Q3=3(x+y/3)^2+14*y^2/3. Conjugation permits y>=0.
    # The true |x| bound is <sqrt(limit), a safe integer loop bound.
    for y in range(isqrt(3 * limit // 14) + 1):
        for x in range(-isqrt(limit), isqrt(limit) + 1):
            record(3 * x * x + 2 * x * y + 5 * y * y, 3)
    totals = Counter()
    for p in primes:
        if p in (2, 7):
            continue
        if splitting(p) == -1:
            assert labels[p] == 0, (p, labels[p])
            totals["inert"] += 1
        else:
            assert labels[p] in (1, 2, 3), p
            assert psi8(p) == (1 if labels[p] in (1, 2) else -1), p
            totals[{1: "principal_split", 2: "class2_split", 3: "odd_class_split"}[labels[p]]] += 1
    return labels, dict(sorted(totals.items()))


def direct_primitive_counts(limit):
    counts = array("I", [0]) * (limit + 1)
    for y in range(isqrt(limit // 14) + 1):
        for x in range(isqrt(limit - 14 * y * y) + 1):
            if gcd(x, y) == 1:
                n = x * x + 14 * y * y
                if n:
                    counts[n] += (1 if x == 0 else 2) * (1 if y == 0 else 2)
    return counts


def arithmetic_count(n, labels, spf):
    """Return (R*, individually locally admissible, r, t, s)."""
    a, exps = n, {}
    while a > 1:
        p, e = spf[a], 0
        while a % p == 0:
            a //= p
            e += 1
        exps[p] = e
    e, f = exps.pop(2, 0), exps.pop(7, 0)
    if e > 1 or f > 1:
        return 0, False, 0, 0, 0
    r, t, s = 0, (e + f) % 2, 0
    for p, exponent in exps.items():
        label = labels[p]
        if label == 0:
            return 0, False, 0, 0, 0
        s += 1
        if label == 2:
            t ^= exponent % 2
        elif label == 3:
            if exponent % 2:
                r += 1
            else:
                t ^= (exponent // 2) % 2
    locally = r % 2 == 0
    expected = (
        (1 << s) if r and r % 2 == 0
        else (1 << (s + 1)) if r == 0 and t == 0
        else 0
    )
    return expected, locally, r, t, s


def local_values(p, h):
    modulus = p ** h
    return {
        (x * x + 14 * y * y) % modulus
        for x in range(modulus) for y in range(modulus)
        if x % p or y % p
    }


def local_checks():
    rows = []
    for p, max_h in ((2, 9), (7, 3), (3, 2), (5, 2), (11, 1), (13, 1)):
        for h in range(1, max_h + 1):
            modulus = p ** h
            values = local_values(p, h)
            sums = {(a + b) % modulus for a in values for b in values}
            if p == 2:
                assert sums == set(range(modulus)), (p, h)
                if h >= 4:
                    expected = {
                        n for n in range(modulus)
                        if n % 8 in (1, 7)
                        or (n % 2 == 0 and (n // 2) % 8 in (1, 7))
                    }
                    assert values == expected, (p, h)
            elif p == 7:
                assert sums == {
                    n for n in range(modulus) if h == 1 or n % 49
                }, (p, h)
                if h >= 2:
                    qr = {1, 2, 4}
                    expected = {
                        n for n in range(modulus)
                        if n % 7 in qr
                        or (n % 7 == 0 and (n // 7) % 7 in qr)
                    }
                    assert values == expected, (p, h)
            else:
                assert sums == set(range(modulus)), (p, h)
                if splitting(p) == -1:
                    assert values == {n for n in range(modulus) if n % p}, (p, h)
                else:
                    assert values == set(range(modulus)), (p, h)
            rows.append({
                "prime": p, "depth": h, "modulus": modulus,
                "primitive_value_residues": len(values), "selfsum_residues": len(sums),
            })
    return rows


def signed_orientation_checks():
    cases = 0
    for s in range(8):
        # Entries are already-exponentiated ideal classes in Z/4.
        for classes in product(range(4), repeat=s):
            for prefix in (0, 2):
                distribution = [0, 0, 0, 0]
                distribution[prefix] = 1
                for c in classes:
                    nxt = [0, 0, 0, 0]
                    for a, multiplicity in enumerate(distribution):
                        # Equal classes still correspond to two distinct ideals.
                        nxt[(a + c) % 4] += multiplicity
                        nxt[(a - c) % 4] += multiplicity
                    distribution = nxt
                r = sum(c % 2 for c in classes)
                t = (prefix // 2 + sum(c == 2 for c in classes)) % 2
                formula = (
                    (1 << s) if r and r % 2 == 0
                    else (1 << (s + 1)) if r == 0 and t == 0
                    else 0
                )
                assert 2 * distribution[0] == formula, (classes, prefix, distribution)
                cases += 1
    return cases


def absorber_checks(flags, labels):
    pool = (3, 5, 13, 19)
    witnesses = ((1, 0), (0, 1), (-2, 1), (-1, 2))
    for q, (x, y) in zip(pool, witnesses):
        assert flags[q] and labels[q] == 3 and psi8(q) == -1
        assert 3 * x * x + 2 * x * y + 5 * y * y == q
    assert {(a + b) % 4 for a in (1, 3) for b in (1, 3)} == {0, 2}
    templates, profile_checks = [], 0
    for mask in range(16):
        dividing = [q for i, q in enumerate(pool) if mask >> i & 1]
        nondividing = [q for q in pool if q not in dividing]
        left, right = set(), set()
        if len(dividing) >= 2:
            left.update(dividing[:2])
            right.update(dividing[:2])
        elif dividing:
            left.update((dividing[0], nondividing[0]))
            right.update((dividing[0], nondividing[1]))
        else:
            left.update(pool[:2])
            right.update(pool[2:])
        assert len(left) == len(right) == 2
        for q in pool:
            a, b = int(q in left), int(q in right)
            if a and b:
                assert q in dividing
            elif a or b:
                assert q not in dividing
            q2 = q * q
            allowed_a = [
                z for z in range(q2)
                if (z % q != 0 if a == 0 else z % q == 0 and z != 0)
            ]
            allowed_b = {
                z for z in range(q2)
                if (z % q != 0 if b == 0 else z % q == 0 and z != 0)
            }
            for target in range(q2):
                if (target % q == 0) != (q in dividing):
                    continue
                assert any((target - z) % q2 in allowed_b for z in allowed_a), (mask, q, target)
                profile_checks += 1
        templates.append({
            "divisibility_mask": mask,
            "left_exact_one": sorted(left), "right_exact_one": sorted(right),
        })
    return {
        "pool": list(pool), "form_witnesses": [list(v) for v in witnesses],
        "templates": templates, "exact_profile_target_residues_checked": profile_checks,
    }


def run(limit):
    assert limit >= 100
    flags, primes, spf = prime_tables(limit)
    labels, prime_summary = reduced_prime_labels(limit, flags, primes)
    coordinates = direct_primitive_counts(limit)
    counts, failures = Counter(), []
    for n in range(1, limit + 1):
        predicted, locally, r, t, s = arithmetic_count(n, labels, spf)
        assert predicted == coordinates[n], (n, predicted, coordinates[n], r, t, s)
        if predicted:
            counts["principal_primitive_values"] += 1
            counts["primitive_coordinate_pairs"] += predicted
        if locally:
            counts["individually_locally_admissible_values"] += 1
            if not predicted:
                counts["individual_local_global_failures"] += 1
                if len(failures) < 30:
                    failures.append(n)
    controls = []

    def reject(name, n, wrong_count):
        correct = arithmetic_count(n, labels, spf)[0]
        assert correct != wrong_count and coordinates[n] == correct
        controls.append({
            "name": name, "integer": n,
            "rejected_count": wrong_count, "actual_count": correct,
        })

    reject("omit_nonprincipal_ramified_class_2", 2, 2)
    reject("omit_nonprincipal_ramified_class_7", 7, 2)
    reject("omit_even_oddclass_power_contribution", 9, 4)
    reject("omit_genus_parity", 3, 2)
    reject("double_active_orientation_count", 15, 8)
    reject("omit_ramified_cancellation", 18, 0)
    reject("reduce_prime_exponent_modulo_2_only", 81, 0)
    for n in (2, 7, 9):
        expected, locally, *_ = arithmetic_count(n, labels, spf)
        assert locally and expected == 0
    family, a = [], 0
    while 3 ** (4 * a + 2) <= limit:
        n = 3 ** (4 * a + 2)
        expected, locally, *_ = arithmetic_count(n, labels, spf)
        assert locally and expected == 0
        family.append(n)
        a += 1
    return {
        "schema": "primitive-norm-minus56-finite-check-v1",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": {
            "field_discriminant": -56, "form": [1, 0, 14],
            "coordinate_and_class_formula_integer_range": [1, limit],
            "unrestricted_additive_theorem_proved_by_this_check": False,
            "interpretation": "Finite coordinate/class-count comparison and residue/profile checks; infinite arithmetic and all-depth local statements require separate mathematical proofs.",
        },
        "prime_classification": prime_summary,
        "integer_comparison": dict(sorted(counts.items())),
        "first_individual_local_global_failures": failures,
        "infinite_failure_family_tested_terms": family,
        "local_residue_checks": local_checks(),
        "signed_orientation_class_sequences_checked": signed_orientation_checks(),
        "absorber_checks": absorber_checks(flags, labels),
        "deliberately_corrupted_formulas_rejected": controls,
        "status": "pass",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.limit)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    print(json.dumps({
        "status": result["status"], "limit": args.limit,
        "integer_comparison": result["integer_comparison"],
        "source_sha256": result["source_sha256"],
        "output": str(args.output) if args.output else None,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
