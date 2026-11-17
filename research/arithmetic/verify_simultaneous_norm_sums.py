#!/usr/bin/env python3
"""Exact finite check for the prime-1-mod-12 binary conjecture.

Index c=12*i+1 and N=12*j+2.  See the accompanying research note for
the square-root sieve proof and the distinction between finite and global claims.
"""
import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np


def membership(limit):
    """Return C-membership for 1,13,... <= limit, using exact integer sieving."""
    size = (limit - 1) // 12 + 1
    good = np.ones(size, dtype=np.bool_)
    root = math.isqrt(limit)
    prime = bytearray(b"\x01") * (root + 1)
    prime[:2] = b"\x00\x00"
    for p in range(2, root + 1):
        if not prime[p]:
            continue
        if p * p <= root:
            prime[p * p::p] = b"\x00" * ((root - p * p) // p + 1)
        if p > 3 and p % 12 != 1:
            good[(-pow(12, -1, p)) % p::p] = False
    return good


def check_controls():
    good = membership(100_000)
    # Independent factorization of every represented integer in this range.
    for i, actual in enumerate(good):
        n, expected, p = 12 * i + 1, True, 2
        while p * p <= n:
            if n % p == 0:
                if p % 12 != 1:
                    expected = False
                while n % p == 0:
                    n //= p
            p += 1
        if n > 1 and n % 12 != 1:
            expected = False
        if bool(actual) != expected:
            raise RuntimeError(f"membership mismatch at {12*i+1}")
    for n in (1, 13, 37, 169, 481):
        if not good[(n - 1) // 12]:
            raise RuntimeError(f"positive control rejected: {n}")
    for n in (25, 49, 121, 145, 289):
        if good[(n - 1) // 12]:
            raise RuntimeError(f"negative control admitted: {n}")
    return {"independent_factorizations": len(good),
            "positive_controls": [1, 13, 37, 169, 481],
            "negative_controls": [25, 49, 121, 145, 289]}


def run(limit):
    if not 2 <= limit <= 2_000_000_000:
        raise ValueError("limit must be between 2 and 2,000,000,000")
    start = time.perf_counter()
    controls = check_controls()
    last = (limit - 2) // 12
    good = membership(12 * last + 1)
    witness = np.full(last + 1, -1, dtype=np.int32)
    pending = np.arange(last + 1, dtype=np.int32)
    failures = []
    distance = 0
    while len(pending):
        smaller = pending // 2 - distance
        in_balance = smaller >= (pending + 2) // 3
        failures.extend(map(int, pending[~in_balance]))
        pending = pending[in_balance]
        smaller = smaller[in_balance]
        hit = good[smaller] & good[pending - smaller]
        witness[pending[hit]] = smaller[hit]
        pending = pending[~hit]
        distance += 1
    failures.sort()
    unbalanced_failures = []
    exception_witnesses = {}
    for j in failures:
        candidates = np.flatnonzero(good[:j // 2 + 1])
        hits = candidates[good[j - candidates]]
        if not len(hits):
            unbalanced_failures.append(12 * j + 2)
        else:
            i = int(hits[0])
            exception_witnesses[str(12 * j + 2)] = [12 * i + 1, 12 * (j - i) + 1]
    # Direct scalar enumeration checks both endpoints and exhausted searches.
    for j in range(min(last, (5000 - 2) // 12) + 1):
        n = 12 * j + 2
        direct = any(good[i] and good[j - i] and 3 * (12 * i + 1) >= n
                     for i in range(j // 2 + 1))
        if direct != bool(witness[j] >= 0):
            raise RuntimeError(f"scalar balanced search disagrees at {n}")
    # Validate every retained witness in chunks, independently of its search path.
    for lo in range(0, last + 1, 1_000_000):
        j = np.arange(lo, min(lo + 1_000_000, last + 1), dtype=np.int64)
        i = witness[lo:lo + len(j)].astype(np.int64)
        present = i >= 0
        if not (np.all(good[i[present]]) and np.all(good[j[present]-i[present]])
                and np.all(3 * i[present] >= j[present])
                and np.all(2 * i[present] <= j[present])):
            raise RuntimeError("invalid balanced witness")
    if int(np.count_nonzero(witness < 0)) != len(failures):
        raise RuntimeError("failure accounting mismatch")
    return {
        "scope": "finite computation only; no unrestricted binary theorem",
        "limit_inclusive": limit,
        "last_target": 12 * last + 2,
        "target_count": last + 1,
        "balanced_exceptions": [12 * j + 2 for j in failures],
        "unbalanced_exceptions": unbalanced_failures,
        "unbalanced_witnesses_for_balanced_exceptions": exception_witnesses,
        "membership_count": int(np.count_nonzero(good)),
        "membership_sha256": hashlib.sha256(good.tobytes()).hexdigest(),
        "witness_encoding": "little-endian int32 index i of a=12*i+1; -1 for balanced exceptions; nearest below half chosen",
        "balanced_witness_sha256": hashlib.sha256(witness.astype("<i4").tobytes()).hexdigest(),
        "search_distances": distance,
        "controls": controls,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.perf_counter() - start, 3),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(args.limit), indent=2) + "\n"
    if args.output:
        args.output.write_text(result, encoding="utf-8")
    print(result, end="")
