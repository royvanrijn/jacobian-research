#!/usr/bin/env python3
"""Exact finite tests for sums of integers supported on primes 1 (mod 3).

Run: python compute_finite_certificate.py --limit 2000000 --output results

The tested targets are N = 6k+2 <= limit.  Index j represents 6j+1.
Numbers at these indices already exclude the primes 2 and 3.  If a number
at an index has an offending prime factor q = 5 (mod 6), then it has an
even number of such factors counted with multiplicity; at least one is
<= sqrt(number).  Thus the sieve below is an exact membership test.

The compressed witness arrays are independently replayed before a receipt
is written.  This proves only the stated finite range, not all targets.
"""

import argparse
import hashlib
import json
import math
import pathlib
import time

import numpy as np


def digest(array):
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def primes_through(n):
    flags = np.ones(n + 1, dtype=np.bool_)
    flags[:2] = False
    for p in range(2, math.isqrt(n) + 1):
        if flags[p]:
            flags[p * p :: p] = False
    return np.flatnonzero(flags)


def allowed_indices(max_index):
    maximum = 6 * max_index + 1
    good = np.ones(max_index + 1, dtype=np.bool_)
    bad_primes = [int(p) for p in primes_through(math.isqrt(maximum)) if p % 3 == 2]
    for p in bad_primes:
        if p == 2:
            continue
        # 5p is the first positive multiple of p congruent to 1 mod 6.
        start = (5 * p - 1) // 6
        good[start::p] = False
    return good, bad_primes


def least_witnesses(good):
    count = len(good)
    missing = np.arange(count, dtype=np.int64)
    witness = np.full(count, -1, dtype=np.int32)
    histogram = []
    for t in np.flatnonzero(good):
        eligible = missing[missing >= 2 * t]
        hits = eligible[good[eligible - t]]
        if len(hits):
            witness[hits] = t
            histogram.append({"summand": int(6 * t + 1), "new_targets": int(len(hits))})
            missing = missing[witness[missing] < 0]
        if not len(missing):
            break
        # Once t exceeds half the largest missing index, no further pair exists.
        if 2 * t > missing[-1]:
            break
    return witness, missing, histogram


def balanced_witnesses(good, chunk_size=1000000):
    count = len(good)
    witness = np.full(count, -1, dtype=np.int32)
    maximum_search_offset = 0
    max_search_target = 2
    for lo in range(0, count, chunk_size):
        missing = np.arange(lo, min(count, lo + chunk_size), dtype=np.int64)
        delta = 0
        while len(missing):
            t = missing // 2 - delta
            eligible_mask = t >= (missing + 2) // 3
            eligible = missing[eligible_mask]
            t = t[eligible_mask]
            if not len(eligible):
                break
            found = good[t] & good[eligible - t]
            hits = eligible[found]
            witness[hits] = t[found]
            if len(hits) and delta > maximum_search_offset:
                maximum_search_offset = delta
                max_search_target = int(6 * hits[-1] + 2)
            missing = eligible[~found]
            delta += 1
    missing = np.flatnonzero(witness < 0)
    return witness, missing, maximum_search_offset, max_search_target


def verify_witnesses(good, witness, balanced):
    present = np.flatnonzero(witness >= 0)
    t = witness[present].astype(np.int64)
    assert np.all(t <= present - t)
    assert np.all(good[t] & good[present - t])
    if balanced:
        assert np.all(t >= (present + 2) // 3)
    for k in np.flatnonzero(witness < 0):
        low = (int(k) + 2) // 3 if balanced else 0
        candidates = np.arange(low, int(k) // 2 + 1)
        assert not np.any(good[candidates] & good[int(k) - candidates])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=2000000)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("finite_certificate"))
    parser.add_argument("--chunk-size", type=int, default=1000000)
    args = parser.parse_args()
    if args.limit < 2:
        parser.error("limit must be at least 2")
    start = time.monotonic()
    max_index = (args.limit - 2) // 6
    good, bad_primes = allowed_indices(max_index)
    unrestricted, misses, histogram = least_witnesses(good)
    balanced, balanced_misses, max_offset, max_offset_target = balanced_witnesses(good, args.chunk_size)
    verify_witnesses(good, unrestricted, False)
    verify_witnesses(good, balanced, True)
    args.output.mkdir(parents=True, exist_ok=True)
    witness_path = args.output / "witnesses.npz"
    np.savez_compressed(witness_path, unrestricted=unrestricted, balanced=balanced)
    receipt = {
        "scope": "Exact finite enumeration only; targets N=6k+2 <= requested_limit",
        "requested_limit": args.limit,
        "largest_target": 6 * max_index + 2,
        "number_of_targets": max_index + 1,
        "allowed_summands_count": int(good.sum()),
        "sieve_prime_limit": math.isqrt(6 * max_index + 1),
        "offending_primes_sieved_count_including_two": len(bad_primes),
        "unrestricted_exceptions": [int(6 * k + 2) for k in misses],
        "largest_least_summand": int(6 * unrestricted.max() + 1),
        "least_summand_histogram": histogram,
        "balanced_condition": "N/3 <= a <= b <= 2N/3",
        "balanced_exceptions": [int(6 * k + 2) for k in balanced_misses],
        "balanced_max_search_offset_in_index_units": max_offset,
        "balanced_max_search_offset_last_target": max_offset_target,
        "all_witnesses_and_all_listed_exceptions_rechecked": True,
        "sha256_allowed_boolean_array": digest(good),
        "sha256_unrestricted_int32_witness_indices": digest(unrestricted),
        "sha256_balanced_int32_witness_indices": digest(balanced),
        "witness_encoding": "N=6k+2; entry t>=0 means a=6t+1, b=6(k-t)+1; -1 means no representation in scope",
        "witness_file": witness_path.name,
        "script_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: v for k, v in receipt.items() if k != "least_summand_histogram"}, indent=2))


if __name__ == "__main__":
    main()
