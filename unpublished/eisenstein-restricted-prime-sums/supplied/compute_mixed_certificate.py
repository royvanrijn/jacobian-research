#!/usr/bin/env python3
"""Exact finite certificate for N=a+3b with a,b supported on primes 1 mod3.

N=6k+4, a=6i+1 and b=6j+1 give i+3j=k.  The integer 1 is allowed.
The default inclusive bound is 2,000,000.  This is only a finite check.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np


def membership(last_index):
    largest = 6 * last_index + 1
    root = math.isqrt(largest)
    primes = np.ones(root + 1, dtype=np.bool_)
    primes[:2] = False
    for p in range(2, math.isqrt(root) + 1):
        if primes[p]:
            primes[p*p::p] = False
    good = np.ones(last_index + 1, dtype=np.bool_)
    for p in np.flatnonzero(primes):
        p = int(p)
        if p % 6 == 5:
            good[(5*p - 1)//6::p] = False
    return good


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=2000000)
    parser.add_argument("--output", type=Path, default=Path("mixed_finite_2000000"))
    args = parser.parse_args()
    if args.limit < 4:
        parser.error("limit must be >=4")
    started = time.monotonic()
    last_index = (args.limit - 4) // 6
    good = membership(last_index)
    missing = np.arange(last_index + 1, dtype=np.int64)
    witnesses = np.full(last_index + 1, -1, dtype=np.int32)
    histogram = []
    for j in np.flatnonzero(good):
        if not len(missing) or 3*j > missing[-1]:
            break
        eligible = missing[missing >= 3*j]
        hits = eligible[good[eligible - 3*j]]
        if len(hits):
            witnesses[hits] = j
            histogram.append({"b": int(6*j+1), "weighted_summand_3b": int(18*j+3), "new_targets": len(hits)})
            missing = missing[witnesses[missing] < 0]
    present = np.flatnonzero(witnesses >= 0)
    js = witnesses[present].astype(np.int64)
    assert np.all(good[js] & good[present-3*js])
    for k in missing:
        choices = np.arange(int(k)//3+1)
        assert not np.any(good[choices] & good[int(k)-3*choices])
    max_j = int(witnesses.max())
    args.output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output / "mixed_witnesses.npz", least_b_index=witnesses)
    receipt = {
        "scope": "Exact finite check only: N=a+3b, a,b have only prime factors1mod3;1 allowed",
        "requested_inclusive_limit": args.limit,
        "largest_target": int(6*last_index+4),
        "number_of_targets": int(last_index+1),
        "exceptions": [int(6*k+4) for k in missing],
        "uncovered_count": len(missing),
        "largest_least_b": 6*max_j+1,
        "largest_least_weighted_summand_3b": 18*max_j+3,
        "largest_least_b_targets": [int(6*k+4) for k in np.flatnonzero(witnesses == max_j)],
        "least_b_histogram": histogram,
        "witness_encoding": "N=6k+4; entry j means b=6j+1 and a=6(k-3j)+1; -1 means no representation",
        "sha256_admissibility_bits": hashlib.sha256(good.tobytes()).hexdigest(),
        "sha256_int32_witness_indices": hashlib.sha256(witnesses.tobytes()).hexdigest(),
        "all_witnesses_and_exceptions_rechecked": True,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic()-started,3),
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__ == "__main__":
    main()
