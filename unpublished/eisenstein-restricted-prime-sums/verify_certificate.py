"""Exact finite certificate for targets n == 2 (mod 6), 2 <= n < 10**8."""
from hashlib import sha256
from math import isqrt
from struct import pack
import json


def membership(last):
    limit = isqrt(6 * last + 1)
    prime = bytearray(b"\1") * (limit + 1)
    prime[:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            prime[p * p::p] = b"\0" * ((limit - p * p) // p + 1)
    good = bytearray(b"\1") * (last + 1)
    for p in range(5, limit + 1, 6):
        if prime[p]:
            first = (5 * p - 1) // 6
            good[first::p] = b"\0" * max(0, (last - first) // p + 1)
    return good


def certificate(good, candidates):
    if not good or candidates != sorted(set(candidates)):
        raise ValueError("Nonempty membership and sorted distinct candidates required")
    for t in candidates:
        if t < 0 or t >= len(good) or not good[t] or 6 * t + 1 > 211:
            raise ValueError("Invalid small-summand index")
        a = 6 * t + 1
        if any(a % d == 0 and d % 3 != 1 for d in range(2, a + 1)):
            raise ValueError("Small summand has a forbidden prime factor")
    counts = dict.fromkeys(candidates, 0)
    witnesses = sha256()
    for k in range(len(good)):
        t = next((t for t in candidates if 2 * t <= k and good[k - t]), None)
        if t is None:
            raise RuntimeError(f"Uncovered target {6 * k + 2}")
        counts[t] += 1
        witnesses.update(pack("<i", t))
    histogram = [{"summand": 6 * t + 1, "new_targets": c}
                 for t, c in counts.items() if c]
    return {"largest_target": 6 * len(good) - 4,
            "number_of_targets": len(good),
            "allowed_summands_count": sum(good),
            "unrestricted_exceptions": [],
            "largest_least_summand": histogram[-1]["summand"],
            "least_summand_histogram": histogram,
            "sha256_allowed_boolean_array": sha256(good).hexdigest(),
            "sha256_unrestricted_int32_witness_indices": witnesses.hexdigest()}


if __name__ == "__main__":
    last = (100_000_000 - 3) // 6
    good = membership(last)
    candidates = [t for t in range(36) if good[t]]
    print(json.dumps(certificate(good, candidates), indent=2))
