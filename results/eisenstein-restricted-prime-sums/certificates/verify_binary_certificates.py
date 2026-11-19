"""Exact binary certificates below10**8: a+b and a+3b, a,b in A."""
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


def certificate(good, weight, cap, candidates=None):
    if not good or weight not in (1, 3) or cap < 1:
        raise ValueError("Nonempty membership, weight1 or3, positive cap required")
    if candidates is None:
        candidates = [t for t in range(min(len(good), (cap - 1) // 6 + 1))
                      if good[t]]
    if candidates != sorted(set(candidates)):
        raise ValueError("Sorted distinct candidate indices required")
    for t in candidates:
        if t < 0 or t >= len(good) or not good[t] or 6 * t + 1 > cap:
            raise ValueError("Invalid small-summand index")
        a = 6 * t + 1
        if any(a % d == 0 and d % 3 != 1 for d in range(2, a + 1)):
            raise ValueError("Small summand has a forbidden prime factor")
    counts = dict.fromkeys(candidates, 0)
    witnesses = sha256()
    order = 2 if weight == 1 else 3
    for k in range(len(good)):
        t = next((t for t in candidates
                  if order * t <= k and good[k - weight * t]), None)
        if t is None:
            raise RuntimeError(f"Uncovered target {6 * k + weight + 1}")
        counts[t] += 1
        witnesses.update(pack("<i", t))
    histogram = [{"summand": 6 * t + 1, "new_targets": c}
                 for t, c in counts.items() if c]
    return {"weight": weight, "largest_target": 6 * len(good) + weight - 5,
            "number_of_targets": len(good), "allowed_count": sum(good),
            "exceptions": [], "largest_least_small_summand": histogram[-1]["summand"],
            "histogram": histogram, "sha256_membership_uint8": sha256(good).hexdigest(),
            "sha256_small_index_int32_le": witnesses.hexdigest()}


if __name__ == "__main__":
    good = membership((100_000_000 - 3) // 6)
    results = {}
    for name, weight, cap in [("unrestricted", 1, 211), ("mixed", 3, 199)]:
        last = (100_000_000 - weight - 2) // 6
        results[name] = certificate(good[:last + 1], weight, cap)
    print(json.dumps(results, indent=2))
