#!/usr/bin/env python3
"""Exact finite support for the discriminant -19 sumset classification.

The written proof supplies the spectral theorem and infinite descent.
This checker derives exact constants and covers every required finite
target through 2**34, in 512 consecutive segments. Standard library only.
"""

if not __debug__:
    raise RuntimeError("This certificate requires active Python assertions; do not use -O or -OO.")

import argparse
import hashlib
import json
import sys
from bisect import bisect_right
from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path


F = Fraction
QR = frozenset(i * i % 19 for i in range(1, 19))
BOUND = 1 << 34
SPAN = 1 << 25
SHORT_BOUND = 4999
ASCII_BITS = bytes([48, 49] + list(range(2, 256)))


def rational(x):
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def prime_flags(limit):
    a = bytearray(b"\x01") * (limit + 1)
    a[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            a[p * p::p] = b"\x00" * ((limit - p * p) // p + 1)
    return a


def primes(limit):
    return [p for p, flag in enumerate(prime_flags(limit)) if flag]


def inert(p):
    return p != 19 and p % 19 not in QR


def sigma_power(p, a):
    return (p ** (a + 1) - 1) // (p - 1)


def cube_root_floor(n):
    lo, hi = 0, 1 << ((n.bit_length() + 2) // 3)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** 3 <= n:
            lo = mid
        else:
            hi = mid
    assert lo ** 3 <= n < (lo + 1) ** 3
    return lo


def square_root_fraction(x):
    x = F(x)
    r = isqrt(x.numerator // x.denominator)
    assert r * r * x.denominator <= x.numerator < (r + 1) ** 2 * x.denominator
    return r


def binary_counts(bound):
    counts = [0] * (bound + 1)
    for y in range(-isqrt(4 * bound // 19), isqrt(4 * bound // 19) + 1):
        s = isqrt(4 * bound - 19 * y * y)
        for x in range(-((s + y) // 2), (s - y) // 2 + 1):
            q = x * x + x * y + 5 * y * y
            assert 0 <= q <= bound
            counts[q] += 1
    return counts


def core_flags(lo, hi, odd_inert_primes):
    """Bit values for odd b in [lo,hi], using the proved M19.21 criterion."""
    assert 1 <= lo <= hi and lo % 2 == hi % 2 == 1
    length = (hi - lo) // 2 + 1
    cycle = bytes(int((lo + 2 * i) % 19 in QR) for i in range(19))
    a = bytearray(cycle * (length // 19) + cycle[:length % 19])
    for p in odd_inert_primes:
        if p * p > hi:
            break
        first = (-lo * ((p + 1) // 2)) % p
        if first < length:
            a[first::p] = b"\x00" * ((length - 1 - first) // p + 1)
    return a


def pack_flags(a):
    # Python's exact base-two parser packs the 0/1 bytes in C, with bit 0
    # corresponding to the first flag. Other byte values cause rejection.
    return int(a[::-1].translate(ASCII_BITS), 2) if a else 0


def required_mask(lo, hi):
    assert 2 <= lo <= hi and lo % 2 == hi % 2 == 0
    length = (hi - lo) // 2 + 1
    cycle = bytes(49 if (lo + 2 * i) % 19 else 48 for i in range(19))
    mask = int((cycle * (length // 19) + cycle[:length % 19])[::-1], 2)
    for n in (4, 20):
        if lo <= n <= hi:
            mask &= ~(1 << ((n - lo) // 2))
    return mask


def demand_cover(covered, required):
    assert required & ~covered == 0, "A required target has no certified pair"


def expect_rejected(name, action):
    try:
        action()
    except (AssertionError, ValueError):
        return name
    raise AssertionError("Negative control accepted: " + name)


def controls(odd_inert_primes):
    direct_bound = 1 << 20
    direct = bytearray((direct_bound + 1) // 2)
    for y in range(-isqrt(4 * direct_bound // 19), isqrt(4 * direct_bound // 19) + 1):
        s = isqrt(4 * direct_bound - 19 * y * y)
        for x in range(-((s + y) // 2), (s - y) // 2 + 1):
            q = x * x + x * y + 5 * y * y
            if q and gcd(x, y) == 1:
                assert q % 2 == 1
                direct[(q - 1) // 2] = 1
    c = core_flags(1, direct_bound - 1, odd_inert_primes)
    expected_p = bytearray(c)
    ramified_slots = len(expected_p[9::19])
    expected_p[9::19] = c[:ramified_slots]
    assert direct == expected_p
    direct_core = bytearray(direct)
    direct_core[9::19] = b"\x00" * ramified_slots
    assert direct_core == c
    assert len(c) == direct_bound // 2

    # A separate trial-factor definition checks every small core flag.
    for i in range(5000):
        n = 2 * i + 1
        rem, allowed = n, True
        for p in range(2, isqrt(n) + 1):
            if rem % p == 0:
                if p % 19 not in QR:
                    allowed = False
                while rem % p == 0:
                    rem //= p
        if rem > 1 and rem % 19 not in QR:
            allowed = False
        assert c[i] == allowed

    short = [2 * i + 1 for i, value in enumerate(c[:2500]) if value]
    assert len(short) == 657 and short[-1] <= SHORT_BOUND
    assert [n for n in range(1, 20, 2) if c[(n - 1) // 2]] == [1, 5, 7, 11, 17]
    packed = pack_flags(c[:500])
    covered = 0
    for a in short:
        if a > 999:
            break
        covered |= packed << ((a - 1) // 2)
    for n in range(2, 1001, 2):
        assert bool(covered & (1 << (n // 2 - 1))) == (n % 19 != 0 and n not in (4, 20))
    primitive_packed = pack_flags(direct[:500])
    primitive_covered = 0
    for i, value in enumerate(direct[:500]):
        if value:
            primitive_covered |= primitive_packed << i
    for n in range(2, 1001, 2):
        assert bool(primitive_covered & (1 << (n // 2 - 1))) == (n % 361 != 0 and n not in (4, 76, 380))
    assert required_mask(2, 1000) == sum(1 << (n // 2 - 1) for n in range(2, 1001, 2)
                                             if n % 19 and n not in (4, 20))

    bad_inert_square = bytearray(c)
    bad_inert_square[4] = 1  # 9 has character +1 but is not in C.
    bad_inert_prime = bytearray(c)
    bad_inert_prime[1] = 1
    rejected = [
        expect_rejected("character-only mask admits inert square 9", lambda: assert_equal(bad_inert_square, direct_core)),
        expect_rejected("false summand 3 would create excluded target 4", lambda: assert_equal(bad_inert_prime, direct_core)),
        expect_rejected("omitted target 2", lambda: demand_cover(covered & ~1, required_mask(2, 1000))),
        expect_rejected("excluded target 4 wrongly required", lambda: demand_cover(covered, required_mask(2, 1000) | 2)),
        expect_rejected("nonbinary membership byte", lambda: pack_flags(bytearray([0, 2, 1]))),
    ]
    bc = binary_counts(38)
    shell = [sum(bc[a] * bc[n - a] for a in range(n + 1)) for n in range(39)]
    assert shell[38] == 4
    # Illegally using the 19-unit bound at n=38 gives
    # 80-(32/3)*sqrt(38)>4. Check this strict contradiction exactly.
    assert (3 * (80 - shell[38])) ** 2 > 32 ** 2 * 38
    rejected.append("19-divisible target 38 violates the unramified spectral formula")
    for x in range(-12, 13):
        for y in range(-12, 13):
            q = x * x + x * y + 5 * y * y
            u, v = -x - 10 * y, 2 * x + y
            assert u * u + u * v + 5 * v * v == 19 * q
            if gcd(x, y) == 1 and q % 19:
                assert gcd(u, v) == 1
    rejected.append(expect_rejected("incorrect ramified matrix", lambda: assert_equal(9 * 9 - 9 + 5, 19 * 5)))
    return {
        "coordinate_bound": direct_bound,
        "all_odd_core_and_primitive_flags_compared": direct_bound // 2,
        "coordinate_primitive_flags_sha256": sha(direct),
        "core_flags_sha256": sha(c),
        "trial_factor_flag_count": 5000,
        "direct_small_sumsets_through": 1000,
        "ramified_matrix_grid": [-12, 12],
        "negative_controls": rejected,
    }, short


def assert_equal(a, b):
    assert a == b


def euler_and_lattice():
    split_product, split_rows = F(1), []
    for p in primes(63):
        if p % 19 in QR:
            value, a = max((F((a + 1) ** 6, p ** a), a) for a in range(11))
            split_product *= value
            split_rows.append({"p": p, "a": a, "factor": rational(value)})
    assert 12 ** 6 < 2 * 11 ** 6
    assert split_product == F(2337302235907620864, 35724352769875) < 7 ** 6
    inert_product, inert_rows = F(3 ** 6), []
    for p in primes(63):
        if p > 2 and inert(p):
            value, a = max((F((a + 1) ** 6 * p ** (5 * a), sigma_power(p, a) ** 6), a)
                           for a in range(5))
            inert_product *= value
            inert_rows.append({"p": p, "a": a, "factor": rational(value)})
    assert 6 ** 6 < 3 * 5 ** 6
    assert 9 ** 6 > 2 * 8 ** 6 and 10 ** 6 < 2 * 9 ** 6
    assert inert_product == F(732525035467020016584703898578215969335561207,
                              2777251580551473425917614639931925200896) < F(81, 10) ** 6
    counts = binary_counts(35)
    cumulative, total = [], 0
    for n in range(1, 36):
        total += counts[n]
        cumulative.append(total)
        assert total <= 2 * n
    assert 43 ** 2 < 100 * 19
    tail_constant = F(440, 301) + F(21, 43)
    assert tail_constant == F(587, 301) < 2
    return {
        "M_constant": 7, "M_power": "1/6", "M_sixth_product": rational(split_product),
        "M_factors": split_rows, "M_exponent_cutoff": 10,
        "E_constant": rational(F(81, 10)), "E_power": "-1/3",
        "E_sixth_product": rational(inert_product), "E_factors": inert_rows,
        "E_odd_exponent_cutoff": 4, "E_dyadic_maximum_at": 6,
        "lattice_finite_bound": 35, "lattice_cumulative_counts": cumulative,
        "lattice_tail_from": 36, "lattice_tail_constant": rational(tail_constant),
    }


def divisor_frontier(split_primes, upper):
    prefix, chosen = 1, []
    for p in split_primes:
        prefix *= p
        chosen.append(p)
        if prefix > upper:
            break
    assert prefix > upper
    records = []

    def visit(index, cap, n, value):
        records.append((n, value))
        if index == len(chosen):
            return
        p = chosen[index]
        for a in range(1, cap + 1):
            n *= p
            if n > upper:
                return
            visit(index + 1, a, n, value * (a + 1))

    visit(0, upper.bit_length(), 1, 1)
    records.sort()
    frontier, best = [], 0
    for n, value in records:
        if value > best:
            frontier.append((n, value))
            best = value
    return frontier, len(records), chosen, prefix


def analytic_certificates():
    limit, scale = 1 << 16, 1 << 64
    ps = primes(limit)
    assert [p for p in ps if p <= 1000] == [p for p in range(2, 1001)
                if all(p % d for d in range(2, isqrt(p) + 1))]
    ip = [p for p in ps if p > 2 and inert(p)]
    sp = [p for p in ps if p % 19 in QR]
    residues = [r for r in range(1, 38, 2) if r % 19 not in QR | {0}]
    assert residues == [3, 13, 15, 21, 27, 29, 31, 33, 37]
    assert sum(residues) == 209 and 9 - F(sum(residues), 38) == F(7, 2)
    b0 = sum((F(2, p * (p - 1)) for p in ip if p <= 1024), F(0))
    b0 += 3 * (F(9, 38 * 1024) + F(9, 1024 ** 2))
    assert b0 < F(9, 25)
    suffix = [0] * (len(ip) + 1)
    for i in range(len(ip) - 1, -1, -1):
        p = ip[i]
        ceil = (scale + p * p - 1) // (p * p)
        assert (ceil - 1) * p * p < scale <= ceil * p * p
        suffix[i] = suffix[i + 1] + ceil
    tail = F(9, 38 * limit) + F(9, limit * limit)
    frontier, states, chosen, prefix = divisor_frontier(sp, 1 << 60)
    keys = [n for n, _ in frontier]
    running = 1
    for n in range(1, 10001):
        rem, value = n, 1
        for p in range(2, isqrt(n) + 1):
            a = 0
            while rem % p == 0:
                rem //= p
                a += 1
            if a and p % 19 in QR:
                value *= a + 1
        if rem > 1 and rem % 19 in QR:
            value *= 2
        running = max(running, value)
        assert running == frontier[bisect_right(keys, n) - 1][1]
    intervals, worst, worst_k = [], F(0), None
    for k in range(34, 60):
        low, high = 1 << k, (1 << (k + 1)) - 1
        d = frontier[bisect_right(keys, high) - 1][1]
        root = cube_root_floor(low)
        e = F(81, 10 * root)
        P = square_root_fraction(F(8 * d, 3) / e)
        assert P < limit
        A = bisect_right(ip, P)
        V = F(suffix[A], scale) + tail
        bound = F(9, 25) + 4 * d * V + e * (F(12, 5) + F(3 * A, 2))
        assert bound < F(4, 5) - F(1, 100)
        if bound > worst:
            worst, worst_k = bound, k
        intervals.append({"k": k, "lower": low, "upper": high, "M_upper": d,
                          "cube_root_floor": root, "E_upper": rational(e), "P": P,
                          "inert_prime_count": A, "tail_upper": rational(V),
                          "criterion_upper": rational(bound), "margin": rational(F(4, 5) - bound)})
    assert worst_k == 34 and intervals[0]["P"] == 494 and intervals[0]["inert_prime_count"] == 47
    core_c1 = F(4 * 7 * 9, 38) + F(81, 10) * F(3, 2) * F(9, 38)
    core_c2 = F(4 * 7 * 9) + F(81, 10) * F(153, 20)
    assert core_c1 == F(7227, 760) and core_c2 == F(62793, 200)
    core_tail = F(9, 25) + core_c1 / 32 + core_c2 / (1 << 20)
    assert core_tail < F(4, 5)
    atomic_c1 = F(4 * 7 * 9, 38) + F(100, 9) * F(3, 2) * F(9, 38)
    atomic_c2 = F(4 * 7 * 9) + F(100, 9) * F(153, 20)
    assert atomic_c1 == F(201, 19) and atomic_c2 == 337
    atomic_tail = F(9, 25) + atomic_c1 / 32 + atomic_c2 / (1 << 15)
    assert atomic_tail < F(4, 5)
    return {
        "inert_residue_classes_mod38": residues, "B0_exact_upper": rational(b0),
        "B0_simple_upper": rational(F(9, 25)), "B0_prime_cutoff": 1024,
        "tail_prime_cutoff": limit, "reciprocal_rounding_scale": scale,
        "divisor_frontier_upper": 1 << 60, "divisor_states": states,
        "split_prime_prefix": chosen, "prefix_product": prefix,
        "frontier": [[n, v] for n, v in frontier],
        "divisor_frontier_trial_factor_range": 10000,
        "intervals": intervals, "worst_interval": worst_k, "worst_upper": rational(worst),
        "inert_core_tail_from": 1 << 60, "inert_core_tail_upper": rational(core_tail),
        "exceptional_quotient_tail_from": 1 << 30,
        "exceptional_quotient_tail_upper": rational(atomic_tail),
        "all_inert_core_onset": BOUND,
    }


def finite_coverage(odd_inert_primes, short):
    segments, previous, total = [], 0, 0
    max_used = 0
    for index, lo in enumerate(range(2, BOUND + 1, SPAN)):
        hi = min(lo + SPAN - 2, BOUND)
        assert lo == previous + 2 and hi % 2 == 0
        b0 = max(1, lo - short[-1])
        c = core_flags(b0, hi - 1, odd_inert_primes)
        membership_hash = sha(c)
        bits = pack_flags(c)
        del c
        required = required_mask(lo, hi)
        covered = 0
        for used, a in enumerate(short, 1):
            shift = (lo - a - b0) // 2
            assert 2 * shift == lo - a - b0
            covered |= bits >> shift if shift >= 0 else bits << (-shift)
            if used % 16 == 0 or used == len(short):
                if required & ~covered == 0:
                    break
        demand_cover(covered, required)
        count = (hi - lo) // 2 + 1
        covered &= (1 << count) - 1
        required_count = required.bit_count()
        total += required_count
        max_used = max(max_used, used)
        size = (count + 7) // 8
        segments.append({"index": index, "lower": lo, "upper": hi,
                         "summand_flags_lower": b0, "summand_flags_upper": hi - 1,
                         "required_targets": required_count, "short_summands_used": used,
                         "membership_sha256": membership_hash,
                         "target_mask_sha256": sha(required.to_bytes(size, "little")),
                         "coverage_sha256": sha(covered.to_bytes(size, "little")),
                         "missing_targets": 0})
        previous = hi
        if (index + 1) % 16 == 0:
            print(f"coverage {index + 1}/512; checked through {hi}", file=sys.stderr, flush=True)
    assert previous == BOUND and len(segments) == 512
    expected = BOUND // 2 - BOUND // 38 - 2
    assert total == expected
    return {"lower": 2, "upper": BOUND, "target_span_per_segment": SPAN,
            "segment_count": len(segments), "eligible_target_count": total,
            "excluded_core_targets": [4, 20], "missing_targets": 0,
            "short_summand_bound": SHORT_BOUND, "short_summands": short,
            "short_summand_count": len(short), "maximum_summands_used": max_used,
            "membership_method": "M19.21: character +1 and deletion of inert primes through sqrt(segment upper)",
            "segments": segments}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ps = primes(isqrt(BOUND))
    ip = [p for p in ps if p > 2 and inert(p)]
    control_result, short = controls(ip)
    euler = euler_and_lattice()
    analytic = analytic_certificates()
    print("local, Euler and all analytic interval checks passed; starting complete finite coverage", file=sys.stderr, flush=True)
    coverage = finite_coverage(ip, short)
    result = {
        "result": "PASS", "coefficient_domain": "exact integers and rational numbers",
        "checker_sha256": sha(Path(__file__).read_bytes()),
        "proof": "PRIMITIVE_NORM_MINUS19_EXACT.md",
        "core_identity": "C+C={n>=2: n even, 19 does not divide n} minus {4,20}",
        "primitive_identity": "P+P={N>=2: N even, 361 does not divide N} minus {4,76,380}",
        "local_and_negative_controls": control_result,
        "euler_and_lattice": euler, "analytic_certificates": analytic,
        "finite_coverage": coverage,
        "scope": "Exact finite certificate. The written proof supplies classical spectral inputs, local interfaces, infinite Euler completeness and exceptional-quotient descent. No formal verification or external review.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"result": "PASS", "eligible_finite_targets": coverage["eligible_target_count"],
                      "segments": coverage["segment_count"], "analytic_intervals": 26,
                      "negative_controls": len(control_result["negative_controls"]),
                      "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
