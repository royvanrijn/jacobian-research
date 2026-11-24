#!/usr/bin/env python3
"""Rebuild the -163 bridge thresholds, replay all witnesses and analyze gap structure."""

import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
from bisect import bisect_right
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import comb, isqrt, prod
from pathlib import Path

from verify_heegner_core_bounds import (
    PRIME_LIMIT,
    ROUNDING,
    ceil_scaled,
    divisor_frontier,
    first_minimum,
    prime_flags,
    sieve_envelopes,
)

if not __debug__:
    raise RuntimeError("Run the completion checker with Python assertions active.")

SOURCE = Path(__file__).resolve().parent
FIRST_BIN, LAST_BIN = 34, 42
PRIME_CAP = 2**25
SCALE = 2**20
TRACKED = (3, 5, 7, 11, 13, 17)


def split_at(p, q):
    return pow(q, (p - 1)//2, p) == 1


def prime_list(bound):
    return [q for q, flag in enumerate(prime_flags(bound)) if flag]


def threshold_input(primes):
    p = 163
    inert = [q for q in primes if q > 2 and pow(q, 81, p) == p - 1]
    split = [q for q in primes if split_at(p, q)]
    assert tuple(inert[:len(TRACKED)]) == TRACKED
    assert (14*Fraction(37, 100))**2*p > 66**2
    assert ((Fraction(7, 2)-3)/6)**2*p > 1
    w = Fraction(1, 27)
    nu = 2*w/(p+1)
    mass, tails = sieve_envelopes(p, inert)
    top = 2**(LAST_BIN+1)-1
    frontier, _ = divisor_frontier(split, top)
    indices = [u for u, _ in frontier]
    odd_frontier, _ = divisor_frontier(inert, top//(2*PRIME_CAP))
    odd_indices = [u for u, _ in odd_frontier]

    harmonic = [0]*(len(inert)+1)
    for i in range(len(inert)-1, -1, -1):
        harmonic[i] = harmonic[i+1] + ceil_scaled(1, inert[i])
    harmonic = [Fraction(x, ROUNDING) for x in harmonic]
    profiles = []
    for vals in product(range(3), repeat=len(TRACKED)):
        b, offset, shell = mass, 1-w, Fraction(1)
        for q, v in zip(TRACKED, vals, strict=True):
            shell *= sum(Fraction(1, q**i) for i in range(v+1))
            if v:
                b -= Fraction(2, q*(q-1))
                offset -= 1
            if v == 2:
                b += Fraction(1, q*q)
                offset += Fraction(1, q)
        profiles.append((b, offset, shell))

    table = [f"{FIRST_BIN} {LAST_BIN} {SCALE} {PRIME_CAP} {len(profiles)}"]
    minima = {}
    classes = (p-1)//2
    for k in range(FIRST_BIN, LAST_BIN+1):
        lower, upper = 2**k, 2**(k+1)-1
        root = isqrt(lower)
        _, divisor_bound = frontier[bisect_right(indices, upper)-1]
        remote_harmonic = classes*(
            Fraction(1, PRIME_LIMIT) + Fraction(7*max(0, (k+2)//2-16), 20*p)
        )
        tail_scores = [
            2*divisor_bound/(3*w)*(
                Fraction(37, 100)*t
                + (Fraction(7, 2)*(h+remote_harmonic)+Fraction(1, 2))/root
                + Fraction(2*classes, lower)
            )
            for t, h in zip(tails, harmonic, strict=True)
        ]
        thresholds = []
        for b, offset, shell in profiles:
            row = []
            for square in (False, True):
                limit = 1-(Fraction(4, p+1) if square else 0)
                error = offset+(1-2*nu if square else 0)

                def bound(i):
                    return (limit-b-tail_scores[i]/shell)/(error+i)

                i = first_minimum(len(TRACKED), len(inert), lambda j: -bound(j))
                # Floor the exact sufficient threshold to retain every possible failure.
                threshold = (3*root*w*bound(i)*SCALE)//1
                assert threshold > 0
                row.append(threshold)
            thresholds.append(row)
            table.append(f"{row[0]} {row[1]}")
        min_threshold = min(map(min, thresholds))
        minima[k] = min_threshold
        # A prime exceeding PRIME_CAP occurs only once. Removing it leaves
        # a tiny odd inert cofactor, whose divisor count has an exact envelope.
        assert PRIME_CAP**2 > upper
        for t in range(1, k+1):
            cofactor_cap = upper//(2**t*(PRIME_CAP+1))
            if cofactor_cap:
                _, d = odd_frontier[bisect_right(odd_indices, cofactor_cap)-1]
                assert 2*(t+3)*d*SCALE < min_threshold, (k, t, d)
    return "\n".join(table)+"\n", minima


def uniform_seven_bound(primes):
    inert = [q for q in primes if pow(q, 3, 7) == 6]
    split = [q for q in primes if split_at(7, q)]
    mass, _ = sieve_envelopes(7, inert)
    assert mass < Fraction(47, 100)
    # The same exponent truncation covers q=2 in the −7 branch.
    assert 12**6 < 2*11**6
    divisor = prod(
        max(Fraction((a+1)**6, q**a) for a in range(11)) for q in split if q < 64
    )
    assert divisor < 8**6
    shell_error = prod(
        max(Fraction((a+1)**6*q**(5*a), ((q**(a+1)-1)//(q-1))**6) for a in range(11))
        for q in primes if q < 64
    )
    assert shell_error < 20**6
    bound = (
        Fraction(47, 100) + 6*8*(Fraction(3, 14*2**10)+Fraction(3, 2**40))
        + 20*(Fraction(3, 14*2**10)+Fraction(7, 2*2**40))
    )
    assert bound < Fraction(1, 2)
    return {"onset": 2**120, "sieve_upper_bound": str(bound), "inert_mass_bound": str(mass)}


def closure_baseline(fibers, required, forbidden, gap_count):
    counts = {(0, 0): 1}
    for root, children in fibers.items():
        forced = sum(child in required for child in children)
        if root in forbidden:
            assert not forced
            continue
        optional = len(children)-forced-sum(child in forbidden for child in children)
        states = [(1+forced+i, 1, comb(optional, i)) for i in range(optional+1)]
        if root not in required and not forced:
            states.append((0, 0, 1))
        updated = defaultdict(int)
        for (size, roots), weight in counts.items():
            for added, root_added, ways in states:
                if size+added <= gap_count:
                    updated[size+added, roots+root_added] += weight*ways
        counts = updated
    return {r: weight for (n, r), weight in counts.items() if n == gap_count}


def mean(distribution):
    return Fraction(sum(k*v for k, v in distribution.items()), sum(distribution.values()))


def analyze_gaps(scans):
    rows = []
    factor_witnesses = []
    total_distribution = {0: 1}
    total_naive_mean = Fraction(0)
    total_observed = 0
    for row in scans:
        p = row["p"]
        gaps = set(row["core_exceptions"])
        largest = max(gaps)
        primes = prime_list(largest)
        split = [q for q in primes if split_at(p, q)]
        blocked = [q for q in primes if not split_at(p, q)]

        def in_core(n):
            return all(n % q for q in blocked)

        local = [n for n in range(2, largest+1, 2) if n % p]
        fibers = defaultdict(list)
        for n in local:
            factors = [q for q in split if n % q == 0]
            assert len(factors) <= 1
            if not factors:
                fibers[n]
                continue
            q, = factors
            assert n % (q*q)
            root = n//q
            fibers[root].append(n)
            if n in gaps:
                assert root in gaps
                summands = [a for a in range(1, n//2+1, 2) if in_core(a)]
                blockers = [next(r for r in blocked if (n-a) % r == 0) for a in summands]
                factor_witnesses.append(dict(p=p, gap=n, inert_part=root, split_prime=q,
                                             smaller_summands=summands,
                                             complementary_forbidden_divisors=blockers))

        # Below twice the first nonunit core value, the only option is 1+(n-1).
        trivial = [n for n in local if n < 2*((p+1)//4)]
        required = {n for n in trivial if not in_core(n-1)} | {largest}
        forbidden = {n for n in trivial if in_core(n-1)}
        assert required <= gaps and not (gaps & forbidden)
        distribution = closure_baseline(fibers, required, forbidden, len(gaps))
        observed = len(gaps & fibers.keys())
        assert distribution.get(observed, 0) > 0

        joined = defaultdict(int)
        for j, a in total_distribution.items():
            for k, b in distribution.items():
                joined[j+k] += a*b
        total_distribution = joined
        total_naive_mean += Fraction(len(gaps)*len(fibers), len(local))
        total_observed += observed
        rows.append(dict(p=p, largest=largest, gaps=len(gaps), admitted=len(local),
                         inert_supported=len(fibers), observed_inert_gaps=observed,
                         closure_baseline_mean=str(mean(distribution))))

    assert len(factor_witnesses) == 10 and total_observed == 90
    tail = Fraction(sum(v for k, v in total_distribution.items() if k >= total_observed),
                    sum(total_distribution.values()))
    return dict(rows=rows, split_bearing_failures=factor_witnesses,
                naive_baseline_mean=str(total_naive_mean),
                closure_baseline_mean=str(mean(total_distribution)),
                closure_baseline_probability_at_least_observed=str(tail))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    args = parser.parse_args()
    primes = prime_list(PRIME_LIMIT)
    table, minima = threshold_input(primes)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp)/"thresholds.txt"
        path.write_text(table, encoding="ascii")
        run = subprocess.run([args.native.resolve(), path], check=True, capture_output=True, text=True)
    coverage = [json.loads(line) for line in run.stdout.splitlines()]
    assert [row["k"] for row in coverage] == list(range(FIRST_BIN, LAST_BIN+1))
    assert all(row["first_failure"] == 0 and row["targets"] > 0 for row in coverage)
    names = (
        "verify_heegner_core_completion.py",
        "verify_minus163_core_bridge.cpp",
        "verify_heegner_core_bounds.py",
        "verify_primitive_norm_minus11.py",
        "heegner_core_sumset_checks.json",
        "sources/jaeschke-strong-pseudoprimes-1993.pdf",
    )
    scans = json.loads((SOURCE/names[4]).read_text())["scans"]
    out = dict(python=platform.python_version(),
               inputs_sha256={name: sha256(SOURCE/name) for name in names},
               thresholds_sha256=hashlib.sha256(table.encode("ascii")).hexdigest(),
               prime_cap=PRIME_CAP, scale=SCALE, threshold_minima=minima, bins=coverage,
               certified_witness_targets=sum(row["targets"] for row in coverage),
               common_tail_minus7=uniform_seven_bound(primes),
               gap_structure=analyze_gaps(scans))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
