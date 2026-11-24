#!/usr/bin/env python3
"""Exact finite identity controls for SFC_CYCLIC_MONOMIAL_BLOCKS.md.

The all-order assertions are proved in the note, not by the finite cases.
Only Python's standard library is required. Failed checks raise explicitly,
including under python -O.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import factorial, gcd, lcm, comb
from pathlib import Path
import json


EXPECTED_NOTE_SHA256 = "43ccca4e636c8cf74b7045a33bad4a410eac258a62a07d00bbb7feb3888f50f6"


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_div_exact(a, b):
    a, b = trim(a), trim(b)
    require(b[-1] == 1, "monic polynomial divisor")
    out = [0] * (len(a) - len(b) + 1)
    while len(a) >= len(b):
        k, c = len(a) - len(b), a[-1]
        out[k] = c
        for j, bj in enumerate(b):
            a[k+j] -= c*bj
        a = trim(a)
    require(all(x == 0 for x in a), "exact cyclotomic polynomial division")
    return trim(out)


@lru_cache(None)
def cyclotomic(n):
    p = [-1] + [0] * (n-1) + [1]
    for d in range(1, n):
        if n % d == 0:
            p = poly_div_exact(p, cyclotomic(d))
    return tuple(p)


class CyclotomicRing:
    def __init__(self, order):
        self.order = order
        self.mod = cyclotomic(order)
        self.degree = len(self.mod)-1
        self.zero = (F(0),) * self.degree
        self.one = (F(1),) + (F(0),) * (self.degree-1)
        self.powers = [self.reduce([0]*i+[1]) for i in range(order)]

    def reduce(self, p):
        p = list(map(F, p))
        while len(p) > self.degree:
            k, c = len(p)-1-self.degree, p[-1]
            if c:
                for j, v in enumerate(self.mod):
                    p[k+j] -= c*v
            p.pop()
        return tuple(p+[F(0)]*(self.degree-len(p)))

    def add(self, a, b):
        return tuple(x+y for x, y in zip(a, b))

    def scale(self, a, c):
        return tuple(c*x for x in a)

    def mul(self, a, b):
        out = [F(0)] * (2*self.degree-1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    out[i+j] += x*y
        return self.reduce(out)


def monomial_series(profile, bound):
    require(profile and all(d >= 0 for d in profile) and any(profile),
            "nonconstant monomial profile")
    out = []
    for n in range(bound+1):
        numerator = 1
        for d in profile:
            numerator *= factorial(d*n)
        out.append(F(numerator, factorial(n)))
    return out


def scalar_mul(a, b, bound):
    out = [F(0)] * (bound+1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:bound+1-i]):
            out[i+j] += x*y
    return out


def kaluza_controls(a):
    bound = len(a)-1
    ratios = [a[i+1]/a[i] for i in range(bound)]
    require(all(ratios[i] <= ratios[i+1] for i in range(bound-1)),
            "log-convex monomial coefficient sequence")
    u = [F(0)] * (bound+1)
    b = [F(0)] * (bound+1)
    for n in range(1, bound+1):
        u[n] = a[n] - sum(u[i]*a[n-i] for i in range(1, n))
        if n >= 2:
            other = sum(u[i]*a[n-1-i]*(ratios[n-1]-ratios[n-1-i])
                        for i in range(1, n))
            require(u[n] == other, "formal Kaluza ratio identity")
        require(u[n] >= 0, "nonnegative reciprocal complement")
        b[n] = a[n] - sum(F(i, n)*b[i]*a[n-i] for i in range(1, n))
        require(b[n] >= a[1]**n/n, "strict logarithm coefficient bound")
    # Independently reconstruct log(phi) from powers of U.
    power = [F(1)]+[F(0)]*bound
    by_powers = [F(0)]*(bound+1)
    for k in range(1, bound+1):
        power = scalar_mul(power, u, bound)
        by_powers = [x+y/k for x, y in zip(by_powers, power)]
    require(b == by_powers, "two independent formal logarithm identities")
    return b


def ring_series_mul(a, b, ring, bound):
    out = [ring.zero]*(bound+1)
    for i, x in enumerate(a):
        if x != ring.zero:
            for j, y in enumerate(b[:bound+1-i]):
                if y != ring.zero:
                    out[i+j] = ring.add(out[i+j], ring.mul(x, y))
    return out


def semigroup_members(generators, bound):
    reachable = [True]+[False]*bound
    for n in range(1, bound+1):
        reachable[n] = any(n >= q and reachable[n-q] for q in generators)
    return reachable


def check_cycles(specs, bound):
    order = lcm(*(q for _, q, _ in specs))
    ring = CyclotomicRing(order)
    direct = [ring.one]+[ring.zero]*bound
    logarithm = [F(0)]*(bound+1)
    for profile, q, rho in specs:
        a = monomial_series(profile, bound)
        b = kaluza_controls(a)
        for n in range(q, bound+1, q):
            logarithm[n] += q*b[n]*rho**n
        for j in range(q):
            factor = [ring.scale(ring.powers[(order//q*j*n) % order],
                                 a[n]*rho**n) for n in range(bound+1)]
            direct = ring_series_mul(direct, factor, ring, bound)
    # exp(logarithm) from E'=E*logarithm'.
    filtered = [F(1)]+[F(0)]*bound
    for n in range(1, bound+1):
        filtered[n] = sum(F(k, n)*logarithm[k]*filtered[n-k]
                          for k in range(1, n+1))
    members = semigroup_members([q for _, q, _ in specs], bound)
    for n in range(bound+1):
        require(direct[n] == ring.scale(ring.one, filtered[n]),
                "direct cyclotomic product equals filtered exponential")
        require((filtered[n] > 0) == members[n], "exact semigroup zero pattern")
    if len(specs) == 1:
        profile, q, rho = specs[0]
        first = monomial_series(profile, 1)[1]
        for n in range(q, bound+1, q):
            require(filtered[n] >= (first*rho)**n/(n//q),
                    "explicit single-cycle nonzero coefficient bound")
    return filtered


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    note = Path(__file__).with_name("SFC_CYCLIC_MONOMIAL_BLOCKS.md")
    note_hash = sha256(note.read_bytes()).hexdigest()
    require(note_hash == EXPECTED_NOTE_SHA256, "proof note hash mismatch")
    bound = 18
    profiles = [(1,), (2,), (1, 1), (2, 1), (3, 2), (0, 2, 1)]
    cases = []
    coefficient_count = 0
    sequences = []
    for profile in profiles:
        for q in range(2, 7):
            specs = [(profile, q, F(1))]
            values = check_cycles(specs, bound)
            sequences.append((q, values))
            cases.append({"profiles": [list(profile)], "cycles": [q],
                          "scales": ["1"], "orders": [0, bound]})
            coefficient_count += bound+1
    for specs in [
        [((2,), 2, F(2)), ((1, 1), 3, F(3))],
        [((1,), 2, F(1, 2)), ((3, 2), 4, F(2, 3))],
        [((2, 1), 3, F(1)), ((3,), 5, F(2))],
        [((1, 1), 1, F(1)), ((2,), 3, F(1, 3))],
    ]:
        values = check_cycles(specs, bound)
        cases.append({"profiles": [list(x[0]) for x in specs],
                      "cycles": [x[1] for x in specs],
                      "scales": [str(x[2]) for x in specs],
                      "orders": [0, bound],
                      "zero_orders": [n for n, v in enumerate(values) if not v]})
        coefficient_count += bound+1
    power_checks = 0
    for q, values in sequences:
        for m in range(1, 6):
            period = q//gcd(q, m)
            require(period <= comb(q+m-1, m), "support count dominates window")
            for n in range(bound//m+1):
                require(bool(values[m*n]) == (n % period == 0),
                        "exact power zero pattern")
                power_checks += 1
    product_checks = 0
    for q, a in sequences[:5]:
        for r, b in sequences[5:10]:
            period = lcm(q, r)
            require(period <= q*r, "product support dominates guaranteed period")
            for n in range(bound+1):
                require(bool(a[n]*b[n]) == (n % period == 0),
                        "independent product zero pattern")
                product_checks += 1
    exponential = [F(1, factorial(n)) for n in range(bound+1)]
    negative_argument = [(-1)**n*x for n, x in enumerate(exponential)]
    control = scalar_mul(exponential, negative_argument, bound)
    require(control == [F(1)]+[F(0)]*bound, "positive-series shortcut control")
    require(exponential[1]**2 > exponential[0]*exponential[2],
            "countercontrol violates log convexity")
    receipt = {
        "status": "PASS",
        "proof_source_sha256": note_hash,
        "scope": "Finite exact identity controls for the all-order cyclic-block theorem.",
        "arithmetic": "fractions and exact cyclotomic quotient rings; no floating point",
        "direct_cyclotomic_cases": cases,
        "coefficient_identities": coefficient_count,
        "power_zero_controls": power_checks,
        "independent_product_controls": product_checks,
        "positive_series_without_log_convexity_shortcut_refuted": True,
        "general_log_convex_interlacing_proved": False,
        "arbitrary_coefficients_or_overlapping_monomials_solved": False,
        "full_SFC_solved": False,
    }
    rendered = json.dumps(receipt, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == "__main__":
    main()
