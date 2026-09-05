#!/usr/bin/env python3
"""Exact support-gap and sharp-family replay for GVC2OC.

The universal two-envelope lemma is a written proof. This script checks its
algebra, a complete small support regression, high operator orders, and the
last nonzero contractions of the sharp family. No moment prefix proves GVC.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts/generated-results/binary-gvc-optimal-cutoff-v1.json"


def constant(d):
    return 1 + (d+1)**2//4


def best_gap(operator, polynomial):
    """Maximize on the breakpoints, retaining a positive endpoint proxy."""
    d = max(map(sum, polynomial))
    slopes = {Fraction(1), Fraction(1, d+1)}
    for support in (operator, polynomial):
        for (a, b), (c, e) in combinations(support, 2):
            if b != e:
                z = Fraction(c-a, b-e)
                if 0 < z < 1:
                    slopes.add(z)

    def gap(z):
        return min(a+b*z for a, b in operator)-max(i+j*z for i, j in polynomial)

    z = max(sorted(slopes), key=gap)
    value = gap(z)
    if value <= 0:
        return None
    assert value >= Fraction(1, constant(d))
    u, v = z.denominator, z.numerator
    integral_gap = min(u*a+v*b for a, b in operator)-max(u*i+v*j for i, j in polynomial)
    assert integral_gap > 0 and Fraction(max(u, v), integral_gap) <= constant(d)
    return {"weight": [u, v], "gap": integral_gap, "normalized_gap": str(value)}


def mixed_terms(e, j, order, q, m):
    """Binomial differentiation of x^q*(x^(e-1)*y^j)^m."""
    x_degree, y_degree = (e-1)*m+q, j*m
    terms = []
    for k in range(m+1):
        dx, dy = e*k, order*(m-k)
        if dx <= x_degree and dy <= y_degree:
            coefficient = (comb(m, k)*factorial(x_degree)//factorial(x_degree-dx)
                           *factorial(y_degree)//factorial(y_degree-dy))
            terms.append((x_degree-dx, y_degree-dy, coefficient))
    return terms


def support_certificate(operator, polynomial):
    """Cover every positive weight by checking both normalized orientations."""
    certificate = best_gap(operator, polynomial)
    if certificate is not None:
        return certificate
    certificate = best_gap([(b, a) for a, b in operator],
                           [(j, i) for i, j in polynomial])
    if certificate is not None:
        certificate["weight"].reverse()
    return certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    A, B, C, E, N = sp.symbols("A B C E N")
    # The normalized inverse gap at an operator breakpoint.
    assert sp.cancel((B+E)/(A*E-B*C) -
                     (1/A+B/N+B*C/(A*N)).subs(N, A*E-B*C)) == 0
    e, j, q = sp.symbols("e j q", integer=True, nonnegative=True)
    order = e*j+1
    assert sp.expand(e*order-((e-1)*order+j*e)) == 1
    assert sp.expand(e*(order-j)*q-((e-1)*order*q+q)) == 0
    assert sp.expand(order*(order*q-(order-j)*q)-j*order*q) == 0

    monomials = [(i, j) for i in range(4) for j in range(4-i)]
    positive = [m for m in monomials if sum(m)]
    operators = [(m,) for m in positive] + list(combinations(positive, 2))
    polynomials = [(m,) for m in monomials] + list(combinations(monomials, 2))
    count = feasible = 0
    for operator in operators:
        for polynomial in polynomials:
            count += 1
            feasible += support_certificate(operator, polynomial) is not None
    assert count == 2475

    families = []
    for d in range(9):
        e = (d+2)//2
        j = d+1-e
        order = e*j+1
        assert order == constant(d)
        certificate = best_gap([(e, 0), (0, order)], [(e-1, j)])
        assert certificate is not None
        assert Fraction(max(certificate["weight"]), certificate["gap"]) == order
        for q in (1, 2):
            m = order*q
            last = mixed_terms(e, j, order, q, m)
            expected = comb(m, j*q)*factorial((e-1)*m+q)*factorial(j*m)
            assert last == [(0, 0, expected)]
            assert not mixed_terms(e, j, order, q, m+1)
            assert all(not mixed_terms(e, j, order, 0, k) for k in range(1, 5))
            families.append({"degree": d, "e": e, "j": j, "operator_degree": order,
                             "multiplier_degree": q, "last_nonzero_power": m,
                             "last_value": str(expected), "certificate": certificate})

    high_orders = []
    for order in (17, 101, 1000000007):
        certificate = best_gap([(2, 0), (0, order)], [(1, 2)])
        assert certificate is not None
        high_orders.append({"operator_degree": order, "polynomial_degree": 3,
                            "certificate": certificate})
    for d in range(1, 101):
        assert max(B*(C+1) for B in range(d+1) for C in range(d+1-B)) == (d+1)**2//4
        assert d+1 <= constant(d)

    data = {
        "format": "binary-gvc-optimal-cutoff-v1", "theorem": "GVC2OC",
        "source_sha256": {Path(__file__).resolve().relative_to(ROOT).as_posix():
                          hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "uniform_coefficient": "C_d=1+floor((d+1)^2/4)",
        "breakpoint_identity": "(B+E)/(A*E-B*C)=1/A+B/N+B*C/(A*N), N=A*E-B*C",
        "small_support_pairs": count, "feasible_pairs": feasible,
        "sharp_family": families, "high_operator_order_controls": high_orders,
        "boundary": "The rational identities and finite regressions are exact; the universal envelope maximum argument and support-certificate necessity remain written proofs. No literature-priority or full formal-verification claim.",
    }
    serialized = json.dumps(data, indent=2)+"\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.read_text() == serialized, "artifact differs; inspect before --write"
    print(f"PASS {count} support pairs, {feasible} feasible, with the uniform gap bound")
    print("PASS symbolic breakpoint and all-degree sharp-family identities")
    print("PASS 18 exact last-nonzero mixed contractions and high-order controls")
    print("PASS " + ("wrote " if args.write else "byte-identical ") + str(OUTPUT))


if __name__ == "__main__":
    main()
