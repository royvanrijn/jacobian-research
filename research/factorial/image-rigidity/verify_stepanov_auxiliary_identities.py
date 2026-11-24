#!/usr/bin/env python3
"""Exact finite guards for STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md.

No unrestricted Strong Factorial claim follows from these finite checks.
"""

import argparse
import json
from math import factorial
from pathlib import Path

import sympy as sp

if not __debug__:
    raise RuntimeError("Assertions must execute; do not run with -O or -OO")

x = sp.Symbol("x")


def same(left, right):
    assert sp.expand(left - right) == 0


def wronskian(polys):
    return sp.det(sp.Matrix([
        [sp.diff(poly, x, i) for poly in polys]
        for i in range(len(polys))
    ]))


def ordinary(m, a, d):
    return sp.Add(*[
        sp.binomial(m, j) * sp.rf(a * m + 1, d * j) * x ** (m - j)
        for j in range(m + 1)
    ])


def marked(poly, n, gamma, d, ell):
    result = poly
    for j in range(ell):
        result = sp.expand((gamma + d * n + j) * result
                           - d * x * sp.diff(result, x))
    return result


def shift_check(p, q):
    m, n = sp.degree(p, x), sp.degree(q, x)
    r = int(m + n)
    polys = [x ** j * p for j in range(n)]
    polys += [x ** j * q for j in range(m)]
    coefficients = sp.Matrix([
        [sp.expand(poly).coeff(x, i) for poly in polys]
        for i in range(r)
    ])
    determinant = coefficients.det()
    resultant = sp.resultant(p, q, x)
    same(determinant, (-1) ** (m * n) * resultant)
    same(wronskian(polys), sp.prod(factorial(i) for i in range(r))
         * determinant)
    # A Taylor jet determinant at a specified point is the same constant.
    jets = sp.Matrix([
        [sp.diff(poly, x, i).subs(x, -2) / factorial(i) for poly in polys]
        for i in range(r)
    ])
    same(jets.det(), determinant)
    if not (p.free_symbols | q.free_symbols) - {x}:
        e = sp.degree(sp.gcd(p, q), x)
        assert coefficients.rank() == r - e
    return resultant


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name(
        'stepanov_auxiliary_identities_verification.json'))
    args = parser.parse_args()
    records = {}

    # Lemma 1: finite guards, independent of the general written CRT proof.
    nodes = [(-1, 2), (2, 3), (5, 1)]
    ranks = []
    for degree in range(9):
        jets = sp.Matrix([
            [sp.binomial(j, k) * sp.Integer(alpha) ** (j - k)
             if j >= k else 0 for j in range(degree + 1)]
            for alpha, multiplicity in nodes for k in range(multiplicity)
        ])
        rank = jets.rank()
        assert rank == min(degree + 1, 6)
        ranks.append(rank)
    records["jet_rank_guard"] = {"nodes_and_multiplicities": nodes,
                                  "degrees": list(range(9)), "ranks": ranks}

    # Lemma 3: fully symbolic small-degree guards also check the sign.
    symbolic_degrees = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 1)]
    for m, n in symbolic_degrees:
        pcs = sp.symbols(f"p0:{m + 1}")
        qcs = sp.symbols(f"q0:{n + 1}")
        p = sum(pcs[j] * x ** j for j in range(m + 1))
        q = sum(qcs[j] * x ** j for j in range(n + 1))
        shift_check(p, q)
    records["symbolic_shift_determinant_degrees"] = symbolic_degrees

    # Lemma 4: generic symbolic quadratics/linear factors, then the shared-root
    # negative control. Both the generic identity and failure at a root matter.
    u, v, w = sp.symbols("u v w")
    p, q = x ** 2 + u * x + v, x + w
    for exponent in (1, 2, 3):
        polys = [p ** (exponent - j) * q ** j
                 for j in range(exponent + 1)]
        expected = (sp.prod(factorial(j) for j in range(exponent + 1))
                    * (p * sp.diff(q, x) - sp.diff(p, x) * q)
                    ** (exponent * (exponent + 1) // 2))
        same(wronskian(polys), expected)
    records["symbolic_product_wronskian_exponents"] = [1, 2, 3]

    p, q = (x + 1) * (x + 2), x + 1
    wron = sp.expand(wronskian([p, q]))
    same(wron, -(x + 1) ** 2)
    assert wron != 0 and wron.subs(x, -1) == 0
    assert p.subs(x, -1) == q.subs(x, -1) == 0
    assert shift_check(p, q) == 0
    records["negative_control"] = {
        "P": str(sp.expand(p)), "Q": str(q), "common_root": -1,
        "nonzero_wronskian": str(wron), "resultant": 0,
    }

    # Lemma 5: a coprime nonsingular certificate and a shared-root example
    # whose singular factor has exactly the necessary root.
    p, q = x + 1, x + 2
    same(sp.diff(p, x), q - p)
    same(sp.diff(q, x), q - p)
    assert sp.gcd(p, q) == 1
    p, q, singular = (x + 1) * (x + 2), x + 1, x + 1
    same(singular * sp.diff(p, x), 2 * p - q)
    same(singular * sp.diff(q, x), q)
    assert sp.gcd(p, q) == singular
    records["derivative_closure_controls"] = {
        "nonsingular": "P=x+1,Q=x+2,s=1; P'=Q-P,Q'=Q-P",
        "singular": "P=(x+1)(x+2),Q=x+1,s=x+1; sP'=2P-Q,sQ'=Q",
    }

    # Actual shape relations, including a=0,1 and currently nontrivial lower
    # exponents. This is a finite algebraic guard, not a claim of new coverage.
    actual_records = []
    for a, d in [(0, 1), (0, 6), (1, 6), (2, 3), (3, 2), (4, 3), (5, 6)]:
        for terminal in (2, 3, 4, 5):
            n = terminal - 1
            gamma = a * n + 1
            p = ordinary(terminal, a, d)
            q = ordinary(n, a, d)
            ta = marked(q, n, gamma, d, a)
            tad = marked(q, n, gamma, d, a + d)
            same(sp.rf(gamma, a) * p, x * ta + tad)
            same(sp.rf(gamma, a) * sp.diff(p, x), terminal * ta)
            same(sp.LC(sp.Poly(p * sp.diff(q, x) - sp.diff(p, x) * q, x)), -1)
            assert sp.degree(p * sp.diff(q, x) - sp.diff(p, x) * q, x) == 2 * terminal - 2
            actual_records.append({"a": a, "d": d, "N": terminal})
    records["actual_shape_identity_cases"] = actual_records

    # Degree-preserving actual shift examples, independent from the operator
    # calculations. Keeping these small avoids an unnecessary large search.
    actual_resultants = []
    for a, d, terminal in [(0, 2, 2), (1, 6, 3), (3, 2, 3), (4, 3, 4)]:
        resultant = shift_check(ordinary(terminal, a, d),
                               ordinary(terminal - 1, a, d))
        actual_resultants.append({"a": a, "d": d, "N": terminal,
                                  "resultant": str(resultant)})
    records["actual_shift_determinant_cases"] = actual_resultants

    receipt = {"status": "PASS", "sympy_version": sp.__version__,
                      "scope": "finite exact guards for written general algebraic proofs; no unrestricted SFC conclusion",
                      "checks": records}
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({"status": "PASS", "output": str(args.output),
                      "check_groups": list(records)}))


if __name__ == "__main__":
    main()
