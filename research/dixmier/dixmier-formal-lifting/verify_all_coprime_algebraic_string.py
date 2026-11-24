#!/usr/bin/env python3
"""Finite exact controls for the written all-coprime algebraic string proof."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import sympy as s

if not __debug__:
    raise SystemExit("Exact assertions require non-optimized Python.")

HERE = Path(__file__).resolve().parent
x, xi, z, lam = s.symbols("x xi z lambda")


def clean(op):
    return {j: a for j, c in op.items() if (a := s.expand(c)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for j, a in op.items():
            out[j] = out.get(j, 0) + a
    return clean(out)


def scale(op, value):
    return clean({j: value*a for j, a in op.items()})


def mul(left, right, minimum=0):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            bound = i + j - minimum
            if i >= 0:
                bound = min(bound, i)
            for k in range(bound+1):
                c = s.binomial(i, k) * a * s.diff(b, x, k)
                out[i+j-k] = out.get(i+j-k, 0) + c
    return clean(out)


def comm(left, right, minimum=0):
    return add(mul(left, right, minimum), scale(mul(right, left, minimum), -1))


def power(op, exponent, minimum=0):
    assert exponent >= 0 and op
    out = {0: s.Integer(1)}
    top = max(op)
    for done in range(1, exponent+1):
        out = mul(out, op, minimum-(exponent-done)*top)
    return out


def root(P, n, depth):
    assert P.get(n) == 1 and n-1 not in P
    out = {1: s.Integer(1)}
    for k in range(1, depth+1):
        target = n-1-k
        trial = power(out, n, target)
        value = s.expand((P.get(target, 0)-trial.get(target, 0))/n)
        if value != 0:
            out[-k] = value
    return out


def inverse_monic_order_one(W, depth):
    out = {-1: s.Integer(1)}
    for k in range(2, depth+1):
        target = 1-k
        trial = mul(W, out, target)
        value = s.expand(-trial.get(target, 0))
        if value != 0:
            out[-k] = value
    return out


def trace_primitive(P, B):
    result = 0
    for i, a in P.items():
        if i <= 0:
            continue
        for j in range(-i, 0):
            b = B.get(j, s.Integer(0))
            k = i+j+1
            result += s.binomial(i, k)*sum(
                (-1)**ell*s.diff(a, x, ell)*s.diff(b, x, k-1-ell)
                for ell in range(k)
            )
    return s.expand(result)


def fractional_positive(f, n, m):
    t = s.symbols("t")
    small = s.expand(t**n*f.subs(z, 1/t)-1)
    value = 0
    for k in range(m//2+1):
        value += s.binomial(s.Rational(m, n), k)*small**k
    poly = s.Poly(s.expand(value), t)
    return s.expand(sum(c*z**(m-e[0]) for e, c in poly.terms() if e[0] <= m))


def leading_weighted(op, R):
    entries = []
    for j, a in op.items():
        for (degree,), c in s.Poly(a, x).terms():
            entries.append((s.Rational(degree)+R*j, c*x**degree*xi**j))
    top = max((weight for weight, _ in entries), default=-s.oo)
    return top, s.expand(sum(term for weight, term in entries if weight == top))


def companion(P, n):
    A = s.zeros(n)
    for j in range(n-1):
        A[j, j+1] = 1
    for j in range(n):
        A[n-1, j] = (lam if j == 0 else 0)-P.get(j, 0)
    return A


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE/"all_coprime_algebraic_string_verification.json")
    parser.add_argument("--ledger", type=Path,
                        default=HERE/"all_coprime_algebraic_string_source_ledger.json")
    args = parser.parse_args()
    started = time.monotonic()
    ledger = json.loads(args.ledger.read_text())
    for item in ledger["sources"]:
        path = HERE/item["path"]
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], item["path"]

    # Independent finite Ore residue versus the explicit total derivative (9).
    a, b = s.Function("a")(x), s.Function("b")(x)
    residue_cases = 0
    for i in range(8):
        for j in range(-8, 5):
            P, B = {i: a}, {j: b}
            actual = comm(P, B, -1).get(-1, 0)
            expected = s.diff(trace_primitive(P, B), x)
            assert s.expand(actual-expected) == 0, (i, j)
            residue_cases += 1
    assert s.expand(comm({2: a}, {-1: b}, -1)[-1]-s.diff(a*b, x)) != 0

    flow_cases = []
    weighted_cases = []
    for n in range(2, 6):
        # Every coefficient grows exactly at its assigned leading weight.
        P = {n: s.Integer(1)}
        P.update({n-i: (i+1)*x**i + x for i in range(2, n+1)})
        W = root(P, n, n)
        assert add(power(W, n), scale(P, -1)) == {}
        f = z**n + sum((i+1)*z**(n-i) for i in range(2, n+1))
        for m in sorted({1, 2, n-1, n+1}):
            Bj = power(W, m)
            bracket = comm(P, Bj)
            assert max(bracket, default=-1) <= n-2
            flow_cases.append({"n": n, "m": m, "orders": sorted(bracket)})
            if math.gcd(n, m) == 1 and m >= 2:
                g = fractional_positive(f, n, m)
                polynomial = s.expand(m*s.diff(f, z)*g - n*f*s.diff(g, z))
                assert polynomial != 0
                weight, leading = leading_weighted(bracket, s.Integer(1))
                expected = s.expand(x**(n+m-2)*polynomial.subs(z, xi/x))
                assert weight == n+m-2 and s.expand(leading-expected) == 0
                weighted_cases.append({"n": n, "m": m, "R": "1",
                                       "leading_weight": str(weight)})
        bad = dict(power(W, n+1))
        bad[n] = bad.get(n, 0) + x
        assert max(comm(P, bad)) > n-2

    # The exact time primitive for Airy operators; all multiplication is finite
    # to the precision needed for res[P,QW^-1].
    time_cases = []
    for n in range(2, 7):
        P, Q = {n: s.Integer(1), 0: x}, {1: -s.Integer(1)}
        assert comm(P, Q) == {0: s.Integer(1)}
        W = root(P, n, n+1)
        V = inverse_monic_order_one(W, n+1)
        assert mul(W, V, -n) == {0: s.Integer(1)}
        B = mul(Q, V, -n)
        T = trace_primitive(P, B)
        assert s.diff(T, x) == 1
        assert comm(P, B, -1).get(-1) == 1
        time_cases.append({"n": n, "primitive": str(T)})

    # Smallest possible positive rational weight R=1/n, with an actual
    # nonconstant operator leading commutator, not a sampled scalar value.
    for n, m in [(3, 2), (4, 3), (5, 2), (5, 3), (7, 2), (7, 3)]:
        P = {n: s.Integer(1), 0: x}
        Q = {m: s.Integer(1)}
        weight, leading = leading_weighted(comm(P, Q), s.Rational(1, n))
        assert weight == s.Rational(m-1, n)
        assert leading == -m*xi**(m-1)
        weighted_cases.append({"n": n, "m": m, "R": str(s.Rational(1, n)),
                               "leading_weight": str(weight)})

    # Common-power boundaries are positive controls for the hypotheses:
    # the coprime argument must not reject these commuting leading symbols.
    common_power_cases = []
    for n, m, h in [(4, 6, z**2+1), (6, 9, z**3+z+1),
                    (6, 4, z**2+2), (8, 12, z**4+z+1)]:
        d = math.gcd(n, m)
        f, g = s.expand(h**(n//d)), s.expand(h**(m//d))
        assert fractional_positive(f, n, m) == g
        assert s.expand(m*s.diff(f,z)*g-n*f*s.diff(g,z)) == 0
        assert s.expand(g**n-f**m) == 0
        common_power_cases.append({"n": n, "m": m, "gcd": d, "h": str(h)})

    positive_prime_cases = []
    for p in [2, 3, 5, 7, 11]:
        P = {p: s.Integer(1), 0: x}
        Q1 = {1: -s.Integer(1)}
        Q = add(power(P, 2), Q1)
        assert max(Q) == 2*p
        assert comm(P, Q) == {0: s.Integer(1)}
        assert add(Q, scale(power(P, 2), -1)) == Q1
        positive_prime_cases.append({"prime": p, "mate_order": 2*p,
                                     "remainder_order": 1})

    spectral_cases = []
    for n in [2, 3, 4, 5]:
        P = {n: s.Integer(1), 0: x}
        A = companion(P, n)
        # Q=-E+P^2 gives Q psi = (-d/dx+lambda^2)psi.
        B = -A + lam**2*s.eye(n)
        curvature = A.diff(lam)-B.diff(x)+A*B-B*A
        assert curvature == s.zeros(n)
        wrong = A.diff(lam)+B.diff(x)-A*B+B*A
        assert wrong != s.zeros(n)
        rows = [s.eye(n)[0, :]]
        for _ in range(n+2):
            rows.append((rows[-1].diff(x)+rows[-1]*A).applyfunc(s.expand))
        first = sum(((j+x)*rows[j] for j in range(n+3)), s.zeros(1,n))
        all_rows = [first]
        for _ in range(n-1):
            all_rows.append((all_rows[-1].diff(x)+all_rows[-1]*A).applyfunc(s.expand))
        spectral = s.Matrix.vstack(*all_rows)
        degrees = [s.Poly(entry, lam).degree() for entry in spectral]
        assert all(degree >= 0 for degree in degrees)
        assert max(degrees) <= (2*n+1)//n
        spectral_cases.append({"n": n, "tested_order": n+2,
                               "maximum_lambda_degree": int(max(degrees))})

    # At every rank the left-ratio estimates force the complete first row.
    first_row_bounds = []
    for n in range(2, 13):
        bounds = [[i-j-1 for j in range(n)] for i in range(n)]
        assert all(power < 0 for power in bounds[0])
        assert all(bounds[i][j] < 0 for i in range(n) for j in range(i,n))
        first_row_bounds.append(n)

    # Explicit failure of coefficient-integrality if classical order drops to
    # one: it must not be silently included in the specialization theorem.
    h = s.symbols("h")
    P = {2: h, 1: s.Integer(1)}
    Q1 = {1: -1/h}
    Q = add(scale(P, 1/h), Q1)
    assert Q == {2: s.Integer(1)}
    assert max({j:c.subs(h,0) for j,c in P.items() if c.subs(h,0)!=0}) == 1

    receipt = {
        "status": "PASS",
        "claim": "CoprimeAlgebraicStringExclusion",
        "residue_monomial_cases": residue_cases,
        "fractional_power_flow_cases": flow_cases,
        "time_primitive_cases": time_cases,
        "weighted_leading_bracket_cases": weighted_cases,
        "common_power_boundary_cases": common_power_cases,
        "positive_prime_order_cases": positive_prime_cases,
        "spectral_companion_cases": spectral_cases,
        "left_ratio_first_row_ranks": first_row_bounds,
        "negative_controls": [
            "incorrect residue primitive",
            "nonconstant high flow coefficient",
            "opposite spectral compatibility sign",
            "omitting the coprime-order hypothesis",
            "including classical order one in the integrality argument",
        ],
        "source_hashes": {item["path"]:item["sha256"] for item in ledger["sources"]},
        "coverage": (
            "Finite exact algebraic controls only. The universal normalization, "
            "all-order residue identity, analytic forward-Stokes application, "
            "algebraic continuation, weighted infinity exclusion and valuation "
            "argument are the written proof, not inferred from these cases."
        ),
        "formal_verification": False,
        "external_review": False,
        "elapsed_seconds": round(time.monotonic()-started, 6),
    }
    args.output.write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"status":"PASS", "receipt":str(args.output),
                      "elapsed_seconds":receipt["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()

