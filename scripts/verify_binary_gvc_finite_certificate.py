#!/usr/bin/env python3
"""Reproduce the finite binary-GVC certificates and adversarial checks.

The universal classification and degree bound are written proofs, not inferred
from these finite tests. --write creates the pinned artifact; the default checks
it byte-for-byte. --input accepts exact rational polynomial term lists and emits
one complete decision/certificate without evaluating a moment prefix.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import comb, factorial
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from jcsearch.binary_gvc_certificate import classify, decode, encode, strict_slope

OUTPUT = ROOT / "artifacts/generated-results/binary-gvc-finite-certificate-v1.json"
x, y = sp.symbols("x y")


def terms(expression) -> dict:
    return {exponent: Fraction(coefficient)
            for exponent, coefficient in sp.Poly(expression, x, y, domain=sp.QQ).terms()
            if coefficient}


def expression(rows):
    return sum(sp.Rational(c) * x**i * y**j for i, j, c in rows)


def differential(symbol, polynomial, m):
    operator = sp.Poly(sp.expand(symbol**m), x, y, domain=sp.QQ)
    p = sp.expand(polynomial**m)
    return sp.expand(sum(c * sp.diff(p, x, i, y, j)
                         for (i, j), c in operator.terms()))


def mixed_differential(symbol, polynomial, multiplier, m):
    operator = sp.Poly(sp.expand(symbol**m), x, y, domain=sp.QQ)
    p = sp.expand(multiplier * polynomial**m)
    return sp.expand(sum(c * sp.diff(p, x, i, y, j)
                         for (i, j), c in operator.terms()))


def independently_check_positive(record):
    """Rebuild the actual coordinate change in SymPy, then check the weights.

    This deliberately does not call the certificate generator's transform or
    slope solver. It checks sufficiency, not the Hall/envelope completeness proof.
    """
    assert record["pure_status"] == "all_powers_zero"
    a = expression(record["input_symbol"])
    p = expression(record["input_polynomial"])
    if record["method"] == "zero input":
        assert a == 0 or p == 0
        if record["input_multiplier"] is not None:
            assert record["mixed_cutoff"] == 1
        return
    chart = record["chart"]
    if chart["kind"] == "identity":
        a1, p1 = a, p
        substitution = {x: x, y: y}
    elif chart["kind"] == "swap":
        substitution = {x: y, y: x}
        a1 = a.subs(substitution, simultaneous=True)
        p1 = p.subs(substitution, simultaneous=True)
    else:
        root = sp.Rational(chart["root"])
        substitution = {x: x, y: y-root*x}
        a1 = a.subs(x, x+root*y)
        p1 = p.subs(substitution, simultaneous=True)
    assert terms(a1) == decode(record["transformed_symbol"])
    assert terms(p1) == decode(record["transformed_polynomial"])
    u, v = record["separation"]["weight"]
    assert type(u) is int and type(v) is int and min(u, v) > 0
    amin = min(u*i+v*j for i, j in terms(a1))
    pmax = max(u*i+v*j for i, j in terms(p1))
    assert amin == record["separation"]["operator_min"]
    assert pmax == record["separation"]["polynomial_max"]
    assert amin-pmax == record["separation"]["gap"] > 0
    R = int(sp.Poly(a, x, y).total_degree())
    d = int(sp.Poly(p, x, y).total_degree())
    assert max(u, v) <= R+d
    if record["input_multiplier"] is not None:
        q = expression(record["input_multiplier"])
        q1 = q.subs(substitution, simultaneous=True)
        assert terms(q1) == decode(record["transformed_multiplier"])
        k = max((u*i+v*j for i, j in terms(q1)), default=0)
        assert k == record["multiplier_weight"]
        assert record["mixed_cutoff"] == k//(amin-pmax)+1
        q_degree = int(sp.Poly(q, x, y).total_degree()) if q else 0
        assert record["degree_only_cutoff"] == (R+d)*q_degree+1
        assert record["mixed_cutoff"] <= record["degree_only_cutoff"]


def run_regressions():
    fixtures = [
        ("zero_symbol", 0, x+y, x),
        ("zero_polynomial", x+y, 0, x),
        ("constant_polynomial", x*y, 3, x**2+y),
        ("constant_symbol_rejected", 1+x, x-y, 1),
        ("ordinary_degree", x*y, x+1, y**3),
        ("equal_degree_repeated_direction", (x+y)**2, (x-y)**2+1, x*y),
        ("equal_degree_simple_direction", x**2-y**2, (x+y)**2, x),
        ("nonhomogeneous", x**2+y**7, x*y**2+y**3, y**6),
        ("rational_rotated", (x-sp.Rational(2,3)*y)**2+y**7,
         x*(y+sp.Rational(2,3)*x)**2+(y+sp.Rational(2,3)*x)**3, x+y),
        ("infinity_direction", y**2+x**7, y*x**2+x**3, y),
        ("hall_but_no_separator", x**2+y**4, x*y**2+y**3, x),
        ("nonrational_factors_rejected", x**2+y**2, x**2-y**2, 1),
    ]
    rows = []
    for label, a, p, q in fixtures:
        record = classify(terms(a), terms(p), terms(q))
        if record["pure_status"] == "all_powers_zero":
            independently_check_positive(record)
            assert all(differential(a, p, m) == 0 for m in range(1, 5))
            m = record["mixed_cutoff"]
            assert mixed_differential(a, p, q, m) == 0
        else:
            failures = [m for m in range(1, 5) if differential(a, p, m) != 0]
            assert failures, "negative fixture needs a separate exact witness"
            record["independent_first_nonzero_pure_power"] = min(failures)
        rows.append({"name": label, "certificate": record})

    # Both these inputs pass the first pure identity and fail at the second.
    for label in ("hall_but_no_separator", "nonrational_factors_rejected"):
        row = next(row for row in rows if row["name"] == label)
        assert row["certificate"]["independent_first_nonzero_pure_power"] == 2

    # Complete small support regression for the arithmetic mediant lemma.
    # It is independent of coefficient choices and of the GVC premise.
    monomials = [(i, j) for i in range(4) for j in range(4-i) if i+j > 0]
    supports = [(m,) for m in monomials] + list(combinations(monomials, 2))
    feasible = 0
    for sa in supports:
        for spoly in supports:
            A, P = dict.fromkeys(sa, Fraction(1)), dict.fromkeys(spoly, Fraction(1))
            result = strict_slope(A, P)
            bound = max(map(sum, sa)) + max(map(sum, spoly))
            brute = [(u, v) for u in range(2, bound+1) for v in range(1, u)
                     if min(u*i+v*j for i, j in sa) > max(u*i+v*j for i, j in spoly)]
            assert bool(brute) == result["feasible"]
            feasible += result["feasible"]

    # Sharp last mixed power for Lambda=dx+dy^(d+1), P=y^d, Q=x^q.
    sharp = []
    for d in range(1, 4):
        for q in range(1, 3):
            m = (d+1)*q
            a, p, Q = x+y**(d+1), y**d, x**q
            record = classify(terms(a), terms(p), terms(Q))
            independently_check_positive(record)
            value = mixed_differential(a, p, Q, m)
            expected = comb(m, q)*factorial(q)*factorial(d*m)
            assert value == expected != 0
            assert mixed_differential(a, p, Q, m+1) == 0
            assert record["mixed_cutoff"] == m+1
            sharp.append({"d": d, "q": q, "last_nonzero_m": m,
                          "last_nonzero_value": str(value), "certificate": record})

    positive = next(row["certificate"] for row in rows if row["name"] == "nonhomogeneous")
    mutant = json.loads(json.dumps(positive))
    mutant["separation"]["weight"] = [1, 1]
    try:
        independently_check_positive(mutant)
    except AssertionError:
        pass
    else:
        raise AssertionError("invalid weight mutation was accepted")
    mutant = json.loads(json.dumps(positive))
    mutant["transformed_polynomial"] = [[0, 0, "1"]]
    try:
        independently_check_positive(mutant)
    except AssertionError:
        pass
    else:
        raise AssertionError("false coordinate-transport mutation was accepted")

    sources = ["jcsearch/binary_gvc_certificate.py",
               "scripts/verify_binary_gvc_finite_certificate.py"]
    return {"format": "binary-gvc-finite-certificate-regressions-v1",
            "theorem": "GVC2SC",
            "assurance": "written universal proof; exact finite algorithm and independent positive-certificate regressions",
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
            "fixtures": rows, "small_support_pairs": len(supports)**2,
            "small_support_feasible": feasible, "sharp_family": sharp,
            "adversarial_mutations_rejected": ["weight", "coordinate transport"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--input", type=Path,
                        help="JSON with symbol, polynomial, and optional multiplier term lists")
    args = parser.parse_args()
    if args.input:
        if args.write:
            parser.error("--write cannot be combined with --input")
        data = json.loads(args.input.read_text())
        q = decode(data["multiplier"]) if "multiplier" in data else None
        result = classify(decode(data["symbol"]), decode(data["polynomial"]), q)
        if result["pure_status"] == "all_powers_zero":
            independently_check_positive(result)
        print(json.dumps(result, indent=2))
        return
    data = run_regressions()
    serialized = json.dumps(data, indent=2) + "\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.read_text() == serialized, "artifact differs; inspect before deliberate --write"
    print("PASS 12 rational-input fixtures and independent coordinate/weight checks")
    print(f"PASS all {data['small_support_pairs']} small support pairs; {data['small_support_feasible']} feasible")
    print("PASS six sharp-cutoff family instances, including the last nonzero mixed power")
    print("PASS invalid weight and coordinate-transport certificates rejected")
    print("PASS " + ("wrote " if args.write else "byte-identical ") + str(OUTPUT))
    print("Universal proof and negative-decision completeness remain written mathematics, not a bounded inference")


if __name__ == "__main__":
    main()
