#!/usr/bin/env python3
"""Replay the exact (6,10) polynomial parametrization gap certificate.

Default replay uses rational polynomial arithmetic only.  The optional
--regenerate-certificate reconstructs the six multipliers with Singular;
the resulting identity is still checked directly before any artifact write.
The normalization and Puiseux interpretation are written proof steps.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import sympy as sp
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts/generated-results/f2-degree-6-10-gap-v1.json"
VARIABLES = ["a", "b", "c", "d", "B", "D", "E"]
ORDERS = [1, 3, 5, 7, 9, 11]


def construct_obstructions():
    R, a, b, c, d, B, D, E = ring(",".join(VARIABLES), sp.QQ)
    u = {2: a, 3: b, 4: c, 5: d}
    powers = {}
    for j in range(-5, 6):
        alpha = sp.QQ(j, 3)
        values = [R.one] + [R.zero] * 21
        # (1+u) f' = alpha*u'*f, for f=(1+u)^alpha.
        for n in range(1, 22):
            values[n] = sum(
                (((alpha + 1) * i - n) * v * values[n-i]
                 for i, v in u.items() if i <= n), R.zero
            ) / n
        powers[j] = values
    coefficients = {5: R.one, 4: B, 3: R.zero, 2: D, 1: E, 0: R.zero}
    rho = {}
    for exponent in range(-1, -12, -1):
        value = sum(
            (v * powers[j][2*j-exponent] for j, v in coefficients.items()
             if 0 <= 2*j-exponent <= 21), R.zero
        )
        if exponent % 2 == 0:
            coefficients[exponent//2] = -value
        else:
            rho[-exponent] = value
    denominators, numerators = zip(*(rho[k].clear_denoms() for k in ORDERS))
    return R, rho, list(denominators), list(numerators), powers


def encode(poly):
    return [[list(m), str(v)] for m, v in sorted(poly.items())]


def decode(R, rows):
    assert rows == sorted(rows, key=lambda row: row[0])
    assert len({tuple(m) for m, _ in rows}) == len(rows)
    assert all(len(m) == 7 and all(type(e) is int and e >= 0 for e in m)
               for m, _ in rows)
    return R.from_dict({tuple(m): sp.QQ(v) for m, v in rows})


def singular_input(numerators):
    lines = ["ring R=0,(a,b,c,d,B,D,E,z),dp;"]
    lines += [f"poly f{k}={str(p).replace('**', '^')};"
              for k, p in zip(ORDERS, numerators)]
    lines += ["ideal I=f1,f3,f5,f7,f9,f11,b*z-1;",
              "option(redSB); matrix L=lift(I,ideal(1)); short=0;",
              'int i; for(i=1;i<=6;i++){print("CERT "+string(L[i,1]));}']
    return "\n".join(lines) + "\n"


def regenerate(R, numerators):
    with tempfile.TemporaryDirectory(prefix="f2-degree-gap-") as temp:
        path = Path(temp) / "certificate.sing"
        path.write_text(singular_input(numerators))
        result = subprocess.run(["Singular", "-q", str(path)], check=True,
                                text=True, capture_output=True, timeout=60)
    R8, *gens = ring(",".join(VARIABLES + ["z"]), sp.QQ)
    local_symbols = {str(g): sp.Symbol(str(g)) for g in gens}
    rows = [R8.from_expr(sp.sympify(line[5:].replace("^", "**"),
                                  locals=local_symbols))
            for line in result.stdout.splitlines() if line.startswith("CERT ")]
    assert len(rows) == 6, result.stdout[-1000:]
    # Substitute z=1/b, discarding the multiplier of b*z-1, then clear b.
    power = max(m[-1]-m[1] for row in rows for m in row)
    multipliers = []
    for row in rows:
        poly = R.zero
        for monomial, coefficient in row.items():
            exponents = list(monomial[:-1])
            exponents[1] += power-monomial[-1]
            poly += R.from_dict({tuple(exponents): coefficient})
        multipliers.append(poly)
    return power, multipliers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--regenerate-certificate", action="store_true")
    args = parser.parse_args()
    R, rho, denominators, numerators, powers = construct_obstructions()
    a, b, c, d, B, D, E = R.gens
    if args.regenerate_certificate:
        power, multipliers = regenerate(R, numerators)
    else:
        existing = json.loads(OUTPUT.read_text())
        assert existing["variables"] == VARIABLES
        power = existing["b_power"]
        multipliers = [decode(R, row) for row in existing["multipliers"]]
    assert power == 12 and len(multipliers) == 6
    residual = sum((x*y for x, y in zip(multipliers, numerators)), R.zero)-b**power
    assert residual == 0, "the literal polynomial identity failed"
    assert residual + numerators[0] != 0, "multiplier mutation was not detected"
    weights = [2, 3, 4, 5, 2, 6, 8]
    for k, multiplier, numerator in zip(ORDERS, multipliers, numerators):
        assert all(sum(w*e for w, e in zip(weights, m)) == 10+k for m in numerator)
        assert all(sum(w*e for w, e in zip(weights, m)) == 36-10-k for m in multiplier)

    # No division by d: the b=0 branch has a literal cubic obstruction.
    branch = rho[5] + 2*a*rho[3]/3 + c*rho[1]/3 + 5*d**3/81
    assert sp.expand(branch.as_expr().subs(sp.Symbol("b"), 0)) == 0
    # When b=d=0 all polynomial parts are even, including the discarded
    # free terms C*p+F, so the pair cannot parametrize its normalization.
    for j in [1, 2, 4, 5]:
        for n in range(1, 2*j+1, 2):
            assert sp.expand(powers[j][n].as_expr().subs(
                {sp.Symbol("b"): 0, sp.Symbol("d"): 0})) == 0

    witness_values = [0, 1, 0, 0, 0, sp.Rational(2, 9), 0]
    substitution = dict(zip(map(sp.Symbol, VARIABLES), witness_values))
    witness_obstructions = [sp.expand(rho[k].as_expr().subs(substitution)) for k in ORDERS]
    assert witness_obstructions == [0, 0, 0, 0, 0, sp.Rational(1, 19683)]
    t, s = sp.symbols("t s")
    p = t**6+t**3
    q = t**10+sp.Rational(5, 3)*t**7+sp.Rational(7, 9)*t**4+sp.Rational(7, 81)*t
    polynomial_part = sum(
        coefficient * sum(powers[j][n].as_expr()*t**(2*j-n) for n in range(2*j+1))
        for j, coefficient in [(5, 1), (4, sp.Symbol("B")),
                                (2, sp.Symbol("D")), (1, sp.Symbol("E"))]
    )
    assert sp.expand(polynomial_part.subs(substitution)-q) == 0
    assert sp.gcd(sp.diff(p, t), sp.diff(q, t)) == 1
    domain = sp.QQ.frac_field(t)
    common = sp.gcd(sp.Poly(p-p.subs(t, s), s, domain=domain),
                    sp.Poly(q-q.subs(t, s), s, domain=domain)).as_expr()
    assert common == s-t

    data = {
        "format": "f2-degree-6-10-gap-v1", "theorem": "PF2D6O1",
        "source_sha256": {Path(__file__).resolve().relative_to(ROOT).as_posix():
                          hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "variables": VARIABLES, "obstruction_orders": ORDERS,
        "obstruction_denominators": [str(x) for x in denominators],
        "obstruction_numerators": [encode(p) for p in numerators],
        "b_power": power, "multipliers": [encode(p) for p in multipliers],
        "identity": "b^12 = sum_i multipliers[i] * obstruction_numerators[i]",
        "multiplier_term_counts": [len(p) for p in multipliers],
        "b_zero_identity": "rho5 + (2*a/3)*rho3 + (c/3)*rho1 = -5*d^3/81",
        "singular_input_sha256": hashlib.sha256(singular_input(numerators).encode()).hexdigest(),
        "sharp_witness": {"p": str(p), "q": str(q),
                          "odd_obstructions": [str(x) for x in witness_obstructions],
                          "normalization_gcd": str(common), "derivative_gcd": "1",
                          "gap": 21, "F2_row": 7},
        "conclusion": "A birational polynomial (6,10) pair has first odd Puiseux gap at most 21; F2 normal row r=9 is impossible.",
        "boundary": "Normalization and Puiseux transport are written proofs. Rows r=5,7 and nonnormal terminal slices remain open. The witness is a target curve, not a Keller map.",
    }
    serialized = json.dumps(data, indent=2) + "\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.read_text() == serialized, "artifact differs; inspect before --write"
    print("PASS exact b^12 ideal-membership certificate (362 multiplier terms)")
    print("PASS complete b=0 branch and even-composition conclusion")
    print("PASS sharp r=7 polynomial normalization witness and gap 21")
    print("PASS " + ("wrote " if args.write else "byte-identical ") + str(OUTPUT))


if __name__ == "__main__":
    main()
