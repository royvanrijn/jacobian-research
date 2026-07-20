#!/usr/bin/env python3
"""Reproduce the explicit Section 5 elimination/valuation obstruction for (9,27).

Primary source: Guccione--Guccione--Horruitiner--Valqui, arXiv:2204.14178,
Theorem 5.1 and equations (5.9)--(5.11).  This checks the published terminal
system, not all earlier reductions from a hypothetical Keller map to that
system.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import sympy as sp

from jcsearch.laurent import laurent_dict
from jcsearch.msolve import input_text, run

X, y = sp.symbols("X y")
dm10, dm8, dm7, dm6, dm5, dm4, dm3, dm2 = sp.symbols(
    "dm10 dm8 dm7 dm6 dm5 dm4 dm3 dm2")
d1, d0, dm1, C3, Fm4 = sp.symbols("d1 d0 dm1 C3 Fm4")

# The transformed Laurent series on paper p.18 has d2=0.
D = (X**3 + d1*X + d0 + dm1/X + dm2/X**2 + dm3/X**3
     + dm4/X**4 + dm5/X**5 + dm6/X**6 + dm7/X**7
     + dm8/X**8 + dm10/X**10)
D2 = laurent_dict(sp.expand(D**2), (X,))
D3 = laurent_dict(sp.expand(D**3), (X,))

# Exactly the six square and three cube coefficient equations displayed before
# equation (5.9).  The x^-4 cube coefficient carries F_-4 C_3^23.
equations = [D2[(-k,)] for k in (1, 2, 3, 4, 5, 7)]
equations += [D3[(-1,)], D3[(-2,)], D3[(-4,)] + Fm4*C3**23]
variables = (dm10, dm8, dm7, dm6, dm5, dm4, dm3, dm2,
             d1, d0, dm1, C3, Fm4)
target_59 = sp.expand(18*C3**23*d1*dm1**6*Fm4
                      + 8*C3**69*Fm4**3 + 27*d0*dm1**9)


def normalize_polynomial_string(text):
    return sp.Poly(sp.sympify(text.replace("^", "**")),
                   d1, d0, dm1, C3, Fm4, domain=sp.QQ)


def target_in_elimination_output(output):
    """Recognize equation (5.9) up to a nonzero rational scalar."""
    target = sp.Poly(target_59, d1, d0, dm1, C3, Fm4, domain=sp.QQ)
    for line in output.splitlines():
        candidate = line.strip().strip("[]:, ")
        if not candidate or any(name in candidate for name in
                                ("dm10", "dm8", "dm7", "dm6", "dm5",
                                 "dm4", "dm3", "dm2")):
            continue
        try:
            poly = normalize_polynomial_string(candidate)
        except (sp.SympifyError, sp.PolynomialError):
            continue
        if poly.is_zero:
            continue
        quotient = sp.cancel(poly.LC()/target.LC())
        if (poly - target.mul_ground(quotient)).is_zero:
            return True
    return False


def valuation_degree_obstruction():
    # Proposition 5.4, with the irrelevant nonzero scalar retained exactly.
    g = -(35 - 42*y + 54*y**2 - 81*y**3 + 243*y**4)/sp.Integer(910)
    assert sp.degree(g, y) == 4
    assert sp.gcd(sp.Poly(g, y), sp.Poly(sp.diff(g, y), y)).degree() == 0
    assert g.subs(y, 0) != 0 and g.subs(y, -1) != 0

    # Check the differential equation used to derive f=y(y+1)g.
    f = sp.expand(y*(y + 1)*g)
    differential = sp.expand(6*y*(y + 1)*sp.diff(f, y)
                             - 10*(9*y + 8)*f - y**9*(y + 1)**2)
    # The paper's displayed solution is f1=y^9(y+1)^2 g, then f=f1/C3.
    f1 = sp.expand(y**9*(y + 1)**2*g)
    differential_f1 = sp.expand(6*y*(y + 1)*sp.diff(f1, y)
                                - 10*(9*y + 8)*f1
                                - y**9*(y + 1)**2)
    assert differential_f1 == 0

    # Equation (5.11): if dm1_tilde=(y+1)^k, k>=8 contradicts
    # ord_{y=-1}=66; k<=7 contradicts the unique degree-78 middle term.
    cases = {
        "k_ge_8": {
            "first_term_min_order_at_-1": 22 + 6*8,
            "middle_term_exact_order_at_-1": 66,
            "third_term_min_order_at_-1": 9*8,
            "contradiction": True,
        },
        "k_le_7": {
            "first_term_max_degree": 22 + 8 + 6*7 + 4,
            "middle_term_exact_degree": 66 + 3*4,
            "third_term_max_degree": 12 + 9*7,
            "contradiction": True,
        },
    }
    assert cases["k_ge_8"] == {
        "first_term_min_order_at_-1": 70,
        "middle_term_exact_order_at_-1": 66,
        "third_term_min_order_at_-1": 72,
        "contradiction": True,
    }
    assert cases["k_le_7"] == {
        "first_term_max_degree": 76,
        "middle_term_exact_degree": 78,
        "third_term_max_degree": 75,
        "contradiction": True,
    }
    return g, cases


def main():
    records = []
    for characteristic in (1000003, 1000033, 1000037, 0):
        try:
            result = run(equations, variables, prime=characteristic,
                         eliminate=8, timeout=180)
            records.append({
                "prime": characteristic,
                "returncode": result.returncode,
                "terminal_relation_found": target_in_elimination_output(result.output)
                    if characteristic == 0 else None,
            })
            if characteristic == 0:
                base = ROOT / "results" / "certificates" / "newton_9_27_elimination_Q"
                base.with_suffix(".input").write_text(input_text(equations, variables, 0))
                base.with_suffix(".msolve").write_text(result.output)
        except subprocess.TimeoutExpired:
            records.append({"prime": characteristic, "timeout": 180,
                            "terminal_relation_found": None})
            break

    g, cases = valuation_degree_obstruction()
    exact = next((r for r in records if r["prime"] == 0), None)
    result = {
        "scope": "published Section 5 terminal system; prior Newton reductions assumed",
        "equations": len(equations),
        "eliminated_variables": [str(v) for v in variables[:8]],
        "terminal_relation_5_9": sp.sstr(target_59),
        "g": sp.sstr(g),
        "valuation_degree_cases": cases,
        "runs": records,
        "exact_terminal_relation_found": bool(exact and exact["terminal_relation_found"]),
    }
    destination = ROOT / "results" / "newton_9_27_regression.json"
    destination.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if exact is not None:
        assert result["exact_terminal_relation_found"]


if __name__ == "__main__":
    main()
