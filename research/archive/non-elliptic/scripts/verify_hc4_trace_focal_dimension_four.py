#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "hc4_trace_focal_dimension_four.json"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)


def canonical_data(n: int):
    N = sp.zeros(n)
    for j in range(1, n):
        N[j - 1, j] = 1
    S = sp.zeros(n)
    for i in range(n):
        S[i, n - 1 - i] = 1
    assert N.T * S == S * N
    T = S * N

    B = sp.zeros(n)
    symbols = {}
    for i in range(n):
        for j in range(i, n):
            symbol = sp.symbols(f"b{i+1}{j+1}")
            symbols[i, j] = symbol
            B[i, j] = B[j, i] = symbol

    J = B * T
    trace = sp.factor(sp.trace(J))
    skew = S * J - (S * J).T
    frobenius = [
        sp.factor(skew[i, j])
        for i in range(n - 1)
        for j in range(i + 1, n - 1)
    ]

    # ell = e_n^*.  Its radical condition on the tangent hyperplane
    # span(e_2^*,...,e_n^*) kills B_{n,j}, j=2,...,n.
    radical_subs = {B[n - 1, j]: 0 for j in range(1, n)}
    return {
        "trace_before": str(trace),
        "trace_after": str(sp.factor(trace.subs(radical_subs))),
        "frob_before": [str(value) for value in frobenius],
        "frob_after": [str(sp.factor(value.subs(radical_subs))) for value in frobenius],
    }


record = {str(n): canonical_data(n) for n in range(3, 7)}
assert record["4"]["trace_after"] == "b33"
assert record["4"]["frob_after"] == ["0", "0", "b33"]
# In dimension five there are genuinely independent survivors after trace zero.
assert "b44" in record["5"]["frob_after"]
assert "-b24 + b33" in record["5"]["frob_after"]

result = {
    "scope": "trace-focal Frobenius comparison for regular nilpotent blocks",
    "status": "dimension-four coincidence verified",
    "dimensions": record,
}
OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
