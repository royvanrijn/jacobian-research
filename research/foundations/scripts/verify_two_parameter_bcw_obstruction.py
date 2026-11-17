#!/usr/bin/env python3
"""Exact obstruction to the no-doubling homogenization of the current K."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import extract_quadratic_cubic
from verify_rank_compressed_bcw_24_route import source_map


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "shared_bcw_33_counterexample.json"


def main() -> None:
    variables, expressions = source_map(json.loads(SOURCE.read_text()))
    quadratic, cubic = extract_quadratic_cubic(expressions, variables)
    jq = sp.Matrix([poly.as_expr() for poly in quadratic]).jacobian(variables)
    jc = sp.Matrix([poly.as_expr() for poly in cubic]).jacobian(variables)
    at_ones = dict.fromkeys(variables, sp.Integer(1))
    jq_value = jq.subs(at_ones)
    jc_value = jc.subs(at_ones)

    # The Keller scaling identity holds on the parabola t=s^2.
    assert (sp.eye(16) + jq_value + jc_value).det() == 1
    assert (sp.eye(16) + 2 * jq_value + 4 * jc_value).det() == 1

    # Independent parameters fail even at the simplest integral point.
    cubic_only = (sp.eye(16) + jc_value).det()
    quadratic_only = (sp.eye(16) + jq_value).det()
    assert cubic_only == -4160
    assert quadratic_only == -78
    assert cubic_only != 1 and quadratic_only != 1

    print("PASS two-parameter obstruction: parabola samples s=1,2 have determinant one")
    print("PASS two-parameter obstruction: at X=(1,...,1), det(I+JC)=-4160")
    print("PASS two-parameter obstruction: current K cannot use the 17-variable shortcut")


if __name__ == "__main__":
    main()
