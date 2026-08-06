#!/usr/bin/env python3
"""Exclude every dense single-shear thickening of the two-row charts.

Fix pivot rows r<s and the internal-gauge quotient U[{r,s},:]=I_2.  The
off-diagonal B support has B[0,r]=B[1,s]=0 and every other entry nonzero.
For either factor column ell and k outside {r,s}, add the single shear
U[k,ell]=a and localize at a and every supported B entry.  The resulting
coefficient matrix C=U*B has exact rank two.  Its only diagonal entry is
C[k,k]=a*B[ell,k], so mu_1=k!*(4-k)!*a*B[ell,k] is nonzero on the torus.

This is a direct factor-chart calculation after quotienting the internal
GL_2 gauge.  No determinantal ideal or moment elimination is used.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations
import json
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_single_shear.json"
)


def chart_record(r: int, s: int, factor_column: int, shear_row: int) -> dict[str, object]:
    a = sp.Symbol("a")
    b = {
        (0, column): sp.Symbol(f"b0{column}")
        for column in range(5)
        if column != r
    }
    b.update(
        {
            (1, column): sp.Symbol(f"b1{column}")
            for column in range(5)
            if column != s
        }
    )

    u = sp.zeros(5, 2)
    u[r, 0] = 1
    u[s, 1] = 1
    u[shear_row, factor_column] = a
    factor = sp.zeros(2, 5)
    for position, coefficient in b.items():
        factor[position] = coefficient
    coefficient_matrix = u * factor

    pivot_minor = sp.expand(
        factor[0, r] * factor[1, s] - factor[0, s] * factor[1, r]
    )
    expected_minor = -b[(0, s)] * b[(1, r)]
    assert sp.expand(pivot_minor - expected_minor) == 0

    diagonal = tuple(sp.expand(coefficient_matrix[index, index]) for index in range(5))
    expected_diagonal = [sp.Integer(0)] * 5
    expected_diagonal[shear_row] = a * b[(factor_column, shear_row)]
    assert diagonal == tuple(expected_diagonal)

    first_moment = sp.expand(
        sum(
            factorial(index)
            * factorial(4 - index)
            * coefficient_matrix[index, index]
            for index in range(5)
        )
    )
    expected_moment = sp.expand(
        factorial(shear_row)
        * factorial(4 - shear_row)
        * a
        * b[(factor_column, shear_row)]
    )
    assert sp.expand(first_moment - expected_moment) == 0

    localizer = sp.expand(a * sp.prod(b.values()))
    assert a in localizer.free_symbols
    assert b[(factor_column, shear_row)] in localizer.free_symbols
    return {
        "pivot_rows": [r, s],
        "factor_column": factor_column,
        "shear_row": shear_row,
        "gauge": f"U[{r},:]=(1,0), U[{s},:]=(0,1)",
        "shear": f"U[{shear_row},{factor_column}]=a",
        "B_zero_positions": [[0, r], [1, s]],
        "B_torus_coordinate_count": len(b),
        "exact_rank_two_minor": str(pivot_minor),
        "only_diagonal_coefficient": str(diagonal[shear_row]),
        "mu_1": str(first_moment),
        "coefficient_torus_localizer_sha256": sha256(
            str(localizer).encode()
        ).hexdigest(),
        "conclusion": "mu_1 is nonzero on the dense single-shear torus",
    }


def main() -> None:
    records = []
    for r, s in combinations(range(5), 2):
        for factor_column in range(2):
            for shear_row in range(5):
                if shear_row not in {r, s}:
                    records.append(chart_record(r, s, factor_column, shear_row))
    assert len(records) == 60

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-single-shear-v1",
        "field": "characteristic zero",
        "chart_count": len(records),
        "parameterization": "C=U*B with U[pivot_rows,:]=I_2",
        "internal_GL2_gauge_removed": True,
        "B_support": "B[0,r]=B[1,s]=0; all other B entries nonzero",
        "rank_certificate": "the (r,s)-column minor of B is -B[0,s]*B[1,r]",
        "moment_certificate": "mu_1=k!*(4-k)!*a*B[ell,k] is nonzero",
        "global_conclusion": (
            "all sixty dense genuinely non-coordinate single-shear factor "
            "tori are SIC-safe"
        ),
        "scope": (
            "one nonpivot entry of U and the dense off-diagonal B torus; "
            "multiple shears and coefficient boundaries remain open"
        ),
        "charts": records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS internal GL2 gauge removed on all sixty single-shear charts")
    print("PASS every chart has exact coefficient rank two")
    print("PASS mu_1 is a nonzero localized monomial on every chart")


if __name__ == "__main__":
    main()
