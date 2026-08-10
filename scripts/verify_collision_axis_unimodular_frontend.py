#!/usr/bin/env python3
"""Verify the collision-axis unimodular frontend and census pruning counts."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.low_degree_census import Support  # noqa: E402


SUPPORT_LEDGER = (
    ROOT
    / "artifacts/generated-results/global_low_degree_census_02_supports.json"
)


def elementary_matrix(
    size: int, row: int, column: int, entry: sp.Expr
) -> sp.Matrix:
    matrix = sp.eye(size)
    matrix[row, column] = entry
    return matrix


def verify_minimal_elementary_length(t: sp.Symbol) -> None:
    # A one-factor word always leaves the first entry of its first column equal
    # to one, so its first-column integral cannot vanish.
    entry = sp.Symbol("u")
    for row in range(3):
        for column in range(3):
            if row == column:
                continue
            column_one = elementary_matrix(3, row, column, entry)[:, 0]
            assert column_one[0] == 1

    # Two factors suffice over every characteristic-zero field.
    f = 24 * t
    g = t - sp.Rational(3, 2) * t**2
    matrix = elementary_matrix(3, 0, 1, f) * elementary_matrix(3, 1, 0, g)
    assert sp.expand(matrix.det()) == 1
    assert matrix.subs(t, 0) == sp.eye(3)
    assert all(
        sp.integrate(entry, (t, 0, 1)) == 0 for entry in matrix[:, 0]
    )

    # Even a relative elementary operation can leave the normalized moment
    # fiber, so that fiber is not an E_3(k[t], (t))-space.
    feasible_column = sp.Matrix([1 - 2 * t, 2 * t - 3 * t**2, 0])
    transformed = elementary_matrix(3, 1, 0, t) * feasible_column
    assert transformed.subs(t, 0) == sp.Matrix([1, 0, 0])
    assert sp.integrate(transformed[1], (t, 0, 1)) == -sp.Rational(1, 6)


def verify_three_occurrence_frontier(t: sp.Symbol) -> sp.Matrix:
    # The restriction h(t)=F(t e_1) has exactly three nonlinear occurrences.
    restricted_map = sp.Matrix([t - t**2, t**2 - t**3, 0])
    column = restricted_map.diff(t)
    nonlinear_occurrences = sum(
        sum(monomial[0] >= 2 for monomial, _coefficient in sp.Poly(entry, t).terms())
        for entry in restricted_map
        if entry != 0
    )
    assert nonlinear_occurrences == 3
    assert restricted_map.subs(t, 0) == sp.zeros(3, 1)
    assert restricted_map.subs(t, 1) == sp.zeros(3, 1)
    assert column.subs(t, 0) == sp.Matrix([1, 0, 0])
    assert sp.gcd(sp.Poly(column[0], t), sp.Poly(column[1], t)) == 1

    # This is an explicit normalized SL_3(k[t])-completion of that column.
    matrix = sp.Matrix(
        [
            [1 - 2 * t, -8 * t, 0],
            [2 * t - 3 * t**2, 1 + 2 * t - 12 * t**2, 0],
            [0, 0, 1],
        ]
    )
    assert matrix[:, 0] == column
    assert matrix.subs(t, 0) == sp.eye(3)
    assert sp.expand(matrix.det()) == 1
    return matrix


def verify_minimal_support_resultant(t: sp.Symbol) -> None:
    # For h=(t-t^p, c(t^q-t^r), 0), failure of unimodularity is
    # equivalent to r^(p-1)=p^(r-q) q^(p-1).  The bounded loop is a
    # regression for the binomial-resultant proof in the canonical note.
    for p in range(2, 13):
        for q in range(2, 13):
            for r in range(q + 1, 13):
                first = sp.Poly(1 - p * t ** (p - 1), t)
                second = sp.Poly(q * t ** (q - 1) - r * t ** (r - 1), t)
                has_common_root = sp.degree(sp.gcd(first, second)) > 0
                binomial_equality = r ** (p - 1) == p ** (r - q) * q ** (p - 1)
                assert has_common_root == binomial_equality

    assert all(
        r ** (p - 1) != p ** (r - q) * q ** (p - 1)
        for p in range(2, 8)
        for q in range(2, 8)
        for r in range(q + 1, 8)
    )
    assert 6**8 == 9**4 * 2**8


def verify_first_jet_lift(t: sp.Symbol, matrix: sp.Matrix) -> None:
    y, z = sp.symbols("y z")
    restricted_map = sp.Matrix([t - t**2, t**2 - t**3, 0])
    first_jet_lift = restricted_map + y * matrix[:, 1] + z * matrix[:, 2]
    jacobian = first_jet_lift.jacobian((t, y, z))
    assert jacobian.subs({y: 0, z: 0}) == matrix
    assert sp.expand(jacobian.det().subs({y: 0, z: 0})) == 1
    # Axis integrability is automatic, but this naive lift is not a global
    # Keller map.  Higher transverse terms must solve the remaining equations.
    assert sp.expand(jacobian.det() - 1) != 0


def pure_axis_counts(support: list[list[list[int]]]) -> tuple[int, int, int]:
    return tuple(
        sum(exponent[1] == 0 and exponent[2] == 0 for exponent in component)
        for component in support
    )  # type: ignore[return-value]


def verify_pinned_census_pruning() -> None:
    ledger = json.loads(SUPPORT_LEDGER.read_text())
    kept_orbits: dict[str, int] = {}
    kept_labelled: dict[str, int] = {}
    for size, rows in ledger["orbits"].items():
        kept = []
        for row in rows:
            counts = pure_axis_counts(row["support"])
            support = Support(
                tuple(
                    tuple(tuple(exponent) for exponent in component)
                    for component in row["support"]
                )
            )
            assert support.collision_axis_counts() == counts
            assert counts[0] >= 1 and counts[1] != 1 and counts[2] != 1
            if support.collision_axis_unimodularity_possible():
                kept.append(row)
        kept_orbits[size] = len(kept)
        kept_labelled[size] = sum(row["orbit_size"] for row in kept)

    assert kept_orbits == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 0,
        "5": 0,
        "6": 450,
    }
    assert kept_labelled == {
        "1": 0,
        "2": 0,
        "3": 0,
        "4": 0,
        "5": 0,
        "6": 900,
    }


def main() -> None:
    t = sp.symbols("t")
    verify_minimal_elementary_length(t)
    matrix = verify_three_occurrence_frontier(t)
    verify_minimal_support_resultant(t)
    verify_first_jet_lift(t, matrix)
    verify_pinned_census_pruning()
    print(
        "PASS collision-axis frontend: support minimum 3; elementary length 2; "
        "degree-seven balanced supports reduce to 0/0/900 in sizes 4/5/6"
    )


if __name__ == "__main__":
    main()
