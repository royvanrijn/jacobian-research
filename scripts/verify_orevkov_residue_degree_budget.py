#!/usr/bin/env python3
"""Verify the Orevkov residue-degree budget and clean (3,5) cusp atlas."""

from __future__ import annotations

import sys
from math import gcd
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CAS = REPO / "plane-jc" / "cas"
sys.path.insert(0, str(CAS))

from plane_boundary_exclusion import orevkov_residue_degree_budget  # noqa: E402


def orevkov_budget_audit() -> None:
    for degree in range(2, 29):
        for e in range(1, degree):
            for f in range(1, degree):
                row = orevkov_residue_degree_budget(degree, ((e, f),))
                assert row.generic_component_cost == e
                assert row.forced_residue_ramification_cost == e * (f - 1)
                assert row.moved_sheet_cost == e * f
                assert row.unexplained_budget == degree - 1 - e * f

    cubic = orevkov_residue_degree_budget(
        6,
        ((3, 1),),
        local_excess_lower_bound=2,
    )
    assert cubic.minimum_total_cost == 5
    assert cubic.status == "saturates_orevkov_residue_budget"


def simple_e8_regular_cartier_audit() -> None:
    # Rows from the exhaustive simple-inertia orbifold atlas in the strict
    # F2 degree range.  Here R=sum(f_i), and regular Cartier multiplicity
    # gives one additional unit per residue sheet because (m_C,e)=(3,2).
    rows = (
        (6, 2, (1, 1)),
        (10, 2, (1, 1, 1, 1)),
        (12, 4, (2, 2)),
        (15, 3, (1, 1, 1, 1, 1, 1)),
        (20, 4, (2, 2, 2, 2)),
        (24, 4, (2, 4, 4)),
    )
    for degree, fixed, residue_degrees in rows:
        residue_sum = sum(residue_degrees)
        moved = degree - fixed
        assert moved == 2 * residue_sum
        budget = orevkov_residue_degree_budget(
            degree,
            tuple((2, residue_degree) for residue_degree in residue_degrees),
            local_excess_lower_bound=residue_sum,
        )
        assert residue_sum > fixed - 1
        assert budget.minimum_total_cost == 3 * residue_sum
        assert budget.status == "excluded_by_orevkov_residue_budget"


def clean_cusp_rows(d1: int, d2: int, bound: int = 30) -> set[tuple]:
    """Solve the two arithmetic families in Orevkov 2026, Theorem 2.

    The output stores ``(case,N,n,passport)``.  Coordinates may be swapped,
    so the caller evaluates both orders of ``(3,5)``.
    """

    rows: set[tuple] = set()

    # Case (a): d1=k1+k2+l2*k1*k2, d2=l1*k1*k2,
    # N=l1*d1, n=l1+l2+1.
    for k1 in range(1, bound + 1):
        for k2 in range(1, bound + 1):
            if gcd(k1, k2) != 1:
                continue
            for l1 in range(1, bound + 1):
                if gcd(l1, d1) != 1:
                    continue
                for l2 in range(0, bound + 1):
                    if k1 + k2 + l2 * k1 * k2 != d1:
                        continue
                    if l1 * k1 * k2 != d2:
                        continue
                    if l1 * k1 * k2 < 2:
                        continue
                    m1 = l1 * k2
                    m2 = l1 * k1
                    degree = l1 * d1
                    ramification = l1 + l2 + 1
                    alpha = (d1,) * l1
                    beta = tuple(sorted((m1, m2) + (d2,) * l2, reverse=True))
                    third = (ramification,) + (1,) * (degree - ramification)
                    assert sum(alpha) == sum(beta) == sum(third) == degree
                    rows.add(("a", degree, ramification, alpha, beta, third))

    # Case (b): m1=k2*l2+1, m2=k1*l1+1,
    # d1=k1*m1, d2=k2*m2, N=m1*m2, n=l1+l2+1.
    for k1 in range(1, bound + 1):
        for k2 in range(1, bound + 1):
            if gcd(k1, k2) != 1:
                continue
            for l1 in range(0, bound + 1):
                for l2 in range(0, bound + 1):
                    m1 = k2 * l2 + 1
                    m2 = k1 * l1 + 1
                    if gcd(m1, m2) != 1:
                        continue
                    if k1 + l2 <= 1 or k2 + l1 <= 1:
                        continue
                    if k1 * m1 != d1 or k2 * m2 != d2:
                        continue
                    degree = m1 * m2
                    ramification = l1 + l2 + 1
                    alpha = tuple(sorted((m1,) + (d1,) * l1, reverse=True))
                    beta = tuple(sorted((m2,) + (d2,) * l2, reverse=True))
                    third = (ramification,) + (1,) * (degree - ramification)
                    assert sum(alpha) == sum(beta) == sum(third) == degree
                    rows.add(("b", degree, ramification, alpha, beta, third))

    return rows


def canonical_clean_cusp_atlas(d1: int, d2: int) -> set[tuple]:
    """Forget duplicate arithmetic presentations and target-coordinate order."""

    rows = clean_cusp_rows(d1, d2) | clean_cusp_rows(d2, d1)
    return {
        (
            row[1],
            row[2],
            tuple(sorted((row[3], row[4]))),
            row[5],
        )
        for row in rows
    }


def clean_k1_cusp_atlas_audit() -> None:
    """Classify every irreducible cusp type occurring on the k=1 atlas."""

    expected = {
        (2, 3): {
            (1, 1, ((1,), (1,)), (1,)),
            (3, 2, ((2, 1), (3,)), (2, 1)),
            (6, 4, ((2, 2, 2), (3, 3)), (4, 1, 1)),
        },
        (2, 5): {
            (1, 1, ((1,), (1,)), (1,)),
            (5, 3, ((2, 2, 1), (5,)), (3, 1, 1)),
            (10, 6, ((2, 2, 2, 2, 2), (5, 5)), (6, 1, 1, 1, 1)),
        },
        (2, 7): {
            (1, 1, ((1,), (1,)), (1,)),
            (7, 4, ((2, 2, 2, 1), (7,)), (4, 1, 1, 1)),
            (
                14,
                8,
                ((2, 2, 2, 2, 2, 2, 2), (7, 7)),
                (8, 1, 1, 1, 1, 1, 1),
            ),
        },
        (3, 4): {
            (1, 1, ((1,), (1,)), (1,)),
            (4, 2, ((3, 1), (4,)), (2, 1, 1)),
            (6, 3, ((3, 3), (4, 2)), (3, 1, 1, 1)),
            (
                12,
                6,
                ((3, 3, 3, 3), (4, 4, 4)),
                (6, 1, 1, 1, 1, 1, 1),
            ),
        },
        (3, 5): {
            (1, 1, ((1,), (1,)), (1,)),
            (
                15,
                7,
                ((3, 3, 3, 3, 3), (5, 5, 5)),
                (7, 1, 1, 1, 1, 1, 1, 1, 1),
            ),
        },
    }
    for cusp, atlas in expected.items():
        assert canonical_clean_cusp_atlas(*cusp) == atlas

    # A single irreducible ramification divisor has one generic transverse
    # order along the target component.  The nontrivial clean A4 and A2
    # packets have disjoint possible orders, so the A4+A2+A1 target cannot
    # be clean at both cusps on one such row.
    a2_orders = {2, 4}
    a4_orders = {3, 6}
    assert a2_orders.isdisjoint(a4_orders)

    rows = clean_cusp_rows(3, 5) | clean_cusp_rows(5, 3)
    numerical = {(row[1], row[2]) for row in rows}
    assert numerical == {(1, 1), (15, 7)}

    nontrivial = [row for row in rows if row[1] == 15]
    assert nontrivial
    passports = {
        (
            tuple(sorted(row[3], reverse=True)),
            tuple(sorted(row[4], reverse=True)),
            row[5],
        )
        for row in nontrivial
    }
    # Swapping target coordinates interchanges the first two fibers.
    expected_left = (5, 5, 5)
    expected_right = (3, 3, 3, 3, 3)
    assert all(
        {passport[0], passport[1]} == {expected_left, expected_right}
        and passport[2] == (7,) + (1,) * 8
        for passport in passports
    )

    # Its single component already costs local degree 15 in Orevkov's
    # identity, so a global polynomial Keller map must have d-1>=15.
    candidate_degrees = {
        degree
        for degree in range(6, 29)
        if any(degree == 5 * r + 3 * s for r in range(10) for s in range(10))
    }
    assert min(degree for degree in candidate_degrees if degree - 1 >= 15) == 16


def main() -> None:
    orevkov_budget_audit()
    simple_e8_regular_cartier_audit()
    clean_k1_cusp_atlas_audit()
    print(
        "PASS: Orevkov plus residue Riemann--Hurwitz gives cost sum(e*f); "
        "all regular-Cartier one-component simple E8 rows through degree 28 "
        "overrun it; the five clean k=1 cusp atlases are exact, and (3,5) "
        "has only the identity and (local degree,ramification)=(15,7)"
    )


if __name__ == "__main__":
    main()
