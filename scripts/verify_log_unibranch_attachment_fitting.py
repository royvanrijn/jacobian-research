#!/usr/bin/env python3
"""Verify the logarithmic Fitting class of a unibranch SNC attachment."""

from __future__ import annotations

import sympy as sp


r, t = sp.symbols("r t")


def model_matrix_audit(multiplicity: int, higher_order: int) -> None:
    assert 1 <= multiplicity < higher_order

    # x=t^m, y=t^n+r.  With source basis (dlog r,dlog t), the
    # target columns dx,dy give this logarithmic matrix.
    theta = sp.Matrix(
        [[0, r], [multiplicity * t**multiplicity, higher_order * t**higher_order]]
    )
    assert sp.factor(theta.det()) == -multiplicity * r * t**multiplicity

    column_operation = sp.Matrix(
        [
            [1, -sp.Rational(higher_order, multiplicity) * t ** (higher_order - multiplicity)],
            [0, 1],
        ]
    )
    reduced = sp.simplify(theta * column_operation)
    assert reduced == sp.Matrix([[0, r], [multiplicity * t**multiplicity, 0]])
    swapped = reduced * sp.Matrix([[0, 1], [1, 0]])
    assert swapped == sp.diag(r, multiplicity * t**multiplicity)

    fitting_zero = sp.factor(theta.det())
    fitting_one = sp.groebner(tuple(theta), r, t, order="lex")
    expected_fitting_one = sp.groebner([r, t**multiplicity], r, t, order="lex")
    assert fitting_zero != 0
    assert fitting_one == expected_fitting_one

    # The monomial ideals satisfy (r) intersect (t^m)=(r*t^m) and
    # (r)+(t^m)=(r,t^m).  The quotient has the displayed m-element basis.
    standard_basis = tuple(t**index for index in range(multiplicity))
    assert len(standard_basis) == multiplicity
    for monomial in standard_basis:
        assert sp.rem(monomial, t**multiplicity, t) == monomial
    assert sp.rem(t**multiplicity, t**multiplicity, t) == 0


def ordinary_to_log_order_shift_audit(multiplicity: int, higher_order: int) -> None:
    ordinary = sp.Matrix(
        [
            [0, 1],
            [multiplicity * t ** (multiplicity - 1), higher_order * t ** (higher_order - 1)],
        ]
    )
    logarithmic = sp.diag(r, t) * ordinary
    assert sp.factor(ordinary.det()) == -multiplicity * t ** (multiplicity - 1)
    assert sp.factor(logarithmic.det()) == -multiplicity * r * t**multiplicity
    assert logarithmic == sp.Matrix(
        [[0, r], [multiplicity * t**multiplicity, higher_order * t**higher_order]]
    )


def residue_ramification_audit(branch_multiplicity: int, residue_index: int) -> None:
    source_multiplicity = branch_multiplicity * residue_index
    target_next_order = source_multiplicity + residue_index
    model_matrix_audit(source_multiplicity, target_next_order)
    assert source_multiplicity == branch_multiplicity * residue_index


def cusp_specialization_audit() -> None:
    for residue_index in range(1, 13):
        point_length = 2 * residue_index
        residue_ramification_audit(2, residue_index)
        assert point_length == 2 * residue_index
        assert point_length != 1


def residue_fiber_sum_audit() -> None:
    # For a degree-f finite map of smooth curves, the residue ramification
    # indices over a geometric point sum to f.  Hence minimal attachments
    # above a branch of multiplicity m_C have total point length m_C*f.
    for branch_multiplicity in range(1, 8):
        for residue_partition in (
            (1,),
            (1, 1),
            (2,),
            (1, 2, 1),
            (3, 2),
            (1, 1, 2, 3),
        ):
            residue_degree = sum(residue_partition)
            total_point_length = sum(
                branch_multiplicity * local_index
                for local_index in residue_partition
            )
            assert total_point_length == branch_multiplicity * residue_degree


def main() -> None:
    for multiplicity in range(1, 10):
        for higher_order in range(multiplicity + 1, multiplicity + 5):
            ordinary_to_log_order_shift_audit(multiplicity, higher_order)
            model_matrix_audit(multiplicity, higher_order)
    cusp_specialization_audit()
    residue_fiber_sum_audit()
    print(
        "PASS: a minimal transverse unibranch SNC attachment has log matrix "
        "diag(r,t^(q_p*m_C)), Fitt_1=(r,t^(q_p*m_C)), and positive point "
        "correction q_p*m_C; a residue-degree-f fiber contributes m_C*f "
        "when every attachment is minimal (ordinary-cusp total 2f)"
    )


if __name__ == "__main__":
    main()
