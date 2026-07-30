#!/usr/bin/env python3
"""Verify compatible order-5, order-6, and order-7 Ritt sections.

The order-seven tensor presentation is the single source of the finite
tower.  Its quotients by the fifth and sixth powers of the Dickson-base
maximal ideal give the lower levels, so the displayed sections are
compatible by construction rather than independently chosen.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_cellular_extension import (  # noqa: E402
    splitting_solution,
    vector_space_section,
)
from research_degree42_ritt_rotated_tensor_extension import (  # noqa: E402
    exact_rank_certificate,
    matrix_cache,
    rotated_tensor_audit,
)


WORDS = ((2, 3, 7), (3, 2, 7))
ORDERS = (5, 6, 7)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_inverse_limit_sections_q5_q7.json"
)


def matrix_from_strings(entries: list[list[str]]) -> sp.Matrix:
    """Parse a serialized exact rational matrix."""

    return sp.Matrix(
        [[sp.Rational(entry) for entry in row] for row in entries]
    )


def column_space(columns: list[sp.Matrix], ambient: int) -> sp.Matrix:
    """Return a basis matrix for the span of the supplied columns."""

    if not columns:
        return sp.zeros(ambient, 0)
    return sp.Matrix.hstack(*columns).columnspace()


def base_power_subspace(
    tau_action: sp.Matrix,
    zeta_action: sp.Matrix,
    order: int,
) -> sp.Matrix:
    """Return ``(tau,zeta)^order M`` inside the order-seven module."""

    ambient = tau_action.rows
    identity = sp.eye(ambient)
    columns: list[sp.Matrix] = []
    for tau_exponent in range(order + 1):
        operator = (
            tau_action**tau_exponent
            * zeta_action ** (order - tau_exponent)
        )
        columns.extend(operator * identity[:, index] for index in range(ambient))
    basis_columns = column_space(columns, ambient)
    if not basis_columns:
        return sp.zeros(ambient, 0)
    return sp.Matrix.hstack(*basis_columns)


def quotient_map(subspace: sp.Matrix) -> sp.Matrix:
    """Return a rational quotient map whose kernel is ``subspace``."""

    ambient = subspace.rows
    if subspace.cols == 0:
        return sp.eye(ambient)
    rows = subspace.T.nullspace()
    quotient = sp.Matrix.vstack(*(row.T for row in rows))
    assert quotient.rank() == quotient.rows
    assert quotient * subspace == sp.zeros(quotient.rows, subspace.cols)
    return quotient


def induced_map(
    source_quotient: sp.Matrix,
    target_quotient: sp.Matrix,
    ambient_map: sp.Matrix,
) -> sp.Matrix:
    """Descend an ambient map through two displayed quotient maps."""

    source_section = vector_space_section(source_quotient)
    result = target_quotient * ambient_map * source_section
    assert result * source_quotient == target_quotient * ambient_map
    return result


def homogeneous_sections(
    projection: sp.Matrix,
    total_actions: tuple[sp.Matrix, ...],
    spectator_actions: tuple[sp.Matrix, ...],
) -> tuple[sp.Matrix, ...]:
    """Return a basis of the vector space of section differences."""

    total_dimension = projection.cols
    spectator_dimension = projection.rows
    variables = sp.symbols(
        f"x0:{total_dimension * spectator_dimension}"
    )
    correction = sp.Matrix(
        total_dimension, spectator_dimension, variables
    )
    equations = list(projection * correction)
    for total_action, spectator_action in zip(
        total_actions, spectator_actions
    ):
        equations.extend(
            total_action * correction - correction * spectator_action
        )
    coefficient_matrix, _ = sp.linear_eq_to_matrix(equations, variables)
    return tuple(
        sp.Matrix(total_dimension, spectator_dimension, vector)
        for vector in coefficient_matrix.nullspace()
    )


def coordinates_in_basis(
    vector: sp.Matrix,
    basis: tuple[sp.Matrix, ...],
) -> sp.Matrix:
    """Write a flattened matrix in the supplied exact basis."""

    if not basis:
        assert vector.is_zero_matrix
        return sp.zeros(0, 1)
    basis_matrix = sp.Matrix.hstack(
        *(sp.Matrix(list(matrix)) for matrix in basis)
    )
    return basis_matrix.gauss_jordan_solve(sp.Matrix(list(vector)))[0]


def tower_for_word(word: tuple[int, int, int]) -> dict[str, object]:
    """Build the nested section tower from one order-seven presentation."""

    top = rotated_tensor_audit(word, base_order=7, normal_order=1)
    presentation = top["finite_module_presentation"]
    projection_top = matrix_from_strings(presentation["projection"])
    total_actions_top = tuple(
        matrix_from_strings(matrix)
        for matrix in presentation["total_actions"][-2:]
    )
    spectator_actions_top = tuple(
        matrix_from_strings(matrix)
        for matrix in presentation["spectator_actions"][-2:]
    )
    top_section = matrix_from_strings(top["one_section"])
    assert projection_top * top_section == sp.eye(projection_top.rows)
    assert all(
        total_action * top_section == top_section * spectator_action
        for total_action, spectator_action in zip(
            total_actions_top, spectator_actions_top
        )
    )

    levels: dict[int, dict[str, object]] = {}
    for order in ORDERS:
        total_power = base_power_subspace(
            *total_actions_top, order=order
        )
        spectator_power = base_power_subspace(
            *spectator_actions_top, order=order
        )
        total_quotient = quotient_map(total_power)
        spectator_quotient = quotient_map(spectator_power)
        total_actions = tuple(
            induced_map(total_quotient, total_quotient, action)
            for action in total_actions_top
        )
        spectator_actions = tuple(
            induced_map(spectator_quotient, spectator_quotient, action)
            for action in spectator_actions_top
        )
        projection = induced_map(
            total_quotient, spectator_quotient, projection_top
        )
        section = (
            total_quotient
            * top_section
            * vector_space_section(spectator_quotient)
        )
        assert section * spectator_quotient == (
            total_quotient * top_section
        )
        assert projection * section == sp.eye(projection.rows)
        assert all(
            total_action * section == section * spectator_action
            for total_action, spectator_action in zip(
                total_actions, spectator_actions
            )
        )
        splits, _ = splitting_solution(
            projection, total_actions, spectator_actions
        )
        assert splits
        kernel = sp.Matrix.hstack(*projection.nullspace())
        kernel_actions = tuple(
            sp.Matrix.hstack(
                *(
                    kernel.gauss_jordan_solve(
                        action * kernel[:, column]
                    )[0]
                    for column in range(kernel.cols)
                )
            )
            for action in total_actions
        )
        certificate = exact_rank_certificate(
            kernel,
            projection,
            kernel_actions,
            total_actions,
            spectator_actions,
            ("tau", "zeta"),
        )
        homogeneous = homogeneous_sections(
            projection, total_actions, spectator_actions
        )
        levels[order] = {
            "total_quotient": total_quotient,
            "spectator_quotient": spectator_quotient,
            "projection": projection,
            "total_actions": total_actions,
            "spectator_actions": spectator_actions,
            "section": section,
            "homogeneous": homogeneous,
            "dimensions": {
                "kernel": kernel.cols,
                "total": projection.cols,
                "spectator": projection.rows,
            },
            "coboundary_rank": certificate["coboundary_rank"],
            "augmented_rank": certificate["augmented_rank"],
            "section_space_dimension": len(homogeneous),
        }

    transitions = []
    for lower, upper in zip(ORDERS, ORDERS[1:]):
        lower_level = levels[lower]
        upper_level = levels[upper]
        total_transition = induced_map(
            upper_level["total_quotient"],
            lower_level["total_quotient"],
            sp.eye(projection_top.cols),
        )
        spectator_transition = induced_map(
            upper_level["spectator_quotient"],
            lower_level["spectator_quotient"],
            sp.eye(projection_top.rows),
        )
        assert (
            total_transition * upper_level["section"]
            == lower_level["section"] * spectator_transition
        )
        columns = []
        for correction in upper_level["homogeneous"]:
            reduced = (
                total_transition
                * correction
                * vector_space_section(spectator_transition)
            )
            assert reduced * spectator_transition == (
                total_transition * correction
            )
            columns.append(
                coordinates_in_basis(reduced, lower_level["homogeneous"])
            )
        restriction = (
            sp.Matrix.hstack(*columns)
            if columns
            else sp.zeros(len(lower_level["homogeneous"]), 0)
        )
        transitions.append(
            {
                "from_order": upper,
                "to_order": lower,
                "compatible_section_residual_rank": 0,
                "homogeneous_source_dimension": len(
                    upper_level["homogeneous"]
                ),
                "homogeneous_target_dimension": len(
                    lower_level["homogeneous"]
                ),
                "homogeneous_restriction_rank": restriction.rank(),
                "homogeneous_restriction_cokernel_dimension": (
                    len(lower_level["homogeneous"]) - restriction.rank()
                ),
            }
        )

    return {
        "word": word,
        "thick_composite_omission": top["thick_composite_omission"],
        "levels": [
            {
                "base_order": order,
                "dimensions": levels[order]["dimensions"],
                "splits_as_B_module": True,
                "coboundary_rank": levels[order]["coboundary_rank"],
                "augmented_rank": levels[order]["augmented_rank"],
                "section_space_dimension": levels[order][
                    "section_space_dimension"
                ],
            }
            for order in ORDERS
        ],
        "transitions": transitions,
    }


def main() -> None:
    sectors = [tower_for_word(word) for word in WORDS]
    caches = [matrix_cache(word, 7, 1) for word in WORDS]
    output = {
        "schema": "degree42-ritt-inverse-limit-sections-q5-q7.v1",
        "status": "exact compatible finite section tower",
        "base_orders": ORDERS,
        "normal_order": 1,
        "matrix_caches": [
            {
                "path": str(cache.relative_to(ROOT)),
                "sha256": hashlib.sha256(cache.read_bytes()).hexdigest(),
            }
            for cache in caches
        ],
        "sectors": sectors,
        "consequence": (
            "Both rotated first-conormal extensions have sections at base "
            "orders five, six, and seven that commute exactly with "
            "truncation. The finite affine inverse-limit obstruction "
            "through order seven is zero."
        ),
        "theorem_boundary": (
            "This finite tower alone does not prove that the completed "
            "extension splits; an all-order Mittag-Leffler, Ext, or braid-"
            "restriction calculation is still required."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_inverse_limit_sections.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    for sector in sectors:
        omission = sector["thick_composite_omission"]
        dimensions = [
            level["dimensions"] for level in sector["levels"]
        ]
        print(f"PASS: cut-{omission} compatible dimensions {dimensions}")
    print("PASS: both finite inverse-limit obstructions vanish through q7")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
