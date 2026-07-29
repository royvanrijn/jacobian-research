#!/usr/bin/env python3
"""Exact frontier checks toward universal quartic-kernel saturation.

This script does not claim to settle the 24-parameter problem.  It performs
three exact calculations which go beyond the coordinate-axis/plane tests:

1. a universal cotangent-input calculation: the parameter/collision
   bidegrees and six unit pivots reducing the 12-by-31 presentation to a
   cokernel-equivalent 6-by-25 matrix;
2. four dense lines, each involving all 24 kernel-basis tensors, for every
   squarefree cubic-symbol orbit;
3. the full coordinate subspace on the first ten kernel-basis tensors for
   the smooth cubic symbol.

On a line it checks the complete polynomial family over Q[t,x,y,z], not
just sampled fibers.  On the coordinate subspace it works over
Q[p0,...,p9,x,y,z].  Results are written to a generated JSON artifact with
the unresolved universal scope recorded explicitly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Any

import verify_cubic_symbol_double_saturation as cubic_audit


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "universal_cubic_quartic_kernel_saturation_frontier.json"
)

# Every vector has full support.  The last two are deterministic signed
# small-height vectors, retained literally so the certificate is insensitive
# to changes in a random-number generator.
DENSE_DIRECTIONS = (
    tuple([1] * 24),
    tuple(range(1, 25)),
    (
        2,
        2,
        2,
        -1,
        1,
        -1,
        -1,
        -1,
        -1,
        1,
        -1,
        2,
        1,
        1,
        1,
        -1,
        1,
        -2,
        -2,
        2,
        -2,
        2,
        1,
        1,
    ),
    (
        2,
        -1,
        2,
        2,
        1,
        -2,
        2,
        -2,
        -2,
        -2,
        -1,
        2,
        1,
        -1,
        1,
        2,
        1,
        -2,
        2,
        1,
        2,
        1,
        1,
        -2,
    ),
)

EXPECTED_LINE_RESULT = (0, 6, 0, 0, 0)
EXPECTED_SUBSPACE_RESULT = (0, 6, 0, 0, 3)


def combined_tensor(
    coefficients: tuple[int, ...],
) -> dict[tuple[int, int, int], cubic_audit.sp.Expr]:
    """Return the indicated integral combination of the fixed kernel basis."""

    directions = cubic_audit.quartic_kernel_basis_tensors()
    assert len(coefficients) == len(directions) == 24
    tensor = {
        triple: cubic_audit.sp.expand(
            sum(
                coefficient * direction[triple]
                for coefficient, direction in zip(coefficients, directions)
            )
        )
        for triple in directions[0]
    }
    for pair in cubic_audit.itertools.combinations_with_replacement(
        range(3), 2
    ):
        assert (
            cubic_audit.sp.expand(
                cubic_audit.z
                * tensor[tuple(sorted((0, *pair)))]
                - cubic_audit.y
                * tensor[tuple(sorted((1, *pair)))]
                + cubic_audit.x
                * tensor[tuple(sorted((2, *pair)))]
            )
            == 0
        )
    return tensor


def universal_tensor() -> tuple[
    tuple[cubic_audit.sp.Symbol, ...],
    dict[tuple[int, int, int], cubic_audit.sp.Expr],
]:
    """Construct sum u_i psi_i over Q[u_1,...,u_24,x,y,z]."""

    parameters = cubic_audit.sp.symbols("u1:25")
    directions = cubic_audit.quartic_kernel_basis_tensors()
    tensor = {
        triple: cubic_audit.sp.expand(
            sum(
                parameter * direction[triple]
                for parameter, direction in zip(parameters, directions)
            )
        )
        for triple in directions[0]
    }
    for pair in cubic_audit.itertools.combinations_with_replacement(
        range(3), 2
    ):
        assert (
            cubic_audit.sp.expand(
                cubic_audit.z
                * tensor[tuple(sorted((0, *pair)))]
                - cubic_audit.y
                * tensor[tuple(sorted((1, *pair)))]
                + cubic_audit.x
                * tensor[tuple(sorted((2, *pair)))]
            )
            == 0
        )
    return parameters, tensor


def unit_pruned_differential_relations(
    relation_columns: list[list[cubic_audit.sp.Expr]],
) -> tuple[list[list[cubic_audit.sp.Expr]], list[tuple[int, int, str]]]:
    """Remove parameter-independent unit pivots before standard bases.

    Rows are free generators and columns are relations.  Each pass uses
    elementary row and column operations over Q[u_1,...,u_24,x,y,z], then
    deletes the resulting one-by-one identity summand.  Hence the returned
    matrix presents the same cokernel as the original differential matrix.
    """

    matrix = [
        [
            relation_columns[column][row]
            for column in range(len(relation_columns))
        ]
        for row in range(len(relation_columns[0]))
    ]
    pivots: list[tuple[int, int, str]] = []
    while True:
        pivot: tuple[int, int, cubic_audit.sp.Rational] | None = None
        for row, entries in enumerate(matrix):
            for column, entry in enumerate(entries):
                if entry != 0 and not entry.free_symbols:
                    pivot = (
                        row,
                        column,
                        cubic_audit.sp.Rational(entry),
                    )
                    break
            if pivot is not None:
                break
        if pivot is None:
            break

        pivot_row, pivot_column, pivot_value = pivot
        pivots.append((pivot_row, pivot_column, str(pivot_value)))
        pivot_entries = matrix[pivot_row]

        # Clear the pivot column by elementary row operations.
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][pivot_column] == 0:
                continue
            factor = matrix[row][pivot_column] / pivot_value
            matrix[row] = [
                cubic_audit.sp.expand(entry - factor * pivot_entry)
                for entry, pivot_entry in zip(
                    matrix[row], pivot_entries
                )
            ]

        # The pivot column now has only its unit entry.  Clear the pivot
        # row by elementary column operations.
        for column in range(len(matrix[0])):
            if (
                column == pivot_column
                or matrix[pivot_row][column] == 0
            ):
                continue
            factor = matrix[pivot_row][column] / pivot_value
            for row in range(len(matrix)):
                matrix[row][column] = cubic_audit.sp.expand(
                    matrix[row][column]
                    - factor * matrix[row][pivot_column]
                )

        matrix = [
            entries[:pivot_column] + entries[pivot_column + 1 :]
            for row, entries in enumerate(matrix)
            if row != pivot_row
        ]

    columns = [
        [matrix[row][column] for row in range(len(matrix))]
        for column in range(len(matrix[0]))
    ]
    return columns, pivots


def universal_cotangent_input_reduction(
    parameters: tuple[cubic_audit.sp.Symbol, ...],
    tensor: dict[tuple[int, int, int], cubic_audit.sp.Expr],
) -> dict[str, Any]:
    """Certify the universal low-jet bound and the unit-pruned input."""

    cubic = cubic_audit.CUBIC_STRATA["smooth"]
    universal_relations = cubic_audit.differential_relations(cubic, tensor)
    central_relations = cubic_audit.differential_relations(cubic)
    variables = (*parameters, *cubic_audit.BASE_VARIABLES)
    bidegree_counts: dict[tuple[int, int], int] = {}
    nonzero_difference_entries = 0
    for universal_relation, central_relation in zip(
        universal_relations, central_relations
    ):
        for universal_entry, central_entry in zip(
            universal_relation, central_relation
        ):
            difference = cubic_audit.sp.expand(
                universal_entry - central_entry
            )
            if difference == 0:
                continue
            nonzero_difference_entries += 1
            polynomial = cubic_audit.sp.Poly(difference, *variables)
            for monomial, _coefficient in polynomial.terms():
                parameter_degree = sum(monomial[: len(parameters)])
                collision_degree = sum(monomial[len(parameters) :])
                assert 1 <= parameter_degree <= 2
                assert collision_degree >= 3
                bidegree = (parameter_degree, collision_degree)
                bidegree_counts[bidegree] = (
                    bidegree_counts.get(bidegree, 0) + 1
                )

    pruned_columns, pivots = unit_pruned_differential_relations(
        universal_relations
    )
    assert len(universal_relations) == 31
    assert len(universal_relations[0]) == 12
    assert len(pruned_columns) == 25
    assert len(pruned_columns[0]) == 6
    assert len(pivots) == 6
    assert all(value in {"1", "2"} for _, _, value in pivots)
    serialized_columns = [
        [cubic_audit.sp.sstr(entry) for entry in column]
        for column in pruned_columns
    ]
    pruned_sha256 = hashlib.sha256(
        json.dumps(
            serialized_columns,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    return {
        "stratum_used_for_explicit_matrix": "smooth",
        "raw_presentation": {"rows": 12, "columns": 31},
        "unit_pruned_presentation": {"rows": 6, "columns": 25},
        "unit_pivots": [
            {"row": row, "column": column, "value": value}
            for row, column, value in pivots
        ],
        "nonzero_parameter_difference_entries": (
            nonzero_difference_entries
        ),
        "parameter_collision_bidegree_term_counts": [
            {
                "parameter_degree": parameter_degree,
                "collision_degree": collision_degree,
                "term_count": count,
            }
            for (parameter_degree, collision_degree), count in sorted(
                bidegree_counts.items()
            )
        ],
        "consequence": (
            "the universal and central cotangent presentations agree "
            "modulo (x,y,z)^3"
        ),
        "unit_pruned_matrix_sha256": pruned_sha256,
    }


def audit_line(
    task: tuple[str, int],
) -> dict[str, Any]:
    """Audit one dense line for one squarefree cubic symbol."""

    stratum, direction_index = task
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    result = cubic_audit.run_singular_family(
        cubic_audit.CUBIC_STRATA[stratum],
        combined_tensor(DENSE_DIRECTIONS[direction_index]),
        timeout=600,
    )
    assert result == EXPECTED_LINE_RESULT, (task, result)
    return {
        "stratum": stratum,
        "direction_index": direction_index,
        "result": {
            "cotangent_saturation_generators": result[0],
            "ext2_multiplicity": result[1],
            "parameter_torsion_generators": result[2],
            "collision_axis_radical_difference": result[3],
            "central_presentation_difference": result[4],
        },
    }


def singular_version() -> str:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-v"],
        input="",
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.splitlines()[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--skip-lines", action="store_true")
    parser.add_argument("--skip-ten-space", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    assert all(
        len(coefficients) == 24 and all(coefficients)
        for coefficients in DENSE_DIRECTIONS
    )
    universal_parameters, universal = universal_tensor()
    tensor_components = [
        {
            "triple": list(triple),
            "polynomial": cubic_audit.sp.sstr(universal[triple]),
        }
        for triple in sorted(universal)
    ]
    tensor_sha256 = hashlib.sha256(
        json.dumps(
            tensor_components,
            separators=(",", ":"),
            sort_keys=True,
        ).encode()
    ).hexdigest()
    cotangent_input_reduction = universal_cotangent_input_reduction(
        universal_parameters, universal
    )
    print(
        "PASS universal cotangent input reduction:",
        cotangent_input_reduction["unit_pruned_presentation"],
        cotangent_input_reduction["consequence"],
    )

    previous_artifact: dict[str, Any] = {}
    if args.output.exists():
        previous_artifact = json.loads(
            args.output.read_text(encoding="utf-8")
        )

    line_results: list[dict[str, Any]] = []
    if not args.skip_lines:
        tasks = [
            (stratum, direction_index)
            for stratum in sorted(cubic_audit.SQUAREFREE_STRATA)
            for direction_index in range(len(DENSE_DIRECTIONS))
        ]
        context = multiprocessing.get_context("fork")
        with context.Pool(processes=args.workers) as pool:
            for result in pool.imap_unordered(audit_line, tasks):
                line_results.append(result)
                print(
                    "PASS dense line:",
                    result["stratum"],
                    result["direction_index"],
                    result["result"],
                    flush=True,
                )
        line_results.sort(
            key=lambda item: (item["stratum"], item["direction_index"])
        )
    else:
        line_results = previous_artifact.get("dense_line_results", [])

    ten_space: dict[str, Any] | None = None
    if not args.skip_ten_space:
        directions = cubic_audit.quartic_kernel_basis_tensors()
        result = cubic_audit.run_singular_subspace(
            cubic_audit.CUBIC_STRATA["smooth"],
            directions[:10],
            timeout=900,
        )
        assert result == EXPECTED_SUBSPACE_RESULT, result
        ten_space = {
            "stratum": "smooth",
            "basis_indices": list(range(10)),
            "result": {
                "cotangent_saturation_generators": result[0],
                "ext2_multiplicity": result[1],
                "parameter_space_radical_difference": result[2],
                "central_pruned_presentation_difference": result[3],
                "pruned_presentation_rank": result[4],
            },
        }
        print("PASS smooth coordinate ten-space:", ten_space["result"])
    else:
        ten_space = previous_artifact.get(
            "smooth_coordinate_ten_space"
        )

    artifact = {
        "schema": "universal-cubic-quartic-kernel-saturation-frontier-v3",
        "mathematical_status": (
            "exact structural input and finite-subspace computation"
        ),
        "universal_24_parameter_problem": "open",
        "construction": {
            "base_ring": "Q[u1,...,u24,x,y,z]",
            "family": "Phi_h + sum_(i=1)^24 u_i psi_i",
            "parameters": [str(parameter) for parameter in universal_parameters],
            "kernel_basis_dimension": 24,
            "squarefree_strata": sorted(cubic_audit.SQUAREFREE_STRATA),
            "universal_tensor_components": tensor_components,
            "universal_tensor_sha256": tensor_sha256,
        },
        "universal_cotangent_input_reduction": (
            cotangent_input_reduction
        ),
        "dense_directions": [list(item) for item in DENSE_DIRECTIONS],
        "dense_line_results": line_results,
        "smooth_coordinate_ten_space": ten_space,
        "proved_on_each_recorded_family": {
            "relative_cotangent_saturation": True,
            "ext2_relative_multiplicity": 6,
            "ext2_parameter_torsion": False,
            "ext2_support": "parameter space at x=y=z=0",
            "ext2_fitting_consequence": "Fitt_6=(1), Fitt_5=(0)",
        },
        "not_computed": [
            "the full 24-parameter Ext presentation",
            "the full parameter discriminant",
            "all higher-support parameter directions",
            "normal nonhomogeneous lifts",
            "Keller-open compatibility",
        ],
        "software": {
            "python": platform.python_version(),
            "sympy": cubic_audit.sp.__version__,
            "singular": singular_version(),
        },
        "reproduce": (
            ".venv/bin/python "
            "scripts/research_universal_cubic_quartic_kernel_saturation.py"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"PASS wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
