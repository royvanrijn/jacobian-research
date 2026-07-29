#!/usr/bin/env python3
"""Exact frontier checks toward universal quartic-kernel saturation.

This script does not claim to settle the 24-parameter problem.  It performs
two exact calculations which go beyond the coordinate-axis/plane tests:

1. four dense lines, each involving all 24 kernel-basis tensors, for every
   squarefree cubic-symbol orbit;
2. the full coordinate subspace on the first eight kernel-basis tensors for
   the smooth cubic symbol.

On a line it checks the complete polynomial family over Q[t,x,y,z], not
just sampled fibers.  On the coordinate subspace it works over
Q[p0,...,p7,x,y,z].  Results are written to a generated JSON artifact with
the unresolved universal scope recorded explicitly.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import platform
import shutil
import subprocess
import time
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


def audit_line(
    task: tuple[str, int],
) -> dict[str, Any]:
    """Audit one dense line for one squarefree cubic symbol."""

    stratum, direction_index = task
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    started = time.monotonic()
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
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }


def singular_version() -> str:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "--version"],
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.splitlines()[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--skip-lines", action="store_true")
    parser.add_argument("--skip-eight-space", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    assert all(
        len(coefficients) == 24 and all(coefficients)
        for coefficients in DENSE_DIRECTIONS
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

    eight_space: dict[str, Any] | None = None
    if not args.skip_eight_space:
        started = time.monotonic()
        directions = cubic_audit.quartic_kernel_basis_tensors()
        result = cubic_audit.run_singular_subspace(
            cubic_audit.CUBIC_STRATA["smooth"],
            directions[:8],
            timeout=900,
        )
        assert result == EXPECTED_SUBSPACE_RESULT, result
        eight_space = {
            "stratum": "smooth",
            "basis_indices": list(range(8)),
            "result": {
                "cotangent_saturation_generators": result[0],
                "ext2_multiplicity": result[1],
                "parameter_space_radical_difference": result[2],
                "central_pruned_presentation_difference": result[3],
                "pruned_presentation_rank": result[4],
            },
            "elapsed_seconds": round(time.monotonic() - started, 3),
        }
        print("PASS smooth coordinate eight-space:", eight_space["result"])

    artifact = {
        "schema": "universal-cubic-quartic-kernel-saturation-frontier-v1",
        "mathematical_status": "exact finite-subspace computation",
        "universal_24_parameter_problem": "open",
        "construction": {
            "base_ring": "Q[u1,...,u24,x,y,z]",
            "family": "Phi_h + sum_(i=1)^24 u_i psi_i",
            "kernel_basis_dimension": 24,
            "squarefree_strata": sorted(cubic_audit.SQUAREFREE_STRATA),
        },
        "dense_directions": [list(item) for item in DENSE_DIRECTIONS],
        "dense_line_results": line_results,
        "smooth_coordinate_eight_space": eight_space,
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
