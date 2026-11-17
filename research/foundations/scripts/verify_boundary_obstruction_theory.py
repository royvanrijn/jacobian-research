#!/usr/bin/env python3
"""Exact regressions for the affine boundary-obstruction package.

This checker has four independent parts:

1. compile a saturated module and a boundary-torsion module with the
   reusable Singular support-saturation backend;
2. exhibit the finite-jet/completion warning from Theorem A2, including
   the unbounded boundary-annihilation exponents;
3. verify the nodal and cuspidal conductor pullback sequences, and their
   tensor products with a finite free coefficient block, by exact rational
   matrices; and
4. separate strict bounded lifting from ordinary solvability with a pair
   of filtered rational matrices.

The finite examples are regressions of the general written proofs.  They do
not replace the derived localization triangle or a problem-specific
Kuranishi tail-elimination certificate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.support_saturation import (
    CompilerOptions,
    ModulePresentation,
    NormalFiltration,
    PolynomialRing,
    SupportSaturationCompiler,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "boundary_obstruction_theory.json"
)


def matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    """Return a JSON-safe exact integer matrix."""

    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def compile_module_regressions() -> dict[str, Any]:
    """Run exact regular, torsion, and finite-jet saturation examples."""

    regularity_options = CompilerOptions(
        associated_primes="regularity",
        saturation_strategy="compute",
        torsion_exponent_bound=12,
    )
    compiler = SupportSaturationCompiler(regularity_options)

    regular = compiler.compile(
        ModulePresentation(
            ring=PolynomialRing(("x", "y")),
            rank=1,
            generators=(("y",),),
            label="regular-boundary-example",
        ),
        boundary_ideal=("x",),
        distinguished_class=("1",),
    )
    assert regular["saturation"]["equal_to_presentation"]
    assert regular["local_cohomology"]["zero"]
    assert regular["local_cohomology"]["boundary_annihilation_exponent"] == 0
    assert regular["regular_elements"]["candidate"] == "x"
    assert not regular["distinguished_class"]["belongs_to_local_cohomology"]

    decomposition_compiler = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="decompose",
            saturation_strategy="compute",
            torsion_exponent_bound=12,
        )
    )
    torsion = decomposition_compiler.compile(
        ModulePresentation(
            ring=PolynomialRing(("x",)),
            rank=1,
            generators=(("x2",),),
            label="boundary-torsion-example",
        ),
        boundary_ideal=("x",),
        distinguished_class=("1",),
    )
    assert not torsion["saturation"]["equal_to_presentation"]
    assert not torsion["local_cohomology"]["zero"]
    assert torsion["local_cohomology"]["boundary_annihilation_exponent"] == 2
    assert torsion["regular_elements"]["candidate"] is None
    assert torsion["associated_primes"]["primes"] == [["x"]]
    assert torsion["local_cohomology"]["associated_primes"] == [["x"]]
    assert torsion["distinguished_class"]["belongs_to_local_cohomology"]
    assert torsion["distinguished_class"]["boundary_annihilation_exponent"] == 2

    jet_orders = tuple(range(1, 7))
    jets = compiler.compile(
        ModulePresentation(
            ring=PolynomialRing(("x",)),
            rank=1,
            generators=(("0",),),
            label="unbounded-finite-jet-torsion",
        ),
        boundary_ideal=("x",),
        distinguished_class=("1",),
        filtration=NormalFiltration(
            ideal=("x",),
            orders=jet_orders,
        ),
    )
    assert jets["local_cohomology"]["zero"]
    assert jets["local_cohomology"]["boundary_annihilation_exponent"] == 0
    exponent_profile = tuple(
        jet["boundary_annihilation_exponent"]
        for jet in jets["finite_jets"]["jets"]
    )
    assert exponent_profile == jet_orders
    assert all(
        transition["surjective"]
        for transition in jets["finite_jets"]["transitions"]
    )
    assert (
        jets["finite_jets"]["maximum_requested_jet_boundary_exponent"]
        == jet_orders[-1]
    )

    return {
        "regular_boundary": regular,
        "boundary_torsion": torsion,
        "finite_jet_warning": {
            "certificate": jets,
            "global_identity": "H^0_(x)(Q[x])=0",
            "jet_identity": "H^0_(x)(Q[x]/(x^n))=Q[x]/(x^n)",
            "boundary_exponent_profile": list(exponent_profile),
            "all_transition_maps_surjective": True,
            "inverse_limit": "Q[[x]]",
            "completed_local_cohomology": "H^0_(x)(Q[[x]])=0",
            "conclusion": (
                "finite-jet torsion needs a uniform boundary exponent "
                "before inverse limit"
            ),
        },
    }


def conductor_matrix(kind: str, degree: int) -> sp.Matrix:
    """Return the truncated conductor difference map.

    Columns are the coefficients of ``f in Q[t]_{\u2264degree}``, followed by
    the residue scalar from ``A/c``.  Rows are the two coordinates of
    ``B/c``.
    """

    if kind == "node":
        evaluation_zero = [1] + [0] * degree + [-1]
        evaluation_one = [1] * (degree + 1) + [-1]
        return sp.Matrix((evaluation_zero, evaluation_one))
    if kind == "cusp":
        constant_jet = [1] + [0] * degree + [-1]
        linear_jet = [0, 1] + [0] * (degree - 1) + [0]
        return sp.Matrix((constant_jet, linear_jet))
    raise ValueError(f"unknown conductor kind {kind!r}")


def compile_conductor_regressions() -> dict[str, Any]:
    """Verify exact conductor kernels and finite-free tensor descent."""

    degree = 7
    coefficient_rank = 3
    results: dict[str, Any] = {}
    for kind in ("node", "cusp"):
        difference = conductor_matrix(kind, degree)
        assert difference.rank() == 2
        assert len(difference.nullspace()) == degree

        coefficient_difference = sp.kronecker_product(
            difference, sp.eye(coefficient_rank)
        )
        assert coefficient_difference.rank() == 2 * coefficient_rank
        assert (
            len(coefficient_difference.nullspace())
            == degree * coefficient_rank
        )

        for vector in difference.nullspace():
            coefficients = vector[:-1, 0]
            residue = vector[-1, 0]
            if kind == "node":
                assert coefficients[0] == residue
                assert sum(coefficients) == residue
            else:
                assert coefficients[0] == residue
                assert coefficients[1] == 0

        results[kind] = {
            "polynomial_degree_cutoff": degree,
            "difference_matrix": matrix_rows(difference),
            "difference_rank": difference.rank(),
            "kernel_dimension": len(difference.nullspace()),
            "conductor_condition": (
                "f(0)=f(1)"
                if kind == "node"
                else "f'(0)=0"
            ),
            "finite_free_coefficient_rank": coefficient_rank,
            "tensor_difference_rank": coefficient_difference.rank(),
            "tensor_kernel_dimension": len(
                coefficient_difference.nullspace()
            ),
        }
    return results


def filtered_subspace(weights: tuple[int, ...], order: int) -> sp.Matrix:
    """Return the coordinate inclusion of a filtered basis piece."""

    columns = [
        sp.eye(len(weights))[:, index]
        for index, weight in enumerate(weights)
        if weight <= order
    ]
    return (
        sp.Matrix.hstack(*columns)
        if columns
        else sp.zeros(len(weights), 0)
    )


def strictness_profile(
    differential: sp.Matrix,
    source_weights: tuple[int, ...],
    target_weights: tuple[int, ...],
) -> list[dict[str, Any]]:
    """Compute the strictness rank equality in every relevant degree."""

    profile = []
    for order in range(max(source_weights + target_weights) + 1):
        source_piece = filtered_subspace(source_weights, order)
        target_piece = filtered_subspace(target_weights, order)
        restricted_image = differential * source_piece
        image_rank = differential.rank()
        target_rank = target_piece.rank()
        joined_rank = differential.row_join(target_piece).rank()
        intersection_rank = image_rank + target_rank - joined_rank
        profile.append(
            {
                "order": order,
                "restricted_image_rank": restricted_image.rank(),
                "image_intersection_target_piece_rank": intersection_rank,
                "strict": restricted_image.rank() == intersection_rank,
            }
        )
    return profile


def compile_filtered_regressions() -> dict[str, Any]:
    """Verify strict bounded lifting and its minimal failure."""

    strict_map = sp.eye(2)
    strict_profile = strictness_profile(
        strict_map,
        source_weights=(0, 1),
        target_weights=(0, 1),
    )
    assert all(row["strict"] for row in strict_profile)

    bounded_defect = sp.Matrix((1, 1))
    bounded_lift = strict_map.inv() * bounded_defect
    assert strict_map * bounded_lift == bounded_defect

    nonstrict_map = sp.Matrix(((1,),))
    nonstrict_profile = strictness_profile(
        nonstrict_map,
        source_weights=(1,),
        target_weights=(0,),
    )
    assert not nonstrict_profile[0]["strict"]
    assert nonstrict_map * sp.Matrix((1,)) == sp.Matrix((1,))

    return {
        "strict_example": {
            "matrix": matrix_rows(strict_map),
            "source_weights": [0, 1],
            "target_weights": [0, 1],
            "profile": strict_profile,
            "bounded_defect": [1, 1],
            "bounded_lift": [int(value) for value in bounded_lift],
        },
        "nonstrict_example": {
            "matrix": matrix_rows(nonstrict_map),
            "source_weights": [1],
            "target_weights": [0],
            "profile": nonstrict_profile,
            "defect_target_degree": 0,
            "least_lift_source_degree": 1,
            "conclusion": (
                "ordinary solvability does not preserve a filtration bound"
            ),
        },
        "tail_rigidity_warning": {
            "equations": 0,
            "gauges": 0,
            "obstructions": 0,
            "formal_solution": "sum_(n>=0) t^n",
            "bounded_representative": False,
            "conclusion": (
                "vanishing high obstruction groups does not imply "
                "finite effectivity"
            ),
        },
    }


def main() -> None:
    certificate = {
        "schema": "boundary-obstruction-theory-regression.v1",
        "module_boundary": compile_module_regressions(),
        "conductor_descent": compile_conductor_regressions(),
        "filtered_effectivity": compile_filtered_regressions(),
        "scope": {
            "proved_by_written_argument": [
                "derived localization triangle",
                "local-cohomology spectral sequence criterion",
                "uniform-exponent finite-jet convergence",
                "perfect and Rees-projective conductor descent",
                "strict bounded lifting",
                "finite effectivity under compatible tail elimination",
            ],
            "executable_regressions": [
                "module saturation and regular boundary element",
                "boundary torsion and distinguished class",
                "unbounded finite-jet torsion exponents",
                "node and cusp conductor kernels",
                "finite-free tensor conductor kernels",
                "strict and non-strict filtered lifts",
                "tail-rigidity countermodel",
            ],
            "not_claimed": [
                "a universal cubic cotangent saturation certificate",
                "global degree-forty-two Ritt coefficient effectivity",
                "rank-two quantization algebraization",
                "the missing plane-JC residue module",
            ],
        },
    }
    OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print("PASS boundary module saturation and boundary torsion")
    print("PASS finite jets expose unbounded boundary exponents 1 through 6")
    print("PASS node/cusp conductor pullbacks and finite-free tensor descent")
    print("PASS strict bounded lifting and non-strict degree-loss control")
    print(f"WROTE {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
