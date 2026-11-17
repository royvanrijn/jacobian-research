#!/usr/bin/env python3
"""Verify the reusable finite-module Postnikov cellular formalism."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.hessian_ritt_cellular import (  # noqa: E402
    EquivariantModuleSurjection,
    FiniteModuleRepresentation,
    PostnikovModuleTower,
    postnikov_braid_totalization,
)
from jcsearch.ritt_complex import MoveType, symmetric_braid_complex  # noqa: E402


CONORMAL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_conormal_transitivity.json"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "cellular_postnikov_transitivity.json"
)


def matrix(entries: list[list[str]]) -> sp.Matrix:
    """Parse an exact string matrix from a generated artifact."""

    return sp.Matrix(
        [[sp.Rational(entry) for entry in row] for row in entries]
    )


def zero_module(
    variable_names: tuple[str, ...],
    dimension: int,
    name: str,
) -> FiniteModuleRepresentation:
    """Return a module with trivial coordinate actions."""

    return FiniteModuleRepresentation(
        variable_names,
        tuple(sp.zeros(dimension, dimension) for _ in variable_names),
        name,
    )


def synthetic_multiflag_audit() -> dict[str, object]:
    """Check that arbitrary-length towers extract exact kernel layers."""

    variables = ("x",)
    modules = tuple(
        zero_module(variables, dimension, f"synthetic-N{index}")
        for index, dimension in enumerate((3, 2, 1, 0))
    )
    projections = (
        sp.Matrix([[1, 0, 0], [0, 1, 0]]),
        sp.Matrix([[1, 0]]),
        sp.zeros(0, 1),
    )
    maps = tuple(
        EquivariantModuleSurjection(
            modules[index],
            modules[index + 1],
            projection,
            f"synthetic-stage-{index}",
        )
        for index, projection in enumerate(projections)
    )
    tower = PostnikovModuleTower(maps, "synthetic three-layer flag")
    assert tower.module_dimensions == (3, 2, 1, 0)
    assert tower.layer_dimensions == (1, 1, 1)
    assert tower.split_profile == (True, True, True)
    return {
        "module_dimensions": tower.module_dimensions,
        "layer_dimensions": tower.layer_dimensions,
        "split_profile": tower.split_profile,
    }


def degree_30_audit() -> dict[str, object]:
    """Model the sector-only degeneration ``I_boundary=K``."""

    variables = ("tau", "zeta")
    sector = zero_module(variables, 1, "degree-30 sector conormal fiber")
    terminal = FiniteModuleRepresentation.zero(
        variables, "degree-30 zero spectator"
    )
    map_ = EquivariantModuleSurjection(
        sector,
        terminal,
        sp.zeros(0, 1),
        "degree-30 terminal projection",
    )
    tower = PostnikovModuleTower((map_,), "degree-30 sector-only tower")
    braid = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    totalization = postnikov_braid_totalization(
        braid,
        base_dimension=2,
        tower=tower,
        layer_names=("sector",),
        name="degree-30 Postnikov associated graded",
    )
    assert tower.module_dimensions == (1, 0)
    assert tower.layer_dimensions == (1,)
    assert totalization.complex.cohomology_dimensions == (2, 1, 0)
    return {
        "module_dimensions": tower.module_dimensions,
        "layer_dimensions": tower.layer_dimensions,
        "split_profile": tower.split_profile,
        "cellular_cohomology": totalization.complex.cohomology_dimensions,
    }


def degree_42_audit() -> dict[str, object]:
    """Replay the actual base-square conormal representation as a tower."""

    artifact = json.loads(CONORMAL_ARTIFACT.read_text())
    calculation = artifact["calculation"]
    presentation = calculation["finite_module_presentation"]
    variables = tuple(presentation["variables"])
    total = FiniteModuleRepresentation(
        variables,
        tuple(matrix(action) for action in presentation["total_actions"]),
        "degree-42 total conormal modulo base square",
    )
    spectator = FiniteModuleRepresentation(
        variables,
        tuple(
            matrix(action)
            for action in presentation["spectator_actions"]
        ),
        "degree-42 spectator conormal modulo base square",
    )
    terminal = FiniteModuleRepresentation.zero(
        variables, "degree-42 terminal zero module"
    )
    conormal_projection = EquivariantModuleSurjection(
        total,
        spectator,
        matrix(presentation["projection"]),
        "degree-42 conormal projection",
    )
    terminal_projection = EquivariantModuleSurjection(
        spectator,
        terminal,
        sp.zeros(0, spectator.dimension),
        "degree-42 terminal projection",
    )
    tower = PostnikovModuleTower(
        (conormal_projection, terminal_projection),
        "degree-42 two-layer conormal tower",
    )
    assert tower.module_dimensions == (6, 2, 0)
    assert tower.layer_dimensions == (4, 2)
    assert tower.split_profile == (False, True)

    braid = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    totalization = postnikov_braid_totalization(
        braid,
        base_dimension=3,
        tower=tower,
        layer_names=("effective-sector", "spectator"),
        name="degree-42 base-square Postnikov associated graded",
    )
    assert totalization.complex.cohomology_dimensions == (3, 6, 0)
    return {
        "module_dimensions": tower.module_dimensions,
        "layer_dimensions": tower.layer_dimensions,
        "split_profile": tower.split_profile,
        "cellular_cohomology": totalization.complex.cohomology_dimensions,
        "interpretation": (
            "the four-dimensional first kernel is the effective sector "
            "after non-flat base change; HRCELL5 identifies the missing "
            "two source dimensions as a Tor image"
        ),
    }


def main() -> None:
    synthetic = synthetic_multiflag_audit()
    degree_30 = degree_30_audit()
    degree_42 = degree_42_audit()
    output = {
        "schema": "cellular-postnikov-transitivity.v1",
        "status": "exact finite-module tower formalism",
        "synthetic_multiflag": synthetic,
        "degree_30_sector_only": degree_30,
        "degree_42_base_square": degree_42,
        "conormal_artifact": str(CONORMAL_ARTIFACT.relative_to(ROOT)),
        "conormal_artifact_sha256": hashlib.sha256(
            CONORMAL_ARTIFACT.read_bytes()
        ).hexdigest(),
        "proved_by_checker": [
            "commuting coordinate actions define finite modules",
            "every declared projection is surjective and equivariant",
            "kernel actions are restricted exactly",
            "compatible sections are solved by exact rational linear algebra",
            "arbitrary tower lengths feed their kernel layers into the cellular totalization",
            "the degree-30 and actual degree-42 finite representations have the asserted tower profiles",
        ],
        "theorem_boundary": (
            "the derived multi-flag theorem is proved in the canonical "
            "note; this checker validates its finite-module and cellular "
            "linear-algebra realization, not the all-degree Hessian--Ritt "
            "homotopy-limit comparison"
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_cellular_postnikov_transitivity.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: arbitrary finite Postnikov towers extract exact layers")
    print("PASS: compatible splitting is solved by exact rational algebra")
    print("PASS: degree 30 is a one-layer sector-only tower")
    print("PASS: degree 42 is a non-split two-layer conormal tower")
    print("PASS: tower layers feed the cellular totalization")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
