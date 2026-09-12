#!/usr/bin/env python3
"""Verify the cellular Hessian--Ritt cotangent prototype.

The calculation has two levels.

1. Exact rational matrices verify the vertex/move/cell totalizations.
2. Previously certified conormal and jet dimensions decorate the relative
   half-path blocks in degrees 30 and 42.

The degree-42 sum is the associated graded of the sector/spectator flag.  It
does not assert a splitting of the completed cotangent transitivity triangle.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.hessian_ritt_cellular import (  # noqa: E402
    braid_totalization,
    ritt_cellular_coboundaries,
)
from jcsearch.ritt_complex import (  # noqa: E402
    MoveType,
    degree_thirty_braid_decorations,
    permutation_ritt_complex,
    symmetric_braid_complex,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_ritt_cellular_cotangent_prototype.json"
)


def summary(model) -> dict[str, object]:
    return {
        "chain_dimensions": model.complex.dimensions,
        "differential_ranks": model.complex.ranks,
        "cohomology_dimensions": model.complex.cohomology_dimensions,
        "vertex_module_summands": len(model.degree_zero_labels),
        "move_module_summands": len(model.degree_one_labels),
        "two_cell_module_summands": len(model.degree_two_labels),
    }


def audit_cellular_incidence() -> dict[str, object]:
    braid = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    braid_model = ritt_cellular_coboundaries(braid)
    assert braid_model.complex.cohomology_dimensions == (1, 0, 0)

    permutohedron_boundary = permutation_ritt_complex(
        (2, 3, 5, 7), MoveType.CHEBYSHEV
    )
    boundary_model = ritt_cellular_coboundaries(permutohedron_boundary)
    relation_counts = {
        relation: sum(
            cell.relation == relation
            for cell in permutohedron_boundary.two_cells
        )
        for relation in ("commuting", "braid")
    }
    assert relation_counts == {"commuting": 6, "braid": 8}
    assert boundary_model.complex.cohomology_dimensions == (1, 0, 1)
    return {
        "degree_30_filled_braid": summary(braid_model),
        "four_factor_coxeter_two_skeleton": {
            **summary(boundary_model),
            "cell_counts": relation_counts,
            "interpretation": (
                "H2=1 is the top class of the permutohedron boundary; "
                "a three-cell, not another cotangent relation, kills it"
            ),
        },
    }


def audit_degree_30() -> dict[str, object]:
    braid = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    sectors = []
    for decoration in degree_thirty_braid_decorations():
        model = braid_totalization(
            braid,
            base_dimension=2,
            defect_dimensions=(1,),
            defect_names=(
                f"conormal-cut-{decoration.composite_omission}",
            ),
            name=(
                "degree-30 cellular cotangent, omitted cut "
                f"{decoration.composite_omission}"
            ),
        )
        assert model.complex.cohomology_dimensions == (2, 1, 0)
        slice_ = decoration.transverse_slice
        sectors.append(
            {
                "composite_omission": decoration.composite_omission,
                "prime_omission": decoration.prime_omission,
                "linear_totalization": summary(model),
                "cellular_H1": {
                    "conormal_fiber_dimension": 1,
                    "support": "z=0",
                    "annihilator_power": decoration.conductor_power,
                },
                "cellular_H2_dimension": 0,
                "completed_comparison": {
                    "path_nilpotence_index": decoration.nilpotence_index,
                    "transverse_slice_exponents": slice_.exponents,
                    "transverse_slice_hilbert_vector": slice_.hilbert_vector,
                    "transverse_slice_length": slice_.length,
                    "point_cotangent_homology_ranks": (
                        slice_.point_cotangent_homology_ranks
                    ),
                    "augmentation_kernel_length": (
                        slice_.augmentation_ideal_length
                    ),
                },
            }
        )
    return {
        "degree": 30,
        "reduced_component": "Dickson A2, with power boundary z=0",
        "base_tangent_dimension": 2,
        "sectors": sectors,
        "conclusion": (
            "all three associated-graded cellular cotangent complexes have "
            "H1 dimension one and H2 zero, while their completed transverse "
            "algebras and z-thicknesses differ"
        ),
    }


def audit_degree_42() -> dict[str, object]:
    braid = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    # Exact completed local lengths from HR42C.  Differences from the
    # Dickson column and between the two path columns give the two filtered
    # coefficient-module dimensions.
    jets = (
        {
            "order": 1,
            "base_dimension": 1,
            "sector_dimension": 0,
            "spectator_dimension": 0,
        },
        {
            "order": 2,
            "base_dimension": 3,
            "sector_dimension": 1,
            "spectator_dimension": 1,
        },
        {
            "order": 3,
            "base_dimension": 6,
            "sector_dimension": 5,
            "spectator_dimension": 3,
        },
        {
            "order": 4,
            "base_dimension": 10,
            "sector_dimension": 13,
            "spectator_dimension": 6,
        },
    )
    jet_results = []
    for jet in jets:
        model = braid_totalization(
            braid,
            base_dimension=jet["base_dimension"],
            defect_dimensions=(
                jet["sector_dimension"],
                jet["spectator_dimension"],
            ),
            defect_names=("sector-z8", "spectator-z"),
            name=f"degree-42 order-{jet['order']} associated graded",
        )
        expected = (
            jet["base_dimension"],
            jet["sector_dimension"] + jet["spectator_dimension"],
            0,
        )
        assert model.complex.cohomology_dimensions == expected
        jet_results.append({**jet, "totalization": summary(model)})

    conormal = braid_totalization(
        braid,
        base_dimension=2,
        defect_dimensions=(1, 1),
        defect_names=("sector-z8", "spectator-z"),
        name="degree-42 first conormal associated graded",
    )
    assert conormal.complex.cohomology_dimensions == (2, 2, 0)
    return {
        "degree": 42,
        "word": "2 o 7 o 3",
        "ideal_flag": "I_6 < I_7=I_boundary < K",
        "first_conormal_totalization": summary(conormal),
        "cellular_H1": {
            "sector": {
                "fiber_dimension": 1,
                "minimal_annihilator": "z^8",
            },
            "spectator": {
                "fiber_dimension": 1,
                "minimal_annihilator": "z",
            },
        },
        "cellular_H2_dimension": 0,
        "completed_jet_totalizations": jet_results,
        "extension_status": (
            "the direct sum is the associated graded of the exact ideal "
            "flag; verify_degree42_ritt_cellular_extension.py proves that "
            "the nested order-three and order-four module extensions are "
            "non-split, and verify_degree42_ritt_tensor_extension.py proves "
            "that the completed ideal-module extension is non-split; "
            "verify_degree42_ritt_conormal_transitivity.py proves that the "
            "first-Postnikov conormal projection is non-split and hence "
            "that the cotangent transitivity connecting morphism is "
            "nonzero; verify_degree42_ritt_postnikov_overlap.py proves "
            "that the completed sector-to-total conormal map is injective "
            "and separates its truncated base-change Tor"
        ),
    }


def main() -> None:
    result = {
        "schema": "hessian-ritt-cellular-cotangent-prototype.v1",
        "status": "exact linear and finite-jet prototype",
        "cellular_incidence": audit_cellular_incidence(),
        "degree_30": audit_degree_30(),
        "degree_42": audit_degree_42(),
        "theorem_boundary": {
            "proved_by_this_checker": [
                "all displayed rational matrices form chain complexes",
                "the stated H1 and H2 dimensions follow by exact rank",
                "commuting-square and braid-cell orientations satisfy d1*d0=0",
                "the certified degree-30 and degree-42 dimensions assemble consistently",
            ],
            "not_proved": [
                "the coefficient system is quasi-isomorphic to the full derived Hessian intersection",
                "the individual higher cotangent homology modules are computed",
                "the filtered H2 obstruction groups vanish beyond the displayed associated graded",
                "an all-degree reusable deformation theorem",
            ],
        },
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_hessian_ritt_cellular_cotangent_prototype.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
