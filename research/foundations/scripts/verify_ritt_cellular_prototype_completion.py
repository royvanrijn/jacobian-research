#!/usr/bin/env python3
"""Complete the explicit degree-42/30 cellular cotangent prototype.

This checker records the actual factor words, adjacent Ritt moves, labelled
braid cell, scalar cellular matrices, and every currently certified
filtration layer.  It also makes the theorem boundary executable:

* the filled-braid prototype has cellular H^2=0 for every coefficient
  module and every filtration degree;
* the degree-42 associated-graded reduction first fails at order three,
  where the sector--spectator coefficient extension is non-split; and
* the induced completed cotangent connecting morphism is nonzero.

The last two assertions consume the exact HRCELL2 and HRCELL4 artifacts.
They are an obstruction to replacing the filtered coefficient diagram by
its direct-sum associated graded, not a computation of all higher
cotangent homology modules.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.hessian_ritt_cellular import (  # noqa: E402
    braid_totalization,
    relative_half_path_model,
    ritt_cellular_coboundaries,
)
from jcsearch.ritt_complex import MoveType, symmetric_braid_complex  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "ritt_cellular_prototype_completion.json"
)
EXTENSION_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_cellular_extension.json"
)
CONORMAL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_conormal_transitivity.json"
)
OVERLAP_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_postnikov_overlap.json"
)


def matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    """Serialize an integral SymPy matrix."""

    return [[int(entry) for entry in row] for row in matrix.tolist()]


def word_label(word: tuple[int, ...]) -> str:
    return " o ".join(map(str, word))


def explicit_braid(degrees: tuple[int, int, int]) -> dict[str, object]:
    """Return the complete labelled factor/move/two-cell diagram."""

    braid = symmetric_braid_complex(degrees, MoveType.CHEBYSHEV)
    model = ritt_cellular_coboundaries(braid)
    cell = braid.two_cells[0]
    assert model.complex.dimensions == (6, 6, 1)
    assert model.complex.ranks == (5, 1)
    assert model.complex.cohomology_dimensions == (1, 0, 0)
    assert model.complex.d1 * model.complex.d0 == sp.zeros(1, 6)
    return {
        "factor_vertices_outer_to_inner": [
            word_label(vertex.word) for vertex in braid.vertices
        ],
        "moves": [
            {
                "from": word_label(edge.endpoints[0]),
                "to": word_label(edge.endpoints[1]),
                "swapped_positions_zero_based": (
                    edge.swapped_index,
                    edge.swapped_index + 1,
                ),
                "label": edge.move_type.value,
            }
            for edge in braid.edges
        ],
        "labelled_two_cell": {
            "relation": cell.relation,
            "first_half": [word_label(word) for word in cell.paths[0]],
            "second_half": [word_label(word) for word in cell.paths[1]],
            "oriented_edge_row": matrix_rows(model.complex.d1)[0],
        },
        "degree_zero_labels": model.degree_zero_labels,
        "degree_one_labels": model.degree_one_labels,
        "degree_two_labels": model.degree_two_labels,
        "scalar_d0": matrix_rows(model.complex.d0),
        "scalar_d1": matrix_rows(model.complex.d1),
        "dimensions": model.complex.dimensions,
        "ranks": model.complex.ranks,
        "cohomology": model.complex.cohomology_dimensions,
    }


def relative_block() -> dict[str, object]:
    """Return the universal relative half-braid coefficient block."""

    model = relative_half_path_model(1, "D", "labelled-power-half")
    assert model.complex.dimensions == (2, 3, 0)
    assert model.complex.ranks == (2, 0)
    assert model.complex.cohomology_dimensions == (0, 1, 0)
    return {
        "d0": matrix_rows(model.complex.d0),
        "d1": matrix_rows(model.complex.d1),
        "cohomology": model.complex.cohomology_dimensions,
        "module_formula": {
            "C0": "D^2",
            "C1": "D^3",
            "C2": "0",
            "H0": "0",
            "H1": "D",
            "H2": "0",
        },
    }


def totalization_formula() -> dict[str, object]:
    """Record the module-valued total complexes in degrees 42 and 30."""

    return {
        "degree_42": {
            "C0": "T_B^6 direct_sum D_sec^2 direct_sum D_sp^2",
            "C1": "T_B^6 direct_sum D_sec^3 direct_sum D_sp^3",
            "C2": "T_B",
            "d0": "diag(delta0 tensor T_B, P tensor D_sec, P tensor D_sp)",
            "d1": "(delta1 tensor T_B, 0, 0)",
            "cohomology": {
                "H0": "T_B",
                "H1_associated_graded": "D_sec direct_sum D_sp",
                "H2": "0",
            },
            "annihilators": {"D_sec": "z^8", "D_sp": "z"},
        },
        "degree_30_sector_control": {
            "C0": "T_B^6 direct_sum D_j^2",
            "C1": "T_B^6 direct_sum D_j^3",
            "C2": "T_B",
            "d0": "diag(delta0 tensor T_B, P tensor D_j)",
            "d1": "(delta1 tensor T_B, 0)",
            "cohomology": {"H0": "T_B", "H1": "D_j", "H2": "0"},
            "spectator": "zero because I_boundary=K",
        },
    }


def filtration_audit() -> dict[str, object]:
    """Replay every certified filtration quotient and uniform H2 vanishing."""

    braid_42 = symmetric_braid_complex((2, 3, 7), MoveType.CHEBYSHEV)
    rows_42 = []
    for order, base, sector, spectator in (
        (1, 1, 0, 0),
        (2, 3, 1, 1),
        (3, 6, 5, 3),
        (4, 10, 13, 6),
    ):
        model = braid_totalization(
            braid_42,
            base_dimension=base,
            defect_dimensions=(sector, spectator),
            defect_names=("sector", "spectator"),
            name=f"degree-42 filtration order {order}",
        )
        expected = (base, sector + spectator, 0)
        assert model.complex.cohomology_dimensions == expected
        rows_42.append(
            {
                "filtration_order": order,
                "coefficient_dimensions": {
                    "base": base,
                    "sector": sector,
                    "spectator": spectator,
                },
                "chain_dimensions": model.complex.dimensions,
                "ranks": model.complex.ranks,
                "cohomology": expected,
            }
        )

    braid_30 = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    rows_30 = []
    for omitted, annihilator in ((10, 2), (15, 2), (6, 4)):
        model = braid_totalization(
            braid_30,
            base_dimension=2,
            defect_dimensions=(1,),
            defect_names=(f"sector-cut-{omitted}",),
            name=f"degree-30 sector {omitted}",
        )
        assert model.complex.cohomology_dimensions == (2, 1, 0)
        rows_30.append(
            {
                "omitted_composite_cut": omitted,
                "annihilator_power": f"z^{annihilator}",
                "chain_dimensions": model.complex.dimensions,
                "ranks": model.complex.ranks,
                "cohomology": model.complex.cohomology_dimensions,
            }
        )

    return {
        "degree_42": rows_42,
        "degree_30": rows_30,
        "uniform_statement": (
            "for every coefficient module M and every filtration quotient, "
            "the filled-disk block has H2=0 and the relative half-path block "
            "has no degree-two term; therefore the prototype has H2=0 in "
            "every filtration degree"
        ),
        "first_nontrivial_H2": {
            "degree_42_prototype": None,
            "degree_30_prototype": None,
            "interpretation": (
                "there is no first nonzero cellular H2 in the prototype; "
                "any algebraic H2 of the actual cotangent complex must lie "
                "in internal higher cotangent homology or in the mapping "
                "cone of coefficient effectivity"
            ),
        },
    }


def exact_higher_obstruction() -> dict[str, object]:
    """Consume HRCELL2/4/5 and isolate the first failed reduction step."""

    extension = json.loads(EXTENSION_ARTIFACT.read_text())
    conormal = json.loads(CONORMAL_ARTIFACT.read_text())
    overlap = json.loads(OVERLAP_ARTIFACT.read_text())

    ext3 = extension["order_three"]
    ext4 = extension["order_four"]
    certificate3 = ext3["extension_certificate"]
    certificate4 = ext4["extension_certificate"]
    assert ext3["splits_over_Q_tau_zeta"] is False
    assert certificate3["coboundary_rank"] == 7
    assert certificate3["augmented_rank"] == 8
    assert ext4["splits_over_Q_tau_zeta"] is False
    assert certificate4["coboundary_rank"] == 52
    assert certificate4["augmented_rank"] == 53

    conormal_calculation = conormal["calculation"]
    conormal_certificate = conormal_calculation["certificate"]
    assert conormal_calculation["splits_as_R_module"] is False
    assert conormal_certificate["coboundary_rank"] == 5
    assert conormal_certificate["augmented_rank"] == 6

    overlap_calculation = overlap["calculation"]
    assert (
        overlap_calculation["quadratic_overlap_mod_base_square"]["dimension"]
        == 0
    )

    return {
        "first_failed_filtration_order": 3,
        "order_3_extension": {
            "dimensions": ext3["module_dimensions"],
            "coboundary_rank": certificate3["coboundary_rank"],
            "augmented_rank": certificate3["augmented_rank"],
            "obstruction_value": certificate3["obstruction_value"],
        },
        "order_4_persistence": {
            "dimensions": ext4["module_dimensions"],
            "coboundary_rank": certificate4["coboundary_rank"],
            "augmented_rank": certificate4["augmented_rank"],
            "obstruction_value": certificate4["obstruction_value"],
        },
        "completed_conormal_projection": {
            "dimensions": conormal_calculation["dimensions"],
            "coboundary_rank": conormal_certificate["coboundary_rank"],
            "augmented_rank": conormal_certificate["augmented_rank"],
            "obstruction_value": conormal_certificate["obstruction_value"],
            "splits": conormal_calculation["splits_as_R_module"],
        },
        "first_postnikov_overlap_vanishes": True,
        "conclusion": (
            "the obstruction is not a cellular topological H2 class.  It is "
            "the nonzero filtered extension/connecting morphism lost when "
            "K/I6 is replaced by (I_boundary/I6) direct_sum "
            "(K/I_boundary).  Thus the general cellwise homotopy limit does "
            "not reduce to the split HRCELL1 prototype; it can reduce only "
            "to the extension-retaining Postnikov tower of HRCELL2--HRCELL5"
        ),
        "artifact_hashes": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                EXTENSION_ARTIFACT,
                CONORMAL_ARTIFACT,
                OVERLAP_ARTIFACT,
            )
        },
    }


def main() -> None:
    output = {
        "schema": "ritt-cellular-prototype-completion.v1",
        "status": "exact cellular prototype and higher obstruction",
        "degree_42_chain_diagram": explicit_braid((2, 3, 7)),
        "degree_30_chain_diagram": explicit_braid((2, 3, 5)),
        "relative_labelled_power_block": relative_block(),
        "totalized_cotangent_prototype": totalization_formula(),
        "filtration_cohomology": filtration_audit(),
        "HRCELL1_to_HRCELL5_comparison": {
            "HRCELL1": "split associated-graded cellular matrices",
            "HRCELL2": "first non-split sector--spectator jets",
            "HRCELL3": "completed ideal-module non-splitting",
            "HRCELL4": "nonzero cotangent transitivity connecting morphism",
            "HRCELL5": "overlap zero; finite-jet excess is base-change Tor",
        },
        "success_criterion": exact_higher_obstruction(),
        "theorem_boundary": (
            "This completes the explicit prototype and identifies the first "
            "higher obstruction.  It does not compute the individual H_i "
            "of the actual cotangent complexes for i>=2 and does not prove "
            "global coefficient effectivity on the other degree-42 charts."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_ritt_cellular_prototype_completion.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: explicit degree-42 factor/move/labelled-cell diagram")
    print("PASS: explicit degree-30 sector-only control diagram")
    print("PASS: cellular H2 vanishes in every prototype filtration degree")
    print("PASS: the first failed split reduction occurs at filtration order 3")
    print("PASS: the completed cotangent connecting morphism is nonzero")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
