#!/usr/bin/env python3
"""Exact nodal persistence of the non-Cartier different through order six.

The nodal gauge cokernel is Q[y,z](-3).  Hence, after removing formal gauge
directions, its complete quartic and quintic pieces have bases

    y*eta, z*eta

and

    y^2*eta, y*z*eta, z^2*eta.

The complete sextic quotient is

    y^3*eta, y^2*z*eta, y*z^2*eta, z^3*eta.

On the resulting nine-parameter polynomial tensor family this checker
computes the intrinsic support module and the collision Nakayama quotient
of J=Ann_B(Omega).  Strict weighted-Rees packets commute J with every
geometric parameter specialization.  The conclusion is exact through the
complete sextic normal-form quotient; order seven and higher are not covered.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import verify_cubic_formal_gauge_cokernel_atlas as atlas  # noqa: E402
import verify_cubic_symbol_double_saturation as cubic  # noqa: E402
from verify_nodal_cubic_formal_slice import coefficient_column  # noqa: E402
import verify_universal_cubic_cotangent_saturation as smooth  # noqa: E402


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "nodal_sextic_different_persistence.json"
)


def tensor_record(
    tensor: dict[tuple[int, int, int], sp.Expr],
) -> list[str]:
    return [sp.sstr(tensor[triple]) for triple in smooth.TRIPLES]


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def compatible_degree_dimension(
    compatibility: sp.Matrix, degree: int
) -> tuple[int, int]:
    """Return dimensions of the homogeneous compatible space and target."""

    input_monomials = cubic.homogeneous_monomials(degree)
    columns = []
    for component in range(len(smooth.TRIPLES)):
        for monomial in input_monomials:
            vector = sp.zeros(len(smooth.TRIPLES), 1)
            vector[component] = monomial
            columns.append(
                coefficient_column(
                    compatibility * vector,
                    degree + 1,
                )
            )
    constraint_matrix = sp.Matrix.hstack(*columns)
    rank = constraint_matrix.rank()
    return constraint_matrix.cols - rank, constraint_matrix.rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="intentionally replace the pinned generated artifact",
    )
    args = parser.parse_args()
    cubic.FACTOR_SINGULAR_EXPRESSIONS = False
    compatibility = smooth.compatibility_matrix()
    nodal_tensor = atlas.symbol_tensor(cubic.CUBIC_STRATA["nodal"])
    gauge = smooth.gauge_matrix(nodal_tensor)
    eta = atlas.symbol_tensor(cubic.Z**3)

    assert (
        compatibility * eta
    ).applyfunc(sp.expand) == sp.zeros(6, 1)

    quartic_basis = cubic.quartic_kernel_basis_tensors()
    quartic_tensors = (quartic_basis[0], quartic_basis[1])
    quartic_action = sp.Matrix.hstack(
        *[
            coefficient_column(gauge[:, column] * variable, 4)
            for column in range(gauge.cols)
            for variable in cubic.BASE_VARIABLES
        ]
    )
    quartic_slice = sp.Matrix.hstack(
        *[
            coefficient_column(
                sp.Matrix(
                    [tensor[triple] for triple in smooth.TRIPLES]
                ),
                4,
            )
            for tensor in quartic_tensors
        ]
    )
    assert quartic_action.rank() == 22
    assert quartic_action.row_join(quartic_slice).rank() == 24

    quintic_monomials = (cubic.y**2, cubic.y * cubic.z, cubic.z**2)
    quintic_tensors = tuple(
        {
            triple: sp.expand(monomial * eta[row])
            for row, triple in enumerate(smooth.TRIPLES)
        }
        for monomial in quintic_monomials
    )
    compatible_quintic_dimension, constraint_target_dimension = (
        compatible_degree_dimension(compatibility, 5)
    )
    quintic_action = sp.Matrix.hstack(
        *[
            coefficient_column(gauge[:, column] * monomial, 5)
            for column in range(gauge.cols)
            for monomial in cubic.homogeneous_monomials(2)
        ]
    )
    quintic_slice = sp.Matrix.hstack(
        *[
            coefficient_column(monomial * eta, 5)
            for monomial in quintic_monomials
        ]
    )
    assert compatible_quintic_dimension == 42
    assert constraint_target_dimension == 168
    assert quintic_action.rank() == 39
    assert quintic_action.row_join(quintic_slice).rank() == 42

    sextic_monomials = (
        cubic.y**3,
        cubic.y**2 * cubic.z,
        cubic.y * cubic.z**2,
        cubic.z**3,
    )
    sextic_tensors = tuple(
        {
            triple: sp.expand(monomial * eta[row])
            for row, triple in enumerate(smooth.TRIPLES)
        }
        for monomial in sextic_monomials
    )
    compatible_sextic_dimension, sextic_target_dimension = (
        compatible_degree_dimension(compatibility, 6)
    )
    sextic_action = sp.Matrix.hstack(
        *[
            coefficient_column(gauge[:, column] * monomial, 6)
            for column in range(gauge.cols)
            for monomial in cubic.homogeneous_monomials(3)
        ]
    )
    sextic_slice = sp.Matrix.hstack(
        *[
            coefficient_column(monomial * eta, 6)
            for monomial in sextic_monomials
        ]
    )
    assert compatible_sextic_dimension == 64
    assert sextic_target_dimension == 216
    assert sextic_action.rank() == 60
    assert sextic_action.row_join(sextic_slice).rank() == 64

    tensors = quartic_tensors + quintic_tensors + sextic_tensors
    computation = cubic.run_singular_subspace_certificate(
        cubic.CUBIC_STRATA["nodal"], tensors, timeout=1800
    )
    expected_computation = {
        "parameter_count": 9,
        "cotangent_saturation_generators": 0,
        "support_module_dimension": 11,
        "support_ext3_vector_dimension": 0,
        "support_ext2_dimension": 9,
        "support_ext2_multiplicity": 6,
        "support_ext2_parameter_axis_radical_difference": 0,
        "support_ext2_central_pruned_presentation_difference": 0,
        "support_ext2_pruned_presentation_rank": 3,
        "support_ext2_collision_square_action_generators": 0,
        "different_generator_module_dimension": 9,
        "different_generator_module_multiplicity": 6,
        "different_generator_parameter_axis_radical_difference": 0,
        "different_generator_central_pruned_presentation_difference": 0,
        "different_generator_pruned_presentation_rank": 6,
    }
    assert computation == expected_computation

    base_change = cubic.run_singular_rees_base_change_certificate(
        cubic.CUBIC_STRATA["nodal"], tensors, timeout=1800
    )
    expected_base_change = {
        "parameter_count": 9,
        "cotangent_rees_torsion_generators": 0,
        "cotangent_initial_presentation_difference": 0,
        "annihilator_cokernel_rees_torsion_generators": 0,
        "annihilator_cokernel_initial_presentation_difference": 0,
    }
    assert base_change == expected_base_change

    tensor_records = [tensor_record(tensor) for tensor in tensors]
    artifact = {
        "schema": "nodal-sextic-different-persistence.v1",
        "case": "nodal-complete-sextic-normal-form-quotient",
        "mathematical_scope": (
            "Exact characteristic-zero calculation on the nine-parameter "
            "nodal normal-form tensor family consisting of the complete "
            "two-dimensional quartic, three-dimensional quintic, and "
            "four-dimensional sextic gauge cokernel pieces. Strict Rees "
            "packets make the annihilator "
            "intrinsic after every geometric specialization. Every fiber "
            "has saturated cotangent module, the same length-six C1 defect, "
            "and a six-dimensional collision Nakayama quotient J/nJ, so "
            "the Kahler different is not Cartier. Gauge directions are "
            "removed only through sextic order; corrections of order seven "
            "and above, normality, and Keller-open compatibility remain "
            "outside the claim."
        ),
        "basis_conventions": {
            "tensor_component_order": [
                list(triple) for triple in smooth.TRIPLES
            ],
            "quartic_basis_indices_zero_based": [0, 1],
            "quartic_quotient_basis": ["y*eta", "z*eta"],
            "quintic_quotient_basis": [
                "y^2*eta",
                "y*z*eta",
                "z^2*eta",
            ],
            "sextic_quotient_basis": [
                "y^3*eta",
                "y^2*z*eta",
                "y*z^2*eta",
                "z^3*eta",
            ],
            "eta_source_cubic": "Z^3",
        },
        "graded_dimensions": {
            "quartic_compatible": 24,
            "quartic_gauge_image": 22,
            "quartic_quotient": 2,
            "quintic_compatible": compatible_quintic_dimension,
            "quintic_gauge_image": 39,
            "quintic_quotient": 3,
            "sextic_compatible": compatible_sextic_dimension,
            "sextic_gauge_image": 60,
            "sextic_quotient": 4,
        },
        "representative_tensor_components": tensor_records,
        "representative_sha256": canonical_sha256(tensor_records),
        "exact_computation": computation,
        "fiberwise_base_change_certificate": base_change,
        "proved": [
            (
                "the displayed quintic and sextic directions complement "
                "the complete degree-five and degree-six gauge images"
            ),
            "C2 passes on every geometric parameter fiber",
            (
                "C1 fails on every geometric parameter fiber by the "
                "constant multiplicity-six support-hull defect"
            ),
            (
                "the intrinsic Kahler different has six minimal local "
                "generators at every collision and is not Cartier"
            ),
        ],
        "not_proved": [
            (
                "persistence after compatible corrections of order seven "
                "and above"
            ),
            "normality or Keller-open compatibility",
        ],
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_nodal_sextic_different_persistence.py"
        ),
    }
    if args.refresh:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        pinned = json.loads(OUTPUT.read_text(encoding="utf-8"))
        assert artifact == pinned, "stale nodal sextic persistence artifact"
    print("PASS: nodal quintic compatible/gauge dimensions are 42=39+3")
    print("PASS: nodal sextic compatible/gauge dimensions are 64=60+4")
    print("PASS: strict Rees packets commute the intrinsic different with fibers")
    print("PASS: every order-six normal-form fiber has dim(J/nJ)=6")
    print("PASS: the Kahler different is non-Cartier through sextic order")
    action = "wrote" if args.refresh else "replayed"
    print(f"PASS: {action} {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
