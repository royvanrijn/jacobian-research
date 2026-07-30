#!/usr/bin/env python3
"""Use PGS3 on the HC4 essential-rank-two repeated-root packet.

The result is conditional on an open lower-layer condition at every
repeated root: after eliminating the transverse leading active equation,
one redundant active component is epsilon times a unit.  On that stratum,
the first Segre multiplicity is the total repeated-root excess.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    SingularEssentialGradientNormalSlice,
)


ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_binary_root_partition_segre.json"
)


def quotient_length(
    generators: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> int:
    basis = sp.groebner(generators, *variables, order="grevlex")
    assert basis.is_zero_dimensional
    leading_exponents = tuple(
        polynomial.LM(order=basis.order).exponents
        for polynomial in basis.polys
    )
    bound = (
        max(sum(exponents) for exponents in leading_exponents)
        + sum(max(exponents) for exponents in zip(*leading_exponents))
    )
    standard = 0
    for total_degree in range(bound + 1):
        for first_exponent in range(total_degree + 1):
            exponent = (first_exponent, total_degree - first_exponent)
            if not any(
                all(
                    left >= right
                    for left, right in zip(exponent, leading)
                )
                for leading in leading_exponents
            ):
                standard += 1
    return standard


# A binary quintic with a root of multiplicity e has local form
# x^e*u(x), u(0)!=0.  Euler gives orders e-1 and e for its two partials,
# so the transverse Jacobian algebra is K[x]/(x^(e-1)).
# The lower-layer active-unit condition kills epsilon and leaves length e-1.
epsilon, x = sp.symbols("epsilon x")
local_calibrations: dict[str, dict[str, object]] = {}
for multiplicity in (2, 3, 4):
    jacobian_length = multiplicity - 1
    model = SingularEssentialGradientNormalSlice(
        ambient_dimension=4,
        map_degree=4,
        essential_rank=2,
        singular_locus_dimension=0,
        transverse_jacobian_length=jacobian_length,
    )
    assert model.base_codimension == 2
    assert model.active_truncated_length(
        generic_rank=0,
        torsion_orders=(1,) * jacobian_length,
    ) == jacobian_length
    assert (
        model.order_one_active_segre_contribution
        == jacobian_length
    )

    # On y=1 use h5=x^e*y^(5-e) and h4=x*y^3+y^4.
    # Both active generators are deformed; the radial lower derivative
    # 3*x+4 is a unit at x=0.
    generators = (
        epsilon**4,
        multiplicity * x ** (multiplicity - 1) + epsilon,
        (5 - multiplicity) * x**multiplicity
        + epsilon * (3 * x + 4),
    )
    length = quotient_length(generators, (epsilon, x))
    assert length == jacobian_length
    local_calibrations[f"root_multiplicity_{multiplicity}"] = {
        "transverse_jacobian_length": jacobian_length,
        "active_dvr_profile": {
            "generic_rank": 0,
            "torsion_orders": [1] * jacobian_length,
        },
        "truncated_and_final_length": length,
    }


# Essential binary quintics have at least two distinct roots.  These are all
# singular root partitions of five that still use two essential variables.
root_partitions = (
    (2, 1, 1, 1),
    (3, 1, 1),
    (2, 2, 1),
    (4, 1),
    (3, 2),
)
partition_records: dict[str, dict[str, object]] = {}
for partition in root_partitions:
    repeated_excess = sum(
        multiplicity - 1
        for multiplicity in partition
        if multiplicity >= 2
    )
    assert repeated_excess == 5 - len(partition)
    assert repeated_excess in (1, 2, 3)
    partition_records["+".join(map(str, partition))] = {
        "distinct_roots": len(partition),
        "repeated_root_excess": repeated_excess,
        "forced_sigma2_on_active_unit_stratum": repeated_excess,
        "atlas_rows": {},
    }


atlas_payload = json.loads(ATLAS.read_text())
sigma2_counts: dict[str, dict[str, int]] = {}
for affine_degree in (2, 3):
    atlas_key = f"gradient_degree_4_affine_degree_{affine_degree}"
    codimension_two_rows = [
        row
        for row in atlas_payload["atlases"][atlas_key]
        if row["leading_base_codimension"] == 2
    ]
    counts = {
        str(sigma2): sum(
            row["segre_degrees"][1] == sigma2
            for row in codimension_two_rows
        )
        for sigma2 in range(1, 13)
    }
    expected = (
        [51, 44, 37, 32, 26, 21, 16, 13, 9, 6, 3, 2]
        if affine_degree == 2
        else [50, 43, 36, 31, 25, 20, 15, 12, 8, 5, 3, 1]
    )
    assert list(counts.values()) == expected
    sigma2_counts[f"affine_degree_{affine_degree}"] = counts

    for record in partition_records.values():
        sigma2 = record["forced_sigma2_on_active_unit_stratum"]
        matching_rows = [
            row
            for row in codimension_two_rows
            if row["segre_degrees"][1] == sigma2
        ]
        record["atlas_rows"][f"affine_degree_{affine_degree}"] = len(
            matching_rows
        )


assert partition_records["2+1+1+1"]["atlas_rows"] == {
    "affine_degree_2": 51,
    "affine_degree_3": 50,
}
assert partition_records["3+1+1"]["atlas_rows"] == {
    "affine_degree_2": 44,
    "affine_degree_3": 43,
}
assert partition_records["2+2+1"]["atlas_rows"] == {
    "affine_degree_2": 44,
    "affine_degree_3": 43,
}
assert partition_records["4+1"]["atlas_rows"] == {
    "affine_degree_2": 37,
    "affine_degree_3": 36,
}
assert partition_records["3+2"]["atlas_rows"] == {
    "affine_degree_2": 37,
    "affine_degree_3": 36,
}


payload = {
    "format": "hc4-binary-root-partition-segre-v1",
    "software_assumptions": {
        "python": "exact integer and symbolic polynomial arithmetic",
        "coefficient_field": "characteristic zero",
        "independent_calibration": "Macaulay2 over QQ",
    },
    "scope": (
        "Conditional PGS3 sieve for essential-rank-two singular binary "
        "quintic tops. At every repeated root, assume a redundant active "
        "gradient component has epsilon-order one with unit coefficient "
        "after transverse elimination. Then sigma2 is the total repeated-"
        "root excess. Failure of that active-unit condition is a proper "
        "lower-layer torsion stratum and remains open."
    ),
    "local_theorem": {
        "root_multiplicity": "e",
        "transverse_jacobian_algebra": "K[x]/(x^(e-1))",
        "transverse_jacobian_length": "e-1",
        "active_unit_profile": (
            "rho=0 with e-1 torsion summands, all of order one"
        ),
        "component_segre_contribution": "e-1",
    },
    "local_calibrations": local_calibrations,
    "root_partitions": partition_records,
    "codimension_two_sigma2_counts": sigma2_counts,
    "generic_discriminant_packet": {
        "root_partition": [2, 1, 1, 1],
        "forced_sigma2": 1,
        "rows_for_affine_degree_2": 51,
        "rows_for_affine_degree_3": 50,
        "rows_before_partition_sieve": {
            "affine_degree_2": 260,
            "affine_degree_3": 249,
        },
    },
    "unconditional_rows_excluded": 0,
    "next_exceptional_calculation": (
        "impose the constant-Hessian determinant faces on the failure "
        "of the active-unit condition and compute its higher torsion orders"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: binary root multiplicity e gives transverse length e-1")
print("PASS: active-unit lower layers force the same Segre contribution")
print("PASS: classified all essential singular root partitions of five")
print("PASS: generic double-root tops force sigma_2=1 on the open stratum")
print("PASS: reduced that packet from 260/249 rows to 51/50 rows")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
