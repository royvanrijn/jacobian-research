#!/usr/bin/env python3
"""Verify all-dimensional singular projective-gradient normal slices."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    SingularEssentialGradientNormalSlice,
    truncated_dvr_module_length,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "projective_gradient_singular_slices.json"
)


def quotient_length(
    generators: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> int:
    """Count standard monomials of a zero-dimensional Groebner basis."""

    basis = sp.groebner(
        generators,
        *variables,
        order="grevlex",
    )
    assert basis.is_zero_dimensional
    leading_exponents = tuple(
        polynomial.LM(order=basis.order).exponents
        for polynomial in basis.polys
    )
    search_bound = max(
        sum(exponents) for exponents in leading_exponents
    ) + sum(max(exponents) for exponents in zip(*leading_exponents))
    standard = 0
    for total_degree in range(search_bound + 1):
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


# Dimension-free support and DVR-profile regression ledger.
records: list[dict[str, object]] = []
for ambient_dimension in range(3, 11):
    for map_degree in range(2, 8):
        for essential_rank in range(2, ambient_dimension):
            for singular_dimension in range(essential_rank - 1):
                for jacobian_length in (1, 2, 4):
                    model = SingularEssentialGradientNormalSlice(
                        ambient_dimension=ambient_dimension,
                        map_degree=map_degree,
                        essential_rank=essential_rank,
                        singular_locus_dimension=singular_dimension,
                        transverse_jacobian_length=jacobian_length,
                        component_degree=3,
                    )
                    assert model.joined_support_dimension == (
                        ambient_dimension
                        - essential_rank
                        + singular_dimension
                    )
                    assert model.base_codimension == (
                        essential_rank - singular_dimension
                    )

                    flat_length = model.active_truncated_length(
                        generic_rank=jacobian_length,
                        torsion_orders=(),
                    )
                    order_one_length = model.active_truncated_length(
                        generic_rank=0,
                        torsion_orders=(1,) * jacobian_length,
                    )
                    assert flat_length == map_degree * jacobian_length
                    assert order_one_length == jacobian_length
                    assert model.leading_segre_contribution_bounds(
                        generic_rank=jacobian_length,
                        torsion_orders=(),
                    ) == (
                        3 * jacobian_length,
                        3 * map_degree * jacobian_length,
                    )
                    assert model.leading_segre_contribution_bounds(
                        generic_rank=0,
                        torsion_orders=(1,) * jacobian_length,
                    ) == (
                        3 * jacobian_length,
                        3 * jacobian_length,
                    )
                    assert (
                        model.unit_kernel_gradient_segre_contribution
                        == 3 * jacobian_length
                    )
                    assert (
                        model.order_one_active_segre_contribution
                        == 3 * jacobian_length
                    )
                    records.append(
                        {
                            **asdict(model),
                            "kernel_dimension": model.kernel_dimension,
                            "joined_support_dimension": (
                                model.joined_support_dimension
                            ),
                            "base_codimension": model.base_codimension,
                            "flat_active_truncated_length": flat_length,
                            "order_one_torsion_active_length": (
                                order_one_length
                            ),
                        }
                    )


# A repeated-root binary quintic provides three exact lower-layer profiles.
# On y=1 take h5=x^3*y^2, so the transverse top Jacobian algebra at x=0 is
# B=Q[x]/(x^2), of length mu=2.  With map degree four:
#
#   h4=0       gives R^2,                    length 8 mod epsilon^4;
#   h4=x*y^3   gives R/(epsilon^2)+R/(epsilon), length 3;
#   h4=y^4     gives R/(epsilon)+R/(epsilon),   length 2.
epsilon, x = sp.symbols("epsilon x")
flat_generators = (epsilon**4, 3 * x**2, 2 * x**3)
mixed_generators = (
    epsilon**4,
    3 * x**2 + epsilon,
    2 * x**3 + 3 * epsilon * x,
)
order_one_generators = (
    epsilon**4,
    3 * x**2,
    2 * x**3 + 4 * epsilon,
)
binary_lengths = {
    "h4_zero_flat": quotient_length(
        flat_generators,
        (epsilon, x),
    ),
    "h4_x_y3_mixed_torsion": quotient_length(
        mixed_generators,
        (epsilon, x),
    ),
    "h4_y4_order_one_torsion": quotient_length(
        order_one_generators,
        (epsilon, x),
    ),
}
assert binary_lengths == {
    "h4_zero_flat": 8,
    "h4_x_y3_mixed_torsion": 3,
    "h4_y4_order_one_torsion": 2,
}
assert truncated_dvr_module_length(
    truncation_order=4,
    generic_rank=2,
    torsion_orders=(),
) == 8
assert truncated_dvr_module_length(
    truncation_order=4,
    generic_rank=0,
    torsion_orders=(2, 1),
) == 3
assert truncated_dvr_module_length(
    truncation_order=4,
    generic_rank=0,
    torsion_orders=(1, 1),
) == 2


# HC4PPG6 support packets are exact specializations of the join formula.
hc4_repeated_binary = SingularEssentialGradientNormalSlice(4, 4, 2, 0, 1)
hc4_isolated_ternary = SingularEssentialGradientNormalSlice(4, 4, 3, 0, 1)
hc4_curve_ternary = SingularEssentialGradientNormalSlice(4, 4, 3, 1, 1)
assert hc4_repeated_binary.base_codimension == 2
assert hc4_isolated_ternary.base_codimension == 3
assert hc4_curve_ternary.base_codimension == 2


for bad_call in (
    lambda: SingularEssentialGradientNormalSlice(4, 4, 1, 0, 1),
    lambda: SingularEssentialGradientNormalSlice(4, 4, 3, 2, 1),
    lambda: SingularEssentialGradientNormalSlice(4, 4, 3, 0, 0),
    lambda: hc4_isolated_ternary.active_truncated_length(
        generic_rank=0,
        torsion_orders=(),
    ),
):
    try:
        bad_call()
    except ValueError:
        pass
    else:
        raise AssertionError("invalid singular normal-slice input was accepted")


payload = {
    "format": "projective-gradient-singular-slices-v1",
    "software_assumptions": {
        "python": "exact integer and symbolic polynomial arithmetic",
        "coefficient_field": "characteristic zero",
        "independent_calibration": "Macaulay2 over QQ",
    },
    "theorem": {
        "top_support": (
            "for a singular component C of dimension s and degree d, "
            "the infinity component is Join(P^(n-r-1),C), has dimension "
            "n-r+s, codimension r-s in P^n, and degree d"
        ),
        "transverse_jacobian_length": "mu=length(B_C)",
        "active_dvr_profile": (
            "M=K[[epsilon]]^rho direct_sum "
            "K[[epsilon]]/(epsilon^a_j), with rho+#a_j=mu"
        ),
        "truncated_active_length": (
            "m*rho+sum_j min(m,a_j)"
        ),
        "component_bounds": (
            "d*mu <= leading Segre contribution <= "
            "d*(m*rho+sum_j min(m,a_j))"
        ),
        "unit_kernel_gradient_law": (
            "if a missing kernel component is epsilon times a unit, "
            "the contribution is exactly d*mu"
        ),
        "order_one_active_law": (
            "if rho=0 and all torsion orders are one, "
            "the contribution is exactly d*mu"
        ),
        "flat_filtered_law": (
            "if rho=mu and in(G)=epsilon^q*s, the upper bound is "
            "d*(m*mu-(m-q)*dim(B_C*s))"
        ),
    },
    "regression_range": {
        "ambient_dimensions": [3, 10],
        "map_degrees": [2, 7],
        "all_essential_ranks_and_singular_dimensions": True,
        "records_checked": len(records),
    },
    "binary_quintic_calibration": {
        "top": "h5=x^3*y^2 near the singular point [0:1]",
        "transverse_jacobian_algebra": "QQ[x]/(x^2)",
        "transverse_jacobian_length": 2,
        "profiles": {
            "h4=0": {
                "generic_rank": 2,
                "torsion_orders": [],
                "truncated_length": binary_lengths["h4_zero_flat"],
            },
            "h4=x*y^3": {
                "generic_rank": 0,
                "torsion_orders": [2, 1],
                "truncated_length": (
                    binary_lengths["h4_x_y3_mixed_torsion"]
                ),
            },
            "h4=y^4": {
                "generic_rank": 0,
                "torsion_orders": [1, 1],
                "truncated_length": (
                    binary_lengths["h4_y4_order_one_torsion"]
                ),
            },
        },
        "conclusion": (
            "the same top singularity and transverse Jacobian length do "
            "not determine the lower-layer Segre multiplicity"
        ),
    },
    "hc4_support_specializations": {
        "rank_two_repeated_binary_root": {"base_codimension": 2},
        "rank_three_isolated_ternary_singularity": {
            "base_codimension": 3
        },
        "rank_three_curve_singularity": {"base_codimension": 2},
    },
    "records": records,
    "scope": (
        "The join support and DVR length formula are exact. The profile "
        "(rho,a_j), transverse Jacobian length, and component degree must "
        "still be computed from the actual singular top and lower layers. "
        "No universal numerical Segre vector follows from singular-locus "
        "dimension alone."
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: checked the all-dimensional kernel-vertex/singularity join")
print("PASS: verified exact truncated DVR-module length profiles")
print("PASS: bounded every component between d*mu and its active length")
print("PASS: one binary top realizes flat, mixed, and order-one torsion")
print("PASS: recovered the singular HC4PPG6 support packets")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
