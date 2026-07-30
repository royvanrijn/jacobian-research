#!/usr/bin/env python3
"""Verify the all-dimensional smooth-essential gradient normal slice."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    SmoothEssentialGradientNormalSlice,
    equal_degree_ci_hilbert_function,
    filtered_missing_generator_drop,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "projective_gradient_normal_slices.json"
)


records = []
for ambient_dimension in range(2, 11):
    for map_degree in range(2, 8):
        for essential_rank in range(1, ambient_dimension):
            model = SmoothEssentialGradientNormalSlice(
                ambient_dimension=ambient_dimension,
                map_degree=map_degree,
                essential_rank=essential_rank,
            )
            hilbert_function = model.jacobian_hilbert_function
            assert sum(hilbert_function) == model.jacobian_length
            assert len(hilbert_function) - 1 == (
                model.jacobian_socle_degree
            )
            assert hilbert_function == tuple(reversed(hilbert_function))
            assert model.kernel_dimension == (
                ambient_dimension - essential_rank
            )
            assert model.support_dimension == (
                ambient_dimension - essential_rank - 1
            )
            assert model.base_codimension == essential_rank + 1
            assert model.truncated_active_length == (
                map_degree ** (essential_rank + 1)
            )
            assert model.unit_penultimate_segre_degree == (
                map_degree**essential_rank
            )

            sample_cyclic_dimension = min(
                2,
                model.jacobian_length,
            )
            sample_drop = model.missing_generator_drop(
                epsilon_order=1,
                cyclic_ideal_dimension=sample_cyclic_dimension,
            )
            assert sample_drop == (
                (map_degree - 1) * sample_cyclic_dimension
            )
            records.append(
                {
                    **asdict(model),
                    "kernel_dimension": model.kernel_dimension,
                    "support_dimension": model.support_dimension,
                    "base_codimension": model.base_codimension,
                    "jacobian_hilbert_function": list(hilbert_function),
                    "jacobian_length": model.jacobian_length,
                    "jacobian_socle_degree": (
                        model.jacobian_socle_degree
                    ),
                    "truncated_active_length": (
                        model.truncated_active_length
                    ),
                    "unit_penultimate_segre_degree": (
                        model.unit_penultimate_segre_degree
                    ),
                    "sample_non_socle_drop_at_epsilon_order_1": (
                        sample_drop
                    ),
                }
            )


# Direct helper checks and invalid-input guards.
assert equal_degree_ci_hilbert_function(
    generator_degree=4,
    codimension=3,
) == (1, 3, 6, 10, 12, 12, 10, 6, 3, 1)
assert filtered_missing_generator_drop(
    map_degree=4,
    epsilon_order=1,
    cyclic_ideal_dimension=2,
) == 6
for bad_call in (
    lambda: filtered_missing_generator_drop(
        map_degree=1,
        epsilon_order=1,
        cyclic_ideal_dimension=1,
    ),
    lambda: filtered_missing_generator_drop(
        map_degree=4,
        epsilon_order=4,
        cyclic_ideal_dimension=1,
    ),
    lambda: SmoothEssentialGradientNormalSlice(4, 4, 4),
):
    try:
        bad_call()
    except ValueError:
        pass
    else:
        raise AssertionError("invalid normal-slice input was accepted")


# HC4PPG7 is the isolated-vertex specialization (n,m,r)=(4,4,3).
hc4_rank_three = SmoothEssentialGradientNormalSlice(4, 4, 3)
assert hc4_rank_three.kernel_dimension == 1
assert hc4_rank_three.jacobian_length == 64
assert hc4_rank_three.truncated_active_length == 256
assert hc4_rank_three.jacobian_socle_degree == 9
assert hc4_rank_three.isolated_vertex_affine_degree_lower_bound(
    epsilon_order=1,
    cyclic_ideal_dimension=2,
) == 6


# HC4PPG8 is the unit penultimate-layer specialization (n,m,r)=(4,4,2).
hc4_rank_two = SmoothEssentialGradientNormalSlice(4, 4, 2)
assert hc4_rank_two.kernel_dimension == 2
assert hc4_rank_two.base_codimension == 3
assert hc4_rank_two.jacobian_hilbert_function == (1, 2, 3, 4, 3, 2, 1)
assert hc4_rank_two.jacobian_length == 16
assert hc4_rank_two.truncated_active_length == 64
assert hc4_rank_two.unit_penultimate_segre_degree == 16
try:
    hc4_rank_two.isolated_vertex_affine_degree_lower_bound(
        epsilon_order=1,
        cyclic_ideal_dimension=1,
    )
except ValueError:
    pass
else:
    raise AssertionError("a positive-dimensional vertex was treated as a point")


payload = {
    "format": "projective-gradient-normal-slices-v1",
    "software_assumptions": {
        "python": "dependency-free exact integer arithmetic",
        "coefficient_field": "characteristic zero",
        "independent_calibration": "Macaulay2 over QQ",
    },
    "theorem": {
        "input": (
            "degree-m gradient compactification in P^n; top potential "
            "depends on r<n essential variables and is smooth in P^(r-1)"
        ),
        "kernel_vertex": "P^(n-r-1)",
        "base_codimension": "r+1",
        "jacobian_hilbert_series": "(1+z+...+z^(m-1))^r",
        "jacobian_length": "m^r",
        "jacobian_socle_degree": "r*(m-1)",
        "truncated_active_length": "m^(r+1)",
        "associated_graded": "B tensor k[epsilon]/(epsilon^m)",
        "missing_generator_bound": (
            "in(G)=epsilon^q*s implies length(A*G)>="
            "(m-q)*dim(B*s)"
        ),
        "leading_segre_bound": (
            "sigma_(r+1)<=m^(r+1)-(m-q)*dim(B*s)"
        ),
        "unit_penultimate_law": (
            "if h_m restricts nontrivially to the kernel vertex, "
            "sigma_(r+1)=m^r"
        ),
        "isolated_vertex_law": (
            "if n-r=1, affine degree delta>= "
            "(m-q)*dim(B*s)"
        ),
        "non_socle_corollary": (
            "if n-r=1 and s is nonzero outside Soc(B), "
            "delta>=2*(m-q)"
        ),
    },
    "regression_range": {
        "ambient_dimensions": [2, 10],
        "map_degrees": [2, 7],
        "all_essential_ranks": True,
        "records_checked": len(records),
    },
    "specializations": {
        "HC4PPG7": {
            "parameters": {"n": 4, "m": 4, "r": 3},
            "jacobian_length": 64,
            "truncated_active_length": 256,
            "non_socle_affine_degree_lower_bound": 6,
        },
        "HC4PPG8": {
            "parameters": {"n": 4, "m": 4, "r": 2},
            "jacobian_length": 16,
            "truncated_active_length": 64,
            "unit_penultimate_sigma3": 16,
        },
    },
    "records": records,
    "scope": (
        "The theorem controls the first Segre multiplicity at the generic "
        "kernel vertex and, only for a zero-dimensional vertex, an affine-"
        "degree lower bound. It does not determine later Segre degrees, "
        "singular essential tops, or transport under cotangent/Schur lifts."
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: checked all smooth-essential normal slices for 2<=n<=10")
print("PASS: Jacobian Hilbert functions have length m^r and correct socle")
print("PASS: truncated active lengths are m^(r+1)")
print("PASS: filtered missing-generator bounds are dimension-free")
print("PASS: recovered HC4PPG7 and HC4PPG8 as specializations")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
