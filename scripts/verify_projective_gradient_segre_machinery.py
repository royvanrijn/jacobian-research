#!/usr/bin/env python3
"""Verify the all-dimensional affine-gradient Segre machinery and registry."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (
    ProjectiveGradientSegreRecord,
    SingularEssentialGradientNormalSlice,
    SmoothEssentialGradientNormalSlice,
    affine_gradient_compactification,
    full_polar_map,
    homogeneous_leading_forms,
    integrate_homogeneous_gradient,
    integrability_residuals,
    is_log_concave,
    projective_degrees_from_segre,
    segre_degrees_from_projective,
    total_segre_correction,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "projective_gradient_segre_registry.json"
)


# The triangular transforms are inverse in arbitrary dimension.  The values
# include signs and zeros deliberately; later Segre components need not be
# effective.
for dimension in range(1, 13):
    for map_degree in range(1, 8):
        sigmas = tuple(
            ((-1) ** (index + dimension))
            * (index * index + map_degree * dimension)
            for index in range(1, dimension + 1)
        )
        degrees = projective_degrees_from_segre(map_degree, sigmas)
        assert (
            segre_degrees_from_projective(map_degree, degrees) == sigmas
        )


x0, x, y, t, u = sp.symbols("x0 x y t u")
variables = (x, y, t, u)


# The affine graph and full polar constructions are distinct already for the
# triangular constant-Hessian calibration.
psi2 = x * y + t * u + y**3 / 3
gamma2 = affine_gradient_compactification(psi2, variables, x0)
polar2 = full_polar_map(psi2, variables, x0)
assert gamma2 == (
    x0**2,
    x0 * y,
    x0 * x + y**2,
    x0 * u,
    x0 * t,
)
assert polar2 == (
    x * y + t * u,
    x0 * y,
    x0 * x + y**2,
    x0 * u,
    x0 * t,
)
assert gamma2[0] != polar2[0]


# Integrability is an executable restriction on the leading infinity ideal.
leading = homogeneous_leading_forms(
    tuple(sp.diff(psi2, variable) for variable in variables),
    variables,
)
assert leading == (0, y**2, 0, 0)
assert not any(integrability_residuals(leading, variables))
assert integrate_homogeneous_gradient(leading, variables) == y**3 / 3
try:
    integrate_homogeneous_gradient((x**2, x * y, 0, 0), variables)
except ValueError:
    pass
else:
    raise AssertionError("a nonintegrable infinity tuple was accepted")


calibrations = {}
for map_degree, graph_degrees, polar_degrees in (
    (2, (1, 2, 2, 2, 1), (1, 2, 4, 4, 2)),
    (3, (1, 3, 3, 3, 1), (1, 3, 6, 6, 3)),
):
    graph = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=graph_degrees,
    )
    polar = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=polar_degrees,
    )
    assert graph.affine_degree == 1
    assert is_log_concave(graph.projective_degrees)
    calibrations[f"triangular_r{map_degree}"] = {
        "actual_affine_compactification": asdict(graph),
        "full_polar_comparison": asdict(polar),
        "source": "HC4_PROJECTIVE_POLAR_GEOMETRY.md#1-two-projective-maps",
        "certificate": "scripts/verify_projective_polar_calibrations.m2",
    }

assert calibrations["triangular_r2"][
    "actual_affine_compactification"
]["segre_degrees"] == (0, 2, -6, 15)
assert calibrations["triangular_r3"][
    "actual_affine_compactification"
]["segre_degrees"] == (0, 6, -30, 116)


cotangent_calibrations = {}
for map_degree, plane_degrees, lift_degrees in (
    (2, (1, 2, 1), (1, 2, 3, 2, 1)),
    (3, (1, 3, 1), (1, 3, 5, 3, 1)),
):
    plane = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=plane_degrees,
    )
    lift = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=lift_degrees,
    )
    assert plane.affine_degree == lift.affine_degree == 1
    assert is_log_concave(plane.projective_degrees)
    assert is_log_concave(lift.projective_degrees)
    cotangent_calibrations[f"plane_triangular_r{map_degree}"] = {
        "plane_map": asdict(plane),
        "cotangent_gradient_lift": asdict(lift),
        "potential": f"t*(x+y^{map_degree})+u*y",
        "certificate": (
            "scripts/verify_projective_gradient_segre_families.m2"
        ),
    }

assert cotangent_calibrations["plane_triangular_r2"][
    "cotangent_gradient_lift"
]["segre_degrees"] == (0, 1, 0, -9)
assert cotangent_calibrations["plane_triangular_r3"][
    "cotangent_gradient_lift"
]["segre_degrees"] == (0, 4, -12, 8)


stabilization_calibrations = {}
for map_degree, original_degrees, stable_degrees in (
    (2, (1, 2, 2, 2, 1), (1, 2, 2, 2, 2, 1)),
    (3, (1, 3, 3, 3, 1), (1, 3, 3, 3, 3, 1)),
):
    original = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=original_degrees,
    )
    stable = ProjectiveGradientSegreRecord.from_projective_degrees(
        map_degree=map_degree,
        projective_degrees=stable_degrees,
    )
    assert original.affine_degree == stable.affine_degree == 1
    stabilization_calibrations[f"triangular_r{map_degree}"] = {
        "original": asdict(original),
        "after_one_quadratic_variable": asdict(stable),
        "certificate": (
            "scripts/verify_projective_gradient_segre_families.m2"
        ),
    }


# These family records intentionally distinguish a complete multidegree
# computation from a top-degree-only transport theorem.
top_degree_controls = {
    "plane_quartic_packet_cotangent_target": {
        "status": "conditional top-degree target; no explicit packet",
        "ambient_dimension": 4,
        "map_degree_options": [2, 3],
        "affine_degree": 4,
        "weighted_corrections": {
            "2": total_segre_correction(
                2,
                ambient_dimension=4,
                affine_degree=4,
            ),
            "3": total_segre_correction(
                3,
                ambient_dimension=4,
                affine_degree=4,
            ),
        },
        "source": "HC4_PROJECTIVE_POLAR_GEOMETRY.md#5-cotangent-lifts-of-quartic-plane-packets",
    },
    "meng_yang_doubled_hc6": {
        "status": "proved top degree; individual Segre degrees uncomputed",
        "ambient_dimension": 6,
        "map_degree": 7,
        "affine_degree": 3,
        "weighted_correction": total_segre_correction(
            7,
            ambient_dimension=6,
            affine_degree=3,
        ),
        "source": "HC4_PROJECTIVE_POLAR_GEOMETRY.md#7-meng--yang-control-before-and-after-schur-descent",
    },
    "meng_yang_schur_hc5": {
        "status": "proved top degree; individual Segre degrees uncomputed",
        "ambient_dimension": 5,
        "map_degree": 13,
        "affine_degree": 3,
        "weighted_correction": total_segre_correction(
            13,
            ambient_dimension=5,
            affine_degree=3,
        ),
        "source": "HC4_PROJECTIVE_POLAR_GEOMETRY.md#7-meng--yang-control-before-and-after-schur-descent",
    },
    "homogeneous_cotangent_hn_38": {
        "status": (
            "explicit gradient Keller collision; top generic degree and "
            "individual projective degrees not yet computed"
        ),
        "ambient_dimension": 38,
        "map_degree": 3,
        "affine_degree": None,
        "source": "verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md",
    },
    "rank_reduced_cotangent_hn_44": {
        "status": (
            "explicit gradient Keller collision; top generic degree and "
            "individual projective degrees not yet computed"
        ),
        "ambient_dimension": 44,
        "map_degree": 3,
        "affine_degree": None,
        "source": "extended-geometry/RESTRICTED_MINIMA_FRONTIER.md#4-quartic-hn-consequences",
    },
    "nonhomogeneous_cotangent_hn_40": {
        "status": (
            "explicit gradient Keller collision; top generic degree and "
            "individual projective degrees not yet computed"
        ),
        "ambient_dimension": 40,
        "map_degree": 3,
        "affine_degree": None,
        "source": "extended-geometry/RESTRICTED_MINIMA_FRONTIER.md#4-quartic-hn-consequences",
    },
}
assert top_degree_controls["plane_quartic_packet_cotangent_target"][
    "weighted_corrections"
] == {"2": 12, "3": 77}
assert top_degree_controls["meng_yang_doubled_hc6"][
    "weighted_correction"
] == 117646
assert top_degree_controls["meng_yang_schur_hc5"][
    "weighted_correction"
] == 371290


smooth_normal_slice = SmoothEssentialGradientNormalSlice(
    ambient_dimension=4,
    map_degree=4,
    essential_rank=3,
)
assert smooth_normal_slice.truncated_active_length == 256
assert smooth_normal_slice.isolated_vertex_affine_degree_lower_bound(
    epsilon_order=1,
    cyclic_ideal_dimension=2,
) == 6

singular_normal_slice = SingularEssentialGradientNormalSlice(
    ambient_dimension=4,
    map_degree=4,
    essential_rank=3,
    singular_locus_dimension=0,
    transverse_jacobian_length=1,
)
assert singular_normal_slice.base_codimension == 3


payload = {
    "format": "projective-gradient-segre-registry-v3",
    "conventions": {
        "actual_affine_compactification": "[X0^m:F1^h:...:Fn^h]",
        "full_polar_map_is_separate": True,
        "formula": (
            "g_i=m^i-sum_{k=1}^i binom(i,k)m^(i-k)sigma_k"
        ),
        "top_degree": (
            "g_n equals the affine generic degree when the affine map is "
            "dominant"
        ),
    },
    "all_dimension_checks": {
        "dimensions": [1, 12],
        "map_degrees": [1, 7],
        "transform_pairs_checked": 12 * 7,
        "integrability_reconstruction": "exact Euler reconstruction",
    },
    "smooth_essential_normal_slice_theorem": {
        "status": "proved in all dimensions and map degrees m>=2",
        "parameters": "(ambient dimension n, map degree m, essential rank r)",
        "kernel_vertex": "P^(n-r-1)",
        "base_codimension": "r+1",
        "jacobian_length": "m^r",
        "truncated_active_length": "m^(r+1)",
        "filtered_drop": "(m-q)*dim(B*s)",
        "artifact": (
            "artifacts/generated-results/"
            "projective_gradient_normal_slices.json"
        ),
        "checker": "scripts/verify_projective_gradient_normal_slices.py",
        "independent_checker": (
            "scripts/verify_projective_gradient_normal_slices.m2"
        ),
    },
    "singular_essential_normal_slice_theorem": {
        "status": "proved support and DVR-profile law in all dimensions",
        "parameters": (
            "(ambient dimension n, map degree m, essential rank r, "
            "singular-component dimension s and degree d)"
        ),
        "joined_support": "Join(P^(n-r-1),C)",
        "base_codimension": "r-s",
        "transverse_jacobian_length": "mu",
        "active_dvr_profile": "(generic rank rho, torsion orders a_j)",
        "truncated_active_length": (
            "m*rho+sum_j min(m,a_j)"
        ),
        "artifact": (
            "artifacts/generated-results/"
            "projective_gradient_singular_slices.json"
        ),
        "checker": (
            "scripts/verify_projective_gradient_singular_slices.py"
        ),
        "independent_checker": (
            "scripts/verify_projective_gradient_singular_slices.m2"
        ),
    },
    "complete_calibrations": calibrations,
    "cotangent_calibrations": cotangent_calibrations,
    "quadratic_stabilization_calibrations": stabilization_calibrations,
    "downstream_consumers": {
        "hc4_quintic_infinity_rees_strata": {
            "status": (
                "implemented universal top-gradient and support-codimension "
                "sieve; singular lower-layer normal-cone multiplicities "
                "remain open"
            ),
            "artifact": (
                "artifacts/generated-results/"
                "hc4_quintic_infinity_rees_strata.json"
            ),
            "checker": "scripts/analyze_hc4_quintic_infinity_rees.py",
            "independent_rees_checker": (
                "scripts/verify_hc4_quintic_infinity_rees_strata.m2"
            ),
        },
        "hc4_rank3_vertex_colength": {
            "status": (
                "implemented lower-layer filtered-length obstruction; "
                "both codimension-four affine-degree-two/three rows excluded"
            ),
            "artifact": (
                "artifacts/generated-results/"
                "hc4_rank3_vertex_colength.json"
            ),
            "checker": "scripts/verify_hc4_rank3_vertex_colength.py",
            "independent_checker": (
                "scripts/verify_hc4_rank3_vertex_colength.m2"
            ),
        },
        "hc4_codimension_three_gradient_strata": {
            "status": (
                "implemented rank-two constant-kernel/Schur sieve and "
                "rank-three ordinary-singularity incidence"
            ),
            "artifact": (
                "artifacts/generated-results/"
                "hc4_codim3_gradient_strata.json"
            ),
            "checker": "scripts/verify_hc4_codim3_gradient_strata.py",
            "independent_checker": (
                "scripts/verify_hc4_codim3_gradient_strata.m2"
            ),
        },
        "hc4_binary_root_partition_segre": {
            "status": (
                "implemented PGS3 active-unit sieve on singular essential "
                "binary quintic root partitions"
            ),
            "generic_double_root_packet": {
                "forced_sigma2": 1,
                "rows_for_affine_degrees_2_and_3": [51, 50],
            },
            "artifact": (
                "artifacts/generated-results/"
                "hc4_binary_root_partition_segre.json"
            ),
            "checker": (
                "scripts/verify_hc4_binary_root_partition_segre.py"
            ),
            "independent_checker": (
                "scripts/verify_hc4_binary_root_partition_segre.m2"
            ),
        },
    },
    "top_degree_controls_and_open_records": top_degree_controls,
    "scope": (
        "The transform and compactification constructors are exact in all "
        "dimensions. Complete projective/Segre lists are recorded only "
        "where independently computed. A top generic degree determines "
        "only the weighted aggregate correction, not the individual Segre "
        "degrees. Boundary normalization data and Hessian rank data are "
        "therefore not silently promoted to Segre classes."
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: all-dimensional projective-degree/Segre transforms are inverse")
print("PASS: affine-gradient and full-polar constructors are distinct")
print("PASS: homogeneous integrability reconstructs the infinity potential")
print("PASS: attached the all-dimensional smooth-essential normal slice")
print("PASS: attached the all-dimensional singular-stratum DVR profile")
print("PASS: complete and top-degree-only family records stay separated")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
