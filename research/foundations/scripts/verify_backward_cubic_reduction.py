#!/usr/bin/env python3
"""Verify the backward cubic reduction primitives on MacFarlane's chain."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from audit_macfarlane_g20_dimension_reduction import (  # noqa: E402
    SOURCE_URL,
    build_maps,
)
from jcsearch.backward_cubic import (  # noqa: E402
    collision_compatible_covectors,
    common_image,
    companion_cancellation,
    fixed_linear_covectors,
    lift_point_to_nonzero_companion_slice,
    parametric_companion_cancellation,
    profile_from_cubic_components,
    project_point_to_restriction,
    restrict_fixed_covector,
    surviving_collision_pairs,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "backward_cubic_reduction_calibration.json"
)


def qtext(value: sp.Expr) -> str:
    return str(sp.cancel(value))


def main() -> None:
    # Generic fixed-covector regression in coordinates not adapted to the
    # covector.  In (u,v,s) coordinates the map is the composition of the
    # two elementary shears U=u+v^2, V=v+U^2, with s=x-2y+3z fixed.
    generic_variables = tuple(sp.symbols("generic_x0:3"))
    generic_x, generic_y, generic_z = generic_variables
    generic_level = generic_x - 2 * generic_y + 3 * generic_z
    generic_u_image = generic_x + generic_y**2
    generic_v_image = generic_y + generic_u_image**2
    generic_map = (
        generic_u_image,
        generic_v_image,
        sp.expand(
            (
                generic_level
                - generic_u_image
                + 2 * generic_v_image
            )
            / 3
        ),
    )
    generic_correction = tuple(
        sp.expand(image - variable)
        for image, variable in zip(generic_map, generic_variables)
    )
    generic_covector = sp.Matrix([1, -2, 3])
    assert sp.expand(
        (generic_covector.T * sp.Matrix(generic_correction))[0]
    ) == 0
    assert sp.expand(
        sp.Matrix(generic_map).jacobian(generic_variables).det() - 1
    ) == 0
    generic_fixed = fixed_linear_covectors(
        generic_correction, generic_variables
    )
    assert len(generic_fixed) == 1
    assert sp.Matrix.hstack(
        generic_fixed[0], generic_covector
    ).rank() == 1
    generic_restriction = restrict_fixed_covector(
        generic_map,
        generic_variables,
        generic_covector,
        sp.Integer(5),
        prefix="generic_slice",
    )
    assert sp.expand(
        sp.Matrix(generic_restriction.restricted_map)
        .jacobian(generic_restriction.variables)
        .det()
        - 1
    ) == 0

    # Generic companion regression.  No Keller hypothesis is needed for the
    # identity: the full specialized determinant equals det(DF) exactly.
    companion_x = tuple(sp.symbols("companion_x0:2"))
    companion_y = tuple(sp.symbols("companion_y0:1"))
    companion_q = (
        companion_x[0] * companion_x[1],
        companion_x[1] ** 2,
    )
    companion_c = (
        companion_x[0] ** 3 - companion_x[1] ** 3,
    )
    companion_b = sp.Matrix([[1], [2]])
    generic_companion = companion_cancellation(
        companion_x,
        companion_y,
        companion_q,
        companion_b,
        companion_c,
    )
    determinant_base = sp.expand(
        sp.Matrix(generic_companion.base_map)
        .jacobian(companion_x)
        .det()
    )
    determinant_specialized = sp.expand(
        sp.Matrix(generic_companion.specialized_map)
        .jacobian(companion_x + companion_y)
        .det()
    )
    assert sp.expand(determinant_specialized - determinant_base) == 0
    companion_parameter = sp.Symbol("companion_t")
    generic_parametric = parametric_companion_cancellation(
        companion_x,
        companion_y,
        companion_parameter,
        companion_q,
        companion_b,
        companion_c,
    )
    determinant_parent_family = sp.expand(
        sp.Matrix(generic_parametric.parent_map)
        .jacobian(companion_x + companion_y + (companion_parameter,))
        .det()
    )
    determinant_base_family = sp.expand(
        sp.Matrix(generic_parametric.scaled_base_family)
        .jacobian(companion_x)
        .det()
    )
    assert sp.expand(
        determinant_parent_family - determinant_base_family
    ) == 0
    assert sp.expand(
        sp.Matrix(generic_parametric.special_fiber)
        .jacobian(companion_x + companion_y)
        .det()
        - 1
    ) == 0

    # Collision compatibility is computed in the full fixed-covector space.
    coordinate_covectors = (
        sp.Matrix([1, 0]),
        sp.Matrix([0, 1]),
    )
    compatible_line = collision_compatible_covectors(
        coordinate_covectors,
        [sp.Matrix([0, 0]), sp.Matrix([0, 1])],
    )
    assert len(compatible_line) == 1
    assert sp.Matrix.hstack(
        compatible_line[0], coordinate_covectors[0]
    ).rank() == 1

    data = build_maps()
    x = data["x"]
    w = data["w"]
    tau = data["tau"]
    variables20 = data["variables20"]
    g20 = data["G20"]
    h20 = data["H20"]
    p20 = data["p20"]
    q20 = data["q20"]

    # The sole fixed covector is tau, and the collision lies on tau=1.
    fixed = fixed_linear_covectors(h20, variables20)
    assert len(fixed) == 1
    tau_covector = sp.eye(20)[:, 19]
    assert fixed[0] == tau_covector
    compatible = collision_compatible_covectors(fixed, [p20, q20])
    assert len(compatible) == 1 and compatible[0] == tau_covector

    restriction = restrict_fixed_covector(
        g20, variables20, tau_covector, sp.Integer(1), prefix="backward_m"
    )
    assert restriction.pivot == 19
    p19 = project_point_to_restriction(p20, restriction)
    q19 = project_point_to_restriction(q20, restriction)
    assert surviving_collision_pairs([p19, q19]) == ((0, 1),)
    assert common_image(
        restriction.restricted_map, restriction.variables, [p19, q19]
    ) == p19

    # Rename the generic slice variables back to (x,w) for literal comparison
    # with the companion block.
    rename = dict(zip(restriction.variables, x + w))
    restricted_xw = tuple(
        sp.expand(value.subs(rename, simultaneous=True))
        for value in restriction.restricted_map
    )
    cancellation = companion_cancellation(
        x, w, data["R"], data["B"], data["gamma"]
    )
    assert restricted_xw == cancellation.specialized_map
    assert cancellation.base_map == tuple(data["F13"])

    # The source shear moves both lifted collision points onto y=0.  The
    # remaining base map retains the two-point collision.
    p13 = data["p13"]
    q13 = data["q13"]
    source_p = sp.Matrix(cancellation.source_shear).subs(dict(zip(x + w, p19)))
    source_q = sp.Matrix(cancellation.source_shear).subs(dict(zip(x + w, q19)))
    assert source_p == p13.col_join(sp.zeros(6, 1))
    assert source_q == q13.col_join(sp.zeros(6, 1))
    assert common_image(cancellation.base_map, x, [p13, q13]) == p13

    # The full parent is an isotrivial family over tau != 0.  A second
    # rational slice provides an exact collision distinct from the published
    # tau=1 normalization, while tau=0 is triangular and injective.
    parametric = parametric_companion_cancellation(
        x,
        w,
        tau,
        data["R"],
        data["B"],
        data["gamma"],
    )
    assert parametric.parent_map == tuple(g20)
    second_level = sp.Integer(2)
    p20_level2 = lift_point_to_nonzero_companion_slice(
        p13, x, data["gamma"], second_level
    )
    q20_level2 = lift_point_to_nonzero_companion_slice(
        q13, x, data["gamma"], second_level
    )
    assert surviving_collision_pairs(
        [p20_level2, q20_level2]
    ) == ((0, 1),)
    level2_image = common_image(
        parametric.parent_map,
        x + w + (tau,),
        [p20_level2, q20_level2],
    )
    assert level2_image == sp.Matrix(
        [sp.cancel(value / second_level) for value in p13]
        + [sp.Integer(0)] * 6
        + [second_level]
    )
    collision_parameter = sp.Symbol("collision_a", nonzero=True)
    symbolic_lift = lift_point_to_nonzero_companion_slice(
        p13,
        x,
        data["gamma"],
        collision_parameter,
    )
    assert sp.Matrix(
        [
            sp.expand(collision_parameter * symbolic_lift[index])
            for index in range(13)
        ]
    ) == p13
    gamma_at_p = sp.Matrix(data["gamma"]).subs(dict(zip(x, p13)))
    assert sp.Matrix(
        [
            sp.expand(
                collision_parameter**3 * symbolic_lift[13 + index]
            )
            for index in range(6)
        ]
    ) == gamma_at_p

    profile = profile_from_cubic_components(13, data["C"], x)
    assert profile.cubic_output_rank == 6
    assert profile.homogeneous_dimension == 20
    assert profile.direct_cubic_key == (13, 20, 6)
    assert profile.homogeneous_key == (20, 13, 6)

    # Pair-aware bookkeeping must not require every collision point to survive.
    toy_points = [
        sp.Matrix([0, 0]),
        sp.Matrix([1, 0]),
        sp.Matrix([0, 1]),
    ]
    projection = sp.Matrix([[1, 0]])
    assert surviving_collision_pairs(toy_points, projection) == ((0, 1), (1, 2))

    artifact = {
        "format": "backward-cubic-reduction-calibration-v1",
        "external_source": SOURCE_URL,
        "status": (
            "exact calibration of restriction and companion cancellation; "
            "no 12-variable map is constructed"
        ),
        "generic_theorem_tests": {
            "noncoordinate_fixed_covector": [1, -2, 3],
            "ambient_jacobian_determinant": "1",
            "restricted_jacobian_determinant": "1",
            "companion_determinant_identity": (
                "det D(x+Q+B*y,y-c)=det D(x+Q+B*c)"
            ),
            "parametric_determinant_identity": (
                "det D(V_t)=det D(t^-1*F(t*x))"
            ),
            "special_fiber": "(x,y-c(x)) with determinant 1",
            "collision_compatible_covector_subspace_dimension": 1,
        },
        "chain": [
            {
                "operation": "fixed-covector restriction",
                "source_dimension": 20,
                "target_dimension": 19,
                "covector": "tau",
                "level": "1",
                "collision_pair_survives": True,
            },
            {
                "operation": "stable companion cancellation",
                "source_dimension": 19,
                "target_dimension": 13,
                "cancelled_companion_variables": 6,
                "factorization": "M=A_B o (F13 x I_6) o S_gamma",
                "collision_pair_survives": True,
            },
        ],
        "parametric_family": {
            "relative_factorization": (
                "V=A_(tau^2*B) o (E_tau x I_6) o S_gamma"
            ),
            "nonzero_slice_identity": "E_tau(x)=tau^-1*F13(tau*x)",
            "nonzero_slice_consequence": (
                "every tau!=0 slice is stably equivalent to F13"
            ),
            "special_slice_tau_0": (
                "(x,w-gamma(x)) is a triangular automorphism"
            ),
            "second_exact_collision_level": "2",
            "normalization_policy": (
                "every collision in the parent has tau!=0 and may be "
                "scaled to tau=1"
            ),
            "collision_scheme_over_nonzero_parameter": (
                "stably equivalent to Collision(F13) x G_m x A^6"
            ),
            "collision_arc_weights": {
                "tau": 1,
                "base_x": -1,
                "cubic_companion_w": -3,
                "weighted_boundary_coordinates": "X=tau*x, W=tau^3*w",
            },
        },
        "terminal_profile": {
            "base_dimension": profile.base_dimension,
            "cubic_output_rank": profile.cubic_output_rank,
            "rank_compressed_homogeneous_dimension": profile.homogeneous_dimension,
            "direct_cubic_key": list(profile.direct_cubic_key),
            "homogeneous_key": list(profile.homogeneous_key),
        },
        "search_policy": {
            "direct_degree_three_archive": (
                "order by (base dimension, homogeneous dimension, cubic rank)"
            ),
            "cubic_homogeneous_archive": (
                "order by (homogeneous dimension, base dimension, cubic rank)"
            ),
            "collision_policy": "retain a terminal when any one collision pair survives",
            "pre_specialization_operation": (
                "restrict every collision-compatible fixed linear covector"
            ),
        },
        "first_open_problem": {
            "name": "MacFarlane F13 to a degree-three Keller collision in dimension 12",
            "success_condition": (
                "an exact determinant-preserving upstream schedule with a surviving "
                "collision pair and base dimension at most 12"
            ),
            "secondary_homogeneous_gates": [
                "base dimension 12 with cubic-output rank at most 6",
                "base dimension 13 with cubic-output rank at most 5",
            ],
        },
        "collision_image_F13": [qtext(value) for value in p13],
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()

    print("PASS backward cubic: fixed tau=1 restricts MacFarlane G20 from 20D to 19D")
    print("PASS backward cubic: companion cancellation splits 19D as F13 x I6")
    print("PASS backward cubic: the exact collision pair survives both operations")
    print("PASS backward cubic: every nonzero parent slice scales to F13")
    print("PASS backward cubic: tau=0 is triangular and tau=2 has an exact collision")
    print("PASS backward cubic: direct and homogeneous terminal objectives are distinct")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
