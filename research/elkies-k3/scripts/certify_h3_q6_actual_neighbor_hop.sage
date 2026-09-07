#!/usr/bin/env sage -python
"""Assemble the complete exact certificate for the first H3 q=6 neighbour.

This is the first end-to-end regression target of the elliptic-neighbour
compiler.  It combines the actual resolved H92 RR cover with the exact
genus-one/Jacobian calculation and the Neron--Severi Weyl transport.  The
certificate makes no claim to have produced minimized coordinates for every
rank-three section; it certifies the requested first-hop invariants.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
PREFLIGHT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-compiler-preflight.json"
RR_COVER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-resolved-rr-cover.json"
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
COMPONENTS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-component-sections.json"
COMPONENT_POINTS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-e7-infinity-sections.json"
TRANSPORT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-weyl-section-transport.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-actual-neighbor-hop.json"
PREFLIGHT_SHA256 = "590d981c98b4f4761a47af244d05ecaba767ea8904f9ed6633f06629dc8ac755"
RR_COVER_SHA256 = "928955731227eead4006392cbe13fc5b95a1e2e951acbc2834eea6bdc1bbaacc"
# Replayed with SageMath 10.9 from certify_h92_q6_child_jacobian.sage; the
# preceding c57... pin was stale, while the semantic q6 hop replay against
# this artifact passes every exact child and transport assertion below.
CHILD_SHA256 = "5eb43d9a0d04195e7a6e38ebd337b0e10a3b1a2eb9246a3b02cce4331bcd36ac"
COMPONENTS_SHA256 = "335a9cb6c1060ac170c063f99bb02d4c4357fa2426d37b4dc3efd447ac2b62ad"
COMPONENT_POINTS_SHA256 = "156821384b45fd5e731dce130686b549030d64d450ece78bfb9f9083bbaf3005"
TRANSPORT_SHA256 = "c4b7e38f0ea9fc3f748200ca9923ea3ffe5c0028c979e5f81be6507954d7c822"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--preflight", type=Path, default=PREFLIGHT)
parser.add_argument("--rr-cover", type=Path, default=RR_COVER)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--components", type=Path, default=COMPONENTS)
parser.add_argument("--component-points", type=Path, default=COMPONENT_POINTS)
parser.add_argument("--transport", type=Path, default=TRANSPORT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

expected_hashes = {
    PREFLIGHT: PREFLIGHT_SHA256, RR_COVER: RR_COVER_SHA256, CHILD: CHILD_SHA256,
    COMPONENTS: COMPONENTS_SHA256, COMPONENT_POINTS: COMPONENT_POINTS_SHA256,
    TRANSPORT: TRANSPORT_SHA256,
}
for path, expected in expected_hashes.items():
    selected = {
        PREFLIGHT: args.preflight, RR_COVER: args.rr_cover, CHILD: args.child,
        COMPONENTS: args.components, COMPONENT_POINTS: args.component_points,
        TRANSPORT: args.transport,
    }[path]
    if selected == path:
        assert digest(selected) == expected

preflight = json.loads(args.preflight.read_text())
rr_cover = json.loads(args.rr_cover.read_text())
child = json.loads(args.child.read_text())
components = json.loads(args.components.read_text())
component_points = json.loads(args.component_points.read_text())
transport = json.loads(args.transport.read_text())
assert preflight["status"] == "PASS_EXACT_Q6_ACTUAL_E7_LOCAL_INPUTS"
assert rr_cover["status"] == "PASS_EXACT_Q6_ACTUAL_RESOLVED_RR_COVER"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert components["status"] == "PASS_EXACT_Q6_COMPONENT_SECTION_IDENTITIES"
assert component_points["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert transport["status"] == "PASS_EXACT_Q6_WEYL_SECTION_TRANSPORT"
assert component_points["inputs"]["child_jacobian"]["sha256"] == CHILD_SHA256
assert component_points["source"] == {
    "old_base": "u=infinity (the III* fibre)",
    "candidate_curves": ["old E7_7", "old affine E7"],
    "affine_E7_chart": "Z=0, U=r^2, Y=r^3 in the first x-chart of the E7 resolution",
    "affine_E7_sign": "plus",
    "E7_7_sign": "minus",
    "boundary": "The affine sign is fixed by its explicit E7 chart and the normalized discriminant square root. The complementary infinity point is E7_7 after the E7 resolution graph certificate, not by a Kodaira-symbol inference.",
}
point_signs = {entry["sign"]: entry for entry in component_points["sections"]}
assert set(point_signs) == {"plus", "minus"}
for point in point_signs.values():
    assert point["x_denominator_coefficients_low_to_high"] == ["1"]
    assert point["y_denominator_coefficients_low_to_high"] == ["1"]

neighbour = preflight["neighbour"]
assert neighbour["divisor_squared"] == 0 and neighbour["primitive"]
assert len(neighbour["recorded_weyl_reflections"]) == 22
assert all(value >= 0 for _, value in neighbour["declared_wall_pairings"])

vertical = rr_cover["vertical_condition_matrix"]
assert (vertical["ambient_dimension"], vertical["rows"], vertical["rank"],
        vertical["codimension"], vertical["kernel_dimension"], vertical["h0_D"]) == (10, 8, 8, 8, 2, 2)

assert child["root_data"] == {"rank": 14, "determinant": 3, "type": "E8+E6"}
assert child["mordell_weil_rank_if_rho_19"] == 3
assert components["root_lattice"] == {"rank": 14, "determinant": 3, "type": "E8+E6"}
expected_height = matrix(QQ, [[QQ(8)/3, QQ(1)/3, -1], [QQ(1)/3, QQ(8)/3, 1], [-1, 1, 46]])
height = matrix(QQ, [[QQ(value) for value in row] for row in transport["child"]["height_gram"]])
assert height == expected_height
assert transport["transported_section_mw_projections"]["formulas"] == [
    "4*(-P1)", "(-P1)", "22*(-P1)-P2",
]
assert transport["third_vertical_correction"]["transported_section_new_fiber_degree"] == 1

# Normalize the independently derived layers into the reusable equation-level
# hop contract.  The preflight only proves nonnegativity on its declared wall
# list, and the name retains that scope rather than silently upgrading it to a
# new global nef assertion.
section_degrees = [
    components["curves"]["D_intersections"]["E7_7"],
    components["curves"]["D_intersections"]["affine_E7"],
    transport["third_vertical_correction"]["transported_section_new_fiber_degree"],
]
assert section_degrees == [1, 1, 1]
exec(compile(CORE.read_text(), str(CORE), "exec"))
shioda_tate = certify_shioda_tate_discriminant(
    child["root_data"]["determinant"], height.rows(),
    torsion_order=1, expected_ns_discriminant=948,
)
assert shioda_tate["height_determinant"] == 316
assert shioda_tate["absolute_ns_discriminant"] == 948
third_correction = transport["third_vertical_correction"]
component_transport = certify_component_pairing_transport(
    {
        "old_O": components["curves"]["child_zero"],
        "old_E7_7": components["curves"]["first_section"],
        "old_affine_E7": components["curves"]["second_pre_difference"],
    },
    {
        "old_O": components["curves"]["D_intersections"]["old_O"],
        "old_E7_7": components["curves"]["D_intersections"]["E7_7"],
        "old_affine_E7": components["curves"]["D_intersections"]["affine_E7"],
    },
    third_correction["basis"][:7],
    third_correction["old_E7_component_intersections"],
    {
        "horizontal_degree": third_correction["old_horizontal_new_fiber_degree"],
        "correction_degree": third_correction["correction_new_fiber_degree"],
        "transported_section_degree": third_correction["transported_section_new_fiber_degree"],
    },
    {
        "section_sources": {
            "old_O": "old O",
            "old_E7_7": "old E7_7 exceptional component",
            "old_affine_E7": "old affine E7 component",
        },
        "source_fiber_degrees": {"old_O": 1, "old_E7_7": 1, "old_affine_E7": 1},
        "pairing_basis": ["old_E7_{}".format(index) for index in range(1, 8)],
        "resolved_pairings": {
            "horizontal_part": [0, 0, 0, 0, 0, 0, 0],
            "transported_section": [0, 0, 0, 0, 0, 0, 22],
            "vertical_correction": [0, 0, 0, 0, 0, 0, 22],
        },
        "vertical_correction": {
            "horizontal_degree": 4812, "correction_degree": -4811,
            "transported_section_degree": 1,
        },
    },
)
hop = certify_exact_neighbor_hop(
    {
        "square": neighbour["divisor_squared"],
        "primitive": neighbour["primitive"],
        "old_fiber_degree": neighbour["old_fiber_degree"],
        "nef_on_declared_walls": all(value >= 0 for _, value in neighbour["declared_wall_pairings"]),
        "weyl_reflection_count": len(neighbour["recorded_weyl_reflections"]),
    },
    {
        "complete_resolved_chart_cover": rr_cover["compiler_replay"]["complete_resolved_chart_cover"],
        "ambient_dimension": vertical["ambient_dimension"],
        "condition_rank": vertical["rank"],
        "condition_codimension": vertical["codimension"],
        "kernel_dimension": vertical["kernel_dimension"],
        "h0": vertical["h0_D"],
    },
    {
        "root_lattice": child["root_data"]["type"],
        "root_rank": child["root_data"]["rank"],
        "root_determinant": child["root_data"]["determinant"],
        "mordell_weil_rank": child["mordell_weil_rank_if_rho_19"],
    },
    {
        "height_gram": height.rows(),
        "section_words": transport["transported_section_mw_projections"]["formulas"],
        "section_new_fiber_degrees": section_degrees,
    },
    {
        "rr": {
            "ambient_dimension": 10, "condition_rank": 8,
            "condition_codimension": 8, "kernel_dimension": 2, "h0": 2,
        },
        "child": {
            "root_lattice": "E8+E6", "root_rank": 14,
            "root_determinant": 3, "mordell_weil_rank": 3,
        },
        "height_gram": expected_height.rows(),
        "section_words": ["4*(-P1)", "(-P1)", "22*(-P1)-P2"],
    },
)

payload = {
    "schema": "elkies-k3.h3-q6-actual-neighbor-hop.v1",
    "status": "PASS_EXACT_H3_Q6_ACTUAL_NEIGHBOR_HOP",
    "inputs": {
        "preflight": {"path": str(args.preflight.relative_to(ROOT)), "sha256": digest(args.preflight)},
        "actual_rr_cover": {"path": str(args.rr_cover.relative_to(ROOT)), "sha256": digest(args.rr_cover)},
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "component_sections": {"path": str(args.components.relative_to(ROOT)), "sha256": digest(args.components)},
        "component_points": {"path": str(args.component_points.relative_to(ROOT)), "sha256": digest(args.component_points)},
        "weyl_transport": {"path": str(args.transport.relative_to(ROOT)), "sha256": digest(args.transport)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
    "divisor": hop["divisor"],
    "rr": {
        "ambient_dimension": hop["rr"]["ambient_dimension"],
        "condition_rank": hop["rr"]["condition_rank"],
        "codimension": hop["rr"]["condition_codimension"],
        "kernel_dimension": hop["rr"]["kernel_dimension"], "h0": hop["rr"]["h0"],
    },
    "child": hop["child"],
    "neron_severi": {
        "height_determinant": int(shioda_tate["height_determinant"]),
        "root_determinant": int(shioda_tate["root_determinant"]),
        "torsion_glue_index": int(shioda_tate["torsion_glue_index"]),
        "absolute_discriminant": int(shioda_tate["absolute_ns_discriminant"]),
    },
    "neron_severi_component_transport": component_transport,
    "resolved_component_point_transport": {
        "old_affine_E7": {
            "resolved_chart": component_points["source"]["affine_E7_chart"],
            "binary_quartic_infinity_sign": component_points["source"]["affine_E7_sign"],
        },
        "old_E7_7": {
            "binary_quartic_infinity_sign": component_points["source"]["E7_7_sign"],
        },
        "child_point_signs": sorted(point_signs),
        "coordinate_artifact": {
            "path": str(args.component_points.relative_to(ROOT)),
            "sha256": digest(args.component_points),
        },
    },
    "transported_sections": {
        "source_words": hop["transport"]["section_words"],
        "height_gram": hop["transport"]["height_gram"],
        "new_fiber_degrees": hop["transport"]["section_new_fiber_degrees"],
    },
    "conclusion": "All four first-target q=6 checks are exact and use an actual all-edge H92 resolved E7 module cover: predicted codimension, h0(D)=2, E8+E6/MW3 child, and the transported rank-three Gram matrix. The first two low-height sections are tied to exact child-Jacobian coordinates by the explicit E7 affine chart and its binary-quartic infinity signs. The final certificate also pins three source section degrees, the exact Neron--Severi old-E7 pairing rows of the third correction, and its 4812-4811=1 degree balance. These third-correction pairing rows are not asserted to be a new resolved-chart trace. The torsion-free Shioda--Tate calculation has absolute Neron--Severi discriminant 948.",
    "boundary": "This completes the first q=6 hop certificate. General compilation of arbitrary divisor data and minimized equation-level coordinates for all transported sections remains future compiler work.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H3Q6ACTUALHOP|D2=0|rank=8|codimension=8|h0=2|"
    "child=E8+E6/MW3|gram=PASS|status=PASS_EXACT_H3_Q6_ACTUAL_NEIGHBOR_HOP",
    flush=True,
)
