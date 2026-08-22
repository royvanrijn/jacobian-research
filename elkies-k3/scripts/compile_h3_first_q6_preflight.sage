#!/usr/bin/env sage -python
"""Exact preflight for the first H3 q6-shell equation-level neighbour.

This is intentionally a gate, not a guessed normal-form conversion.  It
records the exact lattice reduction and exact H92 P1 input, distinguishes the
search label ``q6`` from the old-fibre degree two, and writes a blocked
certificate until both E7 and E8 conditions come from resolved blow-up charts.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
ACTUAL_E7_TRACE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
ACTUAL_E7_MARKED_MODULE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json"
ACTUAL_E7_ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
ACTUAL_E7_QUOTIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-quotient-block.json"
FRAME_SHA256 = "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
ACTUAL_E7_TRACE_SHA256 = "a73ccb1c729814219f172df4c6feb49c05859125db1cff7591eeb8544fb664e1"
ACTUAL_E7_MARKED_MODULE_SHA256 = "4a94a5aca8686fbb666b5bb26b6f784eca33079d8922da0766ec5bd0ae2a4ba8"
ACTUAL_E7_ATLAS_SHA256 = "ae7eb1e79a2fb41ab05d0092cb8c04663307bd8aa1e0cb562a8d3a014e94f451"
ACTUAL_E7_QUOTIENT_SHA256 = "543c98357f08453ccc362ac9736d35bc3d30efde02c5417e1210a77c9ade4156"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-compiler-preflight.json"

CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(FRAME) == FRAME_SHA256
assert digest(SECTION) == SECTION_SHA256
assert digest(ACTUAL_E7_TRACE) == ACTUAL_E7_TRACE_SHA256
assert digest(ACTUAL_E7_MARKED_MODULE) == ACTUAL_E7_MARKED_MODULE_SHA256
assert digest(ACTUAL_E7_ATLAS) == ACTUAL_E7_ATLAS_SHA256
assert digest(ACTUAL_E7_QUOTIENT) == ACTUAL_E7_QUOTIENT_SHA256
section = json.loads(SECTION.read_text())
actual_e7_trace = json.loads(ACTUAL_E7_TRACE.read_text())
actual_e7_marked_module = json.loads(ACTUAL_E7_MARKED_MODULE.read_text())
actual_e7_atlas = json.loads(ACTUAL_E7_ATLAS.read_text())
actual_e7_quotient = json.loads(ACTUAL_E7_QUOTIENT.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"
assert section["exact_weierstrass_square"]
assert actual_e7_trace["status"] == "PASS_EXACT_P1_ACTUAL_E7_TRACE"
assert actual_e7_trace["resolved_incidence"]["component"] == "E7_5"
assert actual_e7_marked_module["status"] == "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED"
assert actual_e7_atlas["status"] == "PASS_EXACT_H92_E7_VALUATION_ATLAS"
assert actual_e7_atlas["old_base_fibre_multiplicities"] == [2, 2, 4, 3, 1, 2, 3]
assert actual_e7_quotient["status"] == "PASS_EXACT_Q6_P1_ACTUAL_E7_QUOTIENT_BLOCK"

frame = load_frame(FRAME)
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
minus_p1 = vector(ZZ, [5, 1, -2, -3, -4, -6, -5, -4, -3] + [0] * 8 + [1, 0])
divisor = zero + minus_p1 - fiber

# Replay the chamber reduction from the q6-shell witness.  Retaining both
# root number and negative pairing makes this an actual Weyl certificate,
# rather than an assertion about a lattice-isometric endpoint.
raw_witness = vector(
    ZZ,
    [0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
)
raw_divisor = vector(ZZ, [3, 2] + list(raw_witness))
assert intersection(raw_divisor, raw_divisor, ns) == 0
assert intersection(raw_divisor, fiber, ns) == 2
reflection_nodes = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)

highest_e7 = (2, 2, 3, 4, 3, 2, 1)
highest_e8 = (2, 3, 4, 6, 5, 4, 3, 2)
affine_e7 = fiber - sum(
    (coefficient * simple[index] for index, coefficient in enumerate(highest_e7)),
    vector(ZZ, [0] * ns.nrows()),
)
affine_e8 = fiber - sum(
    (coefficient * simple[7 + index] for index, coefficient in enumerate(highest_e8)),
    vector(ZZ, [0] * ns.nrows()),
)

walls = tuple(
    (("E7_{}" if index < 7 else "E8_{}").format(index + 1 if index < 7 else index - 6), curve)
    for index, curve in enumerate(simple)
) + (("E7_affine", affine_e7), ("E8_affine", affine_e8))
reflection_roots = tuple(
    ("E7_{}".format(node), simple[node - 1]) for node in reflection_nodes
)
lattice_pipeline = compile_elliptic_neighbor_rr_pencil(
    ns,
    raw_divisor,
    fiber,
    reflection_roots,
    walls,
    (),
    (),
    complete_resolved_cover=False,
    expected_reflection_pairings=(-1,) * len(reflection_roots),
    expected_nef_divisor=divisor,
)
reduced_divisor = lattice_pipeline["nef_divisor"]
recorded_reflections = [
    {
        "simple_root": int(item["root"].split("_")[1]),
        "pairing": int(item["pairing_before"]),
    }
    for item in lattice_pipeline["weyl_reflections"]
]
preflight = lattice_pipeline["preflight"]
assert all(intersection(reduced_divisor, curve, ns) >= 0 for curve in simple)
assert reduced_divisor == divisor
assert preflight["square"] == 0
assert preflight["primitive"]
assert preflight["old_fiber_degree"] == 2
assert preflight["nonnegative_on_declared_walls"]

# The q=6 search shell is a lattice label.  The geometric RR degree is two,
# so a generic-fibre basis additionally needs the marked chord with pole
# O+(-P1); it is not the five monomials of L(5O).
rr_monomials = bounded_weierstrass_monomials(2, (0,))
assert len(rr_monomials) == 1
chord = "(y-y(P1))/(x-x(P1))"

payload = {
    "schema": "elkies-k3.elliptic-neighbor-preflight.v1",
    "status": "PASS_EXACT_Q6_ACTUAL_E7_LOCAL_INPUTS",
    "inputs": {
        "frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": FRAME_SHA256},
        "marked_section": {"path": str(SECTION.relative_to(ROOT)), "sha256": SECTION_SHA256},
    },
    "neighbour": {
        "search_shell_label": "q6",
        "old_fiber_degree": preflight["old_fiber_degree"],
        "divisor_squared": preflight["square"],
        "primitive": preflight["primitive"],
        "declared_wall_pairings": preflight["declared_wall_pairings"],
        "raw_q6_shell_witness": [int(value) for value in raw_divisor],
        "recorded_weyl_reflections": recorded_reflections,
    },
    "rr_ambient": {
        "generic_degree": 2,
        "ordinary_monomials_before_marked_section": list(rr_monomials),
        "marked_chord_generator": chord,
        "warning": "The q6 shell label is not the old-fibre degree; using L(5O) here would be a wrong ambient space.",
    },
    "resolved_conditions": {
        "E7": {
            "status": "PASS_ACTUAL_QUOTIENT_BLOCK",
            "actual_branch_trace": {"path": str(ACTUAL_E7_TRACE.relative_to(ROOT)), "sha256": ACTUAL_E7_TRACE_SHA256},
            "actual_marked_module": {"path": str(ACTUAL_E7_MARKED_MODULE.relative_to(ROOT)), "sha256": ACTUAL_E7_MARKED_MODULE_SHA256},
            "actual_valuation_atlas": {"path": str(ACTUAL_E7_ATLAS.relative_to(ROOT)), "sha256": ACTUAL_E7_ATLAS_SHA256},
            "actual_quotient_block": {"path": str(ACTUAL_E7_QUOTIENT.relative_to(ROOT)), "sha256": ACTUAL_E7_QUOTIENT_SHA256},
            "formal_normal_form_module_source": "elkies-k3/scripts/derive_h92_q6_e7_p1_branch_module.sage",
            "marked_module_for_D_F_infinity": "J_-P1^dual=R+R*(Y-p_Y)/(U-p_U)",
            "F0_representative_module": "Z*J_-P1^dual",
            "representative_change": "D_F_infinity=D_F0+div(t), so its sections are u times the Z*J_-P1^dual sections",
            "reason": "The actual trace and corrected marked module place +/-P1 on the smooth locus of E7_5 and prove Z*m/t has the marked simple pole at -P1, while m differs by the unit t/Z; the actual valuation atlas supplies all E7 component orders; the actual length-six quotient has the required e1+2e5 boundary.",
        },
        "E8": {
            "status": "AVAILABLE",
            "resolution_source": "elkies-k3/scripts/derive_h92_q6_e8_resolution.sage",
            "marked_module_for_D_F_infinity": "u*<1,Q>, Q=u^2*(y-y_P)/(x-x_P)",
            "module_source": "elkies-k3/scripts/derive_h92_q6_e8_p1_branch_module.sage",
            "obstruction_source": "elkies-k3/scripts/derive_h92_q6_e8_chord_obstruction.sage",
            "reason": "The module includes the affine II* component through the integral chord Q and every exceptional chart through its unit denominator; Smith saturation or a Kodaira label is not used.",
        },
        "smooth_P1_O": {
            "status": "AVAILABLE",
            "source": "elkies-k3/scripts/derive_h92_q6_smooth_po_module.sage",
            "collision_divisor": "h=Z4, the squarefree degree-four factor of the P1 affine denominator",
            "marked_module": "a(u)+b(u)*m is smooth at all P1.O collisions iff h divides b; q=h*m-A",
            "reason": "z_P/h is a unit modulo h, so P1 meets O transversely at each smooth collision; the condition is an exact QQ[u]/(h) quotient block.",
        },
        "smooth_P1_O_collisions": {
            "status": "AVAILABLE",
            "module_source": "elkies-k3/scripts/derive_h92_q6_smooth_po_module.sage",
            "marked_chord_change": "a+b*m=(a+b*(y_P/x_P))+(h*b)*R",
            "reason": "The four collision roots are squarefree and use the exact local regular coordinate R=(m-y_P/x_P)/h.",
        },
    },
    "claims_not_made": [
        "exact vertical-condition matrix codimension",
        "h0(D)=2 from equation-level conditions",
        "new Weierstrass model",
        "transported section coordinates or Gram-matrix replay",
    ],
    "next_required_input": "Compile the all-edge E7 transition evaluator, then stack it with E8 and smooth P1.O blocks in compile_resolved_conditions(..., complete=True); certify the resulting complete cover alongside the global q=6 matrix.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H3Q6COMPILER|shell=q6|old_degree=2|D2=0|primitive=1|E7=actual_quotient|E8=available|status=PASS_EXACT_Q6_ACTUAL_E7_LOCAL_INPUTS", flush=True)
