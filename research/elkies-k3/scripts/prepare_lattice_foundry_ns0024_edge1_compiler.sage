#!/usr/bin/env sage-python
"""Export the exact abstract-to-equation handoff for NS0024 edge 1.

status: ACTIVE_COMPILER
claim: exact source marking, divisor identity, and resolved-I5 compiler profile
output: artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-edge1-compiler-preparation.json

This is deliberately equation-independent.  It binds q4/orbit1 to the
minimum-pole source basis and proves that the horizontal is P3 and that

    D = O + P3 + 2F - C2 - 2C3 - C4

on the normalized I5 chain.  A finite-field source family can therefore enter
the resolved-RR compiler without repeating a lattice orbit or marking search.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
BASIS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"
EDGE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-r13-route-step01.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-edge1-compiler-preparation.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

source = json.loads(SOURCE.read_text())
basis = json.loads(BASIS.read_text())
search = json.loads(EDGE.read_text())
assert source["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_MW4_SOURCE_AND_NIEMEIER_CERTIFICATE"
assert basis["status"] == "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS"
assert search["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert len(search["neighbors"]) == 1
edge = search["neighbors"][0]
assert (edge["q"], edge["orbit_index"], edge["old_fiber_degree"]) == (4, 1, 2)
assert edge["mw_projection"] == [-1, 0, 0, 0]

frame = matrix(ZZ, source["source"]["root_adapted_gram"])
ns = matrix(ZZ, 19, 19)
ns[0, 1] = ns[1, 0] = 1
ns[2:, 2:] = -frame
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
new_fibre = vector(ZZ, edge["fiber"])

p3_record = next(item for item in basis["basis"] if item["name"] == "P3")
assert p3_record["mw_quotient_coordinates"] == edge["mw_projection"]
assert p3_record["P_dot_O"] == 0
assert p3_record["components_I7_I5_I4"] == [2, 1, 1]
p3_frame = vector(ZZ, p3_record["frame_vector"])
p3_first = ZZ((p3_frame * frame * p3_frame - 2) // 2)
p3 = vector(ZZ, [p3_first, 1] + list(p3_frame))
assert p3 * ns * p3 == -2
assert p3 * ns * old_fibre == 1
assert p3 * ns * old_zero == 0

vertical = new_fibre - old_zero - p3 - 2 * old_fibre
assert vertical[:2] == vector(ZZ, [0, 0])
assert list(vertical[2:]) == [1, 0, 2, 1] + [0] * 13

# Deterministic simple-root chain order from the minimum-basis certificate.
chains = {
    "I7": (7, 12, 10, 9, 8, 11),
    "I5": (1, 0, 2, 3),
    "I4": (5, 4, 6),
}
i5_coefficients = tuple(-vertical[2 + index] for index in chains["I5"])
assert i5_coefficients == (0, -1, -2, -1)
component_pairings = vector(ZZ, new_fibre[2:]) * frame[:, :13]
assert tuple(component_pairings[index] for index in chains["I5"]) == (0, 0, 2, 0)

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-edge1-compiler-preparation.v1",
    "status": "PASS_EXACT_NS0024_EDGE1_COMPILER_PREPARATION",
    "ns_id": "NS0024",
    "source": {
        "root_type": "A3+A4+A6",
        "semistable_fibres": ["I7", "I5", "I4", "8I1"],
        "normalized_supports": {"I7": "0", "I5": "1", "I4": "infinity"},
        "mw_basis": [item["name"] for item in basis["basis"]],
        "mw_basis_pole_profile": [item["P_dot_O"] for item in basis["basis"]],
        "horizontal": {
            "name": "P3",
            "P_dot_O": p3_record["P_dot_O"],
            "components_I7_I5_I4": p3_record["components_I7_I5_I4"],
            "frame_vector": p3_record["frame_vector"],
            "mw_quotient_coordinates": p3_record["mw_quotient_coordinates"],
        },
    },
    "edge": {
        "edge_index": 1,
        "q": 4,
        "orbit_index": 1,
        "old_fibre_degree": 2,
        "new_fibre_in_source_ns": list(map(int, new_fibre)),
        "divisor_identity": "D=O+P3+2F-C2-2C3-C4 on the normalized I5 chain",
        "vertical_component_coefficients_I5_C1_to_C4": list(map(int, i5_coefficients)),
        "new_fibre_pairings_I5_C1_to_C4": [
            int(component_pairings[index]) for index in chains["I5"]
        ],
        "target_root_type": "A1+A2+A4+D5",
        "target_root_rank": 12,
        "target_mw_rank": 5,
        "target_root_determinant": 120,
        "target_semistable_additive_profile": "I1*+I5+I3+I2+7I1",
    },
    "resolved_RR": {
        "ambient_basis": ["1", "t", "t^2", "m"],
        "chord": "m=(y+y(P3))/(x-x(P3))",
        "supported_fibre": "I5@t=1",
        "marked_component": 1,
        "required_component_vanishing": {"2": 1, "3": 2, "4": 1},
        "adapter": "split_multiplicative_toric_chord_condition",
        "expected_dimensions": {"ambient": 4, "condition_rank": 2, "h0": 2},
        "quartic_adapter": "compile_resolved_degree_two_chord_hop",
        "completeness_basis": (
            "The abstract divisor has no vertical correction away from the displayed "
            "I5 subchain; the adapter evaluates every negative component and every "
            "required multiplicity layer on its exact split toric chart."
        ),
    },
    "model_input_contract": {
        "schema": "elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1",
        "required_section": "P3",
        "supported_parameter_counts": [0, 1],
        "coefficient_encoding": "Sage expressions over GF(p) or GF(p)(u), low to high in t",
    },
    "inputs": {
        "paths": [relative(path) for path in (SOURCE, BASIS, EDGE)],
        "sha256": {relative(path): digest(path) for path in (SOURCE, BASIS, EDGE)},
    },
    "proof_boundary": {
        "proved": (
            "Exact identification of q4/orbit1 with P3, its full source-NS divisor "
            "identity, and the complete resolved-I5 condition profile needed by the compiler."
        ),
        "not_proved": (
            "No finite-field MW4 family, resolved RR kernel on an equation, child "
            "Jacobian, characteristic-zero lift, or NS0024 equation is asserted here."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/prepare_lattice_foundry_ns0024_edge1_compiler.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = args.output.resolve()
if args.check:
    if output.read_text() != serialized:
        raise SystemExit("NS0024 edge-1 compiler preparation artifact is stale")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)

print(
    "NS0024EDGE1PREP|q=4|orbit=1|horizontal=P3|I5=0,-1,-2,-1|"
    "ambient=4|expected_rank=2|expected_h0=2|status=PASS",
    flush=True,
)
