#!/usr/bin/env sage -python
"""Certify the low-trace route from exact D12 E5 to the A11 bridge M.

The exact parent section E5 has selected-child Abel--Jacobi vector
``(13,-20,12,-12,0,1)`` and degree 13 over the A11 base.  In the pinned
equation ordering, the shortest useful expression for
``M=(1,0,0,0,0,1)`` is

    M = AJ(E5) - 4*S7 - 2*S14 - 2*S17 + 2*Qminus,

where S7 and S17 are degree-one old identity-shell curves, S14 has degree
three, and Qminus is one of the two exact pointed-opposite markings.  Hence
only L(14O) and L(4O) univariate traces are needed.  This file certifies the
lattice word and source degrees; it does not promote the equation arrow.
"""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

E5 = LOCAL / "q24-d12-missing-e5-section-qq.json"
PARENT_ROUTE = GENERATED / "elkies-k3-h3-q24-d12-degree14-aj-missing-parent-route.json"
BRIDGE = LOCAL / "q24-a11-target-coset-bridge.json"
POINTED = GENERATED / "elkies-k3-h3-a11-pointed-opposite-mw-candidates.json"
MARKING = LOCAL / "q24-a11-equation-marking-orbit64-mod100003.json"
INPUTS = (E5, PARENT_ROUTE, BRIDGE, POINTED, MARKING)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

e5 = json.loads(E5.read_text())
route = json.loads(PARENT_ROUTE.read_text())
bridge = json.loads(BRIDGE.read_text())
pointed = json.loads(POINTED.read_text())
marking = json.loads(MARKING.read_text())
assert e5["status"] == "PASS_EXACT_Q24_D12_MISSING_E5_SECTION_QQ"
assert route["status"] == "PASS_EXACT_Q24_D12_DEGREE14_AJ_MISSING_PARENT_ROUTE"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_MW_PROFILE_ENUMERATION"
assert marking["status"] == "PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003"

representatives = route["target_section_profile"]["representatives_transported_to_A11"]
assert len(representatives) == 1
trace_mw = vector(ZZ, representatives[0]["A11_MW_Abel_Jacobi"])
assert trace_mw == vector(ZZ, (13, -20, 12, -12, 0, 1))
assert representatives[0]["A11_degree"] == 13

shell = [
    vector(ZZ, row)
    for row in bridge["exact_identity_shell"]["MW_vectors_in_equation_order"]
]
assert len(shell) == 18
target = vector(ZZ, bridge["selected_bridge"]["mw"])
assert target == vector(ZZ, (1, 0, 0, 0, 0, 1))

pointed_minus = vector(ZZ, (0, -1, 0, 0, 0, 0))
assert any(vector(ZZ, row["mw"]) == pointed_minus for row in pointed["candidates"])
coefficients = {7: -4, 14: -2, 17: -2}
check = trace_mw + 2 * pointed_minus
for index, coefficient in coefficients.items():
    check += coefficient * shell[index]
assert check == target

equation_degrees = marking["equation_identity_shell_new_fibre_degrees"]
source_degrees = {index: int(equation_degrees[index]) for index in coefficients}
assert source_degrees == {7: 1, 14: 3, 17: 1}

payload = {
    "schema": "elkies-k3.h3-q24-a11-degree13-e5-bridge-route.v1",
    "status": "PASS_EXACT_A11_DEGREE13_E5_BRIDGE_ROUTE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "trace_carrier": {
        "name": "D12_E5",
        "D12_MW_Abel_Jacobi": [0, 0, 0, 0, 1],
        "A11_MW_Abel_Jacobi": [int(value) for value in trace_mw],
        "A11_degree": 13,
        "trace_space": "L(14O)",
        "trace_dimension": 14,
    },
    "bridge_word": {
        "formula": "M=AJ(E5)-4*S7-2*S14-2*S17+2*Qminus",
        "identity_shell_coefficients": {
            str(index): coefficient for index, coefficient in coefficients.items()
        },
        "identity_shell_source_degrees": {
            str(index): degree for index, degree in source_degrees.items()
        },
        "pointed_opposite_MW": [int(value) for value in pointed_minus],
        "pointed_opposite_coefficient": 2,
        "verified_MW_sum": [int(value) for value in check],
    },
    "required_traces": [
        {"curve": "E5", "degree": 13, "space": "L(14O)"},
        {"curve": "S14", "degree": 3, "space": "L(4O)"},
    ],
    "degree_one_inputs": ["S7", "S17"],
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact characteristic-zero section and lattice/degree certificate. "
        "The degree-13 and degree-3 Abel--Jacobi traces and the resulting exact A11 M coordinates remain to be constructed."
    ),
}
output = GENERATED / "elkies-k3-h3-q24-a11-degree13-e5-bridge-route.json"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11E5BRIDGE|carrier_degree=13|aux_degrees=1,3,1|trace_spaces=L14O,L4O|"
    f"target={','.join(map(str, target))}|status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{output.resolve()}", flush=True)
