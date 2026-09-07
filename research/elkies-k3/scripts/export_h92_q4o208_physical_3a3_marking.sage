#!/usr/bin/env sage -python
"""Export the equation-attached physical 3A3 chamber after q4/orbit208."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
EQUATION = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
PARENT = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
CURRENT = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
FRAME = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


route = json.loads(ROUTE.read_text())
equation = json.loads(EQUATION.read_text())
parent = json.loads(PARENT.read_text())
current = json.loads(CURRENT.read_text())
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert equation["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert parent["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

raw_frame = matrix(ZZ, route["selection"]["child_frame"])
raw_gram = block_diagonal_matrix(U2, -raw_frame)
current_frame = load_matrix(ROOT / current["frame_output"])
current_gram = block_diagonal_matrix(U2, -current_frame)

# Start each A3 chain immediately after its identity component in the physical
# I4 cycle.  This fixes the effective simple-root orientation geometrically.
simple_roots = []
missing_components = {}
for fibre_name, fibre_data in equation["physical_fibres"].items():
    cycle = [vector(ZZ, value) for value in fibre_data["components_in_cycle_order"]]
    identity = int(fibre_data["identity_component_index"])
    # Root-adapted scripts use ``-e_i`` for effective component curves, so the
    # basis rows themselves are the negatives of the physical components.
    chain = [-cycle[(identity + offset) % 4] for offset in (1, 2, 3)]
    simple_roots.extend(chain)
    if fibre_data["missing_component_appended"]:
        missing_components[f"{fibre_name}_missing_component"] = cycle[-1]

e = vector(ZZ, [1, 0] + [0] * 17)
f = vector(ZZ, [0, 1] + [0] * 17)
partial = matrix(ZZ, [list(e), list(f)] + [list(root) for root in simple_roots])
assert partial.rank() == 11
assert partial[:2] * raw_gram * partial[:2].transpose() == U2
root_gram = -(partial[2:] * raw_gram * partial[2:].transpose())
assert root_gram.det() == 64

# Smith completion preserves the displayed U and physical simple roots.
smith, left, right = partial.smith_form()
assert smith[:, :11] == matrix.identity(ZZ, 11)
assert not any(smith[:, 11:].list())
completion = right.inverse().change_ring(ZZ)
change = block_diagonal_matrix(left.inverse().change_ring(ZZ), matrix.identity(ZZ, 8))
basis_in_raw = change * completion
assert basis_in_raw[:11] == partial
assert abs(basis_in_raw.det()) == 1
raw_in_basis = basis_in_raw.inverse().change_ring(ZZ)
physical_gram = basis_in_raw * raw_gram * basis_in_raw.transpose()
assert physical_gram[:2, :2] == U2
assert not any(physical_gram[:2, 2:].list())
physical_frame = -physical_gram[2:, 2:]
assert physical_frame[:9, :9] == root_gram
FRAME.write_text(
    "# exact q4/orbit208 equation-physical 3A3 frame\n"
    + "\n".join(" ".join(map(str, row)) for row in physical_frame.rows()) + "\n"
)

current_in_raw = matrix(ZZ, route["transport"]["current_3A3_basis_in_effective_zero_child"])
targets_current = {
    name: vector(ZZ, value)
    for name, value in current["target_fibres_in_root_adapted_hub"].items()
}
targets_physical = {
    name: entries(target * current_in_raw * raw_in_basis)
    for name, target in targets_current.items()
}

parent_curves = {
    name: vector(ZZ, value)
    for name, value in parent["equation_explicit_curves_in_child"].items()
}
# Despite the historical key name, this stored inverse maps parent coordinates
# to the C5-zero child coordinates.
parent_to_raw = matrix(ZZ, route["transport"]["effective_zero_child_to_parent_basis"])
explicit_physical = {
    name: entries(curve * parent_to_raw * raw_in_basis)
    for name, curve in parent_curves.items()
}
explicit_physical.update({
    name: entries(curve * raw_in_basis)
    for name, curve in missing_components.items()
})
assert all(
    vector(ZZ, value) * physical_gram * vector(ZZ, value) == -2
    for value in explicit_physical.values()
)

raw_in_parent = matrix(ZZ, route["transport"]["parent_to_effective_zero_child_basis"])
parent_in_equation = matrix(ZZ, parent["physical_component_chamber_basis_in_equation_A11"])
physical_in_equation = basis_in_raw * raw_in_parent * parent_in_equation
assert abs(physical_in_equation.det()) == 1

inputs = (ROUTE, EQUATION, PARENT, CURRENT)
payload = {
    "schema": "elkies-k3.h3-q4o208-physical-3a3-marking.v1",
    "status": "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING",
    "hub": "q4o208_equation_physical_3A3_C5_zero",
    "zero": "old_A11_component_5",
    "root_data": [9, 36, 64],
    "frame_output": str(FRAME.relative_to(ROOT)),
    "frame_sha256": sha256(FRAME),
    "physical_simple_roots_in_C5_child": [entries(root) for root in simple_roots],
    "physical_3A3_basis_in_C5_child": rows(basis_in_raw),
    "C5_child_basis_in_physical_3A3": rows(raw_in_basis),
    "equation_A11_to_root_adapted_hub_basis": rows(physical_in_equation),
    "target_fibres_in_root_adapted_hub": targets_physical,
    "equation_explicit_curves_in_child": explicit_physical,
    "prefix_operational_score": route["compiler_profile"]["operational_equation_cost_score"],
    "proof_boundary": (
        "Exact equation-attached physical 3A3 chamber.  The nine simple roots are the "
        "nonidentity components of the three exact physical I4 cycles, the completion is "
        "integral unimodular, and all inherited explicit curves, missing I4 components, "
        "marked targets and equation-A11 transport are carried exactly."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"Q4O208PHYS3A3|roots=9|det=64|status={payload['status']}|output={OUTPUT}")
