#!/usr/bin/env sage -python
"""Apply cheap exact marked-nef gates to the A5+A5 q6/q8 shells."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, default=GENERATED / "elkies-k3-h3-a11-o12-q6q8-degree2-all.json")
parser.add_argument("--transport", type=Path, default=GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json")
parser.add_argument("--selected-explicit-zero", action="store_true")
parser.add_argument("--output", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-q6q8-explicit-curve-gate.json")
args = parser.parse_args()
NEIGHBORS = args.neighbors.resolve()
O12_CERT = args.transport.resolve()
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
D12_FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6_SHELL = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
OUTPUT = args.output.resolve()
INPUTS = (NEIGHBORS, O12_CERT, CROSSOVERS, ZERO_MISMATCH, D12_FRAME, Q6_SHELL, IDENTITY, MATCHING)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


neighbors = json.loads(NEIGHBORS.read_text())
o12 = json.loads(O12_CERT.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
q6_shell = json.loads(Q6_SHELL.read_text())
identity = json.loads(IDENTITY.read_text())
matching = json.loads(MATCHING.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

a11_shell = json.loads((LOCAL / "q24-a11-orbit64-q8-all.json").read_text())
a11 = load_matrix(ROOT / a11_shell["frame"])
g_a11 = block_diagonal_matrix(U2, -a11)
parent = load_matrix(ROOT / neighbors["frame"])
g_parent = block_diagonal_matrix(U2, -parent)
a11_to_parent = matrix(
    ZZ,
    o12["selected"]["equation_A11_to_explicit_zero_basis"]
    if args.selected_explicit_zero
    else o12["transport"]["parent_to_child_basis"],
)
assert a11_to_parent * g_a11 * a11_to_parent.transpose() == g_parent

old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
old_simple = [
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(11)
]
# The equation A11 frame is the physical A11 chain, so its affine root has
# coefficient one on every simple component.
old_affine = old_fibre + vector(ZZ, [0, 0] + [1] * 11 + [0] * 6)
a0 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"])
# Reconstruct the 18 exact identity-shell sections and put them in equation order.
d12 = load_matrix(D12_FRAME)
selected_q6 = next(item for item in q6_shell["neighbors"] if int(item["orbit_index"]) == 64)
d12_to_a11 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, selected_q6["child_root_adapted_basis"])
) * matrix(ZZ, selected_q6["neighbor_basis"])
d12_inverse = d12_to_a11.inverse().change_ring(ZZ)
d12_root = d12[:12, :12]
d12_coupling = d12[:12, 12:]
shell = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * d12_coupling.transpose()) * d12_root.inverse()
    section = vector(ZZ, [1, 1] + list(map(ZZ, root_coefficients)) + list(z))
    shell.append(section * d12_inverse)
mapping = matching["matching"]["mappings_abstract_to_equation"][7]
reordered = [None] * 18
for abstract_index, equation_index in enumerate(mapping):
    reordered[equation_index] = shell[abstract_index]
shell = reordered

curves = shell + [old_zero] + old_simple + [old_affine, a0, p24]
names = [f"shell_S{i}" for i in range(18)] + ["old_A11_zero"] + [f"old_A11_component_{i}" for i in range(11)] + ["old_A11_affine", "oldI9_A0", "close_P24"]
assert all(curve * g_a11 * curve == -2 for curve in curves)

target_fibres = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}
counts = {4: 0, 6: 0, 8: 0}
survivors = []
for raw in neighbors["neighbors"]:
    q = int(raw["q"])
    fibre_parent = vector(ZZ, raw["fiber"])
    fibre_a11 = fibre_parent * a11_to_parent
    degrees = [int(curve * g_a11 * fibre_a11) for curve in curves]
    negative = [names[index] for index, value in enumerate(degrees) if value < 0]
    if negative:
        continue
    counts[q] = counts.get(q, 0) + 1
    marked = {
        target: int(fibre_a11 * g_a11 * target_fibre)
        for target, target_fibre in target_fibres.items()
    }
    assert all(value >= 0 for value in marked.values())
    survivors.append({
        "candidate_id": {"q": q, "old_fibre_degree": int(raw["old_fiber_degree"]), "orbit_index": int(raw["orbit_index"])},
        "child": {"ade": raw["child_ade"], "mw_rank": int(raw["child_mw_rank"]), "root_data": raw["child_root_data"]},
        "explicit_curve_degrees": dict(zip(names, degrees)),
        "explicit_degree_zero_curves": [names[index] for index, value in enumerate(degrees) if value == 0],
        "explicit_degree_one_curves": [names[index] for index, value in enumerate(degrees) if value == 1],
        "marked_target_degrees": marked,
        "coordinate_growth_max": int(max(abs(value) for value in fibre_parent)),
        "source_neighbor_record": raw,
    })

survivors.sort(key=lambda item: (
    item["marked_target_degrees"]["pinned_R17"],
    -len(item["explicit_degree_one_curves"]),
    item["coordinate_growth_max"],
))
payload = {
    "schema": "elkies-k3.h3-a5a5-q6q8-explicit-curve-gate.v1",
    "status": "PASS_EXACT_A5A5_Q6Q8_EXPLICIT_CURVE_GATE",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "input_candidate_count": len(neighbors["neighbors"]),
    "survivor_counts_by_q": {str(key): value for key, value in counts.items() if value},
    "survivor_count": len(survivors),
    "survivors": survivors,
    "proof_boundary": "Exact gate against the physical A11 fibre, its zero, A0 and P24. Identity-shell sections and full section walls remain to be checked.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if survivors:
    best = survivors[0]
    print("A5A5GATE|survivors={}|q4={}|q6={}|q8={}|best=q{}o{}|pinned={}|deg0={}|deg1={}|child={}/MW{}|status={}".format(
        len(survivors), counts[4], counts[6], counts[8], best["candidate_id"]["q"], best["candidate_id"]["orbit_index"],
        best["marked_target_degrees"]["pinned_R17"], len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]),
        best["child"]["ade"], best["child"]["mw_rank"], payload["status"]), flush=True)
else:
    print(f"A5A5GATE|survivors=0|q4={counts[4]}|q6={counts[6]}|q8={counts[8]}|status={payload['status']}", flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
