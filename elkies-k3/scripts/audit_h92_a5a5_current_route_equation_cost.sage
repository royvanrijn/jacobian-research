#!/usr/bin/env sage -python
"""Identify the historical 2A5->3A3 edge in the explicit-zero equation shell."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
ZERO_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
O12_CERT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
A5_RANKING = GENERATED / "elkies-k3-h3-a5a5-marked-target-neighbor-ranking.json"
GATE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-explicit-curve-gate.json"
SCORES = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-current-route-equation-cost-audit.json"
INPUTS = (MANIFEST, FINGERPRINT, ZERO_FRAME, O12_CERT, A5_RANKING, GATE, SCORES, CROSSOVERS, ZERO_MISMATCH, MARKING)


manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
zero_frame = json.loads(ZERO_FRAME.read_text())
o12_certificate = json.loads(O12_CERT.read_text())
a5_ranking = json.loads(A5_RANKING.read_text())
gate = json.loads(GATE.read_text())
scores = json.loads(SCORES.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
marking = json.loads(MARKING.read_text())

historical_a11_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
historical_a11_to_a5 = matrix(ZZ, manifest["forward_steps"][2]["transition"])
historical_a5_in_equation = historical_a11_to_a5 * historical_a11_in_equation
route_fibre_historical = vector(ZZ, manifest["forward_steps"][3]["new_fibre_in_parent"])
route_fibre_equation = route_fibre_historical * historical_a5_in_equation

equation_to_explicit = matrix(ZZ, zero_frame["selected"]["equation_A11_to_explicit_zero_basis"])
explicit_to_equation = equation_to_explicit.inverse().change_ring(ZZ)
route_fibre_explicit = route_fibre_equation * explicit_to_equation
assert route_fibre_explicit[1] == 2
q = int(route_fibre_explicit[0] * route_fibre_explicit[1])

frame = matrix(ZZ, zero_frame["selected"]["frame"])
g_explicit = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
assert route_fibre_explicit * g_explicit * route_fibre_explicit == 0
# ``vector(ZZ, existing_vector)`` may return the same mutable Sage object.
# Copy through a Python list so the recorded pre-Weyl representative is not
# silently overwritten by the chamber-reduction loop below.
dominant = vector(ZZ, list(route_fibre_explicit))
for unused in range(10000):
    labels = list(dominant[2:] * frame[:, :10])
    negative = [index for index, value in enumerate(labels) if value < 0]
    if not negative:
        break
    index = negative[0]
    dominant[2 + index] -= labels[index]
else:
    raise RuntimeError("Weyl reduction did not terminate")

matches = [
    item for item in gate["survivors"]
    if vector(ZZ, item["source_neighbor_record"]["fiber"]) == dominant
]
assert not matches and q == 104

# Price the class directly in the explicit zero even though q104 was outside
# the exhaustive low-q shell.
root_rank = 10
root = frame[:root_rank, :root_rank]
root_lattice = IntegralLattice(root)
z = vector(ZZ, route_fibre_explicit[2 + root_rank:])
base = vector(ZZ, [0] * root_rank + list(z))
dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
iterator = root_lattice.enumerate_close_vectors(-dual)
minimum = None
profiles = []
for unused in range(100000):
    shift = vector(ZZ, next(iterator))
    lifted = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
    norm = QQ(lifted * frame * lifted)
    if minimum is None:
        minimum = norm
    elif norm > minimum:
        break
    pole = (norm - 4) / 2
    if pole in ZZ and pole >= 0:
        profiles.append((ZZ(pole), vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))))
profiles.sort(key=lambda item: (item[0], tuple(item[1])))
pole, section = profiles[0]
residual = route_fibre_explicit - vector(ZZ, [-1, 1] + [0] * 17) - section
vertical = vector(ZZ, residual[2:2 + root_rank])
root_edges = [(i, j) for i in range(root_rank) for j in range(i + 1, root_rank) if root[i, j] == -1]
magnitudes = [abs(value) for value in vertical]
layers = 0
previous = 0
for level in sorted(set(value for value in magnitudes if value)):
    active = set(index for index, value in enumerate(magnitudes) if value >= level)
    count = 0
    while active:
        count += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in root_edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    layers += (level - previous) * count
    previous = level
rr = int(2 + 2 * pole + layers)

explicit_to_equation_matrix = matrix(
    ZZ, zero_frame["selected"]["equation_A11_to_explicit_zero_basis"]
).inverse().change_ring(ZZ)
g_equation = explicit_to_equation_matrix * g_explicit * explicit_to_equation_matrix.transpose()
old_simple = [vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)]) for node in range(11)]
old_affine = vector(ZZ, [1, 0] + [1] * 11 + [0] * 6)
a0 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"])
named_curves = old_simple + [old_affine, a0, p24]
named_names = [f"old_A11_component_{node}" for node in range(11)] + ["old_A11_affine", "oldI9_A0", "close_P24"]
named_degrees = {name: int(curve * g_equation * route_fibre_equation) for name, curve in zip(named_names, named_curves)}

# The historical q104 class is not nef in the physical equation chamber.  Use
# the canonical raw representative of its root-Weyl orbit, prove that it
# reduces to the historical dominant vector in the abstract 2A5 chamber, then
# reduce it against the eleven physical old components, both current-I6
# affine components, and the three additional named effective curves.
canonical_raw = vector(ZZ, [
    52, 2, -3, -3, -5, -9, -11, 2, 2, 7, 4, 5, 0, 0, 1, 0, -2, 0, 1,
])
assert canonical_raw * g_explicit * canonical_raw == 0
assert canonical_raw[1] == 2 and canonical_raw[0] * canonical_raw[1] == 104
abstract_reduction = vector(ZZ, list(canonical_raw))
for unused in range(10000):
    labels = list(abstract_reduction[2:] * frame[:, :root_rank])
    negative = [index for index, value in enumerate(labels) if value < 0]
    if not negative:
        break
    index = negative[0]
    abstract_reduction[2 + index] -= labels[index]
else:
    raise RuntimeError("canonical q104 abstract Weyl reduction did not terminate")
assert abstract_reduction == dominant

components = {
    index: vector(ZZ, marking["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
old_fibre_explicit = vector(ZZ, [1, 0] + [0] * 17)
affines = [
    old_fibre_explicit - sum(
        (components[index] for index in chain), vector(ZZ, [0] * 19)
    )
    for chain in marking["physical_2A5"]["chains"]
]
extended_curves = (
    [components[index] for index in range(11)] + affines
    + [old_affine * explicit_to_equation, a0 * explicit_to_equation, p24 * explicit_to_equation]
)
extended_names = (
    [f"old_A11_component_{index}" for index in range(11)]
    + ["first_I6_affine_component", "second_I6_affine_component",
       "old_A11_affine", "oldI9_A0", "close_P24"]
)
assert all(curve * g_explicit * curve == -2 for curve in extended_curves)

def explicit_reflection(curve):
    action = identity_matrix(ZZ, 19) + (g_explicit * curve.column()) * matrix(ZZ, [list(curve)])
    assert action * g_explicit * action.transpose() == g_explicit
    assert action.det() == -1
    return action


repaired_route_explicit = vector(ZZ, list(canonical_raw))
physical_weyl = identity_matrix(ZZ, 19)
physical_reflections = []
for unused in range(10000):
    negative = [
        (index, ZZ(repaired_route_explicit * g_explicit * curve))
        for index, curve in enumerate(extended_curves)
        if repaired_route_explicit * g_explicit * curve < 0
    ]
    if not negative:
        break
    index, pairing = negative[0]
    step = explicit_reflection(extended_curves[index])
    repaired_route_explicit = repaired_route_explicit * step
    physical_weyl = physical_weyl * step
    physical_reflections.append({"curve": extended_names[index], "pairing": int(pairing)})
else:
    raise RuntimeError("physical q104 Weyl reduction did not terminate")

assert len(physical_reflections) == 61
assert list(repaired_route_explicit) == [
    5, 2, 0, -1, -1, -3, -3, 1, 0, -2, 0, -2, 0, 0, 1, 0, -2, 0, 1,
]
assert repaired_route_explicit * g_explicit * repaired_route_explicit == 0
assert repaired_route_explicit[1] == 2
assert repaired_route_explicit[0] * repaired_route_explicit[1] == 10
repaired_extended_degrees = {
    name: int(repaired_route_explicit * g_explicit * curve)
    for name, curve in zip(extended_names, extended_curves)
}
assert [repaired_extended_degrees[f"old_A11_component_{index}"] for index in range(11)] == [
    0, 1, 0, 0, 0, 1, 1, 0, 0, 3, 0,
]
assert [repaired_extended_degrees[name] for name in (
    "first_I6_affine_component", "second_I6_affine_component",
)] == [1, 0]
assert repaired_extended_degrees["old_A11_affine"] == 3
assert repaired_extended_degrees["oldI9_A0"] == 218
assert repaired_extended_degrees["close_P24"] == 2528
assert all(value >= 0 for value in repaired_extended_degrees.values())

repaired_residual = repaired_route_explicit - vector(ZZ, [-1, 1] + [0] * 17) - section
assert all(value == 0 for value in repaired_residual[2 + root_rank:])
repaired_vertical = vector(ZZ, repaired_residual[2:2 + root_rank])
repaired_magnitudes = [abs(value) for value in repaired_vertical]
repaired_layers = 0
previous = 0
for level in sorted(set(value for value in repaired_magnitudes if value)):
    active = set(index for index, value in enumerate(repaired_magnitudes) if value >= level)
    count = 0
    while active:
        count += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in root_edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    repaired_layers += (level - previous) * count
    previous = level
repaired_rr = int(2 + 2 * pole + repaired_layers)
assert repaired_layers == 3 and repaired_rr == 15
target_fibres = {item["target"]: vector(ZZ, item["target_fibre_in_state"]) for item in crossovers["records"] if item["state"] == "equation_A11"}
marked_degrees = {name: int(route_fibre_equation * g_equation * fibre) for name, fibre in target_fibres.items()}
deterministic_transition = matrix(ZZ, o12_certificate["transport"]["parent_to_child_basis"])
deterministic_inverse = deterministic_transition.inverse().change_ring(ZZ)
g_deterministic = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -matrix(ZZ, o12_certificate["child"]["frame"]))
route_deterministic = route_fibre_equation * deterministic_inverse
target_deterministic = target_fibres["pinned_R17"] * deterministic_inverse
deterministic_pinned_degree = int(route_deterministic * g_deterministic * target_deterministic)
assert deterministic_pinned_degree == marked_degrees["pinned_R17"]

payload = {
    "schema": "elkies-k3.h3-a5a5-current-route-equation-cost-audit.v1",
    "status": "PASS_EXACT_CURRENT_ROUTE_PHYSICAL_WEYL_REPAIR_WITHDRAWS_Q104_SCORE",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "historical_edge": {"q": 4, "orbit": 472, "parent": "2A5/MW7", "child": "3A3/MW8"},
    "explicit_zero_edge": {"q": q, "low_q_shell_match": None, "fibre_before_Weyl_reduction": [int(value) for value in route_fibre_explicit], "dominant_fibre": [int(value) for value in dominant]},
    "stored_nonphysical_equation_cost_profile": {"P_dot_O": int(pole), "expected_RR_ambient": rr, "horizontal_section": [int(value) for value in section], "vertical": [int(value) for value in vertical], "vertical_layers": int(layers), "named_explicit_curve_degrees": named_degrees, "marked_target_degrees": marked_degrees, "negative_physical_components": [name for name, value in named_degrees.items() if name.startswith("old_A11_component_") and value < 0], "former_operational_score": 13518, "withdrawn": True},
    "physical_weyl_repair": {"canonical_raw_fibre_in_explicit_zero_coordinates": [int(value) for value in canonical_raw], "abstract_dominant_fibre_in_explicit_zero_coordinates": [int(value) for value in abstract_reduction], "reflection_sequence": physical_reflections, "parent_weyl_isometry_in_explicit_zero_coordinates": [[int(value) for value in row] for row in physical_weyl.rows()], "repaired_fibre_in_explicit_zero_coordinates": [int(value) for value in repaired_route_explicit], "q": 10, "old_fibre_degree": 2, "primitive": gcd([abs(int(value)) for value in g_explicit * repaired_route_explicit]) == 1, "isotropic": True, "P_dot_O": int(pole), "horizontal_section": [int(value) for value in section], "vertical": [int(value) for value in repaired_vertical], "vertical_layers": int(repaired_layers), "expected_RR_ambient": repaired_rr, "named_explicit_curve_degrees": repaired_extended_degrees, "known_effective_curves_nef": True, "operational_score": None},
    "superseded_ranking_bug": {"old_recorded_pinned_degree": a5_ranking["current_route"]["marked_target_degrees"]["pinned_R17"], "old_recorded_fibre_in_deterministic_child": a5_ranking["current_route"]["fibre_in_equation_marked_parent"], "correct_fibre_in_deterministic_child": [int(value) for value in route_deterministic], "correct_invariant_pinned_degree": deterministic_pinned_degree, "cause": "the earlier current-route row is inconsistent with the direct equation-coordinate invariant replay and is superseded; candidate-shell rankings are unaffected"},
    "proof_boundary": "Exact marked physical-chamber rejection of the historical q104 representative and exact 61-reflection Weyl repair to q10, including the parent NS isometry and corrected horizontal/vertical RR profile. This does not compile the q10 equation, construct its child U-splitting, classify the q10 child, or assign a replacement operational score.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("A5A5CURRENT|historical=q4o472|stored=q{}|stored_RR={}|physical=q10|physical_RR={}|PO={}|reflections={}|pinned={}|status={}".format(
    q, rr, repaired_rr, pole, len(physical_reflections), marked_degrees["pinned_R17"], payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
