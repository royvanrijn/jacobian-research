#!/usr/bin/env sage -python
"""Repair the historical 2A5-to-3A3 class in the physical equation chamber.

status: ACTIVE_PROOF
claim: exact physical-nef q10 presentation of the canonical current 3A3 stage
outputs: artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-certificate.json

The historical q104 representative is dominant only in an abstract A5+A5
root basis.  Reduce its canonical pre-Weyl class against the two physical I6
cycles and the other exact equation curves.  The result is a degree-two q10
fibre with the same MW quotient, a 15-dimensional expected RR ambient space,
and a full unimodular landing on the stored current-3A3 Gram.  No Groebner
basis is used.
"""

import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
ZERO_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
OUTPUT = LOCAL / "q24-2a5-direct-physical-q10-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(item) for item in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
marking = json.loads(MARKING.read_text())
zero_frame = json.loads(ZERO_FRAME.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
current_3a3 = json.loads(CURRENT_3A3.read_text())
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert zero_frame["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert zero_mismatch["status"] == "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH"
assert current_3a3["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

frame = matrix(ZZ, zero_frame["selected"]["frame"])
gram = block_diagonal_matrix(U2, -frame)
equation_to_explicit = matrix(
    ZZ, zero_frame["selected"]["equation_A11_to_explicit_zero_basis"]
)
explicit_to_equation = equation_to_explicit.inverse().change_ring(ZZ)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)

# Replay the canonical historical class before any abstract chamber reduction.
historical_a11_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
historical_a5_in_equation = (
    matrix(ZZ, manifest["forward_steps"][2]["transition"])
    * historical_a11_in_equation
)
historical_exit = manifest["forward_steps"][3]
assert historical_exit["parent"] == "2A5/MW7"
assert historical_exit["child"] == "3A3/MW8"
raw_fibre = (
    vector(ZZ, historical_exit["new_fibre_in_parent"])
    * historical_a5_in_equation * explicit_to_equation
)
assert entries(raw_fibre) == [
    52, 2, -3, -3, -5, -9, -11, 2, 2, 7, 4, 5, 0, 0, 1, 0, -2, 0, 1,
]
assert raw_fibre * gram * raw_fibre == 0

# Exact effective curves in the component-9-zero equation marking.
components = {
    index: vector(ZZ, marking["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
chains = marking["physical_2A5"]["chains"]
affine = [
    old_fibre - sum(
        (components[index] for index in chain), vector(ZZ, [0] * 19)
    )
    for chain in chains
]
curves = [components[index] for index in range(11)] + affine + [
    vector(ZZ, marking["physical_2A5"]["child_coordinates"]["old_A11_affine"]),
]
labels = [f"old_A11_component_{index}" for index in range(11)] + [
    "first_I6_affine_component", "second_I6_affine_component", "old_A11_affine",
]
for name in ("oldI9_A0", "close_P24"):
    curve_equation = vector(
        ZZ, zero_mismatch["correct_selected_R3_transport"][name]["child_coordinates"]
    )
    curves.append(curve_equation * explicit_to_equation)
    labels.append(name)
assert all(curve * gram * curve == -2 for curve in curves)


def reflection(curve):
    action = identity_matrix(ZZ, 19) + (gram * curve.column()) * matrix(ZZ, [list(curve)])
    assert action * gram * action.transpose() == gram
    assert action.det() == -1
    return action


repaired = vector(ZZ, list(raw_fibre))
weyl = identity_matrix(ZZ, 19)
reflection_log = []
for unused in range(10000):
    pairings = [ZZ(repaired * gram * curve) for curve in curves]
    negative = next((index for index, value in enumerate(pairings) if value < 0), None)
    if negative is None:
        break
    value = pairings[negative]
    step = reflection(curves[negative])
    repaired = repaired * step
    weyl = weyl * step
    reflection_log.append({"curve": labels[negative], "pairing": int(value)})
else:
    raise ArithmeticError("physical Weyl reduction did not terminate")

assert len(reflection_log) == 61
assert entries(repaired) == [
    5, 2, 0, -1, -1, -3, -3, 1, 0, -2, 0, -2, 0, 0, 1, 0, -2, 0, 1,
]
assert repaired == raw_fibre * weyl
assert repaired * gram * repaired == 0
assert repaired * gram * old_fibre == 2
assert gcd([abs(ZZ(value)) for value in gram * repaired]) == 1
final_pairings = [ZZ(repaired * gram * curve) for curve in curves]
assert all(value >= 0 for value in final_pairings)
assert list(map(int, final_pairings)) == [
    0, 1, 0, 0, 0, 1, 1, 0, 0, 3, 0, 1, 0, 3, 218, 2528,
]

# Complete all-section and finite horizontal-wall nef gates for old degree two.
center = vector(QQ, repaired[2:]) / repaired[1]
closest = vector(ZZ, next(IntegralLattice(frame).enumerate_close_vectors(center)))
distance = (closest - center) * frame * (closest - center)
minimum_section = repaired[1] * (distance - 2) / 2
assert distance == 2 and minimum_section == 0
negative_horizontal_walls = []
degree = ZZ(repaired[1])
w = vector(ZZ, repaired[2:])
for old_degree in range(1, int(degree) + 1):
    m = ZZ(old_degree)
    cross = -degree * m * frame * w.column()
    augmented = block_matrix(ZZ, [
        [degree**2 * frame, cross],
        [cross.transpose(), matrix(ZZ, [[m**2 * (w * frame * w) + 1]])],
    ])
    result = pari(augmented).qfminim(2 * degree**2 - 1)
    for raw in matrix(ZZ, result[2]).transpose().rows():
        if abs(raw[-1]) != 1:
            continue
        value = raw if raw[-1] == 1 else -raw
        x = vector(ZZ, value[:-1])
        x_norm = ZZ(x * frame * x)
        if (x_norm - 2) % (2 * m):
            continue
        k = ZZ((x_norm - 2) // (2 * m))
        intersection = ZZ(
            (w * frame * w // (2 * degree)) * m + degree * k - w * frame * x
        )
        if intersection < 0:
            negative_horizontal_walls.append((int(m), int(intersection), entries(x)))
assert not negative_horizontal_walls

# Root-adapt to the physical nonidentity components and measure the RR profile.
effective_simple = []
effective_simple_names = []
for chain_index, chain in enumerate(chains):
    chain_curves = [components[index] for index in chain]
    for index, curve in zip(chain, chain_curves):
        if old_zero * gram * curve == 0:
            effective_simple.append(curve)
            effective_simple_names.append(f"old_A11_component_{index}")
    effective_simple.append(affine[chain_index])
    effective_simple_names.append(f"I6_affine_{chain_index}")
root_rows = matrix(ZZ, [-curve[2:] for curve in effective_simple])
smith, _, smith_right = root_rows.smith_form()
assert tuple(abs(smith[index, index]) for index in range(10)) == (1,) * 10
completion = smith_right.inverse()
adaptation = root_rows.stack(completion[10:])
assert abs(adaptation.det()) == 1
physical_frame = adaptation * frame * adaptation.transpose()
change = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation)
repaired_physical = repaired * change.inverse().change_ring(ZZ)
root = physical_frame[:10, :10]
base = vector(ZZ, [0] * 10 + list(repaired_physical[12:]))
dual = vector(QQ, base * physical_frame[:, :10]) * root.inverse()
iterator = IntegralLattice(root).enumerate_close_vectors(-dual)
minimum = None
choices = []
for shift in iterator:
    lifted = base + vector(ZZ, list(shift) + [0] * 7)
    norm = QQ(lifted * physical_frame * lifted)
    if minimum is None:
        minimum = norm
    elif norm > minimum:
        break
    pole = (norm - 4) / 2
    if pole in ZZ and pole >= 0:
        choices.append((ZZ(pole), lifted))
pole, lifted = min(choices, key=lambda item: (item[0], tuple(item[1])))
horizontal_physical = vector(ZZ, [pole + 1, 1] + list(lifted))
horizontal = horizontal_physical * change
vertical = repaired_physical - vector(ZZ, [-1, 1] + [0] * 17) - horizontal_physical
assert pole == 5 and all(value == 0 for value in vertical[12:])
assert entries(vertical[2:12]) == [1, 1, 1, 1, 1, 2, 2, 2, 1, 1]
assert entries(horizontal) == [
    6, 1, 1, 0, -1, -2, -3, 0, 0, -2, 0, -2, 0, 0, 1, 0, -2, 0, 1,
]
vertical_layers = 3  # one connected layer on the first A5, two on the second
expected_rr_ambient = int(2 + 2 * pole + vertical_layers)
assert expected_rr_ambient == 15

# Apply the parent Weyl isometry to the canonical historical child basis.
historical_parent_to_child = matrix(ZZ, historical_exit["transition"])
stored_parent_to_child = (
    historical_parent_to_child * historical_a5_in_equation * explicit_to_equation
)
assert vector(ZZ, stored_parent_to_child.row(0)) == raw_fibre
parent_to_child = stored_parent_to_child * weyl
child_to_parent = parent_to_child.inverse().change_ring(ZZ)
assert vector(ZZ, parent_to_child.row(0)) == repaired
assert abs(parent_to_child.det()) == 1
canonical_frame_path = ROOT / current_3a3["frame_output"]
canonical_frame = load_matrix(canonical_frame_path)
canonical_gram = block_diagonal_matrix(U2, -canonical_frame)
assert parent_to_child * gram * parent_to_child.transpose() == canonical_gram
pinned_basis = matrix(ZZ, current_3a3["pinned_R17_basis_in_source"])
assert abs(pinned_basis.det()) == 1

inputs = (
    MANIFEST, FINGERPRINT, MARKING, ZERO_FRAME, ZERO_MISMATCH,
    CURRENT_3A3, canonical_frame_path,
)
payload = {
    "schema": "elkies-k3.q24-2a5-direct-physical-q10.v1",
    "status": "PASS_EXACT_PHYSICAL_NEF_Q10_CURRENT_3A3_PRESENTATION",
    "historical_edge": {"q": 4, "orbit": 472, "equation_zero_q": 104},
    "physical_weyl_repair": {
        "raw_fibre": entries(raw_fibre),
        "reflection_count": len(reflection_log),
        "reflection_log": reflection_log,
        "parent_weyl_isometry": rows(weyl),
        "repaired_fibre": entries(repaired),
        "q": 10,
        "old_fibre_degree": 2,
        "known_effective_curve_pairings": dict(zip(labels, map(int, final_pairings))),
        "all_section_minimum_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "finite_negative_horizontal_walls": [],
    },
    "RR_profile": {
        "P_dot_O": int(pole),
        "horizontal_section": entries(horizontal),
        "physical_simple_roots": effective_simple_names,
        "vertical_coefficients_in_physical_root_adapted_frame": entries(vertical[2:12]),
        "vertical_connected_layers": vertical_layers,
        "expected_RR_ambient": expected_rr_ambient,
    },
    "landing": {
        "root_data": current_3a3["root_data"],
        "MW_rank": 8,
        "parent_to_current_3A3_basis": rows(parent_to_child),
        "current_3A3_to_parent_basis": rows(child_to_parent),
        "forward_determinant": int(parent_to_child.det()),
        "inverse_determinant": int(child_to_parent.det()),
        "canonical_Gram_exact": True,
        "pinned_R17_basis_determinant": int(pinned_basis.det()),
    },
    "proof_boundary": (
        "Exact physical-component, all-section, and finite horizontal-wall nef gates; "
        "exact RR dimension prediction; and a bidirectional unimodular landing on the "
        "canonical current-3A3 Gram. The characteristic-zero H0 basis, quartic, Jacobian, "
        "fibres, and equation marking remain to be compiled."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5DIRECTPHYSICAL|stored_q=104|physical_q=10|degree=2|PO=5|RR=15|"
    "reflections={}|landing_det={}|status={}|output={}".format(
        len(reflection_log), parent_to_child.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
