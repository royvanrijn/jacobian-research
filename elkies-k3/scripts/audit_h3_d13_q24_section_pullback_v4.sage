#!/usr/bin/env sage -python
"""
Exact pullback of the first rank-growing D13 q24 section to the actual
component-nef q8 equation presentation.

This version deliberately avoids qfisom and avoids comparing independently
constructed Weyl-chamber endpoints.

There are three coordinate layers in the stored q8 classifier:

  source H3 NS
      ^
      | q6.neighbor_basis_in_source_ns
  raw q6-child NS
      ^
      | root_mw_basis, then simple_change
  q6 simple-frame NS

The dominant D13 hit was constructed in the q6 SIMPLE frame.
The stored nef representative was constructed directly in the RAW q6 frame.
The classifier itself supplies the exact four simple-root reflections taking
the nef q8 fibre to the dominant D13 hit.  We invert those reflections on the
pinned q24 section, return to raw q6/source-H3 coordinates, and only then replay
the 102 actual E6/E8 component reflections used by the repaired equation-level
q8 construction.

Thus the transport chain is:

 pinned dominant D13 q24 section
   -> q6 simple frame
   -> inverse stored nef-to-dominant reflections
   -> raw q6 frame
   -> source H3 NS
   -> 102 physical E6/E8 reflections
   -> actual component-nef D13 equation presentation.

The script checks every matrix identity and every fibre/section identity before
reporting pullback degrees.  It is a lattice transport audit; it does not yet
construct the rational equation of the q24 section.

Run:
  sage -python ~/Downloads/audit_h3_d13_q24_section_pullback_v4.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector, xgcd
)

REFLECTIONS_Q6 = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

Q24_WITNESS = vector(ZZ, (
    0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
))


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def reflect_row(row, gram, root):
    row = vector(ZZ, row)
    root = vector(ZZ, root)
    assert root * gram * root == -2
    return row + (row * gram * root) * root


def bezout_mate(ns, fibre):
    current = ZZ(0)
    entries = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        entries = [left * item for item in entries]
        entries[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        entries = [-item for item in entries]
    mate = vector(ZZ, entries)
    mate -= (mate * ns * mate // 2) * fibre
    assert mate * ns * mate == 0
    assert mate * ns * fibre == 1
    return mate


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED_D13 = (
    ROOT / "elkies-k3/data/fibrations/"
           "h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)
ORBITS = (
    ROOT / "artifacts/generated-results/"
           "elkies-k3-h3-q6-q8-orbits.json"
)
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts/local/elkies-k3/"
                "d13-q24-section-pullback-v4.json"
)

for path in (FRAME, PINNED_D13, ORBITS):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

source_frame = load_gram(FRAME)
pinned_d13 = load_gram(PINNED_D13)
assert source_frame.nrows() == pinned_d13.nrows() == 17
assert source_frame.det() == pinned_d13.det() == 948

U = matrix(ZZ, ((0, 1), (1, 0)))
source_ns = block_diagonal_matrix(U, -source_frame)
source_F = vector(ZZ, [1, 0] + [0] * 17)
source_O = vector(ZZ, [-1, 1] + [0] * 17)

data = json.loads(ORBITS.read_text())
assert data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"

# ===========================================================================
# 1. Reconstruct and certify the stored q6 coordinate tower.
# ===========================================================================

B6 = matrix(ZZ, data["q6"]["neighbor_basis_in_source_ns"])
assert B6.nrows() == B6.ncols() == 19
assert abs(B6.det()) == 1

G_q6_raw_ns = B6 * source_ns * B6.transpose()
assert G_q6_raw_ns[:2, :2] == U
assert G_q6_raw_ns[:2, 2:] == matrix(ZZ, 2, 17)
assert G_q6_raw_ns[2:, :2] == matrix(ZZ, 17, 2)
q6_raw_frame = -G_q6_raw_ns[2:, 2:]
assert q6_raw_frame.det() == 948

root_mw_basis = matrix(ZZ, data["q6"]["root_mw_basis_in_child"])
assert root_mw_basis.nrows() == root_mw_basis.ncols() == 17
assert abs(root_mw_basis.det()) == 1

root_mw_frame = (
    root_mw_basis * q6_raw_frame * root_mw_basis.transpose()
)
assert root_mw_frame == matrix(ZZ, data["q6"]["root_adapted_gram"])

simple_root_change = matrix(
    ZZ, data["q8"]["simple_root_change_in_root_block"]
)
assert simple_root_change.nrows() == simple_root_change.ncols() == 14
assert abs(simple_root_change.det()) == 1

simple_change = block_diagonal_matrix(
    simple_root_change, identity_matrix(ZZ, 3)
)
simple_to_q6_raw = simple_change * root_mw_basis
assert abs(simple_to_q6_raw.det()) == 1

q6_simple_frame = matrix(ZZ, data["q8"]["simple_frame_gram"])
assert (
    simple_to_q6_raw
    * q6_raw_frame
    * simple_to_q6_raw.transpose()
    == q6_simple_frame
)

Bsimple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), simple_to_q6_raw)
    * B6
)
assert abs(Bsimple.det()) == 1
q6_simple_ns = block_diagonal_matrix(U, -q6_simple_frame)
assert Bsimple * source_ns * Bsimple.transpose() == q6_simple_ns

print(
    "D13Q24COORD|q6_raw=PASS|root_mw=PASS|simple=PASS|"
    "simple_to_source=PASS",
    flush=True,
)

# ===========================================================================
# 2. Pin the dominant D13 basis used by the rootless path.
# ===========================================================================

# Select by exact pinned Gram, not by list order alone.
dominant_matches = [
    hit for hit in data["q8"]["d13_mw4_hits"]
    if matrix(ZZ, hit["d13_root_adapted_gram"]) == pinned_d13
]
assert len(dominant_matches) == 1
dominant_hit = dominant_matches[0]

B8_dom = matrix(ZZ, dominant_hit["neighbor_basis_in_q6_ns"])
raw_d13_dom = matrix(ZZ, dominant_hit["child_frame"])
assert B8_dom.nrows() == B8_dom.ncols() == 19
assert abs(B8_dom.det()) == 1
assert (
    B8_dom * q6_simple_ns * B8_dom.transpose()
    == block_diagonal_matrix(U, -raw_d13_dom)
)

A13_dom = matrix(ZZ, dominant_hit["d13_root_adapted_basis_in_child"])
assert A13_dom.nrows() == A13_dom.ncols() == 17
assert abs(A13_dom.det()) == 1
assert A13_dom * raw_d13_dom * A13_dom.transpose() == pinned_d13

Bdom_simple = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), A13_dom)
    * B8_dom
)
assert (
    Bdom_simple * q6_simple_ns * Bdom_simple.transpose()
    == block_diagonal_matrix(U, -pinned_d13)
)

F_dom_simple = vector(ZZ, Bdom_simple.row(0))
assert F_dom_simple == vector(
    ZZ, [4, 2] + dominant_hit["witness_simple_frame"]
)

print(
    "D13Q24COORD|dominant_neighbor=PASS|"
    "dominant_d13_adapt=PASS|pinned_frame=PASS|"
    f"dominant_mw={','.join(map(str, dominant_hit['mw_projection']))}",
    flush=True,
)

# ===========================================================================
# 3. The exact q24 section in the pinned D13 frame.
# ===========================================================================

G_d13 = block_diagonal_matrix(U, -pinned_d13)
F_pinned = vector(ZZ, [1, 0] + [0] * 17)
O_pinned = vector(ZZ, [-1, 1] + [0] * 17)
P_pinned = vector(ZZ, [23, 1] + list(Q24_WITNESS))
D_pinned = vector(ZZ, [12, 2] + list(Q24_WITNESS))

assert Q24_WITNESS * pinned_d13 * Q24_WITNESS == 48
assert P_pinned * G_d13 * P_pinned == -2
assert P_pinned * G_d13 * F_pinned == 1
assert P_pinned * G_d13 * O_pinned == 22
assert D_pinned * G_d13 * D_pinned == 0
assert D_pinned * G_d13 * F_pinned == 2
assert D_pinned == O_pinned + P_pinned - 10 * F_pinned

root = pinned_d13[:13, :13]
coupling = pinned_d13[:13, 13:]
tail = pinned_d13[13:, 13:]
height = tail - coupling.transpose() * root.inverse() * coupling
z = vector(ZZ, Q24_WITNESS[13:])
assert z == vector(ZZ, (0, -1, 1, 1))
hP = z * height * z
assert hP == 47

# Exact D13 local correction from the pinned section lift.
root_pairing = (
    vector(QQ, Q24_WITNESS)
    * pinned_d13[:, :13]
)
dual = root_pairing * root.inverse()
local_correction = dual * root * dual
assert local_correction == 1
assert (hP + local_correction - 4) / 2 == 22

print(
    "D13Q24SECTION|height=47|D13_correction=1|PdotO=22|"
    "decomposition=D=O+P-10F|status=PASS",
    flush=True,
)

# Map q24 classes into the q6 SIMPLE ambient of the dominant q8 fibration.
P_dom_simple = P_pinned * Bdom_simple
O_dom_simple = O_pinned * Bdom_simple
D_dom_simple = D_pinned * Bdom_simple
assert F_dom_simple == F_pinned * Bdom_simple

assert P_dom_simple * q6_simple_ns * P_dom_simple == -2
assert P_dom_simple * q6_simple_ns * F_dom_simple == 1
assert D_dom_simple == O_dom_simple + P_dom_simple - 10 * F_dom_simple

# ===========================================================================
# 4. Invert the classifier's exact nef -> dominant simple-root reflections.
# ===========================================================================

nef = data["q8"]["nef_representative"]
nef_to_dominant = [
    (int(index), int(label_pairing))
    for index, label_pairing in nef["to_dominant_reflections"]
]
assert nef_to_dominant, "classifier did not store the nef-to-dominant bridge"

# Simple q6 roots are the first 14 positive-frame coordinates; in NS they
# have square -2 because q6_simple_ns = U + (-simple_frame).
simple_ns_roots = tuple(
    vector(ZZ, [0, 0] + [ZZ(i == node) for i in range(17)])
    for node in range(14)
)
for root_class in simple_ns_roots:
    assert root_class * q6_simple_ns * root_class == -2

# First verify the FORWARD bridge on the stored nef fibre.
nef_witness_raw = vector(ZZ, nef["witness_q6_child"])
assert nef_witness_raw * q6_raw_frame * nef_witness_raw == 16
nef_witness_simple = (
    nef_witness_raw
    * root_mw_basis.inverse()
    * simple_change.inverse()
)
F_nef_simple = vector(ZZ, [4, 2] + list(nef_witness_simple))

forward = vector(ZZ, F_nef_simple)
for index1, stored_label in nef_to_dominant:
    node = index1 - 1
    root_class = simple_ns_roots[node]

    # The classifier stored the positive-frame Dynkin label
    # lambda = w * simple_frame[:,node].
    label = vector(ZZ, forward[2:]) * q6_simple_frame[:, node]
    assert label == stored_label

    # NS intersection is the NEGATIVE of that label.
    intersection = forward * q6_simple_ns * root_class
    assert intersection == -stored_label
    forward = reflect_row(forward, q6_simple_ns, root_class)

assert forward == F_dom_simple

# Now invert the SAME isometry on the q24 section, zero and divisor.
P_nef_simple = vector(ZZ, P_dom_simple)
O_nef_simple = vector(ZZ, O_dom_simple)
D_nef_simple = vector(ZZ, D_dom_simple)
F_back = vector(ZZ, F_dom_simple)

for index1, stored_label in reversed(nef_to_dominant):
    root_class = simple_ns_roots[index1 - 1]
    # Reflections are involutions.
    P_nef_simple = reflect_row(P_nef_simple, q6_simple_ns, root_class)
    O_nef_simple = reflect_row(O_nef_simple, q6_simple_ns, root_class)
    D_nef_simple = reflect_row(D_nef_simple, q6_simple_ns, root_class)
    F_back = reflect_row(F_back, q6_simple_ns, root_class)

assert F_back == F_nef_simple
assert P_nef_simple * q6_simple_ns * P_nef_simple == -2
assert P_nef_simple * q6_simple_ns * F_nef_simple == 1
assert D_nef_simple == O_nef_simple + P_nef_simple - 10 * F_nef_simple

# Convert simple q6 coordinates -> raw q6 -> source H3.
P_nef_source = P_nef_simple * Bsimple
O_nef_source = O_nef_simple * Bsimple
D_nef_source = D_nef_simple * Bsimple
F_nef_source = F_nef_simple * Bsimple

assert F_nef_source == vector(ZZ, nef["fiber_source_h3_ns"])
assert P_nef_source * source_ns * P_nef_source == -2
assert P_nef_source * source_ns * F_nef_source == 1
assert D_nef_source == O_nef_source + P_nef_source - 10 * F_nef_source

print(
    "D13Q24BRIDGE|nef_to_dominant_steps={}|forward=PASS|inverse=PASS|"
    "nef_source_match=PASS".format(len(nef_to_dominant)),
    flush=True,
)

# ===========================================================================
# 5. Replay the 102 ACTUAL E6/E8 component reflections used by the equation.
# ===========================================================================

# q6 fibre and transported old zero in source-H3 coordinates, exactly as in
# derive_h92_q6_child_q8_physical_root_target.sage.
q6_F_source = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])

old_simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)

def reflect_old_sequence(value, nodes):
    result = vector(ZZ, value)
    for node1 in nodes:
        root_class = old_simple[node1 - 1]
        result = reflect_row(result, source_ns, root_class)
    return result

q6_O_source = reflect_old_sequence(
    source_O, tuple(reversed(REFLECTIONS_Q6))
)
assert q6_F_source * source_ns * q6_F_source == 0
assert q6_O_source * source_ns * q6_O_source == -2
assert q6_O_source * source_ns * q6_F_source == 1

# Reconstruct the exact physical E6+E8 simple roots from the same complement
# and the same pinned qfminim indices/constants as the physical-target script.
q6_orth = matrix(
    ZZ,
    [
        list(q6_F_source * source_ns),
        list((q6_O_source + q6_F_source) * source_ns),
    ],
).right_kernel_matrix()
physical_child = -(q6_orth * source_ns * q6_orth.transpose())
qf_basis = matrix(
    ZZ, pari(physical_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert qf_basis.rank() == 14

E6_QF_INDICES = (0, 1, 2, 3, 12, 13)
E6_SIMPLE_IN_QF = (
    (-1, -1, -1, -1, 0, 0),
    (0, 0, 0, 1, 0, 0),
    (0, 0, 1, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (1, 0, 0, 0, 0, 1),
    (2, 1, 0, 0, -1, 1),
)
e6_qf = matrix(ZZ, [list(qf_basis[i]) for i in E6_QF_INDICES])
e6_roots = matrix(
    ZZ, [list(vector(ZZ, row) * e6_qf) for row in E6_SIMPLE_IN_QF]
) * q6_orth

physical_roots = qf_basis * q6_orth
e8_roots = physical_roots[4:12, :]

E6_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])

assert -e6_roots * source_ns * e6_roots.transpose() == E6_CARTAN
assert -e8_roots * source_ns * e8_roots.transpose() == E8_CARTAN
assert e6_roots * source_ns * e8_roots.transpose() == matrix(ZZ, 6, 8)

actual_roots = tuple(e6_roots.rows()) + tuple(e8_roots.rows())

F_component = vector(ZZ, F_nef_source)
P_component = vector(ZZ, P_nef_source)
O_component = vector(ZZ, O_nef_source)
D_component = vector(ZZ, D_nef_source)
physical_reflections = []

for unused in range(500):
    pairings = [
        int(F_component * source_ns * root_class)
        for root_class in actual_roots
    ]
    negative = [i for i, value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    root_class = actual_roots[i]
    fibre_pairing = pairings[i]

    F_component = reflect_row(F_component, source_ns, root_class)
    P_component = reflect_row(P_component, source_ns, root_class)
    O_component = reflect_row(O_component, source_ns, root_class)
    D_component = reflect_row(D_component, source_ns, root_class)
    physical_reflections.append((i, fibre_pairing))
else:
    raise RuntimeError("physical q8 component reduction did not terminate")

assert len(physical_reflections) == 102
assert F_component * source_ns * F_component == 0
assert F_component * source_ns * q6_F_source == 2
assert all(
    F_component * source_ns * root_class >= 0
    for root_class in actual_roots
)

assert P_component * source_ns * P_component == -2
assert P_component * source_ns * F_component == 1
assert D_component * source_ns * D_component == 0
assert D_component * source_ns * F_component == 2
assert D_component == O_component + P_component - 10 * F_component
assert O_component * source_ns * O_component == -2
assert O_component * source_ns * F_component == 1

print(
    "D13Q24PHYSICAL|reflections=102|fibre=PASS|section=PASS|"
    "zero=PASS|divisor=PASS",
    flush=True,
)

# ===========================================================================
# 6. Pullback degrees in earlier fibrations.
# ===========================================================================

def invariants_against(label, value, F, O):
    return {
        f"{label}_degree": int(value * source_ns * F),
        f"{label}_zero_intersection": int(value * source_ns * O),
    }

P_q6_degree = int(P_component * source_ns * q6_F_source)
P_q6_O = int(P_component * source_ns * q6_O_source)
P_h3_degree = int(P_component * source_ns * source_F)
P_h3_O = int(P_component * source_ns * source_O)

O8_q6_degree = int(O_component * source_ns * q6_F_source)
O8_h3_degree = int(O_component * source_ns * source_F)

D_q6_degree = int(D_component * source_ns * q6_F_source)
D_q6_O = int(D_component * source_ns * q6_O_source)
D_h3_degree = int(D_component * source_ns * source_F)
D_h3_O = int(D_component * source_ns * source_O)

# Sanity checks from D=O+P-10F8 and F8.Fq6=2.
assert D_q6_degree == O8_q6_degree + P_q6_degree - 20

print(
    "D13Q24PULLBACK|"
    f"P_q6_degree={P_q6_degree}|P_q6_O={P_q6_O}|"
    f"P_h3_degree={P_h3_degree}|P_h3_O={P_h3_O}|"
    f"O8_q6_degree={O8_q6_degree}|O8_h3_degree={O8_h3_degree}|"
    f"D_q6_degree={D_q6_degree}|D_q6_O={D_q6_O}|"
    f"D_h3_degree={D_h3_degree}|D_h3_O={D_h3_O}|"
    f"P_source_vector={','.join(map(str, P_component))}|"
    f"D_source_vector={','.join(map(str, D_component))}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-d13-q24-section-pullback.v4",
    "status": "PASS_EXACT_D13_Q24_SECTION_PULLBACK",
    "coordinate_audit": {
        "q6_raw": True,
        "root_mw": True,
        "q6_simple": True,
        "dominant_neighbor": True,
        "dominant_pinned_d13": True,
        "nef_to_dominant_reflections":
            [[i, p] for i, p in nef_to_dominant],
        "nef_to_dominant_forward_check": True,
        "nef_to_dominant_inverse_check": True,
        "physical_reflection_count": len(physical_reflections),
    },
    "q24_pinned_d13": {
        "mw_coordinates": [0, -1, 1, 1],
        "height": "47",
        "D13_component_correction": "1",
        "P_dot_O": 22,
        "section": list(map(int, P_pinned)),
        "divisor": list(map(int, D_pinned)),
        "decomposition": "D=O+P-10F",
    },
    "actual_component_nef_presentation": {
        "fibre_source_h3_ns": list(map(int, F_component)),
        "zero_source_h3_ns": list(map(int, O_component)),
        "section_source_h3_ns": list(map(int, P_component)),
        "divisor_source_h3_ns": list(map(int, D_component)),
    },
    "pullback_degrees": {
        "P_q6_degree": P_q6_degree,
        "P_q6_zero_intersection": P_q6_O,
        "P_original_H3_degree": P_h3_degree,
        "P_original_H3_zero_intersection": P_h3_O,
        "q8_zero_q6_degree": O8_q6_degree,
        "q8_zero_original_H3_degree": O8_h3_degree,
        "D_q6_degree": D_q6_degree,
        "D_q6_zero_intersection": D_q6_O,
        "D_original_H3_degree": D_h3_degree,
        "D_original_H3_zero_intersection": D_h3_O,
    },
    "boundary": (
        "This certifies the lattice transport of the q24 section/divisor "
        "into the actual component-nef q8 equation presentation. It does "
        "not yet produce rational section coordinates or the D12 equation."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("D13Q24RESULT|status=PASS_EXACT_D13_Q24_SECTION_PULLBACK")
print(f"OUTPUT|{OUTPUT}")
