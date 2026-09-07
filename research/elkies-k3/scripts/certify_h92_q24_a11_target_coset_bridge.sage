#!/usr/bin/env sage -python
"""Certify the cheapest constructed-section bridge to the A11 q8 target.

The exact orbit42 lift supplies eighteen rational identity-shell sections on
the A11 equation.  Under the pinned orbit64 marking they span rank five, but
with index five in the saturated hyperplane whose sixth MW coordinate is zero.
Consequently an arbitrary point in the missing sixth direction need not
generate the q8/orbit12 target together with those exact sections.

This script audits that construction rather than only the abstract MW rank. It
enumerates the complete target coset through pole order five and selects

    M = (1,0,0,0,0,1),  h(M)=47/4,  correction=9/4,  M.O=5.

The desired orbit12 section is then the short exact group-law word

    P12 = M + S6 - 2*S2 - 2*S8,

where the S_i use zero-based equation-shell indices.  The output fixes the
correct lower-pole construction target; it does not yet construct M's
Weierstrass coordinates.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-target-coset-bridge.json",
)
args = parser.parse_args()

IDENTITY_PATH = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING_PATH = LOCAL / "q24-orbit42-identity-halving-qq.json"
ZERO_PATH = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
NEIGHBOURS_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
PARENT_FRAME_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
CHILD_FRAME_PATH = (
    LOCAL
    / "q24-downstream-lift/d12-c10a-zero-q6-frames/"
    "q6-o0064-r11-n132-d12-ad4a027cb197.txt"
)
INPUTS = (
    IDENTITY_PATH,
    MATCHING_PATH,
    ZERO_PATH,
    NEIGHBOURS_PATH,
    PARENT_FRAME_PATH,
    CHILD_FRAME_PATH,
)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


identity = json.loads(IDENTITY_PATH.read_text())
matching = json.loads(MATCHING_PATH.read_text())
zero = json.loads(ZERO_PATH.read_text())
neighbours = json.loads(NEIGHBOURS_PATH.read_text())
assert identity["status"] == "PASS_Q42_ORBIT42_IDENTITY_HALVING_LATTICE_GATE"
assert matching["status"] == "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert len(zero["sections"]) == 18

# Transport the abstract identity shell into the selected equation-side A11
# frame, then reorder it by the certified abstract-to-equation shell mapping.
parent_frame = load_matrix(PARENT_FRAME_PATH)
parent_root = parent_frame[:12, :12]
parent_coupling = parent_frame[:12, 12:]
old_section_classes = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * parent_coupling.transpose()) * parent_root.inverse()
    assert all(value in ZZ for value in root_coefficients)
    old_section_classes.append(
        vector(ZZ, [1, 1] + list(map(ZZ, root_coefficients)) + list(z))
    )

record = next(
    row for row in neighbours["neighbors"] if int(row["orbit_index"]) == 64
)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, record["child_root_adapted_basis"])
) * matrix(ZZ, record["neighbor_basis"])
assert abs(transition.det()) == 1
abstract_mw = [
    vector(ZZ, (section * transition.inverse())[-6:])
    for section in old_section_classes
]
shell_mapping = matching["matching"]["mappings_abstract_to_equation"][7]
equation_mw = [None] * 18
for abstract_index, equation_index in enumerate(shell_mapping):
    equation_mw[equation_index] = abstract_mw[abstract_index]
assert all(item is not None and item[-1] == 0 for item in equation_mw)

known_lattice = matrix(ZZ, [list(item[:5]) for item in equation_mw]).row_module()
known_basis = known_lattice.basis_matrix()
assert known_lattice.rank() == 5
assert abs(known_basis.det()) == 5

# The equation-side construction target selected by the exact fingerprint
# certificate.
target = vector(ZZ, (0, 0, -1, 0, 0, 1))

# Complete MW/correction enumeration through P.O <= 5.  Since
# h+corr=4+2(P.O), every such section has h <= 14.
child_frame = load_matrix(CHILD_FRAME_PATH)
assert child_frame == matrix(ZZ, record["child_root_adapted_frame"])
root_rank = 11
root = child_frame[:root_rank, :root_rank]
coupling = child_frame[:root_rank, root_rank:]
tail = child_frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
denominator = lcm(value.denominator() for value in height.list())
height_integral = matrix(ZZ, denominator * height)
enumeration = pari(height_integral).qfminim(ZZ(14 * denominator), flag=2)
mw_vectors = {tuple([0] * 6)}
for column in matrix(ZZ, enumeration[2]).columns():
    mw_vectors.add(tuple(column))
    mw_vectors.add(tuple(-column))

root_lattice = IntegralLattice(root)


def section_profile(z):
    z = vector(ZZ, z)
    h = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * child_frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 6)
        norm = QQ(lifted * child_frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
    correction = minimum - h
    pole_order = (h + correction - 4) / 2
    class_order = ZZ(1)
    for value in dual:
        class_order = lcm(class_order, ZZ(QQ(value).denominator()))
    return h, correction, pole_order, class_order


coset_rows = []
for values in sorted(mw_vectors):
    z = vector(ZZ, values)
    if z[-1] != 1 or vector(ZZ, (target - z)[:5]) not in known_lattice:
        continue
    h, correction, pole_order, class_order = section_profile(z)
    if pole_order in ZZ and 0 <= pole_order <= 5:
        coset_rows.append(
            {
                "mw": list(map(int, z)),
                "height": str(h),
                "local_correction": str(correction),
                "component_class_order": int(class_order),
                "P_dot_O": int(pole_order),
                "mw_l1": int(sum(abs(value) for value in z)),
            }
        )

coset_rows.sort(
    key=lambda row: (
        row["P_dot_O"],
        QQ(row["height"]),
        row["mw_l1"],
        tuple(row["mw"]),
    )
)
assert coset_rows == [
    {
        "mw": [1, 0, 0, 0, 0, 1],
        "height": "47/4",
        "local_correction": "9/4",
        "component_class_order": 4,
        "P_dot_O": 5,
        "mw_l1": 2,
    },
    {
        "mw": [2, 0, 1, 0, 0, 1],
        "height": "14",
        "local_correction": "0",
        "component_class_order": 1,
        "P_dot_O": 5,
        "mw_l1": 4,
    },
]
selected = vector(ZZ, coset_rows[0]["mw"])

# A short word in the exact equation-shell points.  Verify it literally in
# the selected A11 MW marking.
shell_coefficients = [0] * 18
shell_coefficients[2] = -2
shell_coefficients[6] = 1
shell_coefficients[8] = -2
shell_word = sum(
    (shell_coefficients[index] * equation_mw[index] for index in range(18)),
    vector(ZZ, [0] * 6),
)
assert selected + shell_word == target
assert sum(abs(value) for value in shell_coefficients) == 5

payload = {
    "schema": "elkies-k3.h3-q24-a11-target-coset-bridge.v1",
    "status": "PASS_EXACT_A11_TARGET_COSET_BRIDGE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "equation_frame": {
        "orbit": 64,
        "shell_mapping": 7,
        "MW_height_gram": [[str(value) for value in row] for row in height.rows()],
    },
    "exact_identity_shell": {
        "count": 18,
        "MW_vectors_in_equation_order": [list(map(int, item)) for item in equation_mw],
        "rank": 5,
        "saturated_hyperplane_index": int(abs(known_basis.det())),
        "row_lattice_basis_first_five_coordinates": [
            list(map(int, row)) for row in known_basis.rows()
        ],
    },
    "target": {
        "equation_q8_orbit": 12,
        "MW": list(map(int, target)),
        "height": "13",
        "local_correction": "3",
        "P_dot_O": 6,
    },
    "complete_target_coset_through_P_dot_O_5": coset_rows,
    "selected_bridge": {
        **coset_rows[0],
        "I12_component_up_to_negation": [3, 9],
        "selection_rule": "minimum height, then MW L1, among minimum-pole target-coset sections",
    },
    "group_law_bridge": {
        "formula": "P12=M+S6-2*S2-2*S8",
        "shell_indices_are_zero_based_equation_order": True,
        "shell_coefficients": shell_coefficients,
        "verified_MW_sum": list(map(int, selected + shell_word)),
    },
    "proof_boundary": (
        "This exactly proves the complete minimum-pole target-coset search and "
        "the short group-law reduction to already exact A11 shell sections. It "
        "does not construct characteristic-zero Weierstrass coordinates for M; "
        "without M the q8 Riemann--Roch pencil and 2A5 child remain open."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11TARGETCOSET|shell_rank=5|shell_index=5|minimum_PO=5|"
    "candidates=2|selected_mw=1,0,0,0,0,1|height=47/4|corr=9/4|"
    "formula=M+S6-2S2-2S8|status=PASS_EXACT_A11_TARGET_COSET_BRIDGE",
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
