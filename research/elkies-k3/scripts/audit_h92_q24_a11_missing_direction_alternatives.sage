#!/usr/bin/env sage -python
"""Audit no-Groebner alternatives for the missing A11 MW direction.

The exact orbit42 equation has eighteen identity-shell points and an exact
opposite spinor pair on its D12 parent.  Three of those old curves have degree
one over the A11 base and therefore give exact characteristic-zero points by
Mobius inversion and binary-quartic covariants.  This script decides, purely
in the certified Neron--Severi lattices, whether those extra points can carry
the sixth A11 MW coordinate required by the A11-to-2A5 construction.

It also checks the other construction-compatible q8 orbit.  No section
ansatz, Groebner basis, finite-field search, or Hensel lift is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-missing-direction-alternatives.json",
)
args = parser.parse_args()

PARENT_FRAME_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY_PATH = LOCAL / "q24-orbit42-identity-halving-audit.json"
SPINOR_PATH = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
COVARIANT_PATH = LOCAL / "q24-a11-degree-one-shell-covariants-qq.json"
BRIDGE_PATH = LOCAL / "q24-a11-target-coset-bridge.json"
FINGERPRINT_PATH = LOCAL / "q24-a11-q8-construction-fingerprint.json"
INPUTS = (
    PARENT_FRAME_PATH,
    Q6_PATH,
    IDENTITY_PATH,
    SPINOR_PATH,
    COVARIANT_PATH,
    BRIDGE_PATH,
    FINGERPRINT_PATH,
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


parent = load_matrix(PARENT_FRAME_PATH)
q6 = json.loads(Q6_PATH.read_text())
identity = json.loads(IDENTITY_PATH.read_text())
spinor = json.loads(SPINOR_PATH.read_text())
covariant = json.loads(COVARIANT_PATH.read_text())
bridge = json.loads(BRIDGE_PATH.read_text())
fingerprint = json.loads(FINGERPRINT_PATH.read_text())

assert q6["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert identity["status"] == "PASS_Q42_ORBIT42_IDENTITY_HALVING_LATTICE_GATE"
assert spinor["status"] == "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ"
assert covariant["status"] == "PASS_EXACT_A11_DEGREE_ONE_CURVE_COVARIANTS_QQ"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"
assert fingerprint["status"] == "PASS_Q24_A11_Q8_CONSTRUCTION_FINGERPRINT"

record = next(row for row in q6["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, record["child_root_adapted_basis"])
) * matrix(ZZ, record["neighbor_basis"])
assert abs(transition.det()) == 1
inverse = transition.inverse()

# In the selected frames the last child MW coordinate is literally the fifth
# parent MW coordinate.  Computing this on every ambient basis vector makes
# the statement independent of choices of vertical correction for a section.
source_basis_to_child_sixth = []
for index in range(parent.nrows() + 2):
    basis_vector = vector(ZZ, parent.nrows() + 2)
    basis_vector[index] = 1
    source_basis_to_child_sixth.append(int((basis_vector * inverse)[-1]))
assert source_basis_to_child_sixth == [0] * 18 + [1]

root_rank = 12
root = parent[:root_rank, :root_rank]
coupling = parent[:root_rank, root_rank:]
tail = parent[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling


def fractional_key(values):
    return tuple(QQ(value) - QQ(value).floor() for value in values)


root_inverse = root.inverse()
correction_by_class = {fractional_key(vector(QQ, [0] * root_rank)): QQ(0)}
for index in range(root_rank):
    weight = vector(QQ, root_inverse.row(index))
    key = fractional_key(weight)
    norm = QQ(weight * root * weight)
    if key not in correction_by_class or norm < correction_by_class[key]:
        correction_by_class[key] = norm
assert sorted(correction_by_class.values()) == [QQ(0), QQ(1), QQ(3), QQ(3)]


def correction_for(mw):
    mw = vector(ZZ, mw)
    base = vector(ZZ, [0] * root_rank + list(mw))
    dual = vector(QQ, base * parent[:, :root_rank]) * root_inverse
    return correction_by_class[fractional_key(dual)]


# Reconstruct the complete parent P.O=0 shell abstractly.  Its two
# correction-three vectors are precisely the spinor pair represented by the
# exact characteristic-zero artifact.
scale = ZZ(1)
for value in height.list():
    scale = scale.lcm(ZZ(QQ(value).denominator()))
enumeration = pari((scale * height).change_ring(ZZ)).qfminim(ZZ(4 * scale))
zero_pole_vectors = set()
spinor_vectors = set()
for column in matrix(ZZ, enumeration[2]).columns():
    for sign in (1, -1):
        mw = sign * vector(ZZ, column)
        h = QQ(mw * height * mw)
        correction = correction_for(mw)
        if (h + correction - 4) / 2 == 0:
            key = tuple(map(int, mw))
            zero_pole_vectors.add(key)
            if correction == 3:
                spinor_vectors.add(key)
assert spinor_vectors == {(-1, 0, 0, 0, 0), (1, 0, 0, 0, 0)}

identity_vectors = [
    tuple(map(int, values))
    for values in identity["exact_model_R3_zero"]["identity_vectors"]
]
assert len(identity_vectors) == 18
assert all(values[-1] == 0 for values in identity_vectors)
assert all(values[-1] == 0 for values in spinor_vectors)

degree_one_rows = covariant["points"]
assert [(row["shell_kind"], int(row["equation_shell_index"])) for row in degree_one_rows] == [
    ("identity", 7),
    ("identity", 17),
    ("spinor", 0),
]
assert all(int(row["new_base_degree"]) == 1 for row in degree_one_rows)

selected_bridge = vector(ZZ, bridge["selected_bridge"]["mw"])
target = vector(ZZ, bridge["target"]["MW"])
assert selected_bridge[-1] == target[-1] == 1

compatible = {
    int(row["orbit_index"]): vector(ZZ, row["mw_projection"])
    for row in fingerprint["construction_compatible_orbits"]
}
assert set(compatible) == {12, 2162}
assert compatible[12][-1] == 1
assert compatible[2162][-1] == -1

# The unique primitive parent coordinate carrying the missing quotient has a
# lower parent pole profile, but it is not among the currently explicit
# zero-pole curves.
parent_missing = vector(ZZ, (0, 0, 0, 0, 1))
parent_missing_height = QQ(parent_missing * height * parent_missing)
parent_missing_correction = correction_for(parent_missing)
parent_missing_pole_order = (
    parent_missing_height + parent_missing_correction - 4
) / 2
assert (parent_missing_height, parent_missing_correction, parent_missing_pole_order) == (
    QQ(12),
    QQ(0),
    QQ(4),
)

payload = {
    "schema": "elkies-k3.h3-q24-a11-missing-direction-alternatives.v1",
    "status": "PASS_EXACT_A11_MISSING_DIRECTION_ALTERNATIVES_AUDIT",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "selected_transition": {
        "parent": "D12/MW5",
        "child": "A11/MW6 orbit64",
        "source_basis_to_child_sixth_coordinate": source_basis_to_child_sixth,
        "conclusion": "child_MW_6_equals_parent_MW_5",
    },
    "exact_zero_pole_shell": {
        "identity_count": len(identity_vectors),
        "spinor_vectors": [list(values) for values in sorted(spinor_vectors)],
        "all_parent_fifth_coordinates_zero": True,
        "degree_one_covariant_curves": [
            {
                "kind": row["shell_kind"],
                "index": int(row["equation_shell_index"]),
            }
            for row in degree_one_rows
        ],
        "carries_missing_child_coordinate": False,
    },
    "construction_alternatives": {
        "orbit12_MW": list(map(int, compatible[12])),
        "orbit2162_MW": list(map(int, compatible[2162])),
        "both_require_missing_coordinate_up_to_sign": True,
    },
    "smallest_parent_coordinate_carrier": {
        "D12_MW": list(map(int, parent_missing)),
        "height": str(parent_missing_height),
        "local_correction": str(parent_missing_correction),
        "P_dot_O": int(parent_missing_pole_order),
        "present_in_exact_zero_pole_shell": False,
    },
    "required_child_bridge": {
        "MW": list(map(int, selected_bridge)),
        "P_dot_O": int(bridge["selected_bridge"]["P_dot_O"]),
        "target_MW": list(map(int, target)),
    },
    "proof_boundary": (
        "Exact lattice replay plus the exact characteristic-zero degree-one "
        "covariant certificate. It proves that identity and spinor zero-pole "
        "transports cannot supply the sixth A11 MW coordinate, and that switching "
        "to construction-compatible orbit2162 only reverses its sign. It does not "
        "construct the parent P.O=4 carrier or the child bridge M."
    ),
    "next_required": (
        "Construct either the parent D12 P.O=4 section with MW (0,0,0,0,1) "
        "together with an exact divisor-to-child transport, or directly construct "
        "the A11 P.O=5 bridge M=(1,0,0,0,0,1), preferably by resolved linear "
        "Riemann--Roch rather than a large polynomial Groebner system."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11MISSINGALT|identity=18|spinor=2|degree1=3|"
    "shell_child6=0|orbit12_child6=1|orbit2162_child6=-1|"
    "parent_carrier=0,0,0,0,1|parent_PO=4|"
    f"status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{args.output}", flush=True)
