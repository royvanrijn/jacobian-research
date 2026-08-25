#!/usr/bin/env sage -python
"""Export the exact equation-D13 frame with pinned/suffix fibre markings."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-equation-d13-root-adapted-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-equation-d13-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


source = json.loads(SOURCE.read_text())
assert source["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
frame = matrix(ZZ, source["equation_d13_frame"])
g = block_diagonal_matrix(U2, -frame)
assert frame.dimensions() == (17, 17) and frame.det() == 948
root = frame[:13, :13]
root_minimum = pari(root).qfminim(2)
root_basis = matrix(ZZ, root_minimum[2]).transpose().row_module().basis_matrix()
root_gram = root_basis * root * root_basis.transpose()
root_data = list(map(int, (root_basis.rank(), root_minimum[0], abs(root_gram.det()))))
assert root_data == [13, 312, 4]

identity = identity_matrix(ZZ, 19)
targets = {
    "equation_D13": vector(ZZ, [1, 0] + [0] * 17),
    "pinned_R17": vector(ZZ, source["equation_d13_to_pinned_r17_transition"][0]),
}
transport = None
for index, step in enumerate(source["steps"]):
    local_transport = matrix(ZZ, step["transition"])
    transport = local_transport if transport is None else local_transport * transport
    name = "current_{}_{}".format(index, step["stage"].replace("+", "_"))
    targets[name] = vector(ZZ, transport.row(0))
assert all(value * g * value == 0 for value in targets.values())

FRAME_OUTPUT.write_text(
    "# exact equation-D13 root-adapted frame\n"
    + "\n".join(" ".join(map(str, row)) for row in frame.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-equation-d13-marking.v1",
    "status": "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "hub": "equation_D13",
    "root_data": root_data,
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": rows(identity),
    "root_adapted_hub_to_equation_A11_basis": rows(identity),
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets.items()
    },
    "current_D12_to_pinned_R17_degree": int(
        targets["current_0_D12"] * g * targets["pinned_R17"]
    ),
    "proof_boundary": (
        "Exact marked equation-D13 export from the certified equation-D13-to-pinned-"
        "R17 transport. The identity matrices are the equation-D13 self-marking; their "
        "legacy equation_A11 field names are retained for compatibility with the generic "
        "frontier ranker. No new neighbour is certified by this export alone."
    ),
    "inputs": {
        "paths": [str(SOURCE.relative_to(ROOT))],
        "sha256": {
            str(SOURCE.relative_to(ROOT)): hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "D13MARK|root={}|targets={}|det=1|status={}|output={}".format(
        ",".join(map(str, root_data)), len(targets), payload["status"], OUTPUT
    ),
    flush=True,
)
