#!/usr/bin/env sage -python
"""Export the exact E8+E6/MW3 source marking for the first H3 q8 edge."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
FRAME_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-first-q8-source-frame.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-first-q8-source-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

data = json.loads(SOURCE.read_text())
assert data["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
q6 = data["q6"]
q8 = data["q8"]
frame = matrix(ZZ, q8["simple_frame_gram"])
root_mw = matrix(ZZ, q6["root_mw_basis_in_child"])
simple_root_change = matrix(ZZ, q8["simple_root_change_in_root_block"])
simple_change = block_diagonal_matrix(simple_root_change, identity_matrix(ZZ, 3))
# The independently nef representative is moved by the four recorded source
# root reflections to the first dominant D13 orbit.  Use that chamber-aligned
# class for outgoing-neighbour gates.
hit = q8["d13_mw4_hits"][0]
d13_fibre = vector(ZZ, [4, 2] + list(hit["witness_simple_frame"]))
g = block_diagonal_matrix(U2, -frame)
assert d13_fibre * g * d13_fibre == 0
assert d13_fibre[1] == 2

# Record the full selected dominant D13 child basis in this source coordinate
# system.  Despite its historical field name, this neighbour basis was built
# directly against ``simple_frame_gram``.
neighbor_q6 = matrix(ZZ, hit["neighbor_basis_in_q6_ns"])
d13_adaptation = matrix(ZZ, hit["d13_root_adapted_basis_in_child"])
d13_in_source = (
    block_diagonal_matrix(identity_matrix(ZZ, 2), d13_adaptation)
    * neighbor_q6
)
assert abs(d13_in_source.det()) == 1
assert d13_in_source.row(0) == d13_fibre

FRAME_OUTPUT.write_text("\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n")
payload = {
    "schema": "elkies-k3.h3-first-q8-source-marking.v1",
    "status": "PASS_EXACT_FIRST_Q8_SOURCE_MARKING",
    "hub": "first_q8_E8_plus_E6_MW3",
    "root_data": [14, 312, 3],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_D13_basis_in_root_adapted_hub": [
        [int(x) for x in row] for row in d13_in_source.rows()
    ],
    "root_adapted_hub_basis_in_equation_D13": [
        [int(x) for x in row]
        for row in d13_in_source.inverse().change_ring(ZZ).rows()
    ],
    "target_fibres_in_root_adapted_hub": {
        "equation_D13": [int(x) for x in d13_fibre],
        # The zero-loop search only carries this field as an endpoint label;
        # final pinned identification is composed through the D13 basis above.
        "pinned_R17": [int(x) for x in d13_fibre],
    },
    "proof_boundary": (
        "Exact q6-child frame, primitive nef first-q8 D13 fibre, marked U, and "
        "determinant-one D13 child basis. This does not claim a cheaper route."
    ),
    "inputs": {
        "paths": [str(SOURCE.relative_to(ROOT))],
        "sha256": {str(SOURCE.relative_to(ROOT)): hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIRSTQ8MARK|degree={}|root=14,312,3|det={}|status={}|output={}".format(
        d13_fibre[1], int(d13_in_source.det()), payload["status"], OUTPUT
    ), flush=True,
)
