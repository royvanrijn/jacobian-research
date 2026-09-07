#!/usr/bin/env sage-python
"""Mark q9/orbit1802 and its cheap MW2 child in the current 3A3 frame."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
PINNED = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
Q9 = GENERATED / "elkies-k3-h3-a11-q9d3o1802-equation-marking.json"
MW2 = GENERATED / "elkies-k3-h3-a11-q9d3o1802-q16d4o114440-marking.json"
OUTPUT = GENERATED / "elkies-k3-h3-current_3A3-q9-mw2-marked-frame.json"

source = json.loads(SOURCE.read_text())
pinned = json.loads(PINNED.read_text())
q9 = json.loads(Q9.read_text())
mw2 = json.loads(MW2.read_text())
assert source["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"
assert pinned["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"
assert q9["status"] == "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT"
assert mw2["status"] == "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT"

equation_to_3a3 = matrix(ZZ, source["equation_A11_to_root_adapted_hub_basis"])
q9_to_equation = matrix(ZZ, q9["child_to_equation_A11_basis"])
mw2_to_q9 = matrix(ZZ, mw2["source_to_root_adapted_hub_basis"])
F = vector(ZZ, [1, 0] + [0] * 17)
q9_fibre_in_equation = F * q9_to_equation
mw2_fibre_in_equation = F * mw2_to_q9 * q9_to_equation

# Check the cross-frame intersection in both coordinate systems.  This also
# guards against confusing the current-route A11 basis with equation A11.
q9_fibre_in_3a3 = q9_fibre_in_equation * equation_to_3a3
mw2_fibre_in_3a3 = mw2_fibre_in_equation * equation_to_3a3

payload = dict(source)
targets = dict(source["target_fibres_in_root_adapted_hub"])
targets.update({
    "q9d3o1802": list(map(int, q9_fibre_in_3a3)),
    "q9d3o1802_q16d4o114440_mw2": list(
        map(int, mw2_fibre_in_3a3)
    ),
})
payload.update({
    "schema": "elkies-k3.h3-current-3a3-q9-mw2-marking.v1",
    "hub": "current_3A3_q9_mw2_marked",
    "target_fibres_in_root_adapted_hub": targets,
    "added_target_degrees": {
        name: int(values[1]) for name, values in targets.items()
        if name.startswith("q9d3o1802")
    },
    "proof_boundary": (
        source["proof_boundary"] + " The q9/orbit1802 fibre and its q16/degree-4 "
        "MW2 child are additionally transported through complete determinant-one "
        "basis maps; no crossover edge is asserted."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SOURCE, PINNED, Q9, MW2)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (SOURCE, PINNED, Q9, MW2)
        },
    },
})
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "3A3Q9MARK|q9_degree={}|mw2_degree={}|status={}|output={}".format(
        payload["added_target_degrees"]["q9d3o1802"],
        payload["added_target_degrees"]["q9d3o1802_q16d4o114440_mw2"],
        payload["status"], OUTPUT,
    ), flush=True,
)
