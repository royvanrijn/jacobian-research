#!/usr/bin/env sage-python
"""Add every certified current-suffix fibre to the q9/orbit1802 marking."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_SOURCE = GENERATED / "elkies-k3-h3-a11-q9d3o1802-equation-marking.json"
DEFAULT_SUFFIX = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-a11-q9d3o1802-current-suffix-marking.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--suffix", type=Path, default=DEFAULT_SUFFIX)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

source_path = args.source.resolve()
suffix_path = args.suffix.resolve()
output_path = args.output.resolve()
source = json.loads(source_path.read_text())
suffix = json.loads(suffix_path.read_text())
assert source["status"] == "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT"
assert suffix["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"

# Rows of pinned_in_equation are pinned basis vectors in equation-A11
# coordinates; rows of equation_to_child are equation-A11 basis vectors in
# q9-child coordinates. Row-vector coordinates therefore compose on the right.
pinned_in_equation = matrix(
    ZZ, suffix["equation_A11_to_root_adapted_hub_basis"]
)
equation_to_child = matrix(ZZ, source["equation_A11_to_child_basis"])
suffix_targets = {
    name: vector(ZZ, value) * pinned_in_equation * equation_to_child
    for name, value in suffix["target_fibres_in_root_adapted_hub"].items()
    if name.startswith("current_")
}

payload = dict(source)
targets = dict(source["target_fibres_in_child"])
targets.update({name: list(map(int, value)) for name, value in suffix_targets.items()})
assert targets["current_A11"] == targets["equation_A11"]
assert targets["current_A5A5"] == targets["orbit12"]
assert targets["current_rootless"] == targets["pinned_R17"]
payload.update({
    "schema": "elkies-k3.h3-a11-q9d3o1802-current-suffix-marking.v1",
    "source_hub": "a11_q9d3o1802_explicit_zero_current_suffix_marked",
    "target_fibres_in_child": targets,
    "current_suffix_target_names": sorted(suffix_targets),
    "proof_boundary": (
        source["proof_boundary"] + " Every certified current-route suffix fibre is "
        "additionally transported through the full pinned/equation-A11 and "
        "equation-A11/q9 determinant-one maps; no new edge is asserted."
    ),
    "inputs": {
        "paths": [str(source_path.relative_to(ROOT)), str(suffix_path.relative_to(ROOT))],
        "sha256": {
            str(source_path.relative_to(ROOT)): hashlib.sha256(source_path.read_bytes()).hexdigest(),
            str(suffix_path.relative_to(ROOT)): hashlib.sha256(suffix_path.read_bytes()).hexdigest(),
        },
    },
})
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q9SUFFIXMARK|added={}|A11={}|A5A5={}|rootless={}|status={}|output={}".format(
        len(suffix_targets), int(targets["current_A11"] == targets["equation_A11"]),
        int(targets["current_A5A5"] == targets["orbit12"]),
        int(targets["current_rootless"] == targets["pinned_R17"]),
        payload["status"], output_path,
    ),
    flush=True,
)
