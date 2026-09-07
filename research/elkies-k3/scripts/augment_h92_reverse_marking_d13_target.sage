#!/usr/bin/env sage -python
"""Add the exact equation-D13 fibre to an equation-A11-marked reverse hub.

The source marking and the canonical equation-D13 marking both use the same
equation-A11 coordinate system.  This keeps reverse searches reproducible and
lets the marked frontier rank against D13 directly, instead of using ADE/MW
types or the much less informative orbit12 degree as a proxy.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--marking", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--d13-route",
    type=Path,
    default=LOCAL / "q24-equation-d13-to-pinned-r17.json",
)
parser.add_argument(
    "--crossovers",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json",
)
args = parser.parse_args()

source_path = args.marking.resolve()
d13_path = args.d13_route.resolve()
crossovers_path = args.crossovers.resolve()
output_path = args.output.resolve()
source = json.loads(source_path.read_text())
d13 = json.loads(d13_path.read_text())
crossovers = json.loads(crossovers_path.read_text())
assert source["status"] == "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
assert d13["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"

# Row-basis convention: convert the D13 fibre to the canonical pinned basis,
# then from that basis to the certified equation-A11 basis.  The D13 marking's
# legacy equation_A11 matrix fields are self-markings and must not be used for
# this crossover.
d13_to_pinned = matrix(ZZ, d13["equation_d13_to_pinned_r17_transition"])
pinned_in_equation = matrix(ZZ, crossovers["pinned_R17_basis_in_equation_A11"])
d13_in_equation = (
    vector(ZZ, [1] + [0] * 18)
    * d13_to_pinned.inverse().change_ring(ZZ)
    * pinned_in_equation
)
equation_in_source = matrix(ZZ, source["root_adapted_hub_to_equation_A11_basis"])
d13_in_source = d13_in_equation * equation_in_source

frame_path = ROOT / source["frame_output"]
frame = matrix(
    ZZ,
    [
        [ZZ(value) for value in line.split()]
        for line in frame_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ],
)
g = matrix(ZZ, 19, 19)
g[0, 1] = g[1, 0] = 1
g[2:, 2:] = -frame
assert d13_in_source * g * d13_in_source == 0

payload = dict(source)
payload["schema"] = "elkies-k3.h3-reverse-hub-equation-marking-with-d13.v1"
payload["status"] = "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
payload["target_fibres_in_root_adapted_hub"] = dict(
    source["target_fibres_in_root_adapted_hub"]
)
payload["target_fibres_in_root_adapted_hub"]["equation_D13"] = [
    int(value) for value in d13_in_source
]
payload["d13_target_derivation"] = {
    "common_coordinates": "equation_A11 row-basis coordinates",
    "isotropic": True,
    "source_marking": str(source_path.relative_to(ROOT)),
    "source_marking_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    "d13_marking": str(d13_path.relative_to(ROOT)),
    "d13_marking_sha256": hashlib.sha256(d13_path.read_bytes()).hexdigest(),
    "crossover_marking": str(crossovers_path.relative_to(ROOT)),
    "crossover_marking_sha256": hashlib.sha256(crossovers_path.read_bytes()).hexdigest(),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "REVHUBD13|hub={}|d13_degree={}|status={}|output={}".format(
        payload["hub"], int(d13_in_source[1]), payload["status"], output_path
    ),
    flush=True,
)
