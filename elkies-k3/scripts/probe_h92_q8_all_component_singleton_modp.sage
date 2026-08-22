#!/usr/bin/env sage -python
"""Apply exact singleton generic-E7 conditions to a modular q=8 smooth kernel.

The all-component generic-condition compiler identifies coefficients whose
negative leading order on an actual E7 component is unique.  Each such
coefficient must vanish independently of every unresolved leading-residue or
node calculation.  This script restricts those exact coordinate rows to a
declared modular smooth kernel and reports the surviving dimension.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, matrix


ROOT = Path(__file__).resolve().parents[2]
KERNEL = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra7.json"
CONDITIONS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra7.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-singleton-mod-43-extra7.json"


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel", type=Path, default=KERNEL)
parser.add_argument("--conditions", type=Path, default=CONDITIONS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

kernel = json.loads(args.kernel.read_text())
conditions = json.loads(args.conditions.read_text())
assert kernel["status"] == "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK"
assert conditions["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"
assert "kernel_basis_rows" in kernel
assert int(kernel["dimensions"]["kernel"]) == len(kernel["kernel_basis_rows"])
assert len(kernel["ambient_basis"]) == int(kernel["dimensions"]["columns"])
signature = [
    [entry["kind"], int(entry["x_power"]), int(entry["m_power"]),
     int(entry["u_power"]), int(entry["h_power"])]
    for entry in kernel["ambient_basis"]
]
assert conditions["ambient_basis_sha256"] == hashlib.sha256(
    json.dumps(signature, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()
# The ambient-basis lists are authoritative for alignment, including an
# enlarged endpoint file outside the repository.  Equality is exact JSON,
# not a positional assumption about a particular r value.
template_ambient = conditions["inputs"]["endpoint_ambient"]
if template_ambient["path"].startswith("artifacts/"):
    endpoint = ROOT / template_ambient["path"]
    assert endpoint.is_file()
    endpoint_payload = json.loads(endpoint.read_text())
    endpoint_signature = [
        [entry["kind"], int(entry["x_power"]), int(entry["m_power"]),
         int(entry["u_power"]), int(entry["h_power"])]
        for entry in endpoint_payload["ambient_basis"]
    ]
    assert endpoint_signature == signature

prime = int(kernel["prime"])
field = GF(prime)
smooth = matrix(field, [[field(int(value)) for value in row] for row in kernel["kernel_basis_rows"]])
indices = tuple(int(value) for value in conditions["singleton_coordinate_block"]["basis_indices"])
assert all(0 <= index < smooth.ncols() for index in indices)
restriction = matrix(
    field,
    [[smooth[row, index] for row in range(smooth.nrows())] for index in indices],
)
rank = restriction.rank()
assert rank <= smooth.nrows()

payload = {
    "schema": "elkies-k3.h92-q8-all-component-singleton-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_ALL_COMPONENT_SINGLETON_OBSTRUCTION",
    "prime": prime,
    "inputs": {
        "smooth_kernel": path_label(args.kernel),
        "generic_component_conditions": path_label(args.conditions),
    },
    "dimensions": {
        "smooth_kernel": int(smooth.nrows()),
        "singleton_coordinate_rows": len(indices),
        "restriction_rank": int(rank),
        "survivor_dimension": int(smooth.nrows()-rank),
    },
    "singleton_basis_indices": list(indices),
    "conclusion": (
        "Every modular smooth-kernel direction is eliminated by exact unique "
        "negative generic-component E7 leading terms."
        if rank == smooth.nrows() else
        "The displayed singleton generic-component conditions leave the stated "
        "modular survivor dimension; unresolved non-singleton and node conditions "
        "remain."
    ),
    "boundary": (
        "This is a necessary all-component generic-E7 obstruction on one modular "
        "smooth kernel. It does not evaluate non-singleton residue cancellations, "
        "nodes, overlaps, a characteristic-zero kernel, or a q8 pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8ALLSINGLETON|prime={}|smooth_kernel={}|rows={}|rank={}|survivor={}|"
    "status=EXPERIMENTAL_MODULAR_ALL_COMPONENT_SINGLETON_OBSTRUCTION".format(
        prime, smooth.nrows(), len(indices), rank, smooth.nrows()-rank
    ),
    flush=True,
)
