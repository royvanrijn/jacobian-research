#!/usr/bin/env sage -python
"""Turn the safe initial q8 node-leading rows into an exact compiler block.

Only singleton Pareto-minimal bidegrees from the all-node template are used:
regular unit corrections can only increase both local exponents, so they
cannot affect those coefficients.  The condition is normalized by its
nonzero chart unit and then assembled through ``compile_resolved_conditions``.
All later equal/nonminimal node groups remain outside this block.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
TEMPLATE = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-e7-node-principal-bidegrees.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-initial-node-conditions.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--template", type=Path, default=TEMPLATE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
template = json.loads(args.template.read_text())
assert ambient["status"] in (
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT", "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
)
assert template["status"] == "PASS_EXACT_Q8_ALL_E7_NODE_PRINCIPAL_BIDEGREE_TEMPLATE"
assert template["inputs"]["endpoint_ambient"]["sha256"] == digest(args.ambient)
assert template["ambient_dimension"] == ambient["ambient_dimension"]
exec(compile(CORE.read_text(), str(CORE), "exec"))

indices = tuple(int(value) for value in template["exact_initial_unique_basis_indices"])
assert len(indices) == len(set(indices))
assert all(0 <= index < int(ambient["ambient_dimension"]) for index in indices)
ambient_labels = tuple(range(int(ambient["ambient_dimension"])))
block = quotient_condition(
    "six actual E7 node singleton Pareto minima",
    ambient_labels,
    lambda index: tuple(QQ(1) if index == row else QQ(0) for row in indices),
    tuple("normalized node leading coefficient of ambient basis {}".format(index) for index in indices),
    "actual six-node principal-bidegree template; only singleton Pareto-minimal terms",
)
compilation = compile_resolved_conditions(
    ambient_labels, (block,), complete=False, compute_kernel=False
)
assert compilation["rank"] == len(indices)
assert compilation["kernel_dimension"] == len(ambient_labels)-len(indices)

payload = {
    "schema": "elkies-k3.h92-q8-initial-node-conditions.v1",
    "status": "PASS_EXACT_Q8_INITIAL_NODE_CONDITION_BLOCK",
    "inputs": {
        "endpoint_ambient": {"path": path_label(args.ambient), "sha256": digest(args.ambient)},
        "node_template": {"path": path_label(args.template), "sha256": digest(args.template)},
        "compiler_core": {"path": path_label(CORE), "sha256": digest(CORE)},
    },
    "condition_rule": template["exact_initial_constraint_rule"],
    "ambient_basis_indices": list(indices),
    "ambient_basis": [ambient["ambient_basis"][index] for index in indices],
    "matrix": {
        "rows": compilation["condition_rows"], "columns": compilation["ambient_dimension"],
        "rank": compilation["rank"], "kernel_dimension": compilation["kernel_dimension"],
        "normalization": "each nonzero actual leading unit is divided out",
    },
    "boundary": (
        "This contains only the initial singleton Pareto-leading node rows. "
        "It does not evaluate the finite two-variable quotient for the other "
        "node groups, combine smooth/generic/marked/overlap blocks, or certify "
        "h0, a pencil, or a child equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8INITIALNODECONDITIONS|ambient={}|rows={}|rank={}|kernel={}|"
    "status=PASS_EXACT_Q8_INITIAL_NODE_CONDITION_BLOCK".format(
        compilation["ambient_dimension"], compilation["condition_rows"],
        compilation["rank"], compilation["kernel_dimension"],
    ),
    flush=True,
)
