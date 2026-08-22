#!/usr/bin/env sage -python
"""Screen a smooth q=8 kernel by all exact generic E7 leading terms.

The all-component condition template groups ambient monomials having the
same negative order on each *actual* resolved E7 component.  A group whose
restriction to the current candidate space has just one nonzero coefficient
cannot cancel, even when the original ambient group had several monomials.
This gives an exact necessary modular cut.  Unlike the earlier unmarked
screen, it includes E7_5 using its audited exact chord order ord(m)=0 rather
than an inferred order from the singular Weierstrass equation.

This is still only a generic-point screen: groups with two or more live
coefficients require their actual leading-residue relation, and node, marked
branch, and overlap calculations remain outside it.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
KERNEL = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json"
CONDITIONS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-module-mod-43-extra10.json"


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def basis_signature(basis):
    return [
        [entry["kind"], int(entry["x_power"]), int(entry["m_power"]),
         int(entry["u_power"]), int(entry["h_power"])]
        for entry in basis
    ]


def canonical_digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel", type=Path, default=KERNEL)
parser.add_argument("--conditions", type=Path, default=CONDITIONS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

kernel = json.loads(args.kernel.read_text())
conditions = json.loads(args.conditions.read_text())
assert kernel["status"] == "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK"
assert conditions["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"
assert canonical_digest(basis_signature(kernel["ambient_basis"])) == conditions["ambient_basis_sha256"]

finite = GF(int(kernel["prime"]))
smooth = matrix(
    finite,
    [[finite(int(value)) for value in row] for row in kernel["kernel_basis_rows"]],
)
assert smooth.nrows() == int(kernel["dimensions"]["kernel"])
assert smooth.ncols() == len(kernel["ambient_basis"])

# The order is deterministic and records every independent coordinate cut in
# the current candidate basis.  The candidate rows always remain ambient
# coefficient vectors, so a one-label group is an honest linear condition.
candidate = smooth
records = []
for component in conditions["component_conditions"]:
    for group in component["negative_order_groups"]:
        labels = [
            int(label) for label in group["basis_indices"]
            if any(candidate[row, int(label)] for row in range(candidate.nrows()))
        ]
        if len(labels) != 1:
            continue
        label = labels[0]
        coefficient = vector(
            finite, [candidate[row, label] for row in range(candidate.nrows())]
        )
        assert coefficient
        before = candidate.nrows()
        candidate = matrix(
            finite,
            matrix(finite, 1, before, [coefficient]).right_kernel().basis_matrix(),
        ) * candidate
        records.append({
            "component": component["component"],
            "residual_order": int(group["residual_order"]),
            "basis_index": label,
            "basis_label": kernel["ambient_basis"][label],
            "candidate_dimension_before": int(before),
            "reason": "the negative-order group has exactly one live ambient coefficient",
        })
        if not candidate.nrows():
            break
    if not candidate.nrows():
        break

payload = {
    "schema": "elkies-k3.h92-q8-all-component-generic-module-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_ALL_COMPONENT_GENERIC_MODULE_SCREEN",
    "prime": int(kernel["prime"]),
    "inputs": {
        "smooth_kernel": path_label(args.kernel),
        "all_component_conditions": path_label(args.conditions),
    },
    "components_tested": [entry["component"] for entry in conditions["component_conditions"]],
    "smooth_kernel_dimension": int(smooth.nrows()),
    "independent_unique_live_leading_constraints": records,
    "survivor_dimension_after_all_component_generic_conditions": int(candidate.nrows()),
    "boundary": (
        "This is a necessary all-component generic E7 screen modulo one prime. "
        "It leaves every group with multiple live coefficients unresolved and "
        "does not test chart-node membership, the marked branch, overlaps, "
        "characteristic-zero lifting, a pencil, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8ALLCOMPONENTGENERIC|prime={}|smooth_kernel={}|constraints={}|survivor={}|"
    "status=EXPERIMENTAL_MODULAR_ALL_COMPONENT_GENERIC_MODULE_SCREEN".format(
        kernel["prime"], smooth.nrows(), len(records), candidate.nrows()
    ),
    flush=True,
)
