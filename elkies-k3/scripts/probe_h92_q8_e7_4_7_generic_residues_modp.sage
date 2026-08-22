#!/usr/bin/env sage -python
"""Apply exact E7_4/E7_7 generic residue rows to a q=8 modular kernel."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
KERNEL = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json"
RESIDUES = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residues-mod-43-extra10.json"


def signature(basis):
    return [
        [entry["kind"], int(entry["x_power"]), int(entry["m_power"]),
         int(entry["u_power"]), int(entry["h_power"])]
        for entry in basis
    ]


def canonical_digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel", type=Path, default=KERNEL)
parser.add_argument("--residues", type=Path, default=RESIDUES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

kernel = json.loads(args.kernel.read_text())
residues = json.loads(args.residues.read_text())
assert kernel["status"] == "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK"
assert residues["status"] == "PASS_EXACT_Q8_E7_4_7_GENERIC_RESIDUE_ROWS"
assert canonical_digest(signature(kernel["ambient_basis"])) == residues["ambient_basis_sha256"]

finite = GF(int(kernel["prime"]))
smooth = matrix(
    finite,
    [[finite(int(value)) for value in row] for row in kernel["kernel_basis_rows"]],
)
assert smooth.nrows() == int(kernel["dimensions"]["kernel"])

ambient_dimension = smooth.ncols()
candidate = smooth
records = []
for component in residues["components"]:
    for residue_row in component["non_singleton_residue_rows"]:
        row = vector(finite, ambient_dimension)
        for entry in residue_row["entries"]:
            coefficient = QQ(entry["coefficient"])
            numerator = int(coefficient.numerator())
            denominator = int(coefficient.denominator())
            assert denominator % int(kernel["prime"])
            row[int(entry["basis_index"])] = finite(numerator)/finite(denominator)
        functional = vector(finite, [
            sum(candidate[index, column] * row[column] for column in range(ambient_dimension))
            for index in range(candidate.nrows())
        ])
        if not functional:
            continue
        before = candidate.nrows()
        candidate = matrix(
            finite,
            matrix(finite, 1, before, [functional]).right_kernel().basis_matrix(),
        ) * candidate
        records.append({
            "component": component["component"],
            "residual_order": int(residue_row["residual_order"]),
            "component_parameter_power": int(residue_row["component_parameter_power"]),
            "candidate_dimension_before": int(before),
            "entry_count": len(residue_row["entries"]),
        })
        if not candidate.nrows():
            break
    if not candidate.nrows():
        break

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-7-generic-residues-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_E7_4_7_GENERIC_RESIDUE_SCREEN",
    "prime": int(kernel["prime"]),
    "inputs": {"smooth_kernel": path_label(args.kernel), "residue_rows": path_label(args.residues)},
    "smooth_kernel_dimension": int(smooth.nrows()),
    "independent_residue_constraints": records,
    "survivor_dimension_after_E7_4_E7_7_generic_residues": int(candidate.nrows()),
    "boundary": (
        "This applies exact characteristic-zero residue rows after reduction "
        "modulo one prime on E7_4 and E7_7 only. It does not evaluate the other "
        "components, nodes, marked branch, overlaps, a characteristic-zero "
        "kernel, a pencil, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8E747RESIDUESMODP|prime={}|smooth_kernel={}|constraints={}|survivor={}|"
    "status=EXPERIMENTAL_MODULAR_E7_4_7_GENERIC_RESIDUE_SCREEN".format(
        kernel["prime"], smooth.nrows(), len(records), candidate.nrows()
    ),
    flush=True,
)
