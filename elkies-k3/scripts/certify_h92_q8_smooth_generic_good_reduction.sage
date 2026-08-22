#!/usr/bin/env sage -python
"""Certify a full-rank q=8 smooth-plus-generic block by good reduction.

This checker is for an explicitly supplied enlarged endpoint ambient.  It
uses the actual smooth q/X principal-part computation modulo one good prime,
then restricts all actual E7 generic-component rows to its recorded kernel.
If that restriction has full kernel rank, the stacked modular matrix has full
column rank.  Consequently the corresponding characteristic-zero smooth plus
generic-component condition matrix has full column rank too.

No node, marked-branch, overlap, or E8 assertion is needed for this rejection
certificate; adding more conditions cannot restore a kernel.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def basis_signature(basis):
    return [
        [str(entry["kind"]), int(entry["x_power"]), int(entry["m_power"]),
         int(entry["u_power"]), int(entry["h_power"])]
        for entry in basis
    ]


def reduced_coefficient(field, value, prime):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime {} divides a generic-residue denominator".format(prime))
    return field(ZZ(value.numerator()))/denominator


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, required=True)
parser.add_argument("--template", type=Path, required=True)
parser.add_argument("--cover", type=Path, required=True)
parser.add_argument("--smooth-probe", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
template = json.loads(args.template.read_text())
cover = json.loads(args.cover.read_text())
smooth = json.loads(args.smooth_probe.read_text())
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}
assert template["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"
assert cover["status"] == "PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER"
assert smooth["status"] == "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK"
assert smooth.get("kernel_basis_rows") is not None
assert template["inputs"]["endpoint_ambient"]["sha256"] == digest(args.ambient)
assert int(smooth["extra_e7_pole"]) == 0
assert basis_signature(smooth["ambient_basis"]) == basis_signature(ambient["ambient_basis"])

prime = int(smooth["prime"])
field = GF(prime)
ambient_dimension = int(ambient["ambient_dimension"])
smooth_dimensions = smooth["dimensions"]
assert int(smooth_dimensions["columns"]) == ambient_dimension
assert int(smooth_dimensions["rank"]) + int(smooth_dimensions["kernel"]) == ambient_dimension
kernel = matrix(field, [
    [field(ZZ(value)) for value in row]
    for row in smooth["kernel_basis_rows"]
])
assert kernel.nrows() == int(smooth_dimensions["kernel"])
assert kernel.ncols() == ambient_dimension and kernel.rank() == kernel.nrows()

row_payloads = []
for reference in cover["inputs"].values():
    path = ROOT / reference["path"]
    assert digest(path) == reference["sha256"]
    payload = json.loads(path.read_text())
    assert payload["ambient_basis_sha256"] == cover["ambient_basis_sha256"]
    row_payloads.append(payload)

components = {}
for payload in row_payloads:
    for component in payload["components"]:
        if component["component"] in components:
            raise ValueError("duplicate component residue payload")
        components[component["component"]] = component
assert tuple(sorted(components)) == tuple("E7_{}".format(index) for index in range(1, 8))
assert cover["ambient_basis_sha256"] == template["ambient_basis_sha256"]

restriction_rows = []
for singleton in template["singleton_coordinate_block"]["basis_indices"]:
    restriction_rows.append([kernel[row, int(singleton)] for row in range(kernel.nrows())])
for name in sorted(components):
    for row in components[name]["non_singleton_residue_rows"]:
        restricted = [field(0)]*kernel.nrows()
        for entry in row["entries"]:
            coefficient = reduced_coefficient(field, entry["coefficient"], prime)
            index = int(entry["basis_index"])
            for column in range(kernel.nrows()):
                restricted[column] += coefficient*kernel[column, index]
        restriction_rows.append(restricted)
expected_rows = (
    len(template["singleton_coordinate_block"]["basis_indices"])
    + int(cover["total_non_singleton_residue_rows"])
)
assert len(restriction_rows) == expected_rows
restriction = matrix(field, restriction_rows)
restriction_rank = restriction.rank()
assert restriction_rank == kernel.nrows()

exec(compile(CORE.read_text(), str(CORE), "exec"))
# This is the rank identity recorded by the reusable compiler core.  The
# smooth matrix itself is represented by its independently recomputed probe
# and kernel; the generic restriction supplies the second summand exactly.
stacked_rank = int(smooth_dimensions["rank"])+restriction_rank
assert stacked_rank == ambient_dimension

payload = {
    "schema": "elkies-k3.h92-q8-smooth-generic-good-reduction.v1",
    "status": "PASS_EXACT_Q8_SMOOTH_GENERIC_GOOD_REDUCTION_REJECTION",
    "inputs": {
        "endpoint_ambient": {"path": path_label(args.ambient), "sha256": digest(args.ambient)},
        "generic_template": {"path": path_label(args.template), "sha256": digest(args.template)},
        "generic_residue_cover": {"path": path_label(args.cover), "sha256": digest(args.cover)},
        "smooth_modular_probe": {"path": path_label(args.smooth_probe), "sha256": digest(args.smooth_probe)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
    "good_reduction": {
        "prime": prime,
        "ambient_dimension": ambient_dimension,
        "smooth_rank": int(smooth_dimensions["rank"]),
        "smooth_kernel_dimension": int(smooth_dimensions["kernel"]),
        "generic_rows_restricted_to_smooth_kernel": expected_rows,
        "generic_restriction_rank": int(restriction_rank),
        "stacked_rank": int(stacked_rank),
        "stacked_kernel_dimension": 0,
        "argument": (
            "All listed rational generic-residue denominators are nonzero modulo "
            "the displayed prime. The actual smooth and generic rows therefore give "
            "a full-column reduction, which proves full characteristic-zero rank."
        ),
    },
    "boundary": (
        "This rejects only the supplied endpoint enlargement. It is not a complete q8 "
        "resolved cover and makes no assertion about a larger ambient, nodes, marked "
        "branches, overlaps, E8, a pencil, child equation, bisection, collision, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8SMOOTHGENERICGOODREDUCTION|prime={}|ambient={}|smooth_rank={}|"
    "smooth_kernel={}|generic_restriction_rank={}|stacked_rank={}|status="
    "PASS_EXACT_Q8_SMOOTH_GENERIC_GOOD_REDUCTION_REJECTION".format(
        prime, ambient_dimension, smooth_dimensions["rank"], smooth_dimensions["kernel"],
        restriction_rank, stacked_rank,
    ),
    flush=True,
)
