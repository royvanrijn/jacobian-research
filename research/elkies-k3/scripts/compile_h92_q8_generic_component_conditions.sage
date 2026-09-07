#!/usr/bin/env sage -python
"""Compile all exact generic-component q=8 E7 rows into one condition block.

The seven generic residue evaluators already turn every non-singleton
negative-order group into explicit QQ coefficient rows.  This script is their
matrix-level handoff to the elliptic-neighbour compiler.  It also includes the
22 singleton negative-order coordinates, which require no cancellation
calculation.  The result is deliberately an incomplete resolved cover: nodes,
marked branches, overlaps, E8, and smooth collision conditions remain
separate blocks.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, gcd, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
TEMPLATE = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions.json"
COVER = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-generic-e7-residue-cover.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-component-condition-block.json"


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
parser.add_argument("--cover", type=Path, default=COVER)
parser.add_argument(
    "--rank-prime", type=int,
    help="certify full column rank by good reduction at this prime instead of QQ elimination",
)
parser.add_argument(
    "--sparse-exact", action="store_true",
    help="compute rank from the row-normalized sparse integer matrix",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.rank_prime is not None and args.rank_prime <= 1:
    raise ValueError("rank-prime must be a prime number")
if args.rank_prime is not None and args.sparse_exact:
    raise ValueError("rank-prime and sparse-exact are mutually exclusive")

ambient = json.loads(args.ambient.read_text())
template = json.loads(args.template.read_text())
cover = json.loads(args.cover.read_text())
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}
assert template["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"
assert cover["status"] == "PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER"
assert template["inputs"]["endpoint_ambient"]["sha256"] == digest(args.ambient)
assert cover["ambient_basis_sha256"] == template["ambient_basis_sha256"]
if not cover.get("enlarged_ambient", False):
    assert int(ambient["ambient_dimension"]) == 54

# The cover is an index of the three independently derived resolved-chart
# row artifacts.  Read and authenticate each of them through that index.
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
        name = component["component"]
        if name in components:
            raise ValueError("duplicate generic-component residue payload for {}".format(name))
        components[name] = component
assert tuple(sorted(components)) == tuple("E7_{}".format(index) for index in range(1, 8))

exec(compile(CORE.read_text(), str(CORE), "exec"))
ambient_labels = tuple(range(int(ambient["ambient_dimension"])))
singleton_indices = tuple(int(value) for value in template["singleton_coordinate_block"]["basis_indices"])
assert singleton_indices == tuple(sorted(set(singleton_indices)))

residue_rows = []
component_counts = {}
for component_name in sorted(components):
    component = components[component_name]
    rows = component["non_singleton_residue_rows"]
    component_counts[component_name] = len(rows)
    for row_index, row in enumerate(rows):
        entries = {}
        for entry in row["entries"]:
            index = int(entry["basis_index"])
            if index in entries:
                raise ValueError("duplicate basis entry in {} row {}".format(component_name, row_index))
            entries[index] = str(entry["coefficient"])
        if not entries or not any(entries.values()):
            raise ValueError("zero generic residue row in {}".format(component_name))
        coordinate = (
            "parameter_power {}".format(int(row["component_parameter_power"]))
            if "component_parameter_power" in row
            else "monomial {}".format(row["component_monomial"])
        )
        residue_rows.append((
            component_name, int(row["residual_order"]), coordinate, entries,
        ))
assert component_counts == cover["component_row_counts"]
assert len(residue_rows) == int(cover["total_non_singleton_residue_rows"])
if not cover.get("enlarged_ambient", False):
    assert len(residue_rows) == 983

condition_rows = len(singleton_indices)+len(residue_rows)
rank_certificate = None
if args.sparse_exact:
    def normalized_integer_entries(entries):
        """Scale one rational row to a primitive signed integer row exactly."""
        rational_entries = {index: QQ(value) for index, value in entries.items()}
        denominator = ZZ(1)
        for value in rational_entries.values():
            denominator = lcm(denominator, ZZ(value.denominator()))
        integers = {
            index: ZZ(value*denominator)
            for index, value in rational_entries.items()
        }
        content = gcd([abs(value) for value in integers.values()])
        if not content:
            raise ValueError("zero generic residue row")
        integers = {index: value // content for index, value in integers.items()}
        first = integers[min(integers)]
        if first < 0:
            integers = {index: -value for index, value in integers.items()}
        return integers

    sparse_entries = {}
    for row_index, singleton in enumerate(singleton_indices):
        sparse_entries[(row_index, singleton)] = ZZ(1)
    for offset, (_, _, _, entries) in enumerate(residue_rows):
        row_index = len(singleton_indices)+offset
        for index, value in normalized_integer_entries(entries).items():
            sparse_entries[(row_index, index)] = value
    sparse_matrix = matrix(
        ZZ, condition_rows, len(ambient_labels), sparse_entries, sparse=True
    )
    rank = sparse_matrix.rank()
    kernel_dimension = len(ambient_labels)-rank
    matrix_sha256 = hashlib.sha256(json.dumps(
        [[row, column, str(value)] for (row, column), value in sorted(sparse_entries.items())],
        separators=(",", ":"),
    ).encode()).hexdigest()
    rank_certificate = {
        "method": "exact_sparse_integer_rank",
        "normalization": (
            "Each nonzero rational residue row is multiplied by the least common "
            "denominator and divided by its integer content; this preserves row span."
        ),
        "nonzero_entries": len(sparse_entries),
    }
elif args.rank_prime is None:
    singleton_block = quotient_condition(
        "all E7 generic singleton negative orders",
        ambient_labels,
        lambda index: tuple(QQ(1) if index == singleton else QQ(0) for singleton in singleton_indices),
        tuple("ambient coefficient {}".format(index) for index in singleton_indices),
        "actual generic-component valuation template; singleton negative orders",
    )
    residue_block = quotient_condition(
        "all E7 generic non-singleton residue cancellations",
        ambient_labels,
        lambda index: tuple(QQ(row[3][index]) if index in row[3] else QQ(0) for row in residue_rows),
        tuple(
            "{} order {} {}".format(component, order, coordinate)
            for component, order, coordinate, _ in residue_rows
        ),
        "actual all-component normalizations and their exact cleared residue rows",
    )
    compilation = compile_resolved_conditions(
        ambient_labels, (singleton_block, residue_block), complete=False, compute_kernel=False
    )
    assert compilation["condition_rows"] == condition_rows
    assert compilation["rank"] <= compilation["ambient_dimension"]
    assert compilation["kernel_basis"] is None
    rank = compilation["rank"]
    kernel_dimension = compilation["kernel_dimension"]
    matrix_sha256 = hashlib.sha256(str(compilation["condition_matrix"]).encode()).hexdigest()
else:
    finite = GF(args.rank_prime)
    rows = []
    for singleton in singleton_indices:
        rows.append([finite(1) if column == singleton else finite(0) for column in ambient_labels])
    for _, _, _, entries in residue_rows:
        row = [finite(0)]*len(ambient_labels)
        for index, value in entries.items():
            coefficient = QQ(value)
            denominator = finite(ZZ(coefficient.denominator()))
            if not denominator:
                raise ValueError(
                    "rank-prime {} divides a generic-residue denominator".format(args.rank_prime)
                )
            row[index] = finite(ZZ(coefficient.numerator()))/denominator
        rows.append(row)
    reduced_rank = matrix(finite, rows).rank()
    if reduced_rank != len(ambient_labels):
        raise ArithmeticError(
            "good reduction has rank {}, not full column rank {}".format(
                reduced_rank, len(ambient_labels)
            )
        )
    # A rational matrix cannot have rank above its column count, while its
    # full-column good reduction gives the opposite inequality exactly.
    rank = len(ambient_labels)
    kernel_dimension = 0
    matrix_sha256 = None
    rank_certificate = {
        "prime": args.rank_prime,
        "reduced_rank": int(reduced_rank),
        "argument": (
            "The exact matrix has at most its ambient column count as rank; "
            "full column rank after reduction at the displayed good prime proves "
            "the same characteristic-zero rank."
        ),
    }
payload = {
    "schema": "elkies-k3.h92-q8-generic-component-condition-block.v1",
    "status": "PASS_EXACT_Q8_GENERIC_COMPONENT_CONDITION_BLOCK",
    "inputs": {
        "endpoint_ambient": {"path": path_label(args.ambient), "sha256": digest(args.ambient)},
        "generic_condition_template": {"path": path_label(args.template), "sha256": digest(args.template)},
        "generic_residue_cover": {"path": path_label(args.cover), "sha256": digest(args.cover)},
        "compiler_core": {"path": path_label(CORE), "sha256": digest(CORE)},
    },
    "component_row_counts": component_counts,
    "condition_matrix": {
        "ambient_dimension": len(ambient_labels),
        "singleton_rows": len(singleton_indices),
        "non_singleton_residue_rows": len(residue_rows),
        "rows": condition_rows,
        "rank": rank,
        "codimension": rank,
        "kernel_dimension": kernel_dimension,
        "matrix_sha256": matrix_sha256,
        "rank_certificate": rank_certificate,
    },
    "boundary": (
        "This is the actual all-E7-component generic-point condition block only. "
        "It excludes nodes, marked smooth branches, chart overlaps, E8 and smooth "
        "conditions, so it cannot certify q8 h0, a pencil, or a child equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
        "H92Q8GENERICCOMPONENTBLOCK|ambient={}|singleton_rows={}|residue_rows={}|"
        "rank={}|kernel={}|status=PASS_EXACT_Q8_GENERIC_COMPONENT_CONDITION_BLOCK".format(
        len(ambient_labels), len(singleton_indices), len(residue_rows),
        rank, kernel_dimension,
    ),
    flush=True,
)
