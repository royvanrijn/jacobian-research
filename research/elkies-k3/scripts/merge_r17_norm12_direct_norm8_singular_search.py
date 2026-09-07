#!/usr/bin/env python3
"""Merge the exhaustive norm-eight singular-pencil search shards."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_PRIORITY_TABLES = {
    "norm12-orbit-11952": (
        GENERATED
        / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
    ),
    "norm12-orbit-103b2": (
        GENERATED
        / "elkies-k3-r17-norm12-103b2-norm8-pencil-priority-v1.tsv"
    ),
}


def prefix(source_label: str) -> str:
    short_label = source_label.removeprefix("norm12-orbit-")
    return f"elkies-k3-r17-norm12-{short_label}-singular-bisection-search"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def default_inputs(source_label: str) -> list[Path]:
    file_prefix = prefix(source_label)
    initial = [
        GENERATED / f"{file_prefix}-l1le4-v1.json",
        GENERATED / f"{file_prefix}-l1eq5-v1.json",
        GENERATED / f"{file_prefix}-l1eq6-v1.json",
    ]
    return [path for path in initial if path.exists()] + sorted(
        GENERATED.glob(f"{file_prefix}-full-*-v1.json"),
        key=lambda path: int(path.name.split("-full-", 1)[1].split("-", 1)[0]),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, action="append")
parser.add_argument("--source-label", default="norm12-orbit-11952")
parser.add_argument("--priority-table", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

paths = (
    [path.resolve() for path in args.input]
    if args.input
    else default_inputs(args.source_label)
)
if not paths:
    raise ValueError("no singular-search shards found")
priority_table = args.priority_table or DEFAULT_PRIORITY_TABLES.get(args.source_label)
if priority_table is None:
    raise ValueError(
        f"no default priority table for source label {args.source_label}; "
        "pass --priority-table explicitly"
    )
priority_path = priority_table.resolve()
minimum_count_by_parity_mask = {}
with priority_path.open(newline="") as stream:
    for row in csv.DictReader(stream, delimiter="\t"):
        word = tuple(map(int, row["section_basis_w"].split()))
        if len(word) != 17:
            raise ValueError("priority-table section word does not have rank 17")
        parity_mask = sum((entry & 1) << index for index, entry in enumerate(word))
        count = int(row["minimal_unoriented_count"])
        if parity_mask in minimum_count_by_parity_mask:
            raise ValueError(f"duplicate priority-table parity mask: {parity_mask}")
        minimum_count_by_parity_mask[parity_mask] = count
if not minimum_count_by_parity_mask:
    raise ValueError("priority table contains no parity classes")
source_token = args.source_label.removeprefix("norm12-orbit-")
expected_model = (
    "artifacts/generated-results/"
    f"elkies-k3-r17-norm12-orbit{source_token}-direct-fibration-v1.json"
)
expected_model_path = ROOT / expected_model
if not expected_model_path.is_file():
    raise ValueError(f"missing direct model: {expected_model_path}")
expected_model_sha256 = digest(expected_model_path)
records = []
status_histogram = Counter()
minimum_count_histogram = Counter()
even_discriminant_degree_histogram = Counter()
seen_trace_parity_masks = set()
for path in paths:
    payload = json.loads(path.read_text())
    if payload.get("schema") != (
        "elkies-k3.r17-norm12-direct-singular-bisection-search.v1"
    ):
        raise ValueError(f"unexpected schema: {path}")
    inputs = payload.get("inputs", {})
    if expected_model not in inputs:
        raise ValueError(f"shard does not belong to {args.source_label}: {path}")
    if inputs[expected_model] != expected_model_sha256:
        raise ValueError(f"shard direct-model hash is stale: {path}")
    search = payload["search"]
    count = int(search["processed_trace_count"])
    interval = search.get("processed_half_open_range")
    if interval is None:
        if int(search.get("start", 0)) != 0:
            raise ValueError(f"cannot infer interval: {path}")
        interval = [0, count]
    interval = tuple(map(int, interval))
    if interval[1] - interval[0] != count:
        raise ValueError(f"interval/count mismatch: {path}")
    if len(payload["trace_records"]) != count:
        raise ValueError(f"trace-record/count mismatch: {path}")
    trace_indices = [int(record["trace_index"]) for record in payload["trace_records"]]
    legacy_local_indices = list(range(count))
    global_indices = list(range(interval[0], interval[1]))
    if trace_indices not in (legacy_local_indices, global_indices):
        raise ValueError(f"trace indices are neither local nor global sequential: {path}")
    if any(
        int(payload[key]) != 0
        for key in (
            "candidate_count",
            "candidate_collision_count",
            "smooth_atlas_match_count",
        )
    ):
        raise ValueError(f"nonempty candidate output: {path}")
    for trace_record in payload["trace_records"]:
        status_histogram[trace_record["status"]] += 1
        word = tuple(map(int, trace_record["basis_coordinates"]))
        if len(word) != 17:
            raise ValueError(f"trace word does not have rank 17: {path}")
        parity_mask = sum((entry & 1) << index for index, entry in enumerate(word))
        if parity_mask in seen_trace_parity_masks:
            raise ValueError(
                f"duplicate trace parity class across singular-search shards: "
                f"{parity_mask} in {path}"
            )
        seen_trace_parity_masks.add(parity_mask)
        try:
            minimum_count = minimum_count_by_parity_mask[parity_mask]
        except KeyError as error:
            raise ValueError(
                f"trace parity mask is absent from the priority table: {path}"
            ) from error
        declared_minimum_count = trace_record.get(
            "minimum_unoriented_split_member_count"
        )
        if (
            declared_minimum_count is not None
            and int(declared_minimum_count) != minimum_count
        ):
            raise ValueError(
                f"trace/priority minimum-count mismatch: {path}, "
                f"parity={parity_mask}"
            )
        discriminant_degree = int(trace_record["pencil_discriminant_degree"])
        odd_degree = trace_record.get(
            "odd_multiplicity_pencil_discriminant_degree"
        )
        if odd_degree is None:
            profile = trace_record.get("factor_degree_multiplicities")
            if profile is None:
                raise ValueError(f"cannot recover odd discriminant degree: {path}")
            odd_degree = sum(
                int(degree)
                for degree, multiplicity in profile
                if int(multiplicity) % 2
            )
        odd_degree = int(odd_degree)
        if (discriminant_degree - odd_degree) % 2:
            raise ValueError(f"odd square-factor degree: {path}")
        even_half_degree = (discriminant_degree - odd_degree) // 2
        # A parity class with m minimum representatives up to sign has m split
        # members S_x+S_(w-x).  The chosen chord gauge puts one at infinity,
        # where q_infinity=h^2.  The other m-1 finite square members each give
        # at least a double discriminant root.  Equality below proves that they
        # exhaust the entire even-multiplicity part, so no cusp or tangential
        # nonsplit rational member is hidden by taking the polynomial
        # squareclass of the pencil discriminant.
        if even_half_degree != minimum_count - 1:
            raise ValueError(
                "even discriminant part is not exhausted by split members: "
                f"{path}, parity={parity_mask}, half_degree={even_half_degree}, "
                f"minimum_count={minimum_count}"
            )
        minimum_count_histogram[minimum_count] += 1
        even_discriminant_degree_histogram[2 * even_half_degree] += 1
    records.append((interval, path, payload))

records.sort(key=lambda item: item[0])
cursor = 0
for interval, path, _payload in records:
    if interval[0] != cursor:
        raise ValueError(
            f"coverage gap or overlap before {relative(path)}: expected {cursor}, "
            f"found {interval[0]}"
        )
    cursor = interval[1]

expected = int(records[-1][2]["search"]["minimum_norm_eight_translation_classes"])
if len(minimum_count_by_parity_mask) != expected:
    raise ValueError(
        "priority-table class count does not match the singular-search shards: "
        f"{len(minimum_count_by_parity_mask)} != {expected}"
    )
if seen_trace_parity_masks != set(minimum_count_by_parity_mask):
    missing = set(minimum_count_by_parity_mask) - seen_trace_parity_masks
    extra = seen_trace_parity_masks - set(minimum_count_by_parity_mask)
    raise ValueError(
        "trace parity-class coverage differs from the priority table: "
        f"missing={len(missing)}, extra={len(extra)}"
    )
if cursor != expected:
    raise ValueError(f"incomplete coverage: reached {cursor} of {expected}")
if any(
    int(payload["search"]["minimum_norm_eight_translation_classes"]) != expected
    for _interval, _path, payload in records
):
    raise ValueError("inconsistent norm-eight class count")

allowed_statuses = {
    "PASS_EXACT_PENCIL_DISCRIMINANT_FACTORIZATION",
    "PASS_EXACT_ODD_DISCRIMINANT_IRREDUCIBLE",
    "PASS_MODULAR_NO_RATIONAL_SINGULAR_PARAMETER",
    "PASS_MODULAR_FULL_DISCRIMINANT_SPLIT_EXHAUSTION",
}
unexpected_statuses = set(status_histogram) - allowed_statuses
if unexpected_statuses:
    raise ValueError(f"unexpected trace statuses: {sorted(unexpected_statuses)}")

if args.output is None:
    output_path = GENERATED / f"{prefix(args.source_label)}-complete-v1.json"
else:
    output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
result = {
    "schema": "elkies-k3.r17-norm12-direct-norm8-singular-search-merge.v1",
    "status": "PASS_COMPLETE_NO_NONSPLIT_RATIONAL_SINGULAR_MEMBER",
    "source_label": args.source_label,
    "trace_norm": 8,
    "minimum_translation_class_count": expected,
    "covered_half_open_range": [0, cursor],
    "shard_count": len(records),
    "trace_status_histogram": dict(sorted(status_histogram.items())),
    "minimal_unoriented_split_member_count_histogram": dict(
        sorted(minimum_count_histogram.items())
    ),
    "finite_even_discriminant_degree_histogram": dict(
        sorted(even_discriminant_degree_histogram.items())
    ),
    "even_multiplicity_exhaustion": (
        "For every class with m minimum representatives up to sign, the finite "
        "even-multiplicity discriminant degree is exactly 2*(m-1). The m "
        "lattice representatives give m distinct split members; q_infinity=h^2 "
        "accounts for one, so the other m-1 finite split members exhaust the "
        "even part."
    ),
    "candidate_count": 0,
    "proof_boundary": (
        "All minimum norm-eight translation classes are covered. For every "
        "finite pencil, exact factorization or irreducibility, or a projective "
        "finite-field root obstruction, excludes an odd-multiplicity rational "
        "parameter yielding a nonsplit rational normalization. The exact "
        "minimum-representative count exhausts the even-multiplicity part by "
        "split members, including q_infinity=h^2. Thus no cuspidal or tangent "
        "nonsplit member is hidden in the discarded square factor. This does "
        "not address norm-four or norm-six trace systems."
    ),
    "inputs": {
        relative(path): digest(path) for _interval, path, _payload in records
    },
    "reproducing_command": (
        "python3 elkies-k3/scripts/merge_r17_norm12_direct_norm8_singular_search.py "
        f"--source-label {args.source_label} "
        f"--priority-table {relative(priority_path)} "
        + " ".join(f"--input {relative(path)}" for path in paths)
        + f" --output {relative(output_path)}"
    ),
}
result["inputs"][relative(priority_path)] = digest(priority_path)
serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
if args.check:
    if not output_path.is_file() or output_path.read_text() != serialized:
        raise ValueError(f"stored merged certificate differs from replay: {output_path}")
else:
    output_path.write_text(serialized)
print(
    "R17NORM8SINGULARMERGE"
    f"|classes={expected}|shards={len(records)}|candidates=0"
    f"|output={relative(output_path)}|status={result['status']}",
    flush=True,
)
