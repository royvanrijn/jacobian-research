#!/usr/bin/env sage-python
"""Run or merge the target-free A1/MW16 parameter experiment."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import shutil
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/a1_mw16_blind_parameter_experiment_v1.json"
KEY = ROOT / "artifacts/local/elliptic-curves/a1-mw16-parameter-experiment/unblinding-v1.json"
CHECKPOINTS = ROOT / "artifacts/local/elliptic-curves/a1-mw16-parameter-experiment/rows"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/a1_mw16_blind_parameter_experiment_v1.json"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
GENERIC_DIMENSION = 16
OPERATIVE_SCALE = 1_000_000
AUDIT_SCALE = 100_000


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def generic_deepest(rows, legacy):
    gram = tuple(tuple(Fraction(value) for value in row) for row in rows)
    twice_gram = tuple(tuple(int(2*value) for value in row) for row in gram)
    if any(Fraction(twice_gram[i][j], 2) != gram[i][j] for i in range(16) for j in range(16)):
        raise ArithmeticError("generic MW16 Gram is not half-integral")
    oracle = legacy.CosetOracle(twice_gram)
    complete = []
    histogram = Counter()
    for mask in range(1 << 16):
        residue = tuple((mask >> index) & 1 for index in range(16))
        norm, representative, error = oracle.solve(residue)
        if error > 1.0e-6:
            raise ArithmeticError("generic half-lattice CVP failed")
        histogram[norm] += 1
        complete.append((norm, mask, representative))
    maximum = max(row[0] for row in complete)
    deepest = [row for row in complete if row[0] == maximum]
    if maximum != 23 or len(deepest) != 12:
        raise ArithmeticError("expected twelve twice-norm-23 classes")
    return deepest, histogram


def run_row(row, deepest, histogram, legacy, args):
    model = tuple(Fraction(value) for value in row["short_model"])
    generic = tuple((Fraction(point["x"]), Fraction(point["y"])) for point in row["specialized_generic_points"])
    if len(generic) != 16:
        raise ArithmeticError("blind row lost its MW16 basis")
    signatures = legacy.find_mod2_reduction_certificate(model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND)
    if legacy.combined_mod2_rank(signatures, 16) != 16:
        raise ArithmeticError("specialized MW16 is not certified independent")
    specialized_gram, asymmetry = legacy.canonical_height_gram(model, generic)
    rankings = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = legacy.CosetOracle(legacy.rounded_gram(specialized_gram, scale))
        ranked = []
        maximum_error = 0.0
        for unused_norm, mask, generic_representative in deepest:
            residue = tuple((mask >> index) & 1 for index in range(16))
            unused, representative, error = oracle.solve(residue)
            maximum_error = max(maximum_error, error)
            depth = legacy.quadratic_decimal(specialized_gram, representative) / 4
            ranked.append((depth, mask, representative, generic_representative))
        ranked.sort(key=lambda item: (-item[0], item[1]))
        rankings[scale] = (ranked, maximum_error)
    operative = rankings[OPERATIVE_SCALE][0]
    audit = rankings[AUDIT_SCALE][0]
    audit_map = {item[1]: item[2] for item in audit}

    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    discoveries = {}
    covers = []
    for priority, (depth, mask, representative, generic_representative) in enumerate(operative, 1):
        outcome = legacy.engine.run_quartic_search(
            mask=mask,
            representative=representative,
            short_model=model,
            generic_points=generic,
            height_bound=args.height_bound,
            timeout_seconds=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        source = f"initial:priority:{priority}:mask:{mask:#06x}"
        for point in outcome.curve_points:
            discoveries.setdefault(legacy.canonical_point(point), set()).add(source)
        covers.append({
            "priority": priority,
            "mask": mask,
            "generic_representative": list(generic_representative),
            "specialized_representative": list(representative),
            "specialized_depth": str(depth),
            "search": outcome.record,
        })
        print(f"A1MW16BLIND|row={row['row_id']}|chart={priority}/12|status={outcome.record['status']}|points={len(outcome.curve_points)}", flush=True)

    basis, classification = legacy.classify_discovered_group(
        model=model,
        basis=generic,
        discoveries=discoveries,
        relation_chunk_size=args.relation_chunk_size,
        relation_timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("discovered group did not classify exactly")
    quotient_rank = len(basis) - 16
    return {
        "row_id": row["row_id"],
        "status": "PASS_COMPLETE_DEEPEST_MW16_WAVE",
        "generic_rank": 16,
        "generic_mod2_independence_rank": 16,
        "exact_quotient_rank_gain": quotient_rank,
        "certified_rank_lower_bound": 16 + quotient_rank,
        "generic_half_lattice": {
            "complete_class_count": 1 << 16,
            "twice_norm_histogram": {str(key): value for key, value in sorted(histogram.items())},
            "maximum_twice_norm": 23,
            "deepest_class_count": 12,
            "deepest_masks": [item[1] for item in deepest],
        },
        "specialized_ranking": {
            "canonical_height_maximum_asymmetry": str(asymmetry),
            "operative_scale": OPERATIVE_SCALE,
            "audit_scale": AUDIT_SCALE,
            "priority_order_identical": [item[1] for item in operative] == [item[1] for item in audit],
            "representative_disagreement_count": sum(audit_map[item[1]] != item[2] for item in operative),
            "maximum_cvp_distance_error": {str(scale): error for scale, (unused, error) in rankings.items()},
        },
        "cover_records": covers,
        "discoveries": legacy.discovery_records(discoveries),
        "discovered_group_saturation": classification,
        "basis_rank_after": len(basis),
        "wall_seconds": time.monotonic() - started_wall,
        "cpu_seconds": cpu_clock() - started_cpu,
    }


def merge(blind, args):
    records = []
    missing = []
    for row in blind["rows"]:
        path = CHECKPOINTS / f"{row['row_id']}.json"
        if not path.is_file():
            missing.append(row["row_id"])
        else:
            records.append(json.loads(path.read_text()))
    key = json.loads(KEY.read_text())
    parameters = {row["row_id"]: row["parameter"] for row in key["rows"]}
    for record in records:
        record["parameter"] = parameters[record["row_id"]]
    counts = Counter(record["exact_quotient_rank_gain"] for record in records)
    payload = {
        "schema": "elliptic-curves.a1-mw16-blind-parameter-experiment-result.v1",
        "status": "PASS_COMPLETE_TARGET_FREE_PARAMETER_EXPERIMENT" if not missing else "PARTIAL_TARGET_FREE_PARAMETER_EXPERIMENT",
        "family_id": blind["family_id"],
        "scheduled_rows": blind["row_count"],
        "completed_rows": len(records),
        "missing_row_ids": missing,
        "declared_search": {
            "deepest_classes_each_row": 12,
            "height_bound": args.height_bound,
            "timeout_seconds_each_chart": args.timeout_seconds,
            "stack_bytes_each_chart": args.stack_bytes,
        },
        "gain_histogram": {str(key): value for key, value in sorted(counts.items())},
        "maximum_exact_quotient_rank_gain": max((record["exact_quotient_rank_gain"] for record in records), default=None),
        "rows": sorted(records, key=lambda record: record["row_id"]),
        "inputs": {str(path.relative_to(ROOT)): digest(path) for path in (INPUT, Path(__file__))},
        "claim_boundary": (
            "Every nonzero gain is an exact rank lower bound beyond specialized MW16. "
            "A zero is only a bounded deepest-wave miss and gives no rank, Selmer, "
            "saturation, or point-absence upper bound. No known-record target was used."
        ),
    }
    write_json(OUTPUT, payload)
    print(f"A1MW16MERGE|completed={len(records)}/{blind['row_count']}|histogram={dict(sorted(counts.items()))}|status={payload['status']}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=5.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    args = parser.parse_args()
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")
    blind = json.loads(INPUT.read_text())
    if blind.get("status") != "PASS_TARGET_FREE_BLIND_INPUT" or not all(not value for value in blind["blindness"].values()):
        raise ArithmeticError("worker input is not target-free")
    if args.merge:
        merge(blind, args)
        return
    if not 0 <= args.shard_index < args.shard_count:
        raise SystemExit("invalid shard index/count")
    legacy = SourceFileLoader(f"a1_mw16_legacy_{args.shard_index}", str(LEGACY)).load_module()
    legacy.GENERIC_DIMENSION = 16
    deepest, histogram = generic_deepest(blind["rows"][0]["generic_height_gram"], legacy)
    selected = [row for index, row in enumerate(blind["rows"]) if index % args.shard_count == args.shard_index]
    if args.limit is not None:
        selected = selected[:args.limit]
    for row in selected:
        checkpoint = CHECKPOINTS / f"{row['row_id']}.json"
        if checkpoint.is_file() and not args.force:
            print(f"A1MW16BLIND|row={row['row_id']}|status=SKIP_CHECKPOINT", flush=True)
            continue
        result = run_row(row, deepest, histogram, legacy, args)
        write_json(checkpoint, result)
        print(f"A1MW16BLIND|row={row['row_id']}|gain={result['exact_quotient_rank_gain']}|status={result['status']}", flush=True)


if __name__ == "__main__":
    main()
