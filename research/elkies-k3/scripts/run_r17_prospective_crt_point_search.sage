#!/usr/bin/env sage-python
"""Run or merge the frozen identical bounded R17 CRT point searches.

The supervisor is checkpointed by deterministic chunks.  Each candidate runs
in its own process under the frozen 300-second/8-GB envelope.  Timeouts and
backend failures are retained separately from completed bounded misses.  A
point counts only after exact equation verification and a full finite-quotient
independence certificate for MW17 plus the counted extra directions.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import resource
import runpy
import subprocess
import sys
import time
from typing import Any

from sage.all import QQ, ZZ
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"
LOCAL_IMPLEMENTATION = ROOT / "elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage"
DEFAULT_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/r17-prospective-crt-point-search-v1"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-point-search-ledger-v1.json"

SCHEMA = "elkies-k3.r17-prospective-crt-point-search-ledger.v1"
CHUNK_SCHEMA = "elkies-k3.r17-prospective-crt-point-search-chunk.v1"
EXPECTED_CANDIDATE_LIST_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
HEIGHT = 12
MAXR = 32
PP = 0
TIMEOUT_SECONDS = 300
MEMORY_BYTES = 8_000_000_000
RETRIES = 0
CERTIFICATE_PRIME_BOUND = 1000

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from elkies_residual_selmer_gate import monotone_sieve_gate_record  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def load_manifest():
    manifest = json.loads(MANIFEST.read_text())
    if manifest.get("status") != "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS":
        raise ArithmeticError("the prospective cohort is not frozen and unopened")
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_LIST_HASH:
        raise ArithmeticError("the reviewed candidate list changed")
    protocol = manifest["frozen_search_protocol"]
    expected = {
        "search_height": HEIGHT,
        "maxr": MAXR,
        "pp": PP,
        "wall_clock_limit_seconds": TIMEOUT_SECONDS,
        "memory_limit_bytes": MEMORY_BYTES,
        "retries": RETRIES,
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            raise ArithmeticError(f"the frozen search setting {key} changed")
    return manifest


def mwrank_point(curve, triple):
    x_coordinate, y_coordinate, z_coordinate = [ZZ(value) for value in triple]
    if z_coordinate == 0:
        return curve(0)
    return curve([QQ(x_coordinate) / QQ(z_coordinate), QQ(y_coordinate) / QQ(z_coordinate)])


def point_key(point):
    if point.is_zero():
        return ("0", "0")
    x_coordinate, y_coordinate = map(QQ, point[:2])
    return (str(x_coordinate), min(str(y_coordinate), str(-y_coordinate)))


def python_point(point):
    x_coordinate, y_coordinate = map(QQ, point[:2])
    return Fraction(str(x_coordinate)), Fraction(str(y_coordinate))


def certificate_record(signatures, column_count):
    return {
        "prime_bound": CERTIFICATE_PRIME_BOUND,
        "combined_rank": combined_mod2_rank(signatures, column_count),
        "column_count": column_count,
        "signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
    }


def run_single(index: int):
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY_BYTES, MEMORY_BYTES))
    manifest = load_manifest()
    if not (0 <= index < len(manifest["rows"])):
        raise ValueError("candidate index is outside the frozen manifest")
    row = manifest["rows"][index]
    local_implementation = runpy.run_path(str(LOCAL_IMPLEMENTATION))
    family = local_implementation["Family"]()
    parameter = QQ(row["parameter"])
    started = time.monotonic()
    print(f"PHASE build_start sample={row['sample_id']}", flush=True)
    curve, known = family.specialize(parameter)
    search_curve = curve.local_data(2).minimal_model()
    isomorphisms = curve.isomorphisms(search_curve)
    if not isomorphisms:
        raise ArithmeticError("no exact isomorphism to the fixed p=2-minimal search model")
    to_search = isomorphisms[0]
    from_search = ~to_search
    search_ainvs = [ZZ(value) for value in search_curve.a_invariants()]
    if any(QQ(value).denominator() != 1 for value in search_curve.a_invariants()):
        raise ArithmeticError("the fixed search model is not integral")
    print(
        f"PHASE build_done model_bits={max(abs(value).nbits() for value in search_ainvs)}",
        flush=True,
    )
    gate = monotone_sieve_gate_record(
        stages=[
            {
                "stage": "presearch_local_fingerprint_only",
                "residual_upper_bound": None,
                "proof_status": "NO_FINITE_UPPER_BOUND_YET",
            }
        ],
        search_limits={
            "height": HEIGHT,
            "wall_seconds": TIMEOUT_SECONDS,
            "memory_bytes": MEMORY_BYTES,
        },
    )
    if gate["bounded_point_search_authorized"] is not True or gate["theorem_claim_authorized"] is not False:
        raise ArithmeticError("the monotone residual-Selmer gate did not authorize bounded search")

    print("PHASE eclib_init_start", flush=True)
    init_started = time.monotonic()
    mw_curve = mwrank_EllipticCurve(search_ainvs)
    mordell_weil = mwrank_MordellWeil(mw_curve, verbose=False, pp=PP, maxr=MAXR)
    init_seconds = time.monotonic() - init_started
    print(f"PHASE eclib_init_done seconds={init_seconds:.6f}", flush=True)
    print(f"PHASE search_start height={HEIGHT}", flush=True)
    search_started = time.monotonic()
    mordell_weil.search(HEIGHT, verbose=False)
    search_seconds = time.monotonic() - search_started
    raw_points = list(mordell_weil.points())
    print(f"PHASE search_done seconds={search_seconds:.6f} raw={len(raw_points)}", flush=True)

    candidates = {}
    for triple in raw_points:
        search_point = mwrank_point(search_curve, triple)
        if search_point.is_zero():
            continue
        point = from_search(search_point)
        if point.curve() != curve or point not in curve:
            raise ArithmeticError("an eclib point failed exact equation transport")
        candidates.setdefault(point_key(point), point)
    ordered_candidates = sorted(
        candidates.values(),
        key=lambda point: max(
            len(str(abs(QQ(point[0]).numerator()))),
            len(str(QQ(point[0]).denominator())),
            len(str(abs(QQ(point[1]).numerator()))),
            len(str(QQ(point[1]).denominator())),
        ),
    )

    coefficients = [Fraction(0), Fraction(0), Fraction(0), Fraction(str(curve.a4())), Fraction(str(curve.a6()))]
    known_python = [python_point(point) for point in known]
    selected_points = []
    selected_certificates = []
    uncertified_candidates = []
    torsion_certificate_prime = None
    for candidate in ordered_candidates:
        trial = known_python + [python_point(point) for point in selected_points + [candidate]]
        signatures = find_mod2_reduction_certificate(
            coefficients,
            trial,
            prime_bound=CERTIFICATE_PRIME_BOUND,
        )
        rank = combined_mod2_rank(signatures, len(trial))
        if rank == len(trial):
            if torsion_certificate_prime is None:
                torsion_certificate_prime = find_two_torsion_certificate_prime(
                    coefficients, prime_bound=200
                )
            selected_points.append(candidate)
            selected_certificates.append(certificate_record(signatures, len(trial)))
        else:
            uncertified_candidates.append(
                {
                    "point": [str(candidate[0]), str(candidate[1])],
                    "reason": "FINITE_MOD2_CERTIFICATE_DID_NOT_REACH_FULL_COLUMN_RANK",
                    "achieved_rank": rank,
                    "column_count": len(trial),
                }
            )

    extra_count = len(selected_points)
    if extra_count:
        outcome = "CERTIFIED_MW17_ESCAPE"
    elif ordered_candidates:
        outcome = "COMPLETED_UNCERTIFIED_CANDIDATES_NO_COUNTED_ESCAPE"
    else:
        outcome = "BOUNDED_PROTOCOL_NO_ESCAPE_FOUND"
    result = {
        "sample_id": row["sample_id"],
        "manifest_index": index,
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": row["anchor_curve_id"],
        "cohort": row["cohort"],
        "parameter": row["parameter"],
        "status": outcome,
        "failure": None,
        "search_model": {
            "normalization": "exact integral p=2-minimal model returned by the same deterministic Sage local_data path",
            "a_invariants": [str(value) for value in search_ainvs],
            "model_sha256": canonical_hash([str(value) for value in search_ainvs]),
        },
        "limits": {
            "height": HEIGHT,
            "wall_seconds": TIMEOUT_SECONDS,
            "memory_bytes": MEMORY_BYTES,
            "maxr": MAXR,
            "pp": PP,
            "retries": RETRIES,
        },
        "monotone_residual_selmer_gate": gate,
        "timing": {
            "eclib_init_seconds": init_seconds,
            "search_seconds": search_seconds,
            "worker_wall_seconds": time.monotonic() - started,
            "time_to_first_escape_seconds": "NOT_AVAILABLE_FROM_BATCH_ECLIB_SEARCH",
        },
        "raw_mwrank_point_count": len(raw_points),
        "distinct_exact_candidate_count": len(ordered_candidates),
        "certified_independent_extra_directions": extra_count,
        "largest_certified_rank_lower_bound": 17 + extra_count,
        "search_effort_seconds_per_certified_direction": (
            search_seconds / extra_count if extra_count else None
        ),
        "certified_points": [
            {"point": [str(point[0]), str(point[1])], "finite_quotient_certificate": certificate}
            for point, certificate in zip(selected_points, selected_certificates)
        ],
        "rational_two_torsion_exclusion_prime": torsion_certificate_prime,
        "uncertified_candidates": uncertified_candidates,
        "claim_boundary": (
            "A completed miss is bounded-protocol-only. Only directions with a full exact "
            "finite-quotient certificate are counted. No Selmer or rank upper bound follows."
        ),
    }
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def parse_child_result(completed, row, index, elapsed):
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    result_line = next((line for line in reversed(lines) if line.startswith("RESULT_JSON=")), None)
    if completed.returncode == 0 and result_line is not None:
        result = json.loads(result_line[len("RESULT_JSON=") :])
        result["supervisor_wall_seconds"] = elapsed
        return result
    return {
        "sample_id": row["sample_id"],
        "manifest_index": index,
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": row["anchor_curve_id"],
        "cohort": row["cohort"],
        "parameter": row["parameter"],
        "status": "CENSORED_BACKEND_FAILURE",
        "failure": {
            "returncode": completed.returncode,
            "output_tail": lines[-20:],
        },
        "limits": {
            "height": HEIGHT,
            "wall_seconds": TIMEOUT_SECONDS,
            "memory_bytes": MEMORY_BYTES,
            "maxr": MAXR,
            "pp": PP,
            "retries": RETRIES,
        },
        "supervisor_wall_seconds": elapsed,
    }


def timeout_result(exc, row, index, elapsed):
    output = exc.stdout or ""
    if isinstance(output, bytes):
        output = output.decode(errors="replace")
    lines = [line for line in output.splitlines() if line.strip()]
    last_phase = next((line for line in reversed(lines) if line.startswith("PHASE ")), None)
    return {
        "sample_id": row["sample_id"],
        "manifest_index": index,
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": row["anchor_curve_id"],
        "cohort": row["cohort"],
        "parameter": row["parameter"],
        "status": "CENSORED_TIMEOUT",
        "failure": {
            "timeout_seconds": TIMEOUT_SECONDS,
            "last_phase": last_phase,
            "output_tail": lines[-20:],
        },
        "limits": {
            "height": HEIGHT,
            "wall_seconds": TIMEOUT_SECONDS,
            "memory_bytes": MEMORY_BYTES,
            "maxr": MAXR,
            "pp": PP,
            "retries": RETRIES,
        },
        "supervisor_wall_seconds": elapsed,
    }


def write_checkpoint(path: Path, chunk_index: int, chunk_count: int, selected_indices, records):
    document = {
        "schema": CHUNK_SCHEMA,
        "status": (
            "COMPLETE_POINT_SEARCH_CHUNK"
            if len(records) == len(selected_indices)
            else "PARTIAL_CHECKPOINT_POINT_SEARCH_CHUNK"
        ),
        "chunk_index": chunk_index,
        "chunk_count": chunk_count,
        "candidate_list_sha256": EXPECTED_CANDIDATE_LIST_HASH,
        "scheduled_indices": selected_indices,
        "completed_record_count": len(records),
        "records": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run_chunk(chunk_index: int, chunk_count: int, output: Path, limit: int | None):
    manifest = load_manifest()
    selected_indices = [index for index in range(len(manifest["rows"])) if index % chunk_count == chunk_index]
    if limit is not None:
        selected_indices = selected_indices[:limit]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        if (
            old.get("schema") != CHUNK_SCHEMA
            or old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("scheduled_indices") != selected_indices
        ):
            raise ArithmeticError("an existing point-search checkpoint belongs to another chunk")
        records = old["records"]
    completed_ids = {row["sample_id"] for row in records}
    for position, index in enumerate(selected_indices, start=1):
        row = manifest["rows"][index]
        if row["sample_id"] in completed_ids:
            continue
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single-index",
            str(index),
        ]
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=TIMEOUT_SECONDS,
                check=False,
            )
            result = parse_child_result(completed, row, index, time.monotonic() - started)
        except subprocess.TimeoutExpired as exc:
            result = timeout_result(exc, row, index, time.monotonic() - started)
        records.append(result)
        write_checkpoint(output, chunk_index, chunk_count, selected_indices, records)
        print(
            f"R17CRTPOINTS|chunk={chunk_index}/{chunk_count}|"
            f"completed={position}/{len(selected_indices)}|status={result['status']}",
            flush=True,
        )
    write_checkpoint(output, chunk_index, chunk_count, selected_indices, records)


def merge_chunks(chunk_dir: Path, chunk_count: int, output: Path):
    manifest = load_manifest()
    chunks = []
    by_id = {}
    for index in range(chunk_count):
        path = chunk_dir / f"chunk-{index:02d}-of-{chunk_count:02d}.json"
        chunk = json.loads(path.read_text())
        if chunk.get("status") != "COMPLETE_POINT_SEARCH_CHUNK":
            raise ArithmeticError(f"point-search chunk {index} is incomplete")
        chunks.append((path, chunk))
        for row in chunk["records"]:
            if row["sample_id"] in by_id:
                raise ArithmeticError("duplicate point-search sample across chunks")
            by_id[row["sample_id"]] = row
    expected_ids = [row["sample_id"] for row in manifest["rows"]]
    if set(by_id) != set(expected_ids):
        raise ArithmeticError("point-search chunks do not cover the frozen manifest")
    records = [by_id[sample_id] for sample_id in expected_ids]
    status_counts = Counter(row["status"] for row in records)
    cohort_counts = defaultdict(Counter)
    for row in records:
        cohort_counts[row["cohort"]][row["status"]] += 1
    document = {
        "schema": SCHEMA,
        "status": "COMPLETE_FROZEN_BOUNDED_POINT_SEARCH_LEDGER",
        "candidate_list_sha256": EXPECTED_CANDIDATE_LIST_HASH,
        "feature_artifact": {"path": relative(FEATURES), "sha256": digest(FEATURES)},
        "summary": {
            "scheduled_candidates": len(records),
            "status_counts": dict(sorted(status_counts.items())),
            "status_counts_by_cohort": {
                cohort: dict(sorted(counts.items())) for cohort, counts in sorted(cohort_counts.items())
            },
            "certified_escape_rows": sum(row["status"] == "CERTIFIED_MW17_ESCAPE" for row in records),
            "certified_extra_directions": sum(row.get("certified_independent_extra_directions", 0) for row in records),
        },
        "protocol": {
            "engine": "Sage/eclib mwrank_MordellWeil.search",
            "search_model": "deterministic exact integral p=2-minimal model",
            "height": HEIGHT,
            "maxr": MAXR,
            "pp": PP,
            "wall_seconds_including_setup": TIMEOUT_SECONDS,
            "memory_bytes": MEMORY_BYTES,
            "retries": RETRIES,
            "finite_quotient_certificate_prime_bound": CERTIFICATE_PRIME_BOUND,
        },
        "records": records,
        "chunk_provenance": [
            {"path": relative(path), "sha256": digest(path), "records": chunk["completed_record_count"]}
            for path, chunk in chunks
        ],
        "inputs": {relative(MANIFEST): digest(MANIFEST), relative(FEATURES): digest(FEATURES)},
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "commands": [
                f"sage -python elkies-k3/scripts/run_r17_prospective_crt_point_search.sage --chunk-index I --chunk-count {chunk_count}",
                f"sage -python elkies-k3/scripts/run_r17_prospective_crt_point_search.sage --merge --chunk-count {chunk_count}",
            ],
        },
        "claim_boundary": [
            "Only full finite-quotient certificates count as MW17 escapes.",
            "A completed bounded miss is not a rank-17 or exact-rank result.",
            "A timeout or backend failure is censored and is not a bounded miss.",
            "The open monotone sieve authorizes only this bounded search and no theorem claim.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"R17CRTPOINTS|records={len(records)}|status=COMPLETE|output={relative(output)}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=32)
    parser.add_argument("--chunk-dir", type=Path, default=DEFAULT_CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    if args.single_index is not None:
        run_single(args.single_index)
        return
    if args.merge:
        merge_chunks(args.chunk_dir.resolve(), args.chunk_count, (args.output or OUTPUT).resolve())
        return
    if args.chunk_index is None:
        raise SystemExit("chunk mode requires --chunk-index")
    output = args.output or (
        args.chunk_dir / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
    )
    run_chunk(args.chunk_index, args.chunk_count, output.resolve(), args.limit)


if __name__ == "__main__":
    main()
