#!/usr/bin/env sage-python
"""Run the frozen half-lattice detector without reading frozen Q_t values."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
COHORT = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-cohort-v1.json"
FEATURES = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-features-v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-detector-protocol-v1.json"
FEATURE_IMPLEMENTATION = ROOT / "elkies-k3/scripts/run_r17_small_field_class_quotient_features.sage"
BASE_RUNNER = ROOT / "elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage"
DEFAULT_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/r17-small-field-class-quotient-detector-v1"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-detector-ledger-v1.json"

PROTOCOL_STATUS = "FROZEN_AFTER_ALL_Q_BEFORE_ANY_POINT_SEARCH"
CHUNK_SCHEMA = "elkies-k3.r17-small-field-class-quotient-detector-chunk.v1"
SCHEMA = "elkies-k3.r17-small-field-class-quotient-detector-ledger.v1"

feature_implementation = SourceFileLoader(
    "small_field_feature_specialization", str(FEATURE_IMPLEMENTATION)
).load_module()
base_runner = SourceFileLoader("small_field_half_lattice_runner", str(BASE_RUNNER)).load_module()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def load_inputs():
    cohort = json.loads(COHORT.read_text())
    protocol = json.loads(PROTOCOL.read_text())
    if protocol.get("status") != PROTOCOL_STATUS:
        raise ArithmeticError("the post-feature detector protocol is absent or not frozen")
    if protocol.get("candidate_list_sha256") != cohort["commitment"]["candidate_list_sha256"]:
        raise ArithmeticError("the detector protocol names another cohort")
    protocol_body = {
        key: value
        for key, value in protocol.items()
        if key not in {"protocol_definition_sha256", "inputs", "generation"}
    }
    if canonical_hash(protocol_body) != protocol.get("protocol_definition_sha256"):
        raise ArithmeticError("the detector protocol definition hash does not replay")
    if digest(FEATURES) != protocol["phase_boundary"]["feature_artifact_sha256"]:
        raise ArithmeticError("the sealed feature bytes differ from the protocol commitment")
    if protocol["phase_boundary"].get("detector_loads_feature_values") is not False:
        raise ArithmeticError("the protocol no longer keeps Q_t values blind to the detector")
    if any(row.get("outcome_status") != "SEALED_UNTIL_ALL_FEATURES_FREEZE" for row in cohort["rows"]):
        raise ArithmeticError("the base cohort contains a prematurely opened outcome")
    declared_runner = protocol.get("inputs", {}).get(relative(Path(__file__)))
    if declared_runner != digest(Path(__file__)):
        raise ArithmeticError("this detector executable differs from the frozen source")
    if len(protocol["detector_manifest"]) != len(cohort["rows"]):
        raise ArithmeticError("the redacted detector manifest has the wrong row count")
    # Deliberately do not json.load(FEATURES): only its committed bytes are hashed.
    return cohort, protocol


class PublishedR17Family:
    def __init__(self):
        _cohort, self.model, self.sections = feature_implementation.load_inputs()

    def specialize(self, parameter):
        value = Fraction(str(parameter))
        row = {
            "projective_pair": [value.numerator, value.denominator],
        }
        _source, minimal, points, _isomorphism = (
            feature_implementation.specialized_curve_and_points(
                row, self.model, self.sections
            )
        )
        return minimal, points


def runner_row(protocol_row):
    return {
        "sample_id": protocol_row["sample_id"],
        "manifest_index": protocol_row["manifest_index"],
        "match_set_id": protocol_row["sample_id"],
        "anchor_curve_id": None,
        "cohort": protocol_row["family"],
        "parameter": protocol_row["parameter"],
    }


def run_single(index: int):
    _cohort, protocol = load_inputs()
    if not 0 <= index < len(protocol["detector_manifest"]):
        raise ValueError("detector index lies outside the frozen manifest")
    row = protocol["detector_manifest"][index]
    family = PublishedR17Family()
    result = base_runner.run_fibre(runner_row(row), protocol, family)
    result["laboratory_sample_id"] = row["sample_id"]
    result["feature_values_loaded_by_detector"] = False
    result["protocol_definition_sha256"] = protocol["protocol_definition_sha256"]
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def failure_record(row, status, failure, elapsed):
    return {
        "sample_id": row["sample_id"],
        "manifest_index": row["manifest_index"],
        "family": row["family"],
        "parameter": row["parameter"],
        "status": status,
        "failure": failure,
        "analysis_eligible_complete_stage_a": False,
        "feature_values_loaded_by_detector": False,
        "supervisor_wall_seconds": elapsed,
    }


def parse_worker(completed, row, elapsed):
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    result_line = next(
        (line for line in reversed(lines) if line.startswith("RESULT_JSON=")), None
    )
    if completed.returncode == 0 and result_line is not None:
        result = json.loads(result_line[len("RESULT_JSON=") :])
        result["sample_id"] = row["sample_id"]
        result["supervisor_wall_seconds"] = elapsed
        return result
    return failure_record(
        row,
        "CENSORED_FIBRE_WORKER_FAILURE",
        {"returncode": completed.returncode, "output_tail": lines[-40:]},
        elapsed,
    )


def write_chunk(path, chunk_index, chunk_count, indices, records, protocol):
    document = {
        "schema": CHUNK_SCHEMA,
        "status": (
            "COMPLETE_SCHEDULED_DETECTOR_CHUNK"
            if len(records) == len(indices)
            else "PARTIAL_DETECTOR_CHECKPOINT"
        ),
        "candidate_list_sha256": protocol["candidate_list_sha256"],
        "protocol_definition_sha256": protocol["protocol_definition_sha256"],
        "chunk_index": chunk_index,
        "chunk_count": chunk_count,
        "scheduled_indices": indices,
        "completed_record_count": len(records),
        "records": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run_chunk(chunk_index, chunk_count, output, limit):
    _cohort, protocol = load_inputs()
    indices = [
        index for index in range(len(protocol["detector_manifest"]))
        if index % chunk_count == chunk_index
    ]
    if limit is not None:
        indices = indices[:limit]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        if (
            old.get("protocol_definition_sha256") != protocol["protocol_definition_sha256"]
            or old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("scheduled_indices") != indices
        ):
            raise ArithmeticError("the existing detector checkpoint belongs to another schedule")
        records = old["records"]
    completed_ids = {record["sample_id"] for record in records}
    timeout = protocol["fibre_worker_envelope"]["wall_timeout_seconds"]
    for position, index in enumerate(indices, start=1):
        row = protocol["detector_manifest"][index]
        if row["sample_id"] in completed_ids:
            continue
        command = [sys.executable, str(Path(__file__).resolve()), "--single-index", str(index)]
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
                check=False,
            )
            record = parse_worker(completed, row, time.monotonic() - started)
        except subprocess.TimeoutExpired as error:
            output_tail = error.stdout or ""
            if isinstance(output_tail, bytes):
                output_tail = output_tail.decode(errors="replace")
            record = failure_record(
                row,
                "CENSORED_FIBRE_WORKER_TIMEOUT",
                {"output_tail": output_tail.splitlines()[-40:]},
                time.monotonic() - started,
            )
        records.append(record)
        write_chunk(output, chunk_index, chunk_count, indices, records, protocol)
        print(
            f"R17SMALLFIELDDETECTOR|chunk={chunk_index}/{chunk_count}"
            f"|row={position}/{len(indices)}|sample={row['sample_id']}"
            f"|status={record['status']}"
            f"|gain={record.get('stage_a', {}).get('certified_quotient_gain')}",
            flush=True,
        )


def merge_chunks(chunk_dir, chunk_count, output):
    _cohort, protocol = load_inputs()
    records_by_id = {}
    chunks = []
    for chunk_index in range(chunk_count):
        path = chunk_dir / f"chunk-{chunk_index:02d}-of-{chunk_count:02d}.json"
        chunk = json.loads(path.read_text())
        if chunk.get("status") != "COMPLETE_SCHEDULED_DETECTOR_CHUNK":
            raise ArithmeticError(f"detector chunk {chunk_index} is incomplete")
        if chunk.get("protocol_definition_sha256") != protocol["protocol_definition_sha256"]:
            raise ArithmeticError(f"detector chunk {chunk_index} used another protocol")
        for record in chunk["records"]:
            if record["sample_id"] in records_by_id:
                raise ArithmeticError("duplicate detector row across chunks")
            records_by_id[record["sample_id"]] = record
        chunks.append(
            {"path": relative(path), "sha256": digest(path), "record_count": len(chunk["records"])}
        )
    expected_ids = [row["sample_id"] for row in protocol["detector_manifest"]]
    if set(records_by_id) != set(expected_ids):
        raise ArithmeticError("detector chunks do not cover the frozen manifest")
    records = [records_by_id[sample_id] for sample_id in expected_ids]
    document = {
        "schema": SCHEMA,
        "status": "COMPLETE_FROZEN_SMALL_FIELD_DETECTOR_LEDGER",
        "candidate_list_sha256": protocol["candidate_list_sha256"],
        "protocol_definition_sha256": protocol["protocol_definition_sha256"],
        "summary": {
            "scheduled_rows": len(records),
            "status_counts": dict(sorted(Counter(record["status"] for record in records).items())),
            "complete_stage_a_rows": sum(
                record.get("analysis_eligible_complete_stage_a", False) for record in records
            ),
            "stage_a_escape_rows": sum(
                record.get("stage_a", {}).get("certified_quotient_gain", 0) > 0
                for record in records
            ),
            "stage_a_certified_directions": sum(
                record.get("stage_a", {}).get("certified_quotient_gain", 0)
                for record in records
            ),
            "feature_values_loaded_by_detector": False,
        },
        "records": records,
        "chunk_provenance": chunks,
        "inputs": {
            relative(COHORT): digest(COHORT),
            relative(FEATURES): digest(FEATURES),
            relative(PROTOCOL): digest(PROTOCOL),
            relative(BASE_RUNNER): digest(BASE_RUNNER),
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"R17SMALLFIELDDETECTORMERGE|rows={len(records)}"
        f"|escapes={document['summary']['stage_a_escape_rows']}|output={relative(output)}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=16)
    parser.add_argument("--chunk-dir", type=Path, default=DEFAULT_CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    modes = sum(value is not None for value in (args.single_index, args.chunk_index)) + int(args.merge)
    if modes != 1:
        raise SystemExit("choose exactly one of --single-index, --chunk-index, or --merge")
    if args.single_index is not None:
        run_single(args.single_index)
    elif args.chunk_index is not None:
        if not 0 <= args.chunk_index < args.chunk_count:
            raise SystemExit("chunk index is outside chunk count")
        output = args.output or (
            args.chunk_dir / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
        )
        run_chunk(args.chunk_index, args.chunk_count, output.resolve(), args.limit)
    else:
        merge_chunks(
            args.chunk_dir.resolve(), args.chunk_count, (args.output or OUTPUT).resolve()
        )


if __name__ == "__main__":
    main()
