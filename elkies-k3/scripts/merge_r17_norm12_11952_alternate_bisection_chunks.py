#!/usr/bin/env python3
"""Merge disjoint exact alternate-Q80 bisection equation chunks."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ORBIT_TABLE = ROOT / "artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.tsv"
HISTORICAL_FRAME = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
EXPECTED_COUNT = 39147


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chunks", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    loaded = []
    for path in arguments.chunks:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "elkies-k3.bisection-extension-input.v1":
            raise ValueError(f"{path}: unexpected input schema")
        if payload.get("artifact_schema") != "elkies-k3.r17-norm12-11952-alternate-bisections.v1":
            raise ValueError(f"{path}: unexpected artifact schema")
        if payload.get("status") != "PASS_EXACT_ALTERNATE_BISECTION_EQUATION_CHUNK":
            raise ValueError(f"{path}: chunk did not pass exact construction")
        interval = payload.get("interval", {})
        start = int(interval.get("start_zero_based", -1))
        stop = int(interval.get("stop_exclusive", -1))
        if stop - start != len(payload.get("bisections", [])):
            raise ValueError(f"{path}: interval length and record count differ")
        loaded.append((start, stop, path, payload))
    loaded.sort(key=lambda item: item[0])
    cursor = 0
    records = []
    labels: set[str] = set()
    masks: set[int] = set()
    chart_counts = {"finite": 0, "inverted_at_infinity": 0}
    for start, stop, path, payload in loaded:
        if start != cursor:
            raise ValueError(f"{path}: expected interval to start at {cursor}, got {start}")
        for record in payload["bisections"]:
            label = str(record["label"])
            mask = int(record["lattice_orbit_mask"])
            if label in labels or mask in masks:
                raise ValueError(f"{path}: duplicate label or orbit mask")
            labels.add(label)
            masks.add(mask)
            records.append(record)
        for chart, count in payload["construction"]["construction_chart_counts"].items():
            chart_counts[chart] = chart_counts.get(chart, 0) + int(count)
        cursor = stop
    if cursor != EXPECTED_COUNT or len(records) != EXPECTED_COUNT:
        raise ValueError(
            f"chunks cover [0,{cursor}) with {len(records)} records; expected {EXPECTED_COUNT}"
        )

    chunk_manifest = [
        {
            "path": relative(path),
            "sha256": digest(path),
            "start_zero_based": start,
            "stop_exclusive": stop,
        }
        for start, stop, path, _ in loaded
    ]
    result = {
        "schema": "elkies-k3.bisection-extension-input.v1",
        "artifact_schema": "elkies-k3.r17-norm12-11952-alternate-bisections-full.v1",
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "base_parameter": "u",
        "invariant_mw_rank": 17,
        "bisections": records,
        "required_lattice_orbits": {
            "table": relative(ORBIT_TABLE),
            "sha256": digest(ORBIT_TABLE),
            "frame_artifact": relative(HISTORICAL_FRAME),
            "frame_sha256": digest(HISTORICAL_FRAME),
            "vector_key": "alternate_rank17_w",
        },
        "construction": {
            "method": "Proposition F1 exact regular residual chord",
            "record_count": len(records),
            "construction_chart_counts": chart_counts,
            "all_branch_fibres_smooth": True,
            "all_lifted_sections_verified": True,
            "complete_translation_orbit_coverage": True,
        },
        "chunk_manifest": chunk_manifest,
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/merge_r17_norm12_11952_alternate_bisection_chunks.py "
            + " ".join(relative(path) for _, _, path, _ in loaded)
            + f" --output {relative(arguments.output)}"
        ),
        "proof_boundary": (
            "There is exactly one exact equation-level degree-two cover and verified lifted "
            "section for every one of the 39147 section-translation classes in the certified "
            "alternate-Q80 lattice table. Equal-cover claims require the separate exact "
            "squareclass collision checker."
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_text() != serialized:
            raise ValueError("stored merged artifact differs from replay")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "ALTFULLMERGE|status=PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS|"
        f"records={len(records)}|chunks={len(loaded)}|output={relative(arguments.output)}"
    )


if __name__ == "__main__":
    main()
