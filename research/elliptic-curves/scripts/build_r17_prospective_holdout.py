#!/usr/bin/env python3
"""Freeze a label-blind random R17 holdout outside the stratified cohort."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_POPULATION = ROOT / "artifacts/local/elliptic-curves/r17-training-population.jsonl"
DEFAULT_SELECTED = ROOT / "artifacts/local/elliptic-curves/r17-training-selected.jsonl"
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_bisection_gain_ranker_quarantined_replay_v1.json"
)
DEFAULT_OUTPUT = ROOT / "artifacts/local/elliptic-curves/r17-prospective-holdout.jsonl"
DEFAULT_SUMMARY = ROOT / "artifacts/local/elliptic-curves/r17-prospective-holdout-summary.json"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def label_blind_key(parameter: str, salt: str) -> bytes:
    return sha256(f"{salt}|{parameter}".encode()).digest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", type=Path, default=DEFAULT_POPULATION)
    parser.add_argument("--selected", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--count", type=int, default=5_000)
    parser.add_argument("--salt", default="r17-prospective-bisection-holdout-v1")
    args = parser.parse_args()
    if args.count < 1:
        raise SystemExit("--count must be positive")

    selected_parameters = {
        json.loads(line)["parameter"]
        for line in args.selected.read_text().splitlines()
        if line.strip()
    }
    candidates = []
    with args.population.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row["parameter"] in selected_parameters:
                continue
            candidates.append((label_blind_key(row["parameter"], args.salt), row))
    if args.count > len(candidates):
        raise SystemExit("holdout count exceeds the unselected population")
    chosen = [row for _key, row in sorted(candidates, key=lambda item: item[0])[: args.count]]
    chosen.sort(key=lambda row: row["parameter"])
    for row in chosen:
        row["split"] = "prospective_holdout"
        row["selection_lanes"] = ["prospective_random_holdout"]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for row in chosen:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    model = json.loads(args.model.read_text())
    summary = {
        "schema": "elliptic-curves.r17-prospective-bisection-holdout.v1",
        "status": "FROZEN_BEFORE_BISECTION_LABEL_EVALUATION",
        "selection": {
            "kind": "smallest SHA-256 keys among rows outside the original selected cohort",
            "salt": args.salt,
            "count": len(chosen),
            "source_unselected_count": len(candidates),
            "output": str(args.output.resolve()),
            "output_sha256": file_sha256(args.output),
        },
        "commitment": {
            "bisection_label_file_read": False,
            "frozen_model_sha256": model["frozen_model_sha256_before_quarantine_open"],
            "model_artifact_sha256": file_sha256(args.model),
            "evaluation_methods": [
                "frozen learned contrast",
                "weakest-block Nagao",
                "partial conductor-quality proxy",
                "negative log projective height",
            ],
            "budgets": ["1/100", "1/20", "1/10"],
        },
        "inputs": {
            str(args.population.resolve()): file_sha256(args.population),
            str(args.selected.resolve()): file_sha256(args.selected),
            str(args.model.resolve()): file_sha256(args.model),
        },
        "script_sha256": file_sha256(Path(__file__)),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
