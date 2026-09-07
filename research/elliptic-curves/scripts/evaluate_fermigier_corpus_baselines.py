#!/usr/bin/env python3
"""Evaluate frozen cheap-score baselines on the Fermigier labelled corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ECSEARCH = ROOT / "elliptic-curves/ecsearch"
if str(ECSEARCH) not in sys.path:
    sys.path.insert(0, str(ECSEARCH))

from fermigier_baseline_evaluation import evaluate, file_sha256  # noqa: E402


DEFAULT_CONFIG = ROOT / "elliptic-curves/data/fermigier_baseline_rankers_v1.json"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/local/elliptic-curves/fermigier-baseline-evaluation-v1.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    result = evaluate(config, ROOT)
    result["reproducibility"] = {
        "command": "python3 elliptic-curves/scripts/evaluate_fermigier_corpus_baselines.py",
        "config": str(args.config.resolve().relative_to(ROOT)),
        "config_sha256": file_sha256(args.config),
        "evaluator": "elliptic-curves/ecsearch/fermigier_baseline_evaluation.py",
        "evaluator_sha256": file_sha256(
            ROOT / "elliptic-curves/ecsearch/fermigier_baseline_evaluation.py"
        ),
        "entry_point": "elliptic-curves/scripts/evaluate_fermigier_corpus_baselines.py",
        "entry_point_sha256": file_sha256(Path(__file__)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    primary = next(
        row for row in result["rankers"] if row["id"] == "family-residual-staged-composite"
    )
    positions = ",".join(
        f"{key}={value['position_one_based']}"
        for key, value in sorted(primary["positive_positions"].items())
    )
    print(
        "FERMIGIER_BASELINE_EVAL_PASS "
        f"rows={result['corpus']['row_count']} rankers={len(result['rankers'])} "
        f"primary_positions={positions}"
    )


if __name__ == "__main__":
    main()
