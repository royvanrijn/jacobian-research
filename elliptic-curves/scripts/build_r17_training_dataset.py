#!/usr/bin/env python3
"""Build a deterministic, embargoed cheap-feature population for R17."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import log
from pathlib import Path
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
ECSEARCH = ROOT / "elliptic-curves/ecsearch"
K3_SCRIPTS = ROOT / "elkies-k3/scripts"
for directory in (ECSEARCH, K3_SCRIPTS):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from r17_training_data import (  # noqa: E402
    EMBARGOED_PARAMETERS,
    build_cheap_prime_tables,
    conductor_proxy_features,
    cover_character_tables,
    cover_diversity_features,
    deterministic_sample,
    development_lane_memberships,
    nagao_features,
    parameter_text,
    quotient_code_features,
    select_cover_panel,
    split_bucket,
)
from search_h92_q12o5867_rootless_nagao import is_prime, load_family_model  # noqa: E402


DEFAULT_OUTPUT = ROOT / "artifacts/local/elliptic-curves/r17-training-population.jsonl"
DEFAULT_SELECTED = ROOT / "artifacts/local/elliptic-curves/r17-training-selected.jsonl"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections.json"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_jsonl(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=100_000)
    parser.add_argument("--height", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20_260_902)
    parser.add_argument("--split-salt", default="r17-training-v1")
    parser.add_argument("--lane-size", type=int, default=1_000)
    parser.add_argument("--cover-panel-size", type=int, default=128)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--selected-output", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()

    started = perf_counter()
    model = load_family_model(MODEL)
    model_document = json.loads(MODEL.read_text())
    a_coefficients = [int(value) for value in model_document["A_coefficients_low_to_high"]]
    b_coefficients = [int(value) for value in model_document["B_coefficients_low_to_high"]]
    score_primes = [value for value in range(19, 600) if is_prime(value)]
    prime_blocks = [score_primes[index::3] for index in range(3)]
    quotient_primes = [19, 23, 29, 31, 37]
    conductor_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    all_table_primes = sorted(set(score_primes + quotient_primes))
    prime_tables = build_cheap_prime_tables(model, all_table_primes)

    bisection_document = json.loads(BISECTIONS.read_text())
    cover_panel = select_cover_panel(bisection_document["bisections"], args.cover_panel_size)
    cover_primes = [101, 103, 107, 109, 113]
    cover_tables = cover_character_tables(cover_panel, cover_primes)

    parameters = deterministic_sample(
        count=args.count,
        height=args.height,
        seed=args.seed,
        excluded=EMBARGOED_PARAMETERS,
    )
    if EMBARGOED_PARAMETERS.intersection(parameters):
        raise AssertionError("an embargoed parameter entered the development population")

    records = []
    for parameter in parameters:
        quotient = quotient_code_features(parameter, prime_tables, quotient_primes)
        record = {
            "parameter": parameter_text(parameter),
            "projective_pair": list(parameter),
            "height": max(abs(parameter[0]), parameter[1]),
            "split": split_bucket(parameter, args.split_salt),
            "features": {
                "level0": {
                    "numerator": parameter[0],
                    "denominator": parameter[1],
                    "projective_height": max(abs(parameter[0]), parameter[1]),
                },
                "level1_nagao": nagao_features(parameter, prime_tables, prime_blocks),
                "level2_conductor_proxy": conductor_proxy_features(
                    parameter, a_coefficients, b_coefficients, conductor_primes
                ),
                "level2_quotient_code": quotient,
                "level2_cover_diversity": cover_diversity_features(
                    parameter, cover_tables, cover_primes
                ),
            },
        }
        records.append(record)

    training_rows = [row for row in records if row["split"] == "train"]
    code_counts = Counter(
        row["features"]["level2_quotient_code"]["code"] for row in training_rows
    )
    smoothed_population = len(training_rows) + len(code_counts)
    for row in records:
        quotient = row["features"]["level2_quotient_code"]
        frequency = code_counts.get(quotient["code"], 0)
        quotient["training_frequency"] = frequency
        quotient["rarity"] = log_ratio(frequency + 1, smoothed_population)
        quotient["rarity_reference"] = "train split with add-one smoothing"

    lanes = development_lane_memberships(records, args.lane_size, args.split_salt)
    memberships: dict[int, list[str]] = {}
    for lane, indices in lanes.items():
        for index in indices:
            memberships.setdefault(index, []).append(lane)
    selected_rows = []
    for index in sorted(memberships):
        row = dict(records[index])
        row["selection_lanes"] = sorted(memberships[index])
        selected_rows.append(row)

    write_jsonl(args.output, records)
    write_jsonl(args.selected_output, selected_rows)
    summary = {
        "schema": "elliptic-curves.r17-training-data-summary.v1",
        "status": "EXPERIMENTAL_FEATURE_POPULATION_NO_RANK_CLAIM",
        "population": {
            "count": len(records),
            "height_bound": args.height,
            "sampling_seed": args.seed,
            "split_salt": args.split_salt,
            "split_counts": dict(Counter(row["split"] for row in records)),
            "output": str(args.output.resolve()),
            "output_sha256": file_sha256(args.output),
        },
        "embargo": {
            "parameters": [parameter_text(value) for value in sorted(EMBARGOED_PARAMETERS)],
            "all_absent_from_development_population": True,
            "primary_quarantined_replay": ["2456/135", "-9529/5471"],
        },
        "selection": {
            "lane_size": args.lane_size,
            "lane_counts": {name: len(indices) for name, indices in lanes.items()},
            "union_count": len(selected_rows),
            "overlap_count": sum(len(values) > 1 for values in memberships.values()),
            "selected_output": str(args.selected_output.resolve()),
            "selected_output_sha256": file_sha256(args.selected_output),
        },
        "inputs": {
            "model": str(MODEL.resolve()),
            "model_sha256": file_sha256(MODEL),
            "bisection_panel_source": str(BISECTIONS.resolve()),
            "bisection_panel_source_sha256": file_sha256(BISECTIONS),
            "cover_panel_labels_sha256": sha256(
                "\n".join(label for label, _quadratic in cover_panel).encode()
            ).hexdigest(),
        },
        "features": {
            "nagao_prime_blocks": prime_blocks,
            "quotient_code_primes": quotient_primes,
            "conductor_proxy_primes": conductor_primes,
            "cover_character_primes": cover_primes,
            "cover_panel_size": len(cover_panel),
        },
        "proof_boundary": [
            "Every row is heuristic training data, not a rank or Selmer statement.",
            "The conductor feature is a partial local proxy, not an exact conductor.",
            "The quotient code describes E(F_p)/2E(F_p) and E(F_p)/3E(F_p), not a Mordell-Weil quotient.",
            "The cover feature is modular character diversity, not a count of rational split covers.",
            "Point and cover searches remain subject to the repository residual-Selmer gate.",
        ],
        "runtime_seconds": perf_counter() - started,
    }
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


def log_ratio(numerator: int, denominator: int) -> float:
    return log(denominator / numerator)


if __name__ == "__main__":
    main()
