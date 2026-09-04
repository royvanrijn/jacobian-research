#!/usr/bin/env python3
"""Compile the uniform H<=1000 extreme-anchored MW>=18 Nagao campaign."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "artifacts/generated-results"
OUTPUT = ARTIFACTS / "elkies-k3-r17-extreme-anchored-mw18-nagao-h1000-summary-v1.json"
STATUS = "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_CAMPAIGN"
RUN_STATUS = "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_SIEVE"
EXPECTED_RUNS = (
    ("curve-545", "07ca9-orbit-08c1e", "curve545-orbit08c1e"),
    ("curve-545", "07ca9-orbit-1d516", "curve545-orbit1d516"),
    ("curve-531", "08234-orbit-0a9bf", "curve531-orbit0a9bf"),
    ("curve-531", "08234-orbit-12f61", "curve531-orbit12f61"),
    ("curve-531", "08234-orbit-1293d", "curve531-orbit1293d"),
    ("curve-534", "08234-orbit-13d7a", "curve534-orbit13d7a"),
    ("curve-534", "08234-orbit-1a371", "curve534-orbit1a371"),
    ("curve-536", "08234-orbit-19188", "curve536-orbit19188"),
    ("historical-rank28", "orbit-15a68", "historical-rank28-orbit15a68"),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def input_path(stem: str) -> Path:
    return ARTIFACTS / f"elkies-k3-r17-extreme-anchored-mw18-nagao-{stem}-h1000-v1.json"


def build() -> dict:
    records = []
    reference_search = None
    reference_certificate_hash = None
    for anchor_id, label, stem in EXPECTED_RUNS:
        path = input_path(stem)
        document = json.loads(path.read_text())
        if document.get("status") != RUN_STATUS:
            raise ValueError(f"{path.name} does not have the exact bounded-search status")
        cover = document["cover"]
        if (cover["anchor_id"], cover["label"]) != (anchor_id, label):
            raise ValueError(f"{path.name} identifies the wrong anchored cover")
        search = {key: document["search"][key] for key in (
            "numerator_bound", "denominator_bound", "prime_blocks",
            "keep_per_bucket", "height_bucket_width",
        )}
        if reference_search is None:
            reference_search = search
            reference_certificate_hash = cover["certificate_sha256"]
        if search != reference_search:
            raise ValueError("the anchored runs do not use one uniform search protocol")
        if cover["certificate_sha256"] != reference_certificate_hash:
            raise ValueError("the anchored runs do not use one cover certificate")
        finalists = document["finalists"]
        records.append(
            {
                "anchor_id": anchor_id,
                "cover_label": label,
                "source": str(path.relative_to(ROOT)),
                "source_sha256": digest(path),
                "displayed_anchor_jump_over_MW17": cover["displayed_anchor_jump_over_MW17"],
                "stage_one_population_scored": document["stages"][0]["population_scored"],
                "final_survivor_count": document["final_survivor_count"],
                "certified_anchor_survived_heuristic_filter": any(
                    finalist["is_certified_anchor"] for finalist in finalists
                ),
                "top_finalist": finalists[0] if finalists else None,
                "top_five_finalists": finalists[:5],
                "rejected_primes": document["search"]["rejected_primes"],
            }
        )
    records.sort(
        key=lambda record: (
            -(record["top_finalist"] or {}).get("total_score_units_1e12", -10**30),
            record["anchor_id"],
            record["cover_label"],
        )
    )
    return {
        "schema": "elkies-k3.r17-extreme-anchored-mw18-nagao-campaign.v1",
        "status": STATUS,
        "claim": (
            "Uniform bounded Nagao triage of all eight refreshed extreme-anchored "
            "rigid-bisection covers and the historical rank-28 anchored cover."
        ),
        "cover_certificate_sha256": reference_certificate_hash,
        "uniform_search": reference_search,
        "cover_count": len(records),
        "total_stage_one_population_scored": sum(
            record["stage_one_population_scored"] for record in records
        ),
        "total_final_survivor_count": sum(
            record["final_survivor_count"] for record in records
        ),
        "priority_order_by_top_nagao_score": records,
        "proof_boundary": (
            "This is a bounded heuristic ordering. The exact source certificate proves "
            "the covers and nonzero anchors, but no finalist rank or independence claim."
        ),
        "reproducing_command": (
            "run search_r17_extreme_anchored_mw18_nagao.py at H<=1000 for the nine "
            "declared covers, then run this compiler"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise SystemExit("anchored Nagao summary differs from the pinned artifact")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    payload = json.loads(serialized)
    print(
        f"PASS extreme_anchored_mw18_nagao_campaign covers={payload['cover_count']} "
        f"population={payload['total_stage_one_population_scored']} "
        f"survivors={payload['total_final_survivor_count']} output={args.output}"
    )


if __name__ == "__main__":
    main()
