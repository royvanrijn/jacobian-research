#!/usr/bin/env python3
"""Freeze the blinded, hash-ordered detector-v2 Selmer sample.

The source experiment has six lane labels, but the comparison requested for
detector v2 has five cohorts: the two full-cylinder lanes are pooled while
retaining their 356/385 anchor strata.  Every pooled cohort is therefore
sampled symmetrically by anchor.  Stage 1 takes the lexicographically first
pre-existing ``sample_id`` in each anchor stratum (ten fibres total); Stage 2
takes the first three (thirty fibres total).  No arithmetic feature, point
search result, discriminant, conductor, or coefficient size enters the order.

The descent input and unblinding key are deliberately separate files.  A
descent worker needs only the former.  Stage 2 remains committed but is not
authorized until both record controls complete the full pipeline.

status: ACTIVE_PROOF
claim: deterministic blinded Stage-1/Stage-2 sample commitment only
inputs: frozen prospective CRT cohort artifact
outputs: blinded sample and separate unblinding-key artifacts
supersedes/superseded-by: none
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
)
BLINDED_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-quotient-rank-escape-detector-v2-sample-v1.json"
)
KEY_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-quotient-rank-escape-detector-v2-unblinding-key-v1.json"
)

SOURCE_SCHEMA = "elkies-k3.r17-prospective-crt-frozen-cohorts.v1"
SOURCE_STATUS = "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS"
EXPECTED_CANDIDATE_LIST_SHA256 = (
    "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
)
BLINDED_SCHEMA = "elkies-k3.r17-quotient-rank-escape-detector-v2-sample.v1"
KEY_SCHEMA = "elkies-k3.r17-quotient-rank-escape-detector-v2-unblinding-key.v1"
SALT = "r17-quotient-rank-escape-detector-v2-sample-v1"

AGGREGATE_COHORTS = {
    "full_cylinders": {356: "A_356_full", 385: "B_385_full"},
    "matched_ordinary": {356: "C_matched_ordinary", 385: "C_matched_ordinary"},
    "two_only": {356: "D_two_only", 385: "D_two_only"},
    "odd_only": {356: "E_odd_only", 385: "E_odd_only"},
    "random_equal_codimension": {
        356: "F_random_equal_codimension",
        385: "F_random_equal_codimension",
    },
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def blind_id(sample_id: str) -> str:
    return sha256(f"{SALT}|blind|{sample_id}".encode()).hexdigest()[:24]


def build() -> tuple[dict[str, Any], dict[str, Any]]:
    source = json.loads(SOURCE.read_text())
    if source.get("schema") != SOURCE_SCHEMA or source.get("status") != SOURCE_STATUS:
        raise ArithmeticError("the source cohort is not the reviewed frozen commitment")
    if (
        source.get("commitment", {}).get("candidate_list_sha256")
        != EXPECTED_CANDIDATE_LIST_SHA256
    ):
        raise ArithmeticError("the reviewed frozen candidate list changed")

    rows = source.get("rows", [])
    sample_ids = [row.get("sample_id") for row in rows]
    if len(rows) != 2560 or len(set(sample_ids)) != len(rows):
        raise ArithmeticError("the source must contain 2,560 unique hash IDs")
    if any(
        row.get("selection_status") != "FROZEN_VALID_NONSINGULAR_UNOPENED"
        or row.get("outcome_status") != "NOT_OPENED"
        for row in rows
    ):
        raise ArithmeticError("the source rows are not frozen and unopened")

    selected: list[tuple[str, int, str, int, dict[str, Any]]] = []
    for aggregate, lane_by_anchor in AGGREGATE_COHORTS.items():
        for anchor in (356, 385):
            lane = lane_by_anchor[anchor]
            stratum = sorted(
                (
                    row
                    for row in rows
                    if int(row["anchor_curve_id"]) == anchor
                    and row["cohort"] == lane
                ),
                key=lambda row: row["sample_id"],
            )
            if len(stratum) != 256:
                raise ArithmeticError(
                    f"expected 256 rows in {aggregate}/{anchor}, found {len(stratum)}"
                )
            selected.extend(
                (aggregate, anchor, lane, hash_rank, row)
                for hash_rank, row in enumerate(stratum[:3], start=1)
            )

    if len(selected) != 30:
        raise ArithmeticError("Stage 2 must contain exactly thirty fibres")
    if len({row[4]["sample_id"] for row in selected}) != 30:
        raise ArithmeticError("a sampled source row was selected twice")

    blinded_rows = []
    key_rows = []
    for aggregate, anchor, lane, hash_rank, row in selected:
        identifier = blind_id(row["sample_id"])
        blinded_rows.append(
            {
                "blind_id": identifier,
                "chart": row["chart"],
                "parameter": row["parameter"],
                "projective_pair": row["projective_pair"],
                "stage_1_included": hash_rank == 1,
                "stage_2_included": True,
            }
        )
        key_rows.append(
            {
                "blind_id": identifier,
                "source_sample_id": row["sample_id"],
                "source_match_set_id": row["match_set_id"],
                "aggregate_cohort": aggregate,
                "source_lane": lane,
                "anchor_curve_id": anchor,
                "hash_rank_within_anchor_stratum": hash_rank,
                "stage_1_included": hash_rank == 1,
                "stage_2_included": True,
            }
        )

    blinded_rows.sort(key=lambda row: row["blind_id"])
    key_rows.sort(key=lambda row: row["blind_id"])
    if len({row["blind_id"] for row in blinded_rows}) != 30:
        raise ArithmeticError("a blinded ID collision occurred")

    stage1_ids = sorted(
        row["blind_id"] for row in blinded_rows if row["stage_1_included"]
    )
    stage2_ids = sorted(row["blind_id"] for row in blinded_rows)
    if len(stage1_ids) != 10 or len(stage2_ids) != 30:
        raise ArithmeticError("the requested Stage-1/Stage-2 sizes changed")

    key_counts_stage1 = Counter(
        row["aggregate_cohort"]
        for row in key_rows
        if row["stage_1_included"]
    )
    key_counts_stage2 = Counter(row["aggregate_cohort"] for row in key_rows)
    if set(key_counts_stage1.values()) != {2} or set(key_counts_stage2.values()) != {6}:
        raise ArithmeticError("the five aggregate cohorts are not balanced")

    common = {
        "selection_rule": {
            "ordering_key": (
                "lexicographic order of the pre-existing 24-hex sample_id in "
                "each aggregate-cohort/anchor stratum"
            ),
            "stage_1": "hash rank 1 in both anchor strata of each of five cohorts",
            "stage_2": "hash ranks 1,2,3 in both anchor strata of each of five cohorts",
            "arithmetic_results_read": False,
            "coefficient_or_conductor_cost_used": False,
            "stage_2_authorization": (
                "forbidden until complete descents and blind quotient recovery "
                "succeed on both record controls 356 and 385"
            ),
        },
        "commitment": {
            "salt": SALT,
            "stage_1_blind_ids": stage1_ids,
            "stage_2_blind_ids": stage2_ids,
            "stage_1_candidate_count": len(stage1_ids),
            "stage_2_candidate_count": len(stage2_ids),
            "record_controls_in_addition_to_sample": [356, 385],
        },
        "inputs": {
            relative(SOURCE): digest(SOURCE),
            "source_candidate_list_sha256": EXPECTED_CANDIDATE_LIST_SHA256,
        },
    }

    blinded = {
        "schema": BLINDED_SCHEMA,
        "status": "FROZEN_BLINDED_STAGE1_AND_STAGE2_SAMPLE_STAGE2_NOT_AUTHORIZED",
        **deepcopy(common),
        "rows": blinded_rows,
        "claim_boundary": [
            "This file contains no cohort labels, anchor labels, Selmer results, point-search results, conductors, or bad-prime factorizations.",
            "Stage 1 has ten prospective fibres plus the two separately identified record controls.",
            "Stage 2 is a frozen thirty-fibre extension and is not authorized merely because its rows are committed.",
        ],
    }
    blinded["commitment"]["blinded_rows_sha256"] = canonical_hash(blinded_rows)

    key = {
        "schema": KEY_SCHEMA,
        "status": "FROZEN_UNBLINDING_KEY_DO_NOT_READ_DURING_DESCENT",
        **deepcopy(common),
        "blinded_sample_schema": BLINDED_SCHEMA,
        "blinded_rows_sha256": blinded["commitment"]["blinded_rows_sha256"],
        "rows": key_rows,
        "predeclared_reports": [
            "residual-dimension distribution by aggregate cohort",
            "mean, median, and complete value list",
            "uncontrolled bad-prime count",
            "summed and independent local codimensions",
            "maximum single-place suppression",
            "predeclared rank-profile distance from each complete record matrix",
        ],
        "claim_boundary": [
            "The descent runner must not read this file.",
            "Unblinding is permitted only after every attempted row in the active stage has a frozen completed or censored outcome.",
        ],
    }
    key["commitment"]["unblinding_rows_sha256"] = canonical_hash(key_rows)
    return blinded, key


def serialize(document: dict[str, Any]) -> str:
    return json.dumps(document, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BLINDED_OUTPUT)
    parser.add_argument("--unblinding-output", type=Path, default=KEY_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    blinded, key = build()
    outputs = ((args.output.resolve(), blinded), (args.unblinding_output.resolve(), key))
    for path, document in outputs:
        payload = serialize(document)
        if args.check:
            if not path.exists() or path.read_text() != payload:
                raise ArithmeticError(f"stored commitment differs: {path}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload)
    print(
        "R17ESCAPEV2SAMPLE"
        f"|stage1={blinded['commitment']['stage_1_candidate_count']}"
        f"|stage2={blinded['commitment']['stage_2_candidate_count']}"
        "|stage2_authorized=false|status=FROZEN_BLINDED",
        flush=True,
    )


if __name__ == "__main__":
    main()
