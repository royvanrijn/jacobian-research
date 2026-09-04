#!/usr/bin/env sage-python
"""Compute or merge the frozen R17 CRT pre-search arithmetic feature panel.

Chunk workers read only the Phase-1 local certificate and the unopened
Phase-2 manifest.  They compute exact specializations, generic-MW17 local
Kummer presentations, intersection/leave-one-place-out ranks, component
data, the fixed Nagao comparison panel, and a fail-closed monotone-sieve
record.  They never call a Mordell--Weil point search.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys
from typing import Any

from sage.all import QQ


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
PHASE1 = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-local-stability-v1.json"
LOCAL_IMPLEMENTATION = ROOT / "elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage"
DEFAULT_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/r17-prospective-crt-features-v1"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"

SCHEMA = "elkies-k3.r17-prospective-crt-arithmetic-features.v1"
CHUNK_SCHEMA = "elkies-k3.r17-prospective-crt-arithmetic-features-chunk.v1"
EXPECTED_CANDIDATE_LIST_HASH = "5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb"
PLACES = (2, 13, 37, 53, 67, 71)
SEARCH_LIMITS = {
    "height": 12,
    "wall_seconds": 300,
    "memory_bytes": 8_000_000_000,
}

sys.path.insert(0, str(ROOT / "elkies-k3/scripts"))
from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    Candidate,
    DEFAULT_PRIME_BLOCKS,
    FamilyModel,
    SCORE_SCALE,
    build_residue_tables,
    local_symbol_record,
    projective_index,
    score_block,
)

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from elkies_residual_selmer_gate import monotone_sieve_gate_record  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def write_json_document(path: Path, document):
    serialized = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as compressed:
                compressed.write(serialized)
    else:
        path.write_bytes(serialized)


def load_inputs():
    manifest = json.loads(MANIFEST.read_text())
    phase1 = json.loads(PHASE1.read_text())
    if manifest.get("status") != "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS":
        raise ArithmeticError("the prospective cohort is not frozen and unopened")
    if manifest["commitment"]["candidate_list_sha256"] != EXPECTED_CANDIDATE_LIST_HASH:
        raise ArithmeticError("the reviewed frozen candidate list changed")
    if any(row["outcome_status"] != "NOT_OPENED" for row in manifest["rows"]):
        raise ArithmeticError("a point-search outcome was opened before feature computation")
    if phase1.get("frozen_cylinder_definition_sha256") != manifest["commitment"][
        "phase1_frozen_cylinder_definition_sha256"
    ]:
        raise ArithmeticError("the feature job is not bound to the frozen Phase-1 cylinders")
    return manifest, phase1


def local_rows(fingerprint):
    local = fingerprint["comparison_payload"]["known_mw17_localization"]
    if fingerprint["rational_prime"] == 2:
        rows = local["canonical_source_order_relation_rows"]
    else:
        rows = fingerprint["diagnostic_raw_odd_squareclass_rows"]
    if len(rows) != 17:
        raise ArithmeticError("a local MW17 presentation does not have seventeen source rows")
    return rows


def localization_panel(local_implementation, fingerprints):
    by_prime = {int(row["rational_prime"]): local_rows(row) for row in fingerprints}
    cumulative = []
    columns = [[] for _ in range(17)]
    for prime in PLACES:
        rows = by_prime[prime]
        columns = [left + right for left, right in zip(columns, rows)]
        rank = local_implementation["f2_rank"](columns)
        cumulative.append(
            {
                "through_prime": prime,
                "stacked_localization_rank": rank,
                "simultaneous_source_kernel_dimension": 17 - rank,
            }
        )
    full_rank = local_implementation["f2_rank"](columns)
    delete_one = []
    for omitted in PLACES:
        kept = [[] for _ in range(17)]
        for prime in PLACES:
            if prime != omitted:
                kept = [left + right for left, right in zip(kept, by_prime[prime])]
        rank = local_implementation["f2_rank"](kept)
        delete_one.append(
            {
                "deleted_prime": omitted,
                "stacked_localization_rank": rank,
                "rank_drop": full_rank - rank,
                "simultaneous_source_kernel_dimension": 17 - rank,
            }
        )
    return {
        "matrix_convention": (
            "Rows are the seventeen fixed generic sections; columns concatenate each "
            "place's exact local squareclass presentation. Ranks and source kernels are "
            "invariant under local target-coordinate changes."
        ),
        "full_stacked_localization_rank": full_rank,
        "full_simultaneous_source_kernel_dimension": 17 - full_rank,
        "cumulative_in_predeclared_place_order": cumulative,
        "leave_one_place_out": delete_one,
    }


def nagao_context(local_implementation):
    family = local_implementation["Family"]()
    model = FamilyModel(
        source=LINEAGE_PLACEHOLDER,
        source_sha256=digest(LOCAL_IMPLEMENTATION),
        a_coefficients=tuple(Fraction(value) for value in family.a_coefficients),
        b_coefficients=tuple(Fraction(value) for value in family.b_coefficients),
        a_degree=8,
        b_degree=12,
        coordinate="native integer parameter of norm12-orbit-074d9",
        coefficient_source_keys=("A_coefficients_low_to_high", "B_coefficients_low_to_high"),
    )
    tables, rejected = build_residue_tables(model, DEFAULT_PRIME_BLOCKS)
    return model, tables, rejected


# FamilyModel records provenance only; the real input hashes are stored in the
# chunk/final artifacts below.
LINEAGE_PLACEHOLDER = MANIFEST


def nagao_features(parameter: int, tables):
    candidate = Candidate(numerator=parameter, denominator=1, height=abs(parameter))
    inverse_cache = {}
    for block in tables:
        candidate = score_block(candidate, block, inverse_cache)
    symbols = []
    for block_number, block in enumerate(tables, start=1):
        for prime, table in block.items():
            symbol = table[projective_index(parameter, 1, prime)]
            record = local_symbol_record(symbol)
            record.update({"prime": prime, "block": block_number})
            symbols.append(record)
    return {
        "status": "HEURISTIC_COMPARISON_FEATURE_ONLY",
        "score_scale": SCORE_SCALE,
        "prime_blocks": [list(block) for block in DEFAULT_PRIME_BLOCKS],
        "block_score_units_1e12": list(candidate.block_score_units),
        "total_score_units_1e12": candidate.total_score_units,
        "good_prime_count": candidate.good_primes,
        "bad_reduction_prime_count": candidate.bad_primes,
        "local_symbols": symbols,
    }


def feature_row(row, local_implementation, family, target_hashes, tables):
    parameter = QQ(row["parameter"])
    result = {
        "sample_id": row["sample_id"],
        "match_set_id": row["match_set_id"],
        "anchor_curve_id": int(row["anchor_curve_id"]),
        "cohort": row["cohort"],
        "parameter": row["parameter"],
        "status": None,
        "failure": None,
    }
    try:
        curve, points = family.specialize(parameter)
        fingerprints = []
        intended_hashes = target_hashes[str(row["anchor_curve_id"])]
        for prime in PLACES:
            fingerprint = local_implementation["local_fingerprint"](curve, points, prime)
            target_hash = intended_hashes.get(str(prime))
            fingerprint["anchor_target_comparison"] = (
                {
                    "status": "EXACT_MATCH" if fingerprint["comparison_sha256"] == target_hash else "EXACT_MISMATCH",
                    "target_comparison_sha256": target_hash,
                }
                if target_hash is not None
                else {"status": "NOT_AN_ANCHOR_TARGET_PLACE", "target_comparison_sha256": None}
            )
            fingerprints.append(fingerprint)
        intended = [
            local["anchor_target_comparison"]["status"] == "EXACT_MATCH"
            for local in fingerprints
            if local["anchor_target_comparison"]["status"] != "NOT_AN_ANCHOR_TARGET_PLACE"
        ]
        result.update(
            {
                "status": "PASS_EXACT_PRESEARCH_ARITHMETIC_PANEL",
                "exact_specialization": {
                    "canonical_short_ainvs": [local_implementation["rational_text"](value) for value in curve.a_invariants()],
                    "discriminant": local_implementation["rational_text"](curve.discriminant()),
                    "generic_mw17_section_count_verified": len(points),
                },
                "local_fingerprints": fingerprints,
                "anchor_fingerprint_survival": {
                    "matched_intended_place_count": sum(intended),
                    "intended_place_count": len(intended),
                    "all_intended_places_match": all(intended),
                },
                "localization_intersections": localization_panel(local_implementation, fingerprints),
                "nagao_comparison": nagao_features(int(parameter), tables),
                "monotone_residual_selmer_gate": monotone_sieve_gate_record(
                    stages=[
                        {
                            "stage": "presearch_local_fingerprint_only",
                            "residual_upper_bound": None,
                            "proof_status": "NO_FINITE_UPPER_BOUND_YET",
                        }
                    ],
                    search_limits=SEARCH_LIMITS,
                ),
                "selmer_measurement_status": {
                    "complete_two_selmer_group": "NOT_COMPUTED",
                    "proved_residual_upper_bound": None,
                    "bounded_search_authorization_only": True,
                },
            }
        )
    except Exception as exc:
        result["status"] = "FAIL_EXACT_PRESEARCH_ARITHMETIC_PANEL"
        result["failure"] = {"type": type(exc).__name__, "message": str(exc)}
    return result


def run_chunk(chunk_index: int, chunk_count: int, output: Path):
    if not (0 <= chunk_index < chunk_count):
        raise ValueError("chunk index is outside the declared chunk count")
    manifest, phase1 = load_inputs()
    local_implementation = runpy.run_path(str(LOCAL_IMPLEMENTATION))
    family = local_implementation["Family"]()
    target_hashes = {
        anchor: {
            str(local["rational_prime"]): local["comparison_sha256"]
            for local in record["local_fingerprints"]
        }
        for anchor, record in phase1["anchors"].items()
    }
    _model, tables, rejected = nagao_context(local_implementation)
    selected = [
        row for index, row in enumerate(manifest["rows"]) if index % chunk_count == chunk_index
    ]
    feature_rows = []
    for ordinal, row in enumerate(selected, start=1):
        feature_rows.append(feature_row(row, local_implementation, family, target_hashes, tables))
        if ordinal % 25 == 0:
            print(
                f"R17CRTFEATURES|chunk={chunk_index}/{chunk_count}|completed={ordinal}/{len(selected)}",
                flush=True,
            )
    document = {
        "schema": CHUNK_SCHEMA,
        "status": "COMPLETE_PRESEARCH_FEATURE_CHUNK",
        "chunk_index": chunk_index,
        "chunk_count": chunk_count,
        "candidate_list_sha256": EXPECTED_CANDIDATE_LIST_HASH,
        "selected_row_count": len(selected),
        "rejected_nagao_model_primes": list(rejected),
        "rows": feature_rows,
        "inputs": {relative(MANIFEST): digest(MANIFEST), relative(PHASE1): digest(PHASE1)},
    }
    write_json_document(output, document)
    print(
        f"R17CRTFEATURES|chunk={chunk_index}/{chunk_count}|rows={len(selected)}|status=COMPLETE",
        flush=True,
    )


def merge_chunks(chunk_dir: Path, chunk_count: int, output: Path):
    manifest, phase1 = load_inputs()
    chunks = []
    for index in range(chunk_count):
        path = chunk_dir / f"chunk-{index:02d}-of-{chunk_count:02d}.json"
        document = json.loads(path.read_text())
        if document.get("schema") != CHUNK_SCHEMA or document.get("status") != "COMPLETE_PRESEARCH_FEATURE_CHUNK":
            raise ArithmeticError(f"feature chunk {index} is incomplete")
        if document["chunk_index"] != index or document["chunk_count"] != chunk_count:
            raise ArithmeticError("feature chunk coordinates changed")
        if document["candidate_list_sha256"] != EXPECTED_CANDIDATE_LIST_HASH:
            raise ArithmeticError("a feature chunk belongs to another candidate list")
        chunks.append((path, document))
    by_id = {}
    for _path, chunk in chunks:
        for row in chunk["rows"]:
            if row["sample_id"] in by_id:
                raise ArithmeticError("duplicate sample across feature chunks")
            by_id[row["sample_id"]] = row
    expected_ids = [row["sample_id"] for row in manifest["rows"]]
    if set(by_id) != set(expected_ids):
        raise ArithmeticError("feature chunks do not cover the frozen candidate list exactly")
    rows = [by_id[sample_id] for sample_id in expected_ids]
    status_counts = Counter(row["status"] for row in rows)
    survival = defaultdict(Counter)
    for row in rows:
        if row["status"] == "PASS_EXACT_PRESEARCH_ARITHMETIC_PANEL":
            survival[row["cohort"]][
                str(row["anchor_fingerprint_survival"]["matched_intended_place_count"])
            ] += 1
        else:
            survival[row["cohort"]]["FAILED_PANEL"] += 1
    document = {
        "schema": SCHEMA,
        "status": (
            "PASS_COMPLETE_EXACT_PRESEARCH_ARITHMETIC_PANEL"
            if status_counts == {"PASS_EXACT_PRESEARCH_ARITHMETIC_PANEL": len(rows)}
            else "COMPLETE_WITH_RETAINED_PRESEARCH_FAILURES"
        ),
        "candidate_list_sha256": EXPECTED_CANDIDATE_LIST_HASH,
        "phase1_frozen_cylinder_definition_sha256": phase1["frozen_cylinder_definition_sha256"],
        "summary": {
            "scheduled_candidates": len(rows),
            "status_counts": dict(sorted(status_counts.items())),
            "anchor_fingerprint_matched_place_count_by_cohort": {
                cohort: dict(sorted(counts.items())) for cohort, counts in sorted(survival.items())
            },
            "complete_two_selmer_groups_computed": 0,
            "finite_proved_residual_upper_bounds_computed": 0,
            "bounded_search_authorizations": sum(
                row.get("monotone_residual_selmer_gate", {}).get("bounded_point_search_authorized", False)
                for row in rows
            ),
        },
        "rows": rows,
        "chunk_provenance": [
            {"path": relative(path), "sha256": digest(path), "rows": chunk["selected_row_count"]}
            for path, chunk in chunks
        ],
        "inputs": {relative(MANIFEST): digest(MANIFEST), relative(PHASE1): digest(PHASE1)},
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "commands": [
                f"sage -python elkies-k3/scripts/run_r17_prospective_crt_arithmetic_features.sage --chunk-index I --chunk-count {chunk_count}",
                f"sage -python elkies-k3/scripts/run_r17_prospective_crt_arithmetic_features.sage --merge --chunk-count {chunk_count}",
            ],
        },
        "claim_boundary": [
            "Every row is a pre-search arithmetic measurement on a frozen candidate.",
            "The local Kummer presentations concern the known generic MW17 subgroup, not the complete global Selmer group.",
            "No finite residual-Selmer upper bound is claimed when the monotone sequence contains only NO_FINITE_UPPER_BOUND_YET.",
            "Nagao values are comparison features and did not select or rebalance any row.",
            "The explicit Hilbert/Tate field remains deferred uniformly; it was not used to tune a cylinder.",
        ],
    }
    write_json_document(output, document)
    print(
        f"R17CRTFEATURES|rows={len(rows)}|status={document['status']}|output={relative(output)}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=8)
    parser.add_argument("--chunk-dir", type=Path, default=DEFAULT_CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    if args.merge:
        merge_chunks(args.chunk_dir.resolve(), args.chunk_count, (args.output or OUTPUT).resolve())
        return
    if args.chunk_index is None:
        raise SystemExit("chunk mode requires --chunk-index")
    output = args.output or (
        args.chunk_dir / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
    )
    run_chunk(args.chunk_index, args.chunk_count, output.resolve())


if __name__ == "__main__":
    main()
