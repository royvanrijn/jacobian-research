#!/usr/bin/env python3
"""Relocate provenance strings after the 2026-08-24 elliptic cleanup.

This is a one-time, deterministic metadata migration.  It does not recompute
any curve, point, rank, conductor, or bounded search.  The pre-migration bytes
are preserved under ``artifacts/snapshots/pre-cleanup-2026-08-24``.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
ACTIVE = ROOT / "artifacts/generated-results/elliptic-curves"
ARCHIVE = ROOT / "archive/elliptic-curves"
SNAPSHOT = ARCHIVE / "artifacts/snapshots/pre-cleanup-2026-08-24"
MANIFEST = ARCHIVE / "MANIFEST.tsv"

RENAMED_ACTIVE = {
    "icarm_7fff_zip_281_282_285_286.json":
        "icarm_7fff_zip_public_source_281_282_285_286.json",
    "icarm_7fff_zip_sequence_analysis.json":
        "icarm_7fff_zip_independence_analysis_v1.json",
    "newfamily_rank_gain_batch_v1.json":
        "newfamily_exact_subgroup_rank_gain_batch_v1.json",
    "newfamily_rank13_exact_bounds_v1.json":
        "newfamily_pari_exact_rank13_batch_v1.json",
    "newfamily_rank14_t83_6_exact_v1.json":
        "newfamily_rank14_t83_6_pari_exact_rank_v1.json",
}

INTERMEDIATE_SOURCE_HASHES = {
    "8b55b0c3708e21e56795ce6f4c393a92139012db9ff7d3f75d2972d1c20c7781":
        "elliptic-curves/cas/search_fermigier_exceptional_pair_simultaneous_h200000.py",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()


def load_archive_mappings() -> dict[str, str]:
    mappings: dict[str, str] = {}
    rows = MANIFEST.read_text().splitlines()[1:]
    for row in rows:
        old, new, _digest, _kind = row.split("\t")
        mappings[old] = new
    return mappings


def active_artifact_mappings() -> dict[str, str]:
    mappings: dict[str, str] = {}
    for artifact in ACTIVE.iterdir():
        if artifact.suffix == ".json" or artifact.name.endswith(".json.gz"):
            mappings[f"artifacts/generated-results/{artifact.name}"] = str(
                artifact.relative_to(ROOT)
            )
    for old_name, new_name in RENAMED_ACTIVE.items():
        mappings[f"artifacts/generated-results/elliptic-curves/{old_name}"] = (
            f"artifacts/generated-results/elliptic-curves/{new_name}"
        )
    return mappings


def active_source_hash_mappings() -> dict[str, str]:
    result = subprocess.run(
        [
            "git",
            "ls-tree",
            "-r",
            "--name-only",
            "HEAD",
            "elliptic-curves/cas",
            "elliptic-curves/scripts",
            "elliptic-curves/tools",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    mappings: dict[str, str] = {}
    for relative in result.stdout.splitlines():
        current = ROOT / relative
        if not current.is_file():
            continue
        original = subprocess.run(
            ["git", "show", f"HEAD:{relative}"],
            cwd=ROOT,
            check=True,
            capture_output=True,
        ).stdout
        old = sha256_bytes(original)
        new = sha256_bytes(current.read_bytes())
        if old != new:
            mappings[old] = new
    for old, relative in INTERMEDIATE_SOURCE_HASHES.items():
        mappings[old] = sha256_bytes((ROOT / relative).read_bytes())
    return mappings


def replace_strings(
    value: Any,
    path_mappings: dict[str, str],
    hash_mappings: dict[str, str],
) -> Any:
    if isinstance(value, dict):
        return {
            replace_strings(key, path_mappings, hash_mappings): replace_strings(
                child, path_mappings, hash_mappings
            )
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [
            replace_strings(child, path_mappings, hash_mappings)
            for child in value
        ]
    if isinstance(value, str):
        if value in hash_mappings:
            return hash_mappings[value]
        if "artifacts/generated-results" not in value and "elliptic-curves/" not in value:
            return value
        result = value
        for old, new in path_mappings.items():
            result = result.replace(str(ROOT / old), str(ROOT / new))
            if old.startswith("artifacts/"):
                result = re.sub(
                    rf"(?<!archive/elliptic-curves/){re.escape(old)}",
                    new,
                    result,
                )
            elif old.startswith("elliptic-curves/"):
                result = re.sub(
                    rf"(?<!archive/){re.escape(old)}",
                    new,
                    result,
                )
        return result
    return value


def preserve_pair_snapshot(data: dict[str, Any]) -> None:
    snapshot = data.get("prior_parameter_snapshot")
    if not isinstance(snapshot, dict) or not isinstance(snapshot.get("sources"), dict):
        return
    sources = snapshot["sources"]
    rewritten: dict[str, Any] = {}
    for raw_path, record in sources.items():
        name = Path(raw_path).name
        frozen = SNAPSHOT / name
        if frozen.is_file():
            rewritten[str(frozen.relative_to(ROOT))] = record
        else:
            rewritten[raw_path] = record
    snapshot["sources"] = rewritten


def stable_projection(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: stable_projection(child)
            for key, child in value.items()
            if key not in {
                "generated_at_utc",
                "wall_seconds",
                "search_wall_seconds",
                "pari_milliseconds",
            }
        }
    if isinstance(value, list):
        return [stable_projection(child) for child in value]
    return value


def common_digest_mode(data: Any) -> bool | None:
    if not isinstance(data, dict):
        return None
    expected = data.get("result_sha256")
    if not isinstance(expected, str):
        return None
    for drop_generated_at in (False, True):
        candidate = copy.deepcopy(data)
        candidate.pop("result_sha256", None)
        if drop_generated_at:
            candidate.pop("generated_at_utc", None)
        if sha256_bytes(stable_json(candidate)) == expected:
            return drop_generated_at
    return None


def refresh_result_digest(
    name: str, data: Any, common_mode: bool | None
) -> None:
    if not isinstance(data, dict):
        return
    if common_mode is not None:
        stable = copy.deepcopy(data)
        stable.pop("result_sha256", None)
        if common_mode:
            stable.pop("generated_at_utc", None)
        data["result_sha256"] = sha256_bytes(stable_json(stable))
        return
    if name == "elliptic_fermigier_exceptional_transport.json":
        stable = {
            key: data[key]
            for key in (
                "anchors",
                "exceptional_quotients",
                "transport",
                "fiber_products",
                "outcome",
                "sources",
            )
        }
        data["result_sha256"] = sha256_bytes(stable_json(stable))
    elif name == "elliptic_fermigier_exceptional_quotient_ball.json":
        stable = {
            "anchors": data["anchors"],
            "direction_balls": data["direction_balls"],
            "affine_transport": data["affine_transport"],
            "outcome": data["outcome"],
            "sources": {
                key: value
                for key, value in data["sources"].items()
                if key != "script_sha256"
            },
        }
        data["result_sha256"] = sha256_bytes(stable_json(stable))
    elif name == "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json":
        stable = {
            "schema_version": data["schema_version"],
            "source": data["source"],
            "search_box": data["search_box"],
            "local_sieve": data["local_sieve"],
            "prior_parameter_snapshot": data["prior_parameter_snapshot"],
            "directions": data["directions"],
            "individual_beyond_anchor_incidences": data[
                "individual_beyond_anchor_incidences"
            ],
            "pair_results": data["pair_results"],
            "third_parameter_certifications": data[
                "third_parameter_certifications"
            ],
            "outcome": data["outcome"],
        }
        data["result_sha256"] = sha256_bytes(
            stable_json(stable_projection(stable))
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--report",
        type=Path,
        default=ARCHIVE / "ARTIFACT_MIGRATION_2026-08-24.tsv",
    )
    args = parser.parse_args()

    path_mappings = load_archive_mappings()
    path_mappings.update(active_artifact_mappings())
    hash_mappings = active_source_hash_mappings()
    rows = []
    for artifact in sorted(ACTIVE.glob("*.json")):
        original_bytes = artifact.read_bytes()
        original = json.loads(original_bytes)
        common_mode = common_digest_mode(original)
        migrated = copy.deepcopy(original)
        if artifact.name == (
            "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json"
        ):
            preserve_pair_snapshot(migrated)
        migrated = replace_strings(migrated, path_mappings, hash_mappings)
        refresh_result_digest(artifact.name, migrated, common_mode)
        if migrated == original:
            continue
        new_bytes = (
            json.dumps(migrated, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        ).encode()
        if new_bytes == original_bytes:
            continue
        rows.append(
            (
                artifact.name,
                sha256_bytes(original_bytes),
                sha256_bytes(new_bytes),
                "metadata paths/hashes only; mathematical payload unchanged",
            )
        )
        if not args.dry_run:
            artifact.write_bytes(new_bytes)

    report = ["filename\told_sha256\tnew_sha256\tscope"]
    report.extend("\t".join(row) for row in rows)
    report_text = "\n".join(report) + "\n"
    if args.dry_run:
        print(report_text, end="")
    else:
        report = args.report if args.report.is_absolute() else ROOT / args.report
        report.write_text(report_text)
        print(f"migrated={len(rows)} report={report.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
