#!/usr/bin/env python3
"""Freeze, replay, and index the second foundry publication cutoff."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import tempfile

import certify_compact_r17_candidates as cert
import refresh_foundry_curve_ledger as v1
from research_curve_refresh import MANIFEST, ROOT
from v3_warm_support import atomic, read, require, sha

ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
SOURCE = ART / "high-rank-foundry-v3"
CONDUCTOR_SOURCE = ART / "foundry_v3_conductors_v1"
CONDUCTOR_LOCAL = ROOT / "artifacts" / "local" / "elliptic-curves" / "foundry-v3-conductors-v1"
CURVES = ART / "foundry_curve_ledger_snapshot_v2.json"
CONDUCTORS = ART / "foundry_v3_conductor_snapshot_v1.json"
CLAIM = "EC-FOUNDRY-V3-LEDGER-CONDUCTORS-SNAPSHOT-20260910"
THRESHOLD = 22


def freeze():
    require(not CURVES.exists() and not CONDUCTORS.exists(),
            "preserve immutable publication snapshots; choose a new version")
    database_path = ROOT / "elliptic-curves" / "data" / "research_curves" / "database.json"
    baseline = [{k: row[k] for k in ("id", "ainvs", "rank_lower_bound")}
                for row in read(database_path)["curves"]]
    # Pin paths before replay. Later source endpoints belong to another cutoff.
    paths = sorted(SOURCE.glob("job-*.json"))
    candidates = []
    for path in paths:
        row = read(path)
        if row.get("status") == "PASS_CERTIFIED_SEARCH" and row.get("rank_lower_bound", 0) >= THRESHOLD:
            candidates.append((str(path.relative_to(ROOT)), row))
    records = v1.select(candidates, json.loads(json.dumps(baseline)))
    curve_snapshot = {
        "status": "FROZEN_FOUNDRY_LEDGER_SELECTION_V2",
        "records": records,
        "baseline": baseline,
        "baseline_database_sha256": sha(database_path),
        "threshold": THRESHOLD,
        "source_paths_considered": len(paths),
        "eligible_endpoints": len(candidates),
        "new_count": sum(e["previous_rank_lower_bound"] is None for e in records.values()),
        "strengthened_count": sum(e["previous_rank_lower_bound"] is not None for e in records.values()),
        "rank_counts": dict(sorted(Counter(
            str(e["result"]["rank_lower_bound"]) for e in records.values()).items())),
        "selection": "All sealed rank-at-least-22 foundry-v3 endpoints at a fixed path cutoff; strongest verified subgroup per rational Q-isomorphism class relative to the 321-curve baseline.",
        "boundary": "Ranks are certified lower bounds. Distinctness is relative to the pinned baseline and selected cohort, not worldwide novelty. No exact-rank or record claim.",
    }
    atomic(CURVES, curve_snapshot, immutable=True)

    conductor_records = {}
    for identifier, entry in records.items():
        source_id = entry["result"]["id"]
        result_path = CONDUCTOR_SOURCE / f"{source_id}.json"
        input_path = CONDUCTOR_LOCAL / "cases" / source_id / "input.json"
        if not result_path.exists() or not input_path.exists():
            continue
        result = read(result_path)
        if result.get("status") != "PASS_INDEPENDENT_CONDUCTOR_REPLAY":
            continue
        source = read(input_path)
        require(result["certificate"]["id"] == source_id == source["id"],
                "conductor source identity differs")
        require(cert.isomorphic(result["certificate"]["curve"], entry["result"]["packet"]["curve"]),
                "conductor curve differs from selected foundry curve")
        conductor_records[identifier] = {
            "input": source,
            "result": result,
            "origin": str(result_path.relative_to(ROOT)),
            "origin_sha256": sha(result_path),
            "input_origin": str(input_path.relative_to(ROOT)),
            "input_sha256": sha(input_path),
            "source_foundry_id": source_id,
        }
    conductor_snapshot = {
        "status": "PASS_FROZEN_FOUNDRY_V3_CONDUCTOR_CUTOFF",
        "records": conductor_records,
        "record_count": len(conductor_records),
        "exact_count": sum(e["result"]["certificate"]["status"] == "EXACT"
                           for e in conductor_records.values()),
        "unknown_count": sum(e["result"]["certificate"]["status"] == "UNKNOWN"
                             for e in conductor_records.values()),
        "curve_snapshot_sha256": sha(CURVES),
        "selection": "Every independently replayed conductor endpoint already available for the fixed selected curve cohort; later queue completions are excluded.",
        "boundary": "Exact only when the remaining cofactor is one. Partial endpoints remain UNKNOWN bounds. No conductor-record claim.",
    }
    atomic(CONDUCTORS, conductor_snapshot, immutable=True)
    print("FOUNDRY_LEDGER_V2_FROZEN", curve_snapshot["new_count"], "new;",
          curve_snapshot["strengthened_count"], "strengthened;",
          curve_snapshot["rank_counts"], flush=True)
    print("FOUNDRY_CONDUCTOR_V1_FROZEN", conductor_snapshot["record_count"], "records;",
          conductor_snapshot["exact_count"], "exact;",
          conductor_snapshot["unknown_count"], "unknown", flush=True)


def check():
    from inventory_conductor_worker_v2 import run as replay_conductor
    from verify_high_rank_foundry_certificate import verify

    curves = read(CURVES)
    require(curves["status"] == "FROZEN_FOUNDRY_LEDGER_SELECTION_V2", "invalid curve snapshot")
    seen = []
    with tempfile.TemporaryDirectory(prefix="foundry-ledger-v2-replay-") as directory:
        folder = Path(directory)
        for identifier, entry in sorted(curves["records"].items()):
            record = entry["result"]
            v1.validate_record(record)
            require(record["rank_lower_bound"] >= curves["threshold"], "below publication threshold")
            require(sha(ROOT / entry["origin"]) == entry["origin_sha256"], "selected source changed")
            require(read(ROOT / entry["origin"]) == record, "embedded source differs")
            path = folder / "certificate.json"
            atomic(path, record)
            verify(path)
            model = record["packet"]["curve"]
            j = v1.bucket(model)
            matches = [row for row in curves["baseline"]
                       if v1.bucket(row["ainvs"]) == j and cert.isomorphic(model, row["ainvs"])]
            if entry["previous_rank_lower_bound"] is None:
                require(not matches, "new curve duplicates baseline")
            else:
                require(len(matches) == 1 and matches[0]["id"] == identifier
                        and matches[0]["rank_lower_bound"] < record["rank_lower_bound"],
                        "invalid subgroup strengthening")
            require(not any(j == old_j and cert.isomorphic(model, old_model)
                            for old_j, old_model in seen), "duplicate selected curve")
            seen.append((j, model))
    require(curves["rank_counts"] == dict(Counter(
        str(e["result"]["rank_lower_bound"]) for e in curves["records"].values())),
        "rank counts differ")

    conductors = read(CONDUCTORS)
    require(conductors["curve_snapshot_sha256"] == sha(CURVES), "conductor cutoff curve binding differs")
    with tempfile.TemporaryDirectory(prefix="foundry-conductor-v1-replay-") as directory:
        folder = Path(directory)
        for identifier, entry in sorted(conductors["records"].items()):
            require(identifier in curves["records"], "conductor curve outside publication cutoff")
            require(sha(ROOT / entry["origin"]) == entry["origin_sha256"], "conductor result changed")
            require(read(ROOT / entry["origin"]) == entry["result"], "embedded conductor differs")
            atomic(folder / "input.json", entry["input"])
            atomic(folder / "certificate.json", entry["result"]["certificate"])
            require(sha(folder / "input.json") == entry["input_sha256"],
                    "embedded conductor input hash differs")
            replay_conductor(folder / "input.json", folder / "certificate.json", True)
    require(conductors["record_count"] == len(conductors["records"]), "conductor count differs")
    require(conductors["exact_count"] == sum(e["result"]["certificate"]["status"] == "EXACT"
                                              for e in conductors["records"].values()),
            "exact conductor count differs")
    require(conductors["unknown_count"] == conductors["record_count"] - conductors["exact_count"],
            "unknown conductor count differs")
    print("PASS_FOUNDRY_LEDGER_V2_REPLAY", len(seen), "curves;",
          curves["new_count"], "new;", curves["strengthened_count"], "strengthened", flush=True)
    print("PASS_FOUNDRY_CONDUCTOR_V1_REPLAY", conductors["record_count"], "records;",
          conductors["exact_count"], "exact", flush=True)


def index():
    curves = read(CURVES)
    conductors = read(CONDUCTORS)
    manifest = read(MANIFEST)
    claims = {row["id"]: row for row in read(ROOT / "MATH_STATUS.json")["entries"]}
    claim = claims[CLAIM]
    curve_rel = str(CURVES.relative_to(ROOT))
    conductor_rel = str(CONDUCTORS.relative_to(ROOT))
    require(claim["state"] == "proved"
            and curve_rel in claim["software_lock"]
            and conductor_rel in claim["software_lock"], "proved snapshot claim required")
    curve_entries = {row["id"]: row for row in manifest["curves"]}
    for identifier, entry in curves["records"].items():
        row = entry["result"]
        curve_entries[identifier] = {
            "id": identifier,
            "family": row["family"],
            "parameter": row["parameter"],
            "rank_lower_bound": row["rank_lower_bound"],
            "kind": "foundry",
            "source": curve_rel,
            "record": identifier,
            "claim": CLAIM,
        }
    conductor_entries = {row["id"]: row for row in manifest["conductors"]}
    for identifier in conductors["records"]:
        conductor_entries[identifier] = {
            "id": identifier,
            "source": conductor_rel,
            "claim": CLAIM,
            "kind": "checkpoint",
        }
    manifest["curves"] = list(curve_entries.values())
    manifest["conductors"] = list(conductor_entries.values())
    manifest["sources"][curve_rel] = sha(CURVES)
    manifest["sources"][conductor_rel] = sha(CONDUCTORS)
    atomic(MANIFEST, manifest)
    print("FOUNDRY_LEDGER_V2_INDEXED", len(curves["records"]), "curves;",
          len(conductors["records"]), "conductors", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["freeze", "check", "index"])
    args = parser.parse_args()
    {"freeze": freeze, "check": check, "index": index}[args.action]()
