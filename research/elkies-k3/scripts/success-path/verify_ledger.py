#!/usr/bin/env python3
"""Verify locked scripts and exact artifacts in the H3 success-path ledger."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
LEDGER = Path(__file__).with_name("ledger.json")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def all_records(ledger: dict) -> list[dict]:
    return ledger["route"] + ledger["shortcut_audits"]


def main() -> None:
    ledger = json.loads(LEDGER.read_text())
    assert ledger["schema"] == "elkies-k3.h3-equation-success-path.v1"
    records = all_records(ledger)
    ids = [record["id"] for record in records]
    assert len(ids) == len(set(ids)), "duplicate success-path stage id"

    script_count = 0
    artifact_count = 0
    for record in records:
        script = record.get("script")
        if script:
            path = ROOT / script["path"]
            assert path.is_file(), f"missing locked script: {path}"
            actual = digest(path)
            assert actual == script["sha256"], (
                f"stale script lock for {record['id']}: "
                f"expected {script['sha256']}, got {actual}"
            )
            script_count += 1

        artifact = record.get("artifact")
        if artifact:
            path = ROOT / artifact["path"]
            assert path.is_file(), f"missing exact artifact: {path}"
            actual = digest(path)
            assert actual == artifact["sha256"], (
                f"stale artifact for {record['id']}: "
                f"expected {artifact['sha256']}, got {actual}"
            )
            payload = json.loads(path.read_text())
            assert payload.get("status") == artifact["status"], (
                f"wrong artifact status for {record['id']}"
            )
            if artifact.get("schema"):
                assert payload.get("schema") == artifact["schema"], (
                    f"wrong artifact schema for {record['id']}"
                )
            artifact_count += 1

    print(
        "SUCCESS_PATH_LEDGER_RESULT|"
        f"route_stages={len(ledger['route'])}|"
        f"shortcut_audits={len(ledger['shortcut_audits'])}|"
        f"scripts={script_count}|artifacts={artifact_count}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()

