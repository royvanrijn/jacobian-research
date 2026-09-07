#!/usr/bin/env python3
"""Audit the active elliptic artifact catalogue and historical archive."""

from __future__ import annotations

import csv
import base64
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
ACTIVE = ROOT / "artifacts/generated-results/elliptic-curves"
CATALOG = ACTIVE / "CATALOG.tsv"
ARCHIVE = ROOT / "archive/elliptic-curves"
ARCHIVE_MANIFEST = ARCHIVE / "MANIFEST.tsv"
SNAPSHOT = ARCHIVE / "artifacts/snapshots/pre-cleanup-2026-08-24"
SNAPSHOT_MANIFEST = SNAPSHOT / "MANIFEST.tsv"
ALLOWED_EVIDENCE = {
    "theorem-certificate",
    "exact-rank-certificate",
    "exact-lower-bound-certificate",
    "conditional-bound",
    "exact-computation",
    "bounded-experiment",
    "source-transcription",
    "reproducibility-fixture",
    "search-plan",
    "partial-reproduction",
}
REFERENCE = re.compile(
    r"^(?:" + re.escape(str(ROOT)) + r"/)?"
    r"(?:artifacts/|archive/|elliptic-curves/)\S+"
    r"\.(?:json(?:\.gz|l)?|py|sage|gp|cpp|m)$"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fail(message: str) -> None:
    raise AssertionError(message)


def load_json(path: Path) -> Any:
    if path.name.endswith(".json.gz"):
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    return json.loads(path.read_text())


def strings(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            if isinstance(key, str):
                yield key
            yield from strings(child)
    elif isinstance(value, list):
        for child in value:
            yield from strings(child)
    elif isinstance(value, str):
        yield value


def status_ids() -> set[str]:
    data = json.loads((ROOT / "MATH_STATUS.json").read_text())
    found: set[str] = set()

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            identifier = value.get("id")
            if isinstance(identifier, str):
                found.add(identifier)
            for child in value.values():
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)

    visit(data)
    return found


def audit_catalog() -> tuple[int, dict[str, Any]]:
    with CATALOG.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    expected_fields = {
        "filename",
        "sha256",
        "evidence",
        "status_ids",
        "description",
    }
    if not rows or set(rows[0]) != expected_fields:
        fail("CATALOG.tsv has the wrong or empty schema")
    filenames = [row["filename"] for row in rows]
    if len(filenames) != len(set(filenames)):
        fail("CATALOG.tsv contains a duplicate filename")
    actual = {
        path.name
        for path in ACTIVE.iterdir()
        if path.name.endswith(".json") or path.name.endswith(".json.gz")
    }
    if set(filenames) != actual:
        fail(
            "catalog/file mismatch: "
            f"missing={sorted(actual - set(filenames))}, "
            f"extra={sorted(set(filenames) - actual)}"
        )
    known_status = status_ids()
    loaded: dict[str, Any] = {}
    for row in rows:
        artifact = ACTIVE / row["filename"]
        if row["evidence"] not in ALLOWED_EVIDENCE:
            fail(f"unknown evidence label for {artifact.name}: {row['evidence']}")
        if sha256(artifact) != row["sha256"]:
            fail(f"catalog hash is stale for {artifact.name}")
        identifiers = [] if row["status_ids"] == "-" else row["status_ids"].split(",")
        unknown = set(identifiers) - known_status
        if unknown:
            fail(f"unknown status id(s) for {artifact.name}: {sorted(unknown)}")
        loaded[artifact.name] = load_json(artifact)
    return len(rows), loaded


def audit_manifest(path: Path, *, archived_paths: bool) -> int:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    filenames = []
    for row in rows:
        target_text = row["archived_path"] if archived_paths else row["filename"]
        target = ROOT / target_text if archived_paths else path.parent / target_text
        filenames.append(target_text)
        if not target.is_file():
            fail(f"manifest target is missing: {target_text}")
        if sha256(target) != row["sha256"]:
            fail(f"manifest hash is stale: {target_text}")
        if archived_paths and (ROOT / row["original_path"]).exists():
            fail(f"archived path was repopulated: {row['original_path']}")
    if len(filenames) != len(set(filenames)):
        fail(f"duplicate target in {path.relative_to(ROOT)}")
    return len(rows)


def audit_references(loaded: dict[str, Any]) -> int:
    # Self-contained sensitivity bundles retain ignored replay inputs by
    # content. Those references remain valid after a clean checkout removes
    # the original local checkpoints; verify the embedded bytes first.
    embedded: set[str] = set()
    for data in loaded.values():
        if not isinstance(data, dict) or data.get("schema") != "elliptic-curves.mw16-sensitivity-evidence.v1":
            continue
        for name, entry in data["files"].items():
            if isinstance(entry.get("text"), str):
                content = entry["text"].encode()
            elif isinstance(entry.get("base64"), str):
                content = base64.b64decode(entry["base64"], validate=True)
            else:
                fail(f"missing embedded content: {name}")
            if hashlib.sha256(content).hexdigest() != entry["sha256"]:
                fail(f"embedded content checksum mismatch: {name}")
            embedded.add(name)
    checked: set[str] = set()
    for artifact_name, data in loaded.items():
        for value in strings(data):
            if "archive/elliptic-curves/archive/" in value:
                fail(f"duplicated archive prefix in {artifact_name}: {value}")
            if not REFERENCE.fullmatch(value):
                continue
            candidate = Path(value)
            target = candidate if candidate.is_absolute() else ROOT / candidate
            if not target.exists() and str(target.relative_to(ROOT)) not in embedded:
                fail(f"dangling repository reference in {artifact_name}: {value}")
            checked.add(value)
    return len(checked)


def audit_flat_directory() -> None:
    forbidden_prefixes = ("elliptic_", "fermigier_", "icarm_", "newfamily_")
    leftovers = [
        path.name
        for path in (ROOT / "artifacts/generated-results").iterdir()
        if path.is_file() and path.name.startswith(forbidden_prefixes)
    ]
    if leftovers:
        fail(f"elliptic artifacts remain in the flat generated directory: {leftovers}")


def main() -> None:
    try:
        active_count, loaded = audit_catalog()
        archive_count = audit_manifest(ARCHIVE_MANIFEST, archived_paths=True)
        snapshot_count = audit_manifest(SNAPSHOT_MANIFEST, archived_paths=False)
        reference_count = audit_references(loaded)
        audit_flat_directory()
    except (AssertionError, KeyError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ELLIPTIC_ARTIFACT_AUDIT_FAIL: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    print(
        "ELLIPTIC_ARTIFACT_AUDIT_PASS "
        f"active={active_count} archived={archive_count} "
        f"snapshot={snapshot_count} references={reference_count}"
    )


if __name__ == "__main__":
    main()
