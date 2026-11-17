#!/usr/bin/env python3
"""Build a read-only per-entry trace of the mathematical-status registry.

This complements ``render_status.py``.  The renderer remains authoritative for
schema and graph validity; this script exposes the evidence topology needed by
retrospective review without pretending to re-prove mathematical claims.  It
records canonical notes, checkers and their current hashes, generated and
local artifacts, missing locks, dependencies, inverse dependency consumers,
declared consumers, and update/supersession edges for every status row.  It
also preserves the exact claim boundary, proof and assurance types, forbidden
attack classes, and status IDs mentioned in the scope so semantic review does
not collapse to a path-existence audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "MATH_STATUS.json"
HASH_PATTERN = re.compile(r"^sha256:([0-9a-f]{64})$")
BOUNDARY_WORDS = re.compile(
    r"\b(bounded|candidate|cannot|conditional|does not|excludes?|finite|not|"
    r"obstruction|only|open|pending|remain(?:s|ing)?|unknown|withdrawn|without)\b",
    re.IGNORECASE,
)
FINITE_EXHAUSTION_WORDS = re.compile(
    r"\b(bounded|degree|displayed|enumerat(?:e|ed|ion)|finite|height|list|"
    r"modulo|range|scope|source|tested|through|trace|up to|within)\b",
    re.IGNORECASE,
)
STATUS_REFERENCE_TOKEN = re.compile(
    r"(?<![A-Za-z0-9-])[A-Z][A-Z0-9-]{1,}(?![A-Za-z0-9-])"
)
CHECKER_AUDIT_PATTERNS = {
    "early_loop_exit": re.compile(r"^\s*break\b"),
    "beam_or_first_hit": re.compile(r"\b(?:beam[_-]?(?:width|size)|first[_-]?hit)\b", re.I),
    "bounded_or_limit": re.compile(r"\b(?:bounded|limit|max(?:imum)?[_-]?(?:height|norm|support|trials?|candidates?|steps?))\b", re.I),
    "timeout": re.compile(r"\btimeout\b", re.I),
    "saturation": re.compile(r"\bsaturat(?:e|ed|es|ing|ion)\b", re.I),
    "denominator_handling": re.compile(r"\b(?:denominators?|clear_denominators?)\b", re.I),
}


def digest(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def repo_path(value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        return None
    return path


def artifact_class(path: str) -> str:
    if path.startswith("artifacts/generated-results/"):
        return "generated"
    if path.startswith("artifacts/local/"):
        return "local"
    return "source"


def checker_audit_markers(path: Path) -> dict[str, dict[str, Any]]:
    """Return navigation markers only; their presence is not a defect."""

    lines = path.read_text(errors="replace").splitlines()
    markers: dict[str, dict[str, Any]] = {}
    for label, pattern in CHECKER_AUDIT_PATTERNS.items():
        matches = [
            line_number
            for line_number, line in enumerate(lines, 1)
            if pattern.search(line)
        ]
        if matches:
            markers[label] = {
                "count": len(matches),
                "lines": matches[:64],
                "line_list_truncated": len(matches) > 64,
            }
    return markers


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the deterministic full JSON ledger (normally under artifacts/local/)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare --output with the current deterministic ledger",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="return nonzero for broken paths or checker-source hash drift",
    )
    parser.add_argument(
        "--show-findings",
        action="store_true",
        help="print all finding IDs after the compact summary",
    )
    args = parser.parse_args()
    if args.check and args.output is None:
        parser.error("--check requires --output")

    registry = json.loads(INDEX.read_text())
    entries: list[dict[str, Any]] = registry["entries"]
    by_id = {entry["id"]: entry for entry in entries}
    known = set(by_id)
    checker_users: dict[str, list[str]] = defaultdict(list)
    for entry in entries:
        if entry["checker"] is not None:
            checker_users[entry["checker"]].append(entry["id"])
    inverse_dependencies: dict[str, set[str]] = defaultdict(set)
    for entry in entries:
        for dependency in entry["dependencies"]:
            if dependency in known:
                inverse_dependencies[dependency].add(entry["id"])

    finding_ids: dict[str, list[str]] = defaultdict(list)
    rows: list[dict[str, Any]] = []
    for entry in entries:
        entry_id = entry["id"]
        canonical_relative = repo_path(entry["canonical_source"])
        canonical_exists = bool(
            canonical_relative is not None and (ROOT / canonical_relative).is_file()
        )
        canonical_actual_hash = (
            digest(ROOT / canonical_relative)
            if canonical_exists and canonical_relative is not None
            else None
        )
        checker_relative = repo_path(entry["checker"])
        checker_exists = entry["checker"] is None or bool(
            checker_relative is not None and (ROOT / checker_relative).is_file()
        )
        checker_actual_hash = (
            digest(ROOT / checker_relative)
            if entry["checker"] is not None and checker_exists and checker_relative
            else None
        )
        checker_markers = (
            checker_audit_markers(ROOT / checker_relative)
            if entry["checker"] is not None and checker_exists and checker_relative
            else {}
        )
        checker_hash_match = (
            entry["artifact_hash"] is None
            if entry["checker"] is None
            else checker_actual_hash == entry["artifact_hash"]
        )

        locks = []
        for value in entry["software_lock"]:
            relative = repo_path(value)
            exists = bool(relative is not None and (ROOT / relative).is_file())
            locks.append(
                {
                    "path": value,
                    "class": artifact_class(value),
                    "exists": exists,
                }
            )
        generated = [lock["path"] for lock in locks if lock["class"] == "generated"]
        local = [lock["path"] for lock in locks if lock["class"] == "local"]
        missing_generated = [
            lock["path"] for lock in locks
            if lock["class"] == "generated" and not lock["exists"]
        ]
        missing_local = [
            lock["path"] for lock in locks
            if lock["class"] == "local" and not lock["exists"]
        ]
        missing_source = [
            lock["path"] for lock in locks
            if lock["class"] == "source" and not lock["exists"]
        ]

        internal_dependencies = [
            dependency for dependency in entry["dependencies"] if dependency in known
        ]
        external_dependencies = [
            dependency
            for dependency in entry["dependencies"]
            if dependency.startswith("external: ")
        ]
        declared_status_consumers = [
            consumer for consumer in entry["consumers"] if consumer in known
        ]
        document_consumers = [
            consumer for consumer in entry["consumers"] if consumer not in known
        ]
        inferred_consumers = sorted(inverse_dependencies.get(entry_id, set()))
        scope_status_references = sorted(
            {
                token
                for token in STATUS_REFERENCE_TOKEN.findall(entry["scope"])
                if token in known and token != entry_id
            }
        )
        declared_status_edges = set(internal_dependencies)
        for field in (
            "supersedes",
            "replaced_by",
            "closes_problems",
            "narrows_problems",
            "invalidates_assumptions",
        ):
            declared_status_edges.update(entry[field])

        findings: list[str] = []
        if not canonical_exists:
            findings.append("missing_canonical_source")
        if not checker_exists:
            findings.append("missing_checker")
        if not checker_hash_match:
            findings.append("checker_hash_changed")
        if missing_source:
            findings.append("missing_source_lock")
        if missing_generated:
            findings.append("missing_generated_lock")
        if missing_local:
            findings.append("missing_local_lock")
        if entry["state"] == "partial" and not BOUNDARY_WORDS.search(entry["scope"]):
            findings.append("partial_scope_boundary_language_review")
        if (
            entry["state"] == "proved"
            and "exhaustive" in entry["scope"].lower()
            and not FINITE_EXHAUSTION_WORDS.search(entry["scope"])
        ):
            findings.append("unqualified_exhaustive_language_review")
        if (
            entry["kind"] == "open_problem"
            and not entry.get("forbidden_attack_classes")
            and not entry.get("forbidden_attack_review")
        ):
            findings.append("open_problem_without_forbidden_attack_metadata")
        for finding in findings:
            finding_ids[finding].append(entry_id)

        rows.append(
            {
                "id": entry_id,
                "title": entry["title"],
                "scope": entry["scope"],
                "kind": entry["kind"],
                "state": entry["state"],
                "priority": entry["priority"],
                "proof": {
                    "type": entry["proof_type"],
                    "independent_replay": entry["independent_replay"],
                    "formal_verification": entry["formal_verification"],
                    "external_review": entry["external_review"],
                    "external_formal_certificates": entry.get(
                        "external_formal_certificates", []
                    ),
                },
                "canonical_source": {
                    "path": entry["canonical_source"],
                    "exists": canonical_exists,
                    "sha256": canonical_actual_hash,
                },
                "checker": {
                    "path": entry["checker"],
                    "exists": checker_exists,
                    "recorded_sha256": entry["artifact_hash"],
                    "actual_sha256": checker_actual_hash,
                    "hash_matches": checker_hash_match,
                    "shared_status_entries": sorted(
                        checker_users.get(entry["checker"], [])
                    ),
                    "source_audit_markers": checker_markers,
                },
                "software_locks": locks,
                "artifacts": {
                    "generated": generated,
                    "local": local,
                    "missing_generated": missing_generated,
                    "missing_local": missing_local,
                },
                "dependencies": {
                    "internal": internal_dependencies,
                    "external": external_dependencies,
                },
                "scope_status_references": {
                    "all": scope_status_references,
                    "declared_as_edges": sorted(
                        set(scope_status_references) & declared_status_edges
                    ),
                    "mentioned_without_status_edge": sorted(
                        set(scope_status_references) - declared_status_edges
                    ),
                },
                "consumers": {
                    "declared_status": declared_status_consumers,
                    "documents": document_consumers,
                    "inferred_from_dependencies": inferred_consumers,
                    "inferred_not_declared": sorted(
                        set(inferred_consumers) - set(declared_status_consumers)
                    ),
                },
                "updates": {
                    field: entry[field]
                    for field in (
                        "supersedes",
                        "replaced_by",
                        "closes_problems",
                        "narrows_problems",
                        "invalidates_assumptions",
                    )
                },
                "supersedes_notes": entry.get("supersedes_notes", []),
                "forbidden_attack_classes": entry.get(
                    "forbidden_attack_classes", []
                ),
                "forbidden_attack_review": entry.get("forbidden_attack_review"),
                "findings": findings,
            }
        )

    finding_counts = {key: len(value) for key, value in sorted(finding_ids.items())}
    blocking_findings = {
        key: sorted(finding_ids.get(key, []))
        for key in (
            "missing_canonical_source",
            "missing_checker",
            "missing_source_lock",
            "checker_hash_changed",
        )
        if finding_ids.get(key)
    }
    payload = {
        "schema": "discovery-retrospective-registry-audit-v1",
        "authority": "MATH_STATUS.json",
        "authority_sha256": digest(INDEX),
        "entry_count": len(entries),
        "counts": {
            "state": dict(sorted(Counter(entry["state"] for entry in entries).items())),
            "kind": dict(sorted(Counter(entry["kind"] for entry in entries).items())),
            "priority": dict(
                sorted(Counter(entry["priority"] for entry in entries).items())
            ),
            "internal_dependency_edges": sum(
                len(row["dependencies"]["internal"]) for row in rows
            ),
            "external_dependency_edges": sum(
                len(row["dependencies"]["external"]) for row in rows
            ),
            "generated_artifact_locks": sum(
                len(row["artifacts"]["generated"]) for row in rows
            ),
            "local_artifact_locks": sum(
                len(row["artifacts"]["local"]) for row in rows
            ),
            "checkers_also_listed_in_software_lock": sum(
                entry["checker"] is not None
                and entry["checker"] in entry["software_lock"]
                for entry in entries
            ),
            "unique_checker_paths": len(checker_users),
            "shared_checker_paths": sum(
                len(users) > 1 for users in checker_users.values()
            ),
            "entries_with_early_loop_exit_markers": sum(
                "early_loop_exit" in row["checker"]["source_audit_markers"]
                for row in rows
            ),
        },
        "finding_counts": finding_counts,
        "finding_ids": {key: sorted(value) for key, value in sorted(finding_ids.items())},
        "blocking_findings": blocking_findings,
        "entries": rows,
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    if args.output:
        output = args.output if args.output.is_absolute() else ROOT / args.output
        if args.check:
            if not output.is_file() or output.read_text() != encoded:
                raise SystemExit(f"stale discovery registry audit: {output}")
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(encoded)

    print(
        "DISCOVERYREGISTRY|entries={}|proved={}|partial={}|open={}|parked={}|"
        "archived={}|falsified={}|dependency_edges={}|generated_locks={}|"
        "local_locks={}|blocking={}|status=PASS_REPORT".format(
            len(entries),
            payload["counts"]["state"].get("proved", 0),
            payload["counts"]["state"].get("partial", 0),
            payload["counts"]["state"].get("open", 0),
            payload["counts"]["state"].get("parked", 0),
            payload["counts"]["state"].get("archived", 0),
            payload["counts"]["state"].get("falsified", 0),
            payload["counts"]["internal_dependency_edges"],
            payload["counts"]["generated_artifact_locks"],
            payload["counts"]["local_artifact_locks"],
            sum(len(value) for value in blocking_findings.values()),
        )
    )
    for key, count in finding_counts.items():
        print(f"FINDING|{key}|count={count}")
        if args.show_findings:
            print("  " + ", ".join(sorted(finding_ids[key])))

    return 1 if args.strict and blocking_findings else 0


if __name__ == "__main__":
    raise SystemExit(main())
