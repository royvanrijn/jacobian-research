#!/usr/bin/env python3
"""Validate MATH_STATUS.json and render the public mathematical status page."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEX_PATH = ROOT / "MATH_STATUS.json"
STATUS_PATH = ROOT / "STATUS.md"
FIELDS = {
    "id", "kind", "state", "title", "scope", "canonical_source",
    "dependencies", "checker", "proof_type", "independent_replay",
    "formal_verification", "external_review", "artifact_hash",
    "software_lock", "supersedes", "replaced_by", "priority",
}
KINDS = {"theorem", "corollary", "example", "reproduction", "open_problem"}
STATES = {"proved", "partial", "open", "parked", "archived", "falsified"}
PROOF_TYPES = {
    "unclassified",
    "hybrid",
    "formal",
    "reproduction",
    "counterexample",
    "not-applicable",
}
PRIORITIES = {"core", "derived", "reference", "primary", "parked"}
HASH_PATTERN = re.compile(r"sha256:[0-9a-f]{64}")
CORE_ORDER = ["F1", "W1", "S1", "WB1", "C1", "B1", "P1", "M1", "D1", "F2", "R1", "R2", "IA1"]
ACTIVE_OPEN = {
    "OP-CR",
    "OP-LR-REES",
    "OP-LR-II",
    "OP-LR-NE",
    "OP-RITT",
    "OP-SUSP",
}


def load_index() -> dict:
    return json.loads(INDEX_PATH.read_text())


def validate_index(index: dict) -> None:
    assert index.get("schema_version") == 3, "unsupported status schema"
    assert index.get("authority") == "MATH_STATUS.json"
    entries = index.get("entries")
    assert isinstance(entries, list) and entries, "the status registry is empty"

    ids = [item.get("id") for item in entries]
    assert len(ids) == len(set(ids)), "status IDs must be unique"
    known = set(ids)
    for item in entries:
        item_id = item.get("id", "?")
        assert set(item) == FIELDS, f"{item_id}: unexpected schema"
        assert item["id"] and item["title"] and item["scope"]
        assert item["kind"] in KINDS, f"{item_id}: invalid kind"
        assert item["state"] in STATES, f"{item_id}: invalid state"
        assert item["proof_type"] in PROOF_TYPES, f"{item_id}: invalid proof type"
        for field in (
            "independent_replay",
            "formal_verification",
            "external_review",
        ):
            assert isinstance(item[field], bool), f"{item_id}: {field} must be boolean"
        assert item["priority"] in PRIORITIES, f"{item_id}: invalid priority"
        assert isinstance(item["dependencies"], list)
        assert isinstance(item["software_lock"], list)
        assert isinstance(item["supersedes"], list)
        assert isinstance(item["replaced_by"], list)
        if item["kind"] == "open_problem":
            assert item["state"] in {"open", "parked"}
            assert item["proof_type"] == "not-applicable"
            assert not item["independent_replay"]
            assert not item["formal_verification"]
            assert not item["external_review"]
        else:
            assert item["state"] not in {"open", "parked"}
            assert item["proof_type"] != "not-applicable"
        if item["state"] == "falsified":
            assert item["proof_type"] == "counterexample"
        for dependency in item["dependencies"]:
            assert dependency in known or dependency.startswith("external: "), (
                f"{item_id}: unresolved dependency {dependency}"
            )
        for replacement in item["replaced_by"]:
            assert replacement in known, f"{item_id}: unresolved replacement {replacement}"
        assert (ROOT / item["canonical_source"]).is_file(), (
            f"{item_id}: missing canonical source {item['canonical_source']}"
        )
        checker = item["checker"]
        assert checker is None or (ROOT / checker).is_file(), (
            f"{item_id}: missing checker {checker}"
        )
        artifact_hash = item["artifact_hash"]
        if checker is None:
            assert artifact_hash is None, (
                f"{item_id}: artifact hash requires a checker"
            )
            assert not item["independent_replay"], (
                f"{item_id}: independent replay requires a checker"
            )
            assert not item["formal_verification"], (
                f"{item_id}: formal verification requires a checker"
            )
        else:
            assert isinstance(artifact_hash, str) and HASH_PATTERN.fullmatch(
                artifact_hash
            ), f"{item_id}: invalid artifact hash"
            actual_hash = "sha256:" + hashlib.sha256(
                (ROOT / checker).read_bytes()
            ).hexdigest()
            assert artifact_hash == actual_hash, (
                f"{item_id}: stale artifact hash for {checker}"
            )
        for lock in item["software_lock"]:
            assert isinstance(lock, str) and (ROOT / lock).is_file(), (
                f"{item_id}: missing software lock {lock}"
            )

    core = {x["id"] for x in entries if x["priority"] == "core"}
    assert core == set(CORE_ORDER), "the canonical theorem backbone changed"
    active = {x["id"] for x in entries if x["kind"] == "open_problem" and x["state"] == "open"}
    assert active == ACTIVE_OPEN, "the primary continuation queue changed"


def _link(label: str, path: str) -> str:
    return f"[{label}]({path})"


def _items(values: list[str]) -> str:
    if not values:
        return "—"
    return ", ".join(f"`{v}`" if not v.startswith("external: ") else v for v in values)


def _evidence(item: dict) -> str:
    parts = [item["proof_type"]]
    if item["independent_replay"]:
        parts.append("independent replay")
    if item["formal_verification"]:
        parts.append("formal verification")
    if item["external_review"]:
        parts.append("external review")
    if item["artifact_hash"]:
        parts.append(f"`{item['artifact_hash'][:19]}…`")
    if item["software_lock"]:
        parts.append("locks: " + ", ".join(f"`{x}`" for x in item["software_lock"]))
    return "; ".join(parts)


def _table(lines: list[str], entries: list[dict], *, replacements: bool = False) -> None:
    tail = " | Replaced by" if replacements else ""
    lines.extend([
        f"| ID | Result | Scope | Source | Dependencies | Checker | Evidence{tail} |",
        f"|---|---|---|---|---|---|---{'|---' if replacements else ''}|",
    ])
    for item in entries:
        source = _link("source", item["canonical_source"])
        checker = _link("checker", item["checker"]) if item["checker"] else "—"
        extra = f" | {_items(item['replaced_by'])}" if replacements else ""
        lines.append(
            f"| {item['id']} | {item['title']} | {item['scope']} | {source} | "
            f"{_items(item['dependencies'])} | {checker} | {_evidence(item)}{extra} |"
        )
    lines.append("")


def render(index: dict) -> str:
    by_id = {x["id"]: x for x in index["entries"]}
    entries = index["entries"]
    lines = [
        "# Mathematical status",
        "",
        "<!-- Generated by scripts/render_status.py from MATH_STATUS.json; do not edit. -->",
        "",
        "[`MATH_STATUS.json`](MATH_STATUS.json) is the sole status authority. Canonical "
        "sources contain the proofs; this page records their scope, dependency role, "
        "proof classification, and separate assurance signals. A checker establishes "
        "reproducibility, not independent replay, formal verification, or external review.",
        "",
        "## Core theorem chain",
        "",
    ]
    _table(lines, [by_id[i] for i in CORE_ORDER])

    sections = [
        ("Falsified claims", [x for x in entries if x["state"] == "falsified"], True),
        ("Audited high-risk claims", [x for x in entries if x["kind"] == "theorem" and x["state"] == "partial" and x["priority"] == "reference"], False),
        ("Completed reference theorems", [x for x in entries if x["kind"] == "theorem" and x["state"] == "proved" and x["priority"] == "reference"], False),
        ("Derived corollaries", [x for x in entries if x["kind"] == "corollary" and x["state"] in {"proved", "partial"}], False),
        ("Examples and regressions", [x for x in entries if x["kind"] == "example" and x["state"] in {"proved", "partial"}], True),
        ("External reproductions", [x for x in entries if x["kind"] == "reproduction" and x["state"] in {"proved", "partial"}], False),
        ("Active open problems", [x for x in entries if x["kind"] == "open_problem" and x["state"] == "open"], False),
        ("Parked problems", [x for x in entries if x["kind"] == "open_problem" and x["state"] == "parked"], False),
    ]
    for heading, members, replacements in sections:
        lines.extend([f"## {heading}", ""])
        _table(lines, members, replacements=replacements)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="fail if STATUS.md is stale")
    args = parser.parse_args()
    index = load_index()
    validate_index(index)
    rendered = render(index)
    if args.check:
        assert STATUS_PATH.read_text() == rendered, (
            "STATUS.md is stale; run python3 scripts/render_status.py"
        )
    else:
        STATUS_PATH.write_text(rendered)


if __name__ == "__main__":
    main()
