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
REQUIRED_FIELDS = {
    "id", "kind", "state", "title", "scope", "canonical_source",
    "dependencies", "checker", "proof_type", "independent_replay",
    "formal_verification", "external_review", "artifact_hash",
    "software_lock", "supersedes", "closes_problems", "narrows_problems",
    "consumers", "invalidates_assumptions", "replaced_by", "priority",
}
OPTIONAL_FIELDS = {
    "external_formal_certificates",
    "forbidden_attack_classes",
    "supersedes_notes",
}
KINDS = {"theorem", "corollary", "example", "reproduction", "open_problem"}
STATES = {"proved", "partial", "open", "parked", "archived", "falsified"}
PROOF_TYPES = {
    "unclassified",
    "hybrid",
    "formal",
    "reproduction",
    "counterexample",
    "exact_symbolic",
    "exact_modular",
    "not-applicable",
}
PRIORITIES = {"core", "derived", "reference", "primary", "parked"}
HASH_PATTERN = re.compile(r"sha256:[0-9a-f]{64}")
CORE_ORDER = [
    "F1", "W1", "S1", "WB1", "C1", "B1", "P1", "M1", "D1", "F2",
    "R1", "R2", "AMS1", "IA1",
]
ACTIVE_OPEN = {
    "GVC2ALL",
    "OP-CR",
    "OP-CCDM",
    "OP-GVC2-RP",
    "OP-GVC3-MIN",
    "OP-HC4-D5",
    "OP-HC4-SCHUR",
    "OP-KCOMP",
    "OP-KDESC",
    "OP-KMOD",
    "OP-LR-NE",
    "OP-RITT",
    "OP-SIC2-B33",
    "OP-SUSP",
}
GMC2_RETAINED_IDS = {
    "G2F",
    "G2C",
    "G2D",
    "G2E",
    "G2N",
    "G2Q",
    "G2S",
    "G2R",
}


def load_index() -> dict:
    return json.loads(INDEX_PATH.read_text())


def consumer_marker(item: dict) -> str:
    payload = "\0".join((item["state"], item["title"], item["scope"]))
    digest = hashlib.sha256(payload.encode()).hexdigest()[:16]
    return f"<!-- status-consumer: {item['id']} {digest} -->"


def _assert_acyclic(edges: dict[str, list[str]], label: str) -> None:
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node: str) -> None:
        assert node not in visiting, f"{label} contains a cycle through {node}"
        if node in visited:
            return
        visiting.add(node)
        for target in edges[node]:
            visit(target)
        visiting.remove(node)
        visited.add(node)

    for node in edges:
        visit(node)


def validate_index(index: dict) -> None:
    assert index.get("schema_version") == 6, "unsupported status schema"
    assert index.get("authority") == "MATH_STATUS.json"
    entries = index.get("entries")
    assert isinstance(entries, list) and entries, "the status registry is empty"

    ids = [item.get("id") for item in entries]
    assert len(ids) == len(set(ids)), "status IDs must be unique"
    known = set(ids)
    by_id = {item["id"]: item for item in entries}
    update_fields = (
        "supersedes",
        "closes_problems",
        "narrows_problems",
        "consumers",
        "invalidates_assumptions",
    )
    for item in entries:
        item_id = item.get("id", "?")
        assert REQUIRED_FIELDS <= set(item) <= REQUIRED_FIELDS | OPTIONAL_FIELDS, (
            f"{item_id}: unexpected schema"
        )
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
        certificates = item.get("external_formal_certificates", [])
        assert isinstance(certificates, list)
        for certificate in certificates:
            assert set(certificate) == {
                "name", "url", "prover", "scope", "refereed",
            }, f"{item_id}: invalid external formal certificate"
            assert all(
                isinstance(certificate[field], str) and certificate[field]
                for field in ("name", "url", "prover", "scope")
            ), f"{item_id}: incomplete external formal certificate"
            assert certificate["url"].startswith("https://"), (
                f"{item_id}: external formal certificate requires an HTTPS URL"
            )
            assert isinstance(certificate["refereed"], bool)
        if certificates:
            assert item["formal_verification"], (
                f"{item_id}: external formal certificate requires formal verification"
            )
        if item["external_review"]:
            assert any(certificate["refereed"] for certificate in certificates), (
                f"{item_id}: external review requires identified refereed evidence"
            )
        assert item["priority"] in PRIORITIES, f"{item_id}: invalid priority"
        assert isinstance(item["dependencies"], list)
        assert isinstance(item["software_lock"], list)
        assert isinstance(item["replaced_by"], list)
        for field in update_fields:
            assert isinstance(item[field], list), f"{item_id}: {field} must be a list"
            assert len(item[field]) == len(set(item[field])), (
                f"{item_id}: duplicate {field} edge"
            )
            assert all(
                isinstance(target, str) and target for target in item[field]
            ), f"{item_id}: invalid {field} edge"
        supersedes_notes = item.get("supersedes_notes", [])
        assert isinstance(supersedes_notes, list)
        assert all(
            isinstance(note, str) and note for note in supersedes_notes
        ), f"{item_id}: invalid supersedes note"
        forbidden_attacks = item.get("forbidden_attack_classes", [])
        assert isinstance(forbidden_attacks, list)
        if forbidden_attacks:
            assert item["kind"] == "open_problem", (
                f"{item_id}: forbidden attack classes belong only to open problems"
            )
        for attack in forbidden_attacks:
            assert set(attack) == {"attack", "reason", "witnesses"}, (
                f"{item_id}: invalid forbidden attack class"
            )
            assert isinstance(attack["attack"], str) and attack["attack"]
            assert isinstance(attack["reason"], str) and attack["reason"]
            assert isinstance(attack["witnesses"], list) and attack["witnesses"]
            assert len(attack["witnesses"]) == len(set(attack["witnesses"]))
            for witness in attack["witnesses"]:
                assert witness in known, (
                    f"{item_id}: unresolved forbidden-attack witness {witness}"
                )
                assert witness in item["dependencies"], (
                    f"{item_id}: forbidden-attack witness {witness} "
                    "must also be a dependency"
                )
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
        for field in (
            "supersedes",
            "closes_problems",
            "narrows_problems",
            "invalidates_assumptions",
        ):
            for target in item[field]:
                assert target in known, f"{item_id}: unresolved {field} target {target}"
                assert target != item_id, f"{item_id}: self-referential {field} edge"
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

    for item in entries:
        item_id = item["id"]
        for target_id in item["supersedes"]:
            target = by_id[target_id]
            assert item_id in target["replaced_by"], (
                f"{item_id}: supersedes {target_id}, but the target does not name it "
                "in replaced_by"
            )
        for target_id in item["closes_problems"]:
            assert item["state"] == "proved", (
                f"{item_id}: only a proved entry may close a problem"
            )
            target = by_id[target_id]
            assert target["kind"] == "open_problem", (
                f"{item_id}: closes non-problem {target_id}"
            )
            assert target["state"] != "open", (
                f"{item_id}: closed problem {target_id} remains active"
            )
            assert item_id in target["dependencies"], (
                f"{item_id}: closed problem {target_id} does not consume the result"
            )
            assert item_id in target["replaced_by"], (
                f"{item_id}: closed problem {target_id} does not name its replacement"
            )
        for target_id in item["narrows_problems"]:
            assert item["state"] == "proved", (
                f"{item_id}: only a proved entry may narrow a problem"
            )
            target = by_id[target_id]
            assert target["kind"] == "open_problem" and target["state"] == "open", (
                f"{item_id}: narrowed target {target_id} is not an active problem"
            )
            assert item_id in target["dependencies"], (
                f"{item_id}: narrowed problem {target_id} does not consume the result"
            )
            assert item_id in target["scope"], (
                f"{item_id}: narrowed problem {target_id} does not mention the result "
                "in its scope"
            )
        for target_id in item["invalidates_assumptions"]:
            assert item["state"] == "proved", (
                f"{item_id}: only a proved entry may invalidate an assumption"
            )
            target = by_id[target_id]
            assert target["state"] == "falsified", (
                f"{item_id}: invalidated assumption {target_id} is not falsified"
            )
            assert item_id in target["replaced_by"], (
                f"{item_id}: invalidated assumption {target_id} does not name its "
                "replacement"
            )
        for consumer in item["consumers"]:
            if consumer in known:
                target = by_id[consumer]
                assert item_id in target["dependencies"], (
                    f"{item_id}: declared consumer {consumer} does not depend on it"
                )
                assert item_id in target["scope"], (
                    f"{item_id}: declared consumer {consumer} does not acknowledge it "
                    "in its scope"
                )
                continue
            consumer_path = (ROOT / consumer).resolve()
            assert ROOT.resolve() in consumer_path.parents, (
                f"{item_id}: consumer escapes the repository: {consumer}"
            )
            assert consumer_path.is_file(), (
                f"{item_id}: missing document consumer {consumer}"
            )
            marker = consumer_marker(item)
            assert marker in consumer_path.read_text(), (
                f"{item_id}: stale document consumer {consumer}; review it and update "
                f"the marker to {marker}"
            )

    for item in entries:
        for replacement_id in item["replaced_by"]:
            replacement = by_id[replacement_id]
            reciprocal = (
                item["id"] in replacement["supersedes"]
                or item["id"] in replacement["closes_problems"]
                or item["id"] in replacement["invalidates_assumptions"]
            )
            assert reciprocal, (
                f"{item['id']}: replaced_by {replacement_id} lacks a reciprocal "
                "machine-readable update edge"
            )

    _assert_acyclic(
        {item["id"]: item["supersedes"] for item in entries},
        "supersedes graph",
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
    certificates = item.get("external_formal_certificates", [])
    if certificates:
        links = ", ".join(
            _link(certificate["name"], certificate["url"])
            for certificate in certificates
        )
        parts.append(f"external formal certificates: {links}")
    if item["artifact_hash"]:
        parts.append(f"`{item['artifact_hash'][:19]}…`")
    if item["software_lock"]:
        parts.append("locks: " + ", ".join(f"`{x}`" for x in item["software_lock"]))
    updates = []
    for field in (
        "supersedes",
        "closes_problems",
        "narrows_problems",
        "consumers",
        "invalidates_assumptions",
    ):
        if item[field]:
            updates.append(f"{field.replace('_', ' ')} {_items(item[field])}")
    if updates:
        parts.append("updates: " + "; ".join(updates))
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


def _forbidden_attack_sections(lines: list[str], entries: list[dict]) -> None:
    for item in entries:
        attacks = item.get("forbidden_attack_classes", [])
        if not attacks:
            continue
        lines.extend([
            f"### {item['id']}: Forbidden attack classes",
            "",
            "The following routes are obsolete for this programme. The cited results "
            "either remove the proposed defect or provide explicit countermodels:",
            "",
        ])
        for attack in attacks:
            witnesses = _items(attack["witnesses"])
            lines.append(
                f"- **{attack['attack']}.** {attack['reason']} "
                f"Witnesses: {witnesses}."
            )
        lines.append("")


def render(index: dict) -> str:
    by_id = {x["id"]: x for x in index["entries"]}
    entries = index["entries"]
    gmc2_retained = [x for x in entries if x["id"] in GMC2_RETAINED_IDS]
    lines = [
        "# Mathematical status",
        "",
        "<!-- Generated by scripts/render_status.py from MATH_STATUS.json; do not edit. -->",
        "",
        "[`MATH_STATUS.json`](MATH_STATUS.json) is the sole status authority. Canonical "
        "sources contain the proofs; this page records their scope, dependency role, "
        "proof classification, and separate assurance signals. A checker establishes "
        "reproducibility, not independent replay, formal verification, or external review. "
        "External review includes identified refereeing by a formal-proof archive as well "
        "as conventional publication review.",
        "",
        "## Core theorem chain",
        "",
    ]
    _table(lines, [by_id[i] for i in CORE_ORDER])

    sections = [
        ("Falsified claims", [x for x in entries if x["state"] == "falsified"], True),
        ("Audited high-risk claims", [x for x in entries if x["kind"] == "theorem" and x["state"] == "partial" and x["priority"] == "reference"], False),
        ("Completed reference theorems", [x for x in entries if x["kind"] == "theorem" and x["state"] == "proved" and x["priority"] == "reference"], False),
        ("Superseded proof route / retained refinements", gmc2_retained, False),
        ("Derived corollaries", [x for x in entries if x["kind"] == "corollary" and x["state"] in {"proved", "partial"} and x["id"] not in GMC2_RETAINED_IDS], False),
        ("Examples and regressions", [x for x in entries if x["kind"] == "example" and x["state"] in {"proved", "partial"}], True),
        ("External reproductions", [x for x in entries if x["kind"] == "reproduction" and x["state"] in {"proved", "partial"}], False),
        ("Active open problems", [x for x in entries if x["kind"] == "open_problem" and x["state"] == "open"], False),
        ("Parked problems", [x for x in entries if x["kind"] == "open_problem" and x["state"] == "parked"], False),
    ]
    for heading, members, replacements in sections:
        lines.extend([f"## {heading}", ""])
        _table(lines, members, replacements=replacements)
        _forbidden_attack_sections(lines, members)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="fail if STATUS.md is stale")
    parser.add_argument(
        "--consumer-marker",
        metavar="ID",
        help="print the current document-consumer marker for one status entry",
    )
    args = parser.parse_args()
    index = load_index()
    if args.consumer_marker:
        by_id = {item["id"]: item for item in index["entries"]}
        assert args.consumer_marker in by_id, (
            f"unknown status entry {args.consumer_marker}"
        )
        print(consumer_marker(by_id[args.consumer_marker]))
        return
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
