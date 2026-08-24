#!/usr/bin/env python3
"""Validate and render the Elkies--K3 mathematical process ledger.

This is a dependency-free navigation/provenance checker.  It does not verify
the underlying Sage, Magma, or mathematical certificates and it does not
promote any claim beyond MATH_STATUS.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LEDGER = ROOT / "elkies-k3/data/process/elkies_k3_process_ledger.json"
DEFAULT_ATLAS = ROOT / "elkies-k3/ELKIES_K3_PROCESS_ATLAS.md"
LOCAL_ROOT = ROOT / "artifacts/local/elkies-k3"

CHRONOLOGY_BEGIN = "<!-- BEGIN GENERATED PROCESS CHRONOLOGY -->"
CHRONOLOGY_END = "<!-- END GENERATED PROCESS CHRONOLOGY -->"
ROUTES_BEGIN = "<!-- BEGIN GENERATED RANK ROUTES -->"
ROUTES_END = "<!-- END GENERATED RANK ROUTES -->"
MECHANISMS_BEGIN = "<!-- BEGIN GENERATED CANCELLATION MECHANISMS -->"
MECHANISMS_END = "<!-- END GENERATED CANCELLATION MECHANISMS -->"

ADE_TERM = re.compile(r"(?:(\d+))?([ADE])(\d+)")
STATUS_ID = re.compile(r"^[A-Z][A-Z0-9-]+$")


class LedgerError(RuntimeError):
    """A ledger consistency failure."""


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise LedgerError(message)


def unique_index(rows: Iterable[dict[str, Any]], kind: str) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for row in rows:
        row_id = row.get("id")
        require(isinstance(row_id, str) and row_id, f"{kind} has a missing id")
        require(row_id not in index, f"duplicate {kind} id: {row_id}")
        index[row_id] = row
    return index


def ade_rank(ade: str) -> int:
    if ade == "0":
        return 0
    total = 0
    terms = ade.split("+")
    for term in terms:
        match = ADE_TERM.fullmatch(term)
        require(match is not None, f"cannot parse ADE term {term!r} in {ade!r}")
        multiplicity = int(match.group(1) or 1)
        total += multiplicity * int(match.group(3))
    return total


def parse_time(value: str) -> datetime:
    parsed = datetime.fromisoformat(value)
    require(parsed.tzinfo is not None, f"timestamp lacks timezone: {value}")
    return parsed


def status_ids() -> set[str]:
    authority = load_json(ROOT / "MATH_STATUS.json")
    return {entry["id"] for entry in authority["entries"]}


def validate_evidence(reference: str, known_status_ids: set[str]) -> None:
    if reference.startswith("https://") or reference.startswith("git:"):
        return
    if STATUS_ID.fullmatch(reference):
        require(reference in known_status_ids, f"unknown MATH_STATUS id: {reference}")
        return
    require((ROOT / reference).exists(), f"missing evidence path: {reference}")


def validate_ledger(ledger: dict[str, Any]) -> dict[str, Any]:
    require(ledger.get("schema") == "elkies-k3.process-ledger.v1", "wrong ledger schema")
    require(ledger.get("status_authority") == "MATH_STATUS.json", "wrong status authority")

    stages = unique_index(ledger["stages"], "stage")
    transitions = unique_index(ledger["transitions"], "transition")
    events = unique_index(ledger["events"], "event")
    mechanisms = unique_index(ledger["mechanisms"], "mechanism")
    literature = unique_index(ledger["literature"], "literature entry")
    known_status_ids = status_ids()

    for stage in stages.values():
        expected_root_rank = ade_rank(stage["ade"])
        require(
            stage["root_rank"] == expected_root_rank,
            f"{stage['id']}: root rank {stage['root_rank']} != ADE rank {expected_root_rank}",
        )
        budget = stage["picard_rank"] - 2
        require(
            stage["root_rank"] + stage["mw_rank"] == budget,
            f"{stage['id']}: Shioda--Tate balance fails",
        )
        for reference in stage["evidence"]:
            validate_evidence(reference, known_status_ids)

    route_edges: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for transition in transitions.values():
        source = stages.get(transition["from"])
        target = stages.get(transition["to"])
        require(source is not None, f"{transition['id']}: unknown source stage")
        require(target is not None, f"{transition['id']}: unknown target stage")
        if "factorization" in transition:
            left, right = transition["factorization"]
            require(left * right == transition["q"], f"{transition['id']}: q factorization fails")
        require(
            source["picard_rank"] == target["picard_rank"],
            f"{transition['id']}: a neighbor unexpectedly changes Picard rank",
        )
        require(
            source.get("ns_discriminant") == target.get("ns_discriminant")
            or source.get("ns_discriminant") is None
            or target.get("ns_discriminant") is None,
            f"{transition['id']}: a same-surface neighbor changes NS discriminant",
        )
        root_delta = target["root_rank"] - source["root_rank"]
        mw_delta = target["mw_rank"] - source["mw_rank"]
        require(root_delta + mw_delta == 0, f"{transition['id']}: rank exchange fails")
        route_edges[transition["route"]].append(transition)

    event_times = [parse_time(event["observed_at"]) for event in events.values()]
    require(event_times == sorted(event_times), "events are not in chronological order")
    for event in events.values():
        has_globs = bool(event.get("artifact_globs"))
        require(
            (event["time_basis"] == "imported_artifact_mtime") == has_globs,
            f"{event['id']}: artifact mtime events must have globs, and only those events may have globs",
        )
        for reference in event["evidence"]:
            validate_evidence(reference, known_status_ids)

    for mechanism in mechanisms.values():
        for reference in mechanism["evidence"]:
            validate_evidence(reference, known_status_ids)

    for item in literature.values():
        require(item["url"].startswith("https://"), f"{item['id']}: literature URL is not HTTPS")

    indegree = Counter(transition["to"] for transition in transitions.values())
    repeated_targets = {stage: count for stage, count in indegree.items() if count > 1}
    require(not repeated_targets, f"multiple ledger transitions target one stage: {repeated_targets}")

    return {
        "stage_count": len(stages),
        "transition_count": len(transitions),
        "event_count": len(events),
        "mechanism_count": len(mechanisms),
        "literature_count": len(literature),
        "route_transition_counts": dict(sorted((key, len(value)) for key, value in route_edges.items())),
    }


def markdown_escape(value: Any) -> str:
    text = str(value).replace("\n", " ").replace("|", "\\|")
    return re.sub(r"\s+", " ", text).strip()


def short_date(value: str) -> str:
    return parse_time(value).strftime("%m-%d %H:%M")


def render_chronology(ledger: dict[str, Any]) -> str:
    lines = [
        "| Time | Turn / behavior | How it was found | Mathematical meaning | Reusable rule | Status |",
        "|---|---|---|---|---|---|",
    ]
    for event in ledger["events"]:
        lines.append(
            "| "
            + " | ".join(
                markdown_escape(value)
                for value in (
                    short_date(event["observed_at"]),
                    f"{event['phase']}: {event['behavior']}",
                    event["method"],
                    event["meaning"],
                    event["generalization"],
                    event["status"],
                )
            )
            + " |"
        )
    return "\n".join(lines)


def ordered_route_stages(ledger: dict[str, Any], route: str) -> list[dict[str, Any]]:
    stages = {stage["id"]: stage for stage in ledger["stages"]}
    edges = [edge for edge in ledger["transitions"] if edge["route"] == route]
    require(edges, f"no edges for route {route}")
    sources = {edge["from"] for edge in edges}
    targets = {edge["to"] for edge in edges}
    starts = sources - targets
    require(len(starts) == 1, f"route {route} does not have one start: {starts}")
    by_source = {edge["from"]: edge for edge in edges}
    require(len(by_source) == len(edges), f"route {route} branches; renderer expects a chain")
    result = [stages[starts.pop()]]
    while result[-1]["id"] in by_source:
        result.append(stages[by_source[result[-1]["id"]]["to"]])
    require(len(result) == len(edges) + 1, f"route {route} is disconnected")
    return result


def render_routes(ledger: dict[str, Any]) -> str:
    stage_index = {stage["id"]: stage for stage in ledger["stages"]}
    edge_index = {(edge["route"], edge["from"]): edge for edge in ledger["transitions"]}
    route_labels = {
        "low-q-backtrack": "Low-q reverse backtrack",
        "h3-selected": "H3 selected corridor",
        "q80-generic": "Q80 generic corridor",
        "q80-cm24": "Q80 CM24 specialization shadow",
    }
    lines = [
        "| Route | Stage change | q / old degree | Root rank | MW rank | Rank interpretation | Equation status |",
        "|---|---|---:|---:|---:|---|---|",
    ]
    for route in route_labels:
        stages = ordered_route_stages(ledger, route)
        first = stages[0]
        lines.append(
            f"| {route_labels[route]} | start: {markdown_escape(first['label'])} | - | "
            f"{first['root_rank']} | {first['mw_rank']} | rho={first['picard_rank']}; budget={first['picard_rank'] - 2} | "
            f"{markdown_escape(first['status'])} |"
        )
        for source, target in zip(stages, stages[1:]):
            edge = edge_index[(route, source["id"])]
            degree = edge.get("old_fiber_degree", "-")
            root_delta = target["root_rank"] - source["root_rank"]
            mw_delta = target["mw_rank"] - source["mw_rank"]
            interpretation = f"roots {root_delta:+d}; MW {mw_delta:+d}; sum fixed"
            lines.append(
                f"|  | {markdown_escape(source['label'])} -> {markdown_escape(target['label'])} | "
                f"{edge['q']} / {degree} | {target['root_rank']} | {target['mw_rank']} | "
                f"{interpretation} | {markdown_escape(target['status'])} |"
            )
    return "\n".join(lines)


def render_mechanisms(ledger: dict[str, Any]) -> str:
    lines = [
        "| Mechanism | Diagnostic symptom | Exact response | Meaning | General rule |",
        "|---|---|---|---|---|",
    ]
    for mechanism in ledger["mechanisms"]:
        lines.append(
            "| "
            + " | ".join(
                markdown_escape(value)
                for value in (
                    mechanism["name"],
                    mechanism["symptom"],
                    mechanism["resolution"],
                    mechanism["meaning"],
                    mechanism["generalization"],
                )
            )
            + " |"
        )
    return "\n".join(lines)


def replace_generated_block(document: str, begin: str, end: str, rendered: str) -> str:
    require(document.count(begin) == 1, f"document must contain one {begin}")
    require(document.count(end) == 1, f"document must contain one {end}")
    prefix, remainder = document.split(begin, 1)
    _, suffix = remainder.split(end, 1)
    return f"{prefix}{begin}\n{rendered}\n{end}{suffix}"


def check_document(ledger: dict[str, Any], path: Path) -> None:
    document = path.read_text(encoding="utf-8")
    expected = replace_generated_block(
        document, CHRONOLOGY_BEGIN, CHRONOLOGY_END, render_chronology(ledger)
    )
    expected = replace_generated_block(expected, ROUTES_BEGIN, ROUTES_END, render_routes(ledger))
    expected = replace_generated_block(
        expected,
        MECHANISMS_BEGIN,
        MECHANISMS_END,
        render_mechanisms(ledger),
    )
    require(document == expected, f"generated tables are stale in {path.relative_to(ROOT)}")


def update_document(ledger: dict[str, Any], path: Path) -> None:
    document = path.read_text(encoding="utf-8")
    document = replace_generated_block(
        document, CHRONOLOGY_BEGIN, CHRONOLOGY_END, render_chronology(ledger)
    )
    document = replace_generated_block(document, ROUTES_BEGIN, ROUTES_END, render_routes(ledger))
    document = replace_generated_block(
        document,
        MECHANISMS_BEGIN,
        MECHANISMS_END,
        render_mechanisms(ledger),
    )
    path.write_text(document, encoding="utf-8")


def local_artifact_audit(ledger: dict[str, Any]) -> dict[str, Any]:
    require(LOCAL_ROOT.is_dir(), f"local artifact directory not found: {LOCAL_ROOT}")
    local_stats = [
        (path, path.stat()) for path in LOCAL_ROOT.rglob("*") if path.is_file()
    ]
    require(local_stats, f"local artifact directory is empty: {LOCAL_ROOT}")
    event_rows = []
    all_files: set[Path] = set()
    for event in ledger["events"]:
        patterns = event.get("artifact_globs")
        if not patterns:
            continue
        matches: set[Path] = set()
        for pattern in patterns:
            recursive_pattern = f"{pattern}/*" if pattern.endswith("/**") else pattern
            matches.update(path for path in LOCAL_ROOT.glob(recursive_pattern) if path.is_file())
        require(matches, f"{event['id']}: local artifact globs match no files")
        all_files.update(matches)
        stats = [(path, path.stat()) for path in sorted(matches)]
        earliest = min(stat.st_mtime for _, stat in stats)
        latest = max(stat.st_mtime for _, stat in stats)
        path_time_rows = [
            f"{path.relative_to(LOCAL_ROOT)}\t{stat.st_size}\t{stat.st_mtime_ns}"
            for path, stat in stats
        ]
        event_rows.append(
            {
                "event": event["id"],
                "file_count": len(stats),
                "earliest_mtime": datetime.fromtimestamp(earliest).astimezone().isoformat(),
                "latest_mtime": datetime.fromtimestamp(latest).astimezone().isoformat(),
                "metadata_sha256": hashlib.sha256(
                    ("\n".join(path_time_rows) + "\n").encode()
                ).hexdigest(),
            }
        )
    return {
        "warning": "Imported mtimes establish discovery order only; they are not mathematical certificates.",
        "local_file_count": len(local_stats),
        "local_earliest_mtime": datetime.fromtimestamp(
            min(stat.st_mtime for _, stat in local_stats)
        ).astimezone().isoformat(),
        "local_latest_mtime": datetime.fromtimestamp(
            max(stat.st_mtime for _, stat in local_stats)
        ).astimezone().isoformat(),
        "matched_file_count": len(all_files),
        "events": event_rows,
    }


def build_audit(ledger: dict[str, Any], summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": "elkies-k3.process-ledger-audit.v1",
        "ledger_sha256": hashlib.sha256(canonical_json_bytes(ledger)).hexdigest(),
        "checks": {
            **summary,
            "ade_rank_parser": "PASS",
            "shioda_tate_balance": "PASS",
            "neighbor_rank_exchange": "PASS",
            "transition_endpoints": "PASS",
            "chronology_order": "PASS",
            "evidence_resolution": "PASS",
        },
        "claim_boundary": "Navigation audit only. Underlying theorem and computation status remains in MATH_STATUS.json and cited certificates.",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--check-document", type=Path, nargs="?", const=DEFAULT_ATLAS)
    parser.add_argument("--update-document", type=Path, nargs="?", const=DEFAULT_ATLAS)
    parser.add_argument("--render", choices=("chronology", "routes", "mechanisms", "all"))
    parser.add_argument("--audit-local", action="store_true")
    parser.add_argument("--write-audit", type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_args()
    ledger = load_json(arguments.ledger)
    summary = validate_ledger(ledger)
    if arguments.check_document is not None:
        check_document(ledger, arguments.check_document)
    if arguments.update_document is not None:
        update_document(ledger, arguments.update_document)
    if arguments.render in ("chronology", "all"):
        print(render_chronology(ledger))
    if arguments.render == "all":
        print()
    if arguments.render in ("routes", "all"):
        print(render_routes(ledger))
    if arguments.render == "all":
        print()
    if arguments.render in ("mechanisms", "all"):
        print(render_mechanisms(ledger))
    if arguments.audit_local:
        print(json.dumps(local_artifact_audit(ledger), indent=2, sort_keys=True))
    if arguments.write_audit is not None:
        audit = build_audit(ledger, summary)
        arguments.write_audit.write_bytes(canonical_json_bytes(audit))
    if not any(
        (
            arguments.render,
            arguments.audit_local,
            arguments.write_audit is not None,
            arguments.update_document is not None,
        )
    ):
        print(
            "PASS Elkies--K3 process ledger: "
            f"{summary['stage_count']} stages, {summary['transition_count']} transitions, "
            f"{summary['event_count']} events, {summary['mechanism_count']} mechanisms"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except LedgerError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
