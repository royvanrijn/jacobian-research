#!/usr/bin/env python3
"""Report repository-local SHA-256 edges embedded in generated JSON artifacts.

Generated certificates often pin the exact programs and artifacts used to
produce them.  A changed source hash is not automatically an error: immutable
historical outputs should continue to name the code that actually ran.  This
audit therefore calls an edge ``current`` only when the owning artifact and
the referenced input occur together in the ``software_lock`` of at least one
non-archived ``MATH_STATUS.json`` entry.  All other recognizable edges are
reported as historical/unclassified evidence.

The extractor is deliberately conservative.  It recognizes path-to-hash maps
and unambiguous sibling path/hash fields; it does not guess among several
possible paths or hashes.  Owning checkers remain authoritative for formats
that need richer interpretation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
HASH_PATTERN = re.compile(r"^(?:sha256:)?([0-9a-fA-F]{64})$")
PATH_FIELD_TOKENS = ("artifact", "checker", "file", "input", "path", "script", "source")
HASH_FIELD_TOKENS = ("digest", "hash", "sha256")
KNOWN_TOP_LEVELS = {
    "artifacts",
    "data",
    "elliptic-curves",
    "elkies-k3",
    "formal",
    "jcsearch",
    "lean",
    "papers",
    "schemas",
    "scripts",
}


@dataclass(frozen=True)
class Edge:
    owner: str
    input: str
    expected_sha256: str
    pointer: str


@dataclass(frozen=True)
class Result:
    owner: str
    input: str
    pointer: str
    lock_class: str
    status: str
    expected_sha256: str
    actual_sha256: str | None
    status_entries: tuple[str, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_hash(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    match = HASH_PATTERN.fullmatch(value.strip())
    return match.group(1).lower() if match else None


def normalized_repo_path(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    candidate = value.strip().replace("\\", "/")
    if candidate.startswith("./"):
        candidate = candidate[2:]
    if not candidate or "://" in candidate or candidate.startswith("/"):
        return None
    pure = PurePosixPath(candidate)
    if ".." in pure.parts or len(pure.parts) < 2:
        return None
    if pure.parts[0] not in KNOWN_TOP_LEVELS:
        return None
    return str(pure)


def field_stem(name: str) -> str:
    words = [word for word in re.split(r"[^a-z0-9]+", name.lower()) if word]
    while words and words[-1] in {*PATH_FIELD_TOKENS, *HASH_FIELD_TOKENS}:
        words.pop()
    return "_".join(words)


def child_pointer(pointer: str, key: Any) -> str:
    token = str(key).replace("~", "~0").replace("/", "~1")
    return f"{pointer}/{token}"


def extract_edges(owner: str, payload: Any) -> list[Edge]:
    edges: set[Edge] = set()

    def add(path: Any, digest: Any, pointer: str) -> None:
        normalized_path = normalized_repo_path(path)
        normalized_digest = normalized_hash(digest)
        if normalized_path and normalized_digest:
            edges.add(Edge(owner, normalized_path, normalized_digest, pointer))

    def walk(value: Any, pointer: str) -> None:
        if isinstance(value, dict):
            # The common compact form is {"repository/path": "<sha256>"}.
            for key, item in value.items():
                add(key, item, child_pointer(pointer, key))

            path_fields = {
                str(key): item
                for key, item in value.items()
                if isinstance(item, str)
                and normalized_repo_path(item)
                and any(token in str(key).lower() for token in PATH_FIELD_TOKENS)
            }
            hash_fields = {
                str(key): item
                for key, item in value.items()
                if normalized_hash(item)
                and any(token in str(key).lower() for token in HASH_FIELD_TOKENS)
            }

            # Prefer name-matched pairs such as source/source_sha256.  Generic
            # path+sha256 records are accepted only when the pairing is unique.
            paired: set[tuple[str, str]] = set()
            for path_key, path_value in path_fields.items():
                path_stem = field_stem(path_key)
                matches = [
                    (hash_key, hash_value)
                    for hash_key, hash_value in hash_fields.items()
                    if path_stem and field_stem(hash_key) == path_stem
                ]
                if len(matches) == 1:
                    hash_key, hash_value = matches[0]
                    add(path_value, hash_value, child_pointer(pointer, hash_key))
                    paired.add((path_key, hash_key))
            if len(path_fields) == 1 and len(hash_fields) == 1:
                path_key, path_value = next(iter(path_fields.items()))
                hash_key, hash_value = next(iter(hash_fields.items()))
                if (path_key, hash_key) not in paired:
                    add(path_value, hash_value, child_pointer(pointer, hash_key))

            for key, item in value.items():
                walk(item, child_pointer(pointer, key))
        elif isinstance(value, list):
            for index, item in enumerate(value):
                walk(item, child_pointer(pointer, index))

    walk(payload, "")
    return sorted(edges, key=lambda edge: (edge.input, edge.expected_sha256, edge.pointer))


def tracked_json_artifacts() -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "-z", "artifacts/generated-results"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return sorted(
        path.decode()
        for path in completed.stdout.split(b"\0")
        if path and path.decode().endswith(".json")
    )


def active_lock_entries() -> dict[str, set[str]]:
    registry = json.loads((ROOT / "MATH_STATUS.json").read_text())
    entries_by_path: dict[str, set[str]] = defaultdict(set)
    for entry in registry["entries"]:
        if entry["state"] == "archived":
            continue
        paths: Iterable[Any] = entry.get("software_lock", [])
        for value in paths:
            path = normalized_repo_path(value)
            if path:
                entries_by_path[path].add(entry["id"])
    return entries_by_path


def classify(edges: Iterable[Edge], locks: dict[str, set[str]]) -> list[Result]:
    results: list[Result] = []
    digest_cache: dict[str, str | None] = {}
    for edge in edges:
        common_entries = tuple(sorted(locks.get(edge.owner, set()) & locks.get(edge.input, set())))
        lock_class = "current" if common_entries else "historical/unclassified"
        if edge.input not in digest_cache:
            path = ROOT / edge.input
            digest_cache[edge.input] = sha256(path) if path.is_file() else None
        actual = digest_cache[edge.input]
        status = "missing" if actual is None else ("match" if actual == edge.expected_sha256 else "changed")
        results.append(
            Result(
                owner=edge.owner,
                input=edge.input,
                pointer=edge.pointer,
                lock_class=lock_class,
                status=status,
                expected_sha256=edge.expected_sha256,
                actual_sha256=actual,
                status_entries=common_entries,
            )
        )
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit the full report as JSON")
    parser.add_argument(
        "--include-matches",
        action="store_true",
        help="print matching edges as well as changed or missing edges",
    )
    parser.add_argument(
        "--include-historical",
        action="store_true",
        help="print changed or missing historical/unclassified edges",
    )
    parser.add_argument(
        "--fail-current",
        action="store_true",
        help="return nonzero when a current co-locked edge is changed or missing",
    )
    args = parser.parse_args()

    artifacts = tracked_json_artifacts()
    parse_errors: list[str] = []
    edges: list[Edge] = []
    for owner in artifacts:
        try:
            payload = json.loads((ROOT / owner).read_text())
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            parse_errors.append(f"{owner}: {error}")
            continue
        edges.extend(extract_edges(owner, payload))

    # Collapse repeated encodings of the same dependency inside one artifact.
    unique_edges = {
        (edge.owner, edge.input, edge.expected_sha256): edge
        for edge in edges
    }
    results = classify(unique_edges.values(), active_lock_entries())
    results.sort(key=lambda row: (row.lock_class, row.status, row.owner, row.input))
    counts = Counter((row.lock_class, row.status) for row in results)
    current_failures = [
        row for row in results
        if row.lock_class == "current" and row.status != "match"
    ]

    report = {
        "artifact_count": len(artifacts),
        "edge_count": len(results),
        "counts": {
            lock_class: {
                status: counts[(lock_class, status)]
                for status in ("match", "changed", "missing")
            }
            for lock_class in ("current", "historical/unclassified")
        },
        "parse_errors": parse_errors,
        "results": [asdict(row) for row in results],
    }
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(
            "PROVENANCE|artifacts={}|edges={}|current_match={}|current_changed={}|"
            "current_missing={}|historical_match={}|historical_changed={}|"
            "historical_missing={}|parse_errors={}".format(
                len(artifacts),
                len(results),
                counts[("current", "match")],
                counts[("current", "changed")],
                counts[("current", "missing")],
                counts[("historical/unclassified", "match")],
                counts[("historical/unclassified", "changed")],
                counts[("historical/unclassified", "missing")],
                len(parse_errors),
            )
        )
        for error in parse_errors:
            print(f"PARSE_ERROR|{error}")
        for row in results:
            if row.status == "match" and not args.include_matches:
                continue
            if row.lock_class != "current" and not args.include_historical and not args.include_matches:
                continue
            entries = ",".join(row.status_entries) if row.status_entries else "-"
            actual = row.actual_sha256 or "MISSING"
            print(
                f"{row.status.upper()}|class={row.lock_class}|owner={row.owner}|"
                f"input={row.input}|expected={row.expected_sha256}|actual={actual}|"
                f"entries={entries}|pointer={row.pointer}"
            )

    if parse_errors:
        return 2
    if args.fail_current and current_failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
