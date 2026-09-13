#!/usr/bin/env python3
"""Report where primary EC artifact digests are documented.

The registry is the authoritative location for an ``artifact_hash``.  This
audit distinguishes that registry coverage from a canonical note repeating a
digest prefix, so a source-review pass can select the latter deliberately
instead of bulk-copying hashes into historical prose.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HASH = re.compile(r"sha256:([0-9a-f]{64})$")
STATES = {"proved", "partial"}
REVIEW_PATH = ROOT / "knowledge/canonical_evidence_reviews.json"
REVIEW_FIELDS = {
    "id", "canonical_source", "source_sha256", "evidence_route", "boundary",
}


def rows() -> list[dict[str, str | bool]]:
    status = json.loads((ROOT / "MATH_STATUS.json").read_text())
    result: list[dict[str, str | bool]] = []
    for entry in status["entries"]:
        source = entry.get("canonical_source", "")
        digest = entry.get("artifact_hash")
        if (
            entry.get("state") not in STATES
            or not source.startswith("elliptic-curves/")
            or not digest
        ):
            continue
        match = HASH.fullmatch(digest)
        if not match:
            raise SystemExit(f"invalid artifact hash on {entry['id']}: {digest!r}")
        path = ROOT / source
        if not path.is_file():
            raise SystemExit(f"missing canonical source for {entry['id']}: {source}")
        prefix = match.group(1)[:16]
        source_text = path.read_text(errors="replace").lower()
        generated_locks = [
            lock
            for lock in entry.get("software_lock", [])
            if lock.startswith("artifacts/generated-results/")
        ]
        result.append(
            {
                "id": entry["id"],
                "source": source,
                "digest_prefix": prefix,
                "direct_note_digest": prefix in source_text,
                "registered_generated_artifact": bool(generated_locks),
                "direct_note_artifact_reference": any(
                    lock.lower() in source_text or Path(lock).name.lower() in source_text
                    for lock in generated_locks
                ),
            }
        )
    return result


def validate_source_reviews(coverage: list[dict[str, str | bool]]) -> int:
    """Require a source-hash-bound explanation for every no-artifact row."""
    data = json.loads(REVIEW_PATH.read_text())
    if set(data) != {"schema_version", "reviews"} or data["schema_version"] != 1:
        raise SystemExit("invalid canonical evidence-review schema")
    reviews = data["reviews"]
    if not isinstance(reviews, list):
        raise SystemExit("canonical evidence reviews must be a list")
    by_id: dict[str, dict[str, str]] = {}
    for review in reviews:
        if set(review) != REVIEW_FIELDS or not all(
            isinstance(value, str) and value for value in review.values()
        ):
            raise SystemExit("invalid canonical evidence review")
        if review["id"] in by_id:
            raise SystemExit(f"duplicate canonical evidence review: {review['id']}")
        by_id[review["id"]] = review
    expected = {
        str(row["id"])
        for row in coverage
        if not row["registered_generated_artifact"]
    }
    if set(by_id) != expected:
        missing = sorted(expected - set(by_id))
        extra = sorted(set(by_id) - expected)
        raise SystemExit(
            "canonical evidence reviews do not match no-artifact rows: "
            f"missing {missing}; extra {extra}"
        )
    sources = {str(row["id"]): str(row["source"]) for row in coverage}
    for item_id, review in by_id.items():
        source = sources[item_id]
        if review["canonical_source"] != source:
            raise SystemExit(f"{item_id}: stale canonical evidence-review source")
        actual = hashlib.sha256((ROOT / source).read_bytes()).hexdigest()
        if review["source_sha256"] != actual:
            raise SystemExit(f"{item_id}: stale canonical evidence-review source hash")
    return len(reviews)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit individual source rows")
    args = parser.parse_args()
    coverage = rows()
    reviewed_without_artifact = validate_source_reviews(coverage)
    direct = sum(bool(row["direct_note_digest"]) for row in coverage)
    references = sum(bool(row["direct_note_artifact_reference"]) for row in coverage)
    registered = sum(bool(row["registered_generated_artifact"]) for row in coverage)
    unlinked_registered = sum(
        row["registered_generated_artifact"]
        and not row["direct_note_digest"]
        and not row["direct_note_artifact_reference"]
        for row in coverage
    )
    summary = {
        "registry_artifact_hashes": len(coverage),
        "direct_canonical_note_digests": direct,
        "registered_generated_artifacts": registered,
        "direct_canonical_note_artifact_references": references,
        "canonical_note_digest_not_repeated": len(coverage) - direct,
        "unlinked_registered_generated_artifacts": unlinked_registered,
        "no_registered_generated_artifact": len(coverage) - registered,
        "source_reviews_without_registered_artifact": reviewed_without_artifact,
    }
    if args.json:
        print(json.dumps({"summary": summary, "rows": coverage}, indent=2))
    else:
        print(
            "PASS canonical artifact-digest audit: "
            f"{direct}/{len(coverage)} direct canonical-note prefixes; "
            f"{references}/{registered} registered generated artifacts referenced; "
            f"{len(coverage) - registered} have no registered generated artifact; "
            f"{reviewed_without_artifact}/{len(coverage) - registered} source reviews current; "
            "all primary digests are registry-authoritative"
        )


if __name__ == "__main__":
    main()
