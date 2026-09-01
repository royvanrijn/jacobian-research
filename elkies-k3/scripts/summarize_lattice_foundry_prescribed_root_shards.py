#!/usr/bin/env python3
"""Audit and summarize sharded prescribed-root foundry enumerations."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


SCHEMA = "elkies-k3-lattice-foundry-prescribed-root-shard-summary-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def sorted_counter(counter: Counter) -> dict:
    return {str(key): counter[key] for key in sorted(counter, key=str)}


def build_summary(paths: list[Path]) -> dict:
    if not paths:
        raise ValueError("at least one --input is required")

    inputs: dict[str, str] | None = None
    search_template: dict | None = None
    ambient_owner: dict[str, str] = {}
    occurrences: list[tuple[str, str, int, int, str]] = []
    shards = []

    for path in paths:
        data = load_json(path)
        status = data.get("status", "")
        if not isinstance(status, str) or not status.startswith("PASS_EXACT_"):
            raise ValueError(f"{path}: status is not an exact PASS status")
        search = data.get("search", {})
        if search.get("source_root_rank_range") != [16, 17]:
            raise ValueError(f"{path}: expected source root ranks 16--17")
        if search.get("source_mw_range_for_rho_19") != [0, 1]:
            raise ValueError(f"{path}: expected source MW ranks 0--1")
        if search.get("truncated_by_prefix_limit"):
            raise ValueError(f"{path}: prefix-truncated search is not census input")
        if search.get("source_support_range") != [1, 17] or search.get("all_a_only"):
            raise ValueError(f"{path}: expected the unrestricted support/type window")
        if len(search.get("ns_ids", [])) != 48:
            raise ValueError(f"{path}: expected all 48 foundry NS classes")
        if data.get("skipped_ns_classes"):
            raise ValueError(f"{path}: one or more NS classes were skipped")
        comparable_search = {key: value for key, value in search.items() if key != "ambient_labels"}
        if search_template is None:
            search_template = comparable_search
        elif comparable_search != search_template:
            raise ValueError(f"{path}: search settings differ across shards")
        if inputs is None:
            inputs = data.get("inputs")
        elif data.get("inputs") != inputs:
            raise ValueError(f"{path}: input hashes differ across shards")

        ambients = search.get("ambient_labels", [])
        for ambient in ambients:
            if ambient in ambient_owner:
                raise ValueError(
                    f"ambient {ambient} occurs in both {ambient_owner[ambient]} and {path}"
                )
            ambient_owner[ambient] = str(path)

        shard_mw = Counter()
        shard_types = Counter()
        shard_ns = set()
        for source in data.get("sources", []):
            ns_id = source["ns_id"]
            frame = source["source"]
            mw_rank = frame["mw_rank_for_rho_19"]
            root_rank = frame["root_rank"]
            if mw_rank not in (0, 1) or root_rank not in (16, 17):
                raise ValueError(f"{path}: source outside declared MW0--1 window")
            identity = (ns_id, frame["gram_sha256"], mw_rank, root_rank, frame["root_type"])
            occurrences.append(identity)
            shard_mw[mw_rank] += 1
            shard_types[frame["root_type"]] += 1
            shard_ns.add(ns_id)

        shards.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "ambient_labels": ambients,
                "selected_d5_anchors": data["accounting"]["selected_d5_anchors"],
                "source_occurrences": len(data.get("sources", [])),
                "ns_classes_with_source": len(shard_ns),
                "ns_classes_without_source": sorted(set(search["ns_ids"]) - shard_ns),
                "by_mw_rank": sorted_counter(shard_mw),
                "by_root_type": sorted_counter(shard_types),
            }
        )

    identity_counts = Counter((item[0], item[1]) for item in occurrences)
    identity_metadata: dict[tuple[str, str], tuple[int, int, str]] = {}
    for ns_id, gram_sha256, mw_rank, root_rank, root_type in occurrences:
        key = (ns_id, gram_sha256)
        metadata = (mw_rank, root_rank, root_type)
        if key in identity_metadata and identity_metadata[key] != metadata:
            raise ValueError(f"inconsistent source metadata for {ns_id}/{gram_sha256}")
        identity_metadata[key] = metadata
    distinct = [key + identity_metadata[key] for key in identity_counts]
    distinct_mw = Counter(identity[2] for identity in distinct)
    distinct_types = Counter(identity[4] for identity in distinct)
    occurrence_mw = Counter(identity[2] for identity in occurrences)
    occurrence_types = Counter(identity[4] for identity in occurrences)
    repeated_by_ns = defaultdict(int)
    for identity, count in identity_counts.items():
        if count > 1:
            repeated_by_ns[identity[0]] += count - 1
    all_ns_ids = set(search_template["ns_ids"] if search_template else [])
    covered_ns_ids = {identity[0] for identity in distinct}

    return {
        "schema": SCHEMA,
        "status": "PASS",
        "objective": (
            "Exact audit of the sharded all-NS prescribed-root enumeration for "
            "rank-16/17 source frames (MW0--1 at rho=19) already carrying "
            "same-NS MW15--17 targets."
        ),
        "inputs": inputs,
        "shards": shards,
        "accounting": {
            "ambient_labels": sorted(ambient_owner),
            "ambient_count": len(ambient_owner),
            "selected_d5_anchors": sum(item["selected_d5_anchors"] for item in shards),
            "source_occurrences": len(occurrences),
            "distinct_ns_gram_pairs": len(distinct),
            "repeated_ns_gram_occurrences_across_shards": len(occurrences) - len(distinct),
            "ns_classes_with_source": len(covered_ns_ids),
            "ns_classes_without_source": sorted(all_ns_ids - covered_ns_ids),
            "by_mw_rank_occurrences": sorted_counter(occurrence_mw),
            "by_mw_rank_distinct_ns_gram": sorted_counter(distinct_mw),
            "by_root_type_occurrences": sorted_counter(occurrence_types),
            "by_root_type_distinct_ns_gram": sorted_counter(distinct_types),
            "repeated_occurrences_by_ns": dict(sorted(repeated_by_ns.items())),
        },
        "proof_boundary": [
            "Every accepted source is an exact primitive auxiliary embedding, saturated complement, full norm-two root classification, and genus match in its shard.",
            "The census is complete only for the declared sequential D5-anchor cover, sixth norm bound, and rank-16/17 prescribed-root window.",
            "An equal reduced Gram digest is merged only within the same NS class; unequal reduced Grams are not claimed to be distinct lattice-isometry or J2 classes.",
            "No source in this summary has yet acquired a rational K3 marking, Weierstrass equation, or certified geometric neighbour corridor to a target frame.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    summary = build_summary(args.input)
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"summary differs from {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")

    accounting = summary["accounting"]
    print(
        "FOUNDRYPRESCRIBEDSUMMARY"
        f"|ambients={accounting['ambient_count']}"
        f"|occurrences={accounting['source_occurrences']}"
        f"|distinct={accounting['distinct_ns_gram_pairs']}"
        f"|ns={accounting['ns_classes_with_source']}"
        "|status=PASS"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
