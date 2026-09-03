#!/usr/bin/env python3
"""Certify the compact bounded relative-U obstruction for the first NS0024 edge."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
OUTPUT = GENERATED / "elkies-k3-ns0024-relative-u-first-edge-obstruction-v1.json"
COMPARISON = GENERATED / "elkies-k3-ns0024-completed-frame-comparison-v1.json"
FRAME = GENERATED / "elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt"
SEARCHES = {
    2: (
        GENERATED / "elkies-k3-ns0024-relative-u-degree2-fibre-summary-v1.json",
        list(range(4, 41, 2)),
        18,
        431174,
        429877,
        5,
    ),
    3: (
        GENERATED / "elkies-k3-ns0024-relative-u-degree3-fibre-summary-v1.json",
        [9, 12, 15, 18, 21],
        4,
        13711,
        13704,
        7,
    ),
    4: (
        GENERATED / "elkies-k3-ns0024-relative-u-degree4-fibre-summary-v1.json",
        [16, 20, 24, 28, 32],
        4,
        79701,
        79375,
        8,
    ),
}


def load(path: Path):
    return json.loads(path.read_text())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def certify_search(degree, data, q_values, maximum_t, orbit_total,
                   primitive_total, expected_maximum_mw):
    assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert data["frame"] == relative(FRAME)
    assert data["determinant"] == 950
    assert data["input_root_data"] == [13, 280, 4]
    assert data["input_mw_rank"] == 4
    assert data["summary_only"] is True
    assert data["neighbors"] == []
    summaries = data["summaries"]
    assert [row["q"] for row in summaries] == q_values
    assert all(row["factor_order"] == [row["q"] // degree, degree]
               for row in summaries)
    assert all(row["q"] == degree * (degree + offset)
               for offset, row in enumerate(summaries))
    assert all(row["dominant_orbits_complete"] for row in summaries)
    assert all(row["mw_enumeration_complete"] for row in summaries)
    assert all(row["mw_vector_cap"] is None for row in summaries)
    assert all(not row["search_stopped_early"] for row in summaries)
    assert all(not row["stream_limit_reached"] for row in summaries)
    assert all(row["screened_orbits"] == row["dominant_orbits"]
               for row in summaries)
    assert all(row["primitive_neighbors"] + row["nonprimitive_orbits"]
               == row["dominant_orbits"] for row in summaries)
    assert sum(row["dominant_orbits"] for row in summaries) == orbit_total
    assert sum(row["primitive_neighbors"] for row in summaries) == primitive_total
    maximum_mw = max(
        item["mw_rank"]
        for row in summaries
        for item in row["root_histogram"]
    )
    assert maximum_mw == expected_maximum_mw
    target_profiles = sum(
        item["orbit_count"]
        for row in summaries
        for item in row["root_histogram"]
        if item["root_rank"] == 5
    )
    assert target_profiles == 0
    return {
        "degree": degree,
        "t_values": [0, maximum_t],
        "q_values": q_values,
        "dominant_orbits": orbit_total,
        "primitive_fibres": primitive_total,
        "maximum_child_mw_rank": maximum_mw,
        "root_rank_five_target_profiles": target_profiles,
        "complete_modulo": "source root Weyl group",
    }


def build_payload():
    comparison = load(COMPARISON)
    assert comparison["status"] == "PASS_EXACT_KNOWN_ROUTE_FRAME_COMPARISON"
    assert [(row["target_stage"], len(row["matches"]))
            for row in comparison["matches"]] == [(0, 0), (1, 0), (2, 0), (3, 0)]
    adaptation = comparison["completed_source_root_adaptation"]
    assert adaptation["primitive_root_lattice"] is True
    assert adaptation["root_rank"] == 13
    frame_rows = [
        list(map(int, line.split())) for line in FRAME.read_text().splitlines()
    ]
    assert frame_rows == adaptation["adapted_gram"]

    certified = []
    for degree, (path, q_values, maximum_t, orbits, primitive, maximum_mw) in SEARCHES.items():
        certified.append(
            certify_search(
                degree,
                load(path),
                q_values,
                maximum_t,
                orbits,
                primitive,
                maximum_mw,
            )
        )

    inputs = [FRAME, COMPARISON] + [row[0] for row in SEARCHES.values()]
    return {
        "schema": "elkies-k3.ns0024-relative-u-first-edge-obstruction.v1",
        "status": "PASS_EXACT_BOUNDED_FIRST_EDGE_OBSTRUCTION",
        "source": "completed NS0024 D5+E8/MW4 frame",
        "target": "completed NS0024 3A1+A2/MW12 frame",
        "fibre_relation": "q=d(d+t), with d=F.F' and t=O.F'",
        "known_route_comparison": {
            "completed_stages": 4,
            "isometric_matches_to_known_route": 0,
        },
        "searches": certified,
        "conclusion": (
            "No primitive target fibre with the required root rank five occurs "
            "in the declared (d,t) boxes. Therefore no choice of the second U' "
            "basis vector can give the requested first edge in those boxes."
        ),
        "input_hashes": {relative(path): digest(path) for path in inputs},
        "proof_boundary": {
            "proved": (
                "Complete exact primitive-fibre enumeration modulo the source "
                "root Weyl group for d=2,t=0..18 and d=3,4,t=0..4, including "
                "the zero Mordell-Weil projection, with exact child root ranks."
            ),
            "not_proved": (
                "Nonexistence outside the declared boxes; a global degree or t "
                "bound; a nef/effective zero, horizontal-wall audit, equation, "
                "rational map, or field of definition."
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    rendered = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if args.check:
        assert output.read_text() == rendered
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    print("PASS ns0024 relative-U first-edge obstruction")


if __name__ == "__main__":
    main()
