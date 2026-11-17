#!/usr/bin/env python3
"""Audit the committed plane wild-boundary ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_HASHES = {
    "scripts/verify_wild_boundary_atlas.py":
        "2be79bb21e9b6be71e6649f8c6bdf698fdcfcd06190ae67d9cba9d581f46c67a",
    "scripts/compile_plane_wild_boundary_survivors.py":
        "e887c3a19cd8087497d95298e729b6ce50ebd8b9a631a9c2d7c301c6e3a1f3f5",
    "scripts/verify_plane_wild_boundary_p3_degree7.py":
        "dc347ec72765c98ff99790a4703b92d75196d0ad19ab4ee7a0cbbe0ffecf1aa4",
    "plane-jc/cas/boundary_lattice_prefilter.py":
        "37ea003916a359fcb563315e0327ee663f013e2bc5fcd946cabb646f9f7135a9",
    "artifacts/generated-results/plane_wild_boundary_survivor_atlas.json":
        "442a211cab64d4dc0694f5c956e700e38e1315032bb513862b76ef63f989589f",
    "artifacts/generated-results/plane_wild_boundary_p3_degree7_scan.json":
        "2dfbe8e4e4dbaf23c6e3ad4102d6a4d875dc0087fc7691b5b70dccf18b8355fa",
}


def load(name: str) -> dict[str, object]:
    return json.loads(
        (ROOT / "artifacts" / "generated-results" / name).read_text(
            encoding="utf-8"
        )
    )


def decisions(rows: list[dict[str, object]]) -> dict[str, int]:
    return dict(sorted(Counter(row["decision"] for row in rows).items()))


def audit_survivor_atlas() -> None:
    atlas = load("plane_wild_boundary_survivor_atlas.json")
    assert atlas["format"] == "plane-wild-boundary-survivor-atlas-v2"
    assert "A bounded packet survivor is not a cover and not a Keller map" in (
        atlas["scope"]
    )

    hidden = atlas["hidden_order_rows"]
    prescribed = atlas["prescribed_degree_rows"]
    balanced = atlas["balanced_prescribed_degree_rows"]
    comparisons = atlas["comparison_rows"]
    summary = atlas["summary"]
    assert summary == {
        "balanced_prescribed_degree_decisions": {"obstructed": 23},
        "balanced_prescribed_degree_row_count": 23,
        "comparison_decisions": {
            "local_comparison_only": 19,
            "needs_reconstruction": 6,
        },
        "comparison_row_count": 25,
        "hidden_decisions": {"known_keller": 1, "obstructed": 45},
        "hidden_row_count": 46,
        "odd_hidden_keller_rows": 0,
        "prescribed_degree_decisions": {
            "needs_reconstruction": 20,
            "obstructed": 3,
        },
        "prescribed_degree_row_count": 23,
    }
    assert len(hidden) == 46 and decisions(hidden) == summary["hidden_decisions"]
    assert len(prescribed) == 23
    assert decisions(prescribed) == summary["prescribed_degree_decisions"]
    assert len(balanced) == 23
    assert decisions(balanced) == summary["balanced_prescribed_degree_decisions"]
    assert len(comparisons) == 25
    assert decisions(comparisons) == summary["comparison_decisions"]

    known = [row for row in hidden if row["decision"] == "known_keller"]
    assert len(known) == 1
    assert (known[0]["characteristic"], known[0]["wild_degree"]) == (2, 2)
    assert not any(
        row["decision"] == "known_keller" and row["characteristic"] % 2
        for row in hidden
    )

    support = atlas["balanced_support_certificate"]
    assert support == {
        "admissible_retained_degrees_through_bound": {
            "3": [1, 4, 7, 10],
            "5": [1, 6],
            "7": [1, 8],
        },
        "coefficient_characterization": "A=a0+T*B(T^p)",
        "identity": "A*H_T-A'*H=P^(N-1)*Q*(A-T*A')",
        "rows_checked": 30,
    }
    root_count = atlas["balanced_root_count_certificate"]
    assert root_count["cover_count"] == "q^2+n_q(A)*q"
    assert root_count["boundary_count"] == "q"
    assert root_count["open_count"] == "q^2+(n_q(A)-1)*q"
    assert root_count["rows_checked"] == 8

    support_only = atlas["odd_characteristic_support_only_rows"]
    assert [
        (row["characteristic"], row["cover_degree"])
        for row in support_only
    ] == [(3, 7), (3, 10), (3, 13), (5, 11), (7, 15)]
    assert all(row["former_status"] == "support_only_survivor" for row in support_only)
    assert all(row["decision"] == "obstructed" for row in support_only)
    assert atlas["odd_characteristic_reconstruction_queue"] == []
    # That empty field is only the post-support balanced queue.  Other stored
    # architectures deliberately retain their unresolved reconstruction rows.
    assert decisions(prescribed)["needs_reconstruction"] == 20
    assert decisions(comparisons)["needs_reconstruction"] == 6

    packet_scan = atlas["odd_characteristic_packet_scan"]
    assert packet_scan["bounds"] == {"max_length": 3, "max_multiplicity": 12}
    expected_packets = {
        "3": (240, 186, 54),
        "5": (142, 103, 39),
        "7": (77, 43, 34),
    }
    for prime, (tested, rejected, survivors) in expected_packets.items():
        profile = packet_scan["by_characteristic"][prime]
        assert (
            profile["tested_count"],
            profile["rejected_by_vertical_torsion"],
            profile["survivor_count"],
        ) == (tested, rejected, survivors)
        assert rejected + survivors == tested
        assert all(
            row["status"] == "packet_gate_only"
            for row in profile["minimal_survivors"]
        )
    assert "only packets not rejected" in packet_scan["interpretation"]

    module = atlas["characteristic_zero_boundary_module"]
    assert module["status"] == "coherent_template_not_instantiated"
    assert module["required_order_of_gates"] == [
        "different vanishes away from the omitted boundary",
        "compiled determinant class vanishes in Cl(U) or Pic(U)",
        "finite-support local-cohomology residue vanishes",
        "affineness and polynomial-coordinate reconstruction",
    ]


def audit_degree_seven_scan() -> None:
    scan = load("plane_wild_boundary_p3_degree7_scan.json")
    assert scan["format"] == "plane-wild-boundary-p3-degree7-scan-v2"
    assert (scan["characteristic"], scan["cover_degree"]) == (3, 7)
    assert "No affine-plane source or constant-Jacobian polynomial map" in (
        scan["scope"]
    )
    assert scan["summary"] == {
        "coefficient_rows": 6,
        "geometric_survivors": 0,
        "geometrically_obstructed": 6,
        "normalization_smooth": 6,
        "point_count_obstructed_over_F3": 4,
        "point_count_survivor_ids_over_F3": [[1, 1], [1, 2]],
        "point_count_survivors_over_F27": 0,
        "point_count_survivors_over_F3": 2,
        "relative_different_gate_passed": 6,
    }

    rows = scan["rows"]
    assert [(row["a0"], row["b"]) for row in rows] == [
        (1, 0),
        (1, 1),
        (1, 2),
        (2, 0),
        (2, 1),
        (2, 2),
    ]
    assert all(row["decision"] == "geometrically_obstructed" for row in rows)
    assert all(row["normalization"] == "smooth" for row in rows)
    assert all(
        row["relative_different_support"] == "fierce_boundary_only"
        for row in rows
    )
    prime_survivors = [
        row for row in rows if row["prime_field_decision"] == "point_count_survivor"
    ]
    assert [(row["a0"], row["b"]) for row in prime_survivors] == [(1, 1), (1, 2)]
    for row in prime_survivors:
        assert row["factor_degrees_over_F3"] == [1, 3]
        assert row["decisive_field"] == "F_27"
        assert row["extension_point_counts"] == [
            {
                "boundary": 9,
                "cover": 90,
                "field": "F_9",
                "open": 81,
                "retained_roots": 1,
            },
            {
                "boundary": 27,
                "cover": 837,
                "field": "F_27",
                "open": 810,
                "retained_roots": 4,
            },
        ]


def audit_fail_closed_sources() -> None:
    compiler = (ROOT / "scripts/compile_plane_wild_boundary_survivors.py").read_text(
        encoding="utf-8"
    )
    assert '"status": "packet_gate_only"' in compiler
    assert "A bounded packet survivor is not a cover" in compiler
    assert "if args.write:" in compiler
    assert "assert OUTPUT.read_text() == serialized" in compiler

    degree_seven = (
        ROOT / "scripts/verify_plane_wild_boundary_p3_degree7.py"
    ).read_text(encoding="utf-8")
    assert "Singular is required for the degree-seven scan" in degree_seven
    assert "completed.returncode != 0 or marker not in completed.stdout" in degree_seven
    assert "if args.write:" in degree_seven
    assert "assert OUTPUT.read_text() == serialized" in degree_seven

    atlas_source = (ROOT / "scripts/verify_wild_boundary_atlas.py").read_text(
        encoding="utf-8"
    )
    assert "if args.singular:" in atlas_source
    assert "if args.balanced_singular:" in atlas_source
    assert "if args.thickened_singular:" in atlas_source
    assert "completed.returncode != 0 or marker not in completed.stdout" in atlas_source


def main() -> None:
    for relative_path, expected_hash in EXPECTED_HASHES.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (relative_path, actual_hash, expected_hash)

    audit_survivor_atlas()
    audit_degree_seven_scan()
    audit_fail_closed_sources()
    print(
        "PASS: plane wild-boundary sources and committed ledgers match pinned "
        "hashes, status partitions, point counts, and fail-closed behavior"
    )
    print(
        "BOUNDARY: the empty odd balanced support queue does not erase the "
        "20 prescribed-cover and 6 comparison reconstruction rows; packet "
        "survivors are abstract necessary-gate templates, not covers"
    )
    print(
        "BOUNDARY: the degree-seven result is the complete six-row F_3 "
        "coefficient scan inside one retained-polynomial architecture, not "
        "a theorem about arbitrary plane covers"
    )
    print("NO RECOMPUTATION: no SymPy or Singular calculation was run")


if __name__ == "__main__":
    main()
