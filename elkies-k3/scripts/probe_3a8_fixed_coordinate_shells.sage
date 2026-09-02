#!/usr/bin/env sage-python
"""Scan every eligible residual-class fixed-coordinate shell in N(3A8).

For every nonidentity matrix conjugacy class of the exact order-twelve
residual group whose fixed rank is at least seven, test every rank-seven
coordinate direct summand of a pinned LLL basis of the primitive fixed
lattice. Apply the determinant-5000, discriminant-length-three, MW12--17,
and nontrivial mod-two complement-action gates.

This is exhaustive before residual quotienting for the declared coordinate
languages. Residual canonicalization and ternary/T-NS gates are separate.
"""

from __future__ import annotations

import argparse
import json
import runpy
from pathlib import Path

from sage.all import ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT / "artifacts/generated-results/elkies-k3-3a8-residual-group-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3a8-fixed-coordinate-shell-probe-v1.json"
)
SCAN_COMMON_SOURCE = (
    Path(__file__).resolve().parent / "probe_8a3_fixed_coordinate_shells.sage"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
SCAN_COMMON = runpy.run_path(
    str(SCAN_COMMON_SOURCE), run_name="_fixed_coordinate_scan_common"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")

scan_class = SCAN_COMMON["scan_class"]
digest = COMMON["digest"]


def build(catalog, residual, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == "elkies-k3.3a8-residual-group.v1"
    assert residual["status"] == "PASS_EXACT_3A8_GLUE_AND_RESIDUAL_GROUP"
    assert residual["residual_group"]["order"] == 12
    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == "3A8"
        ),
    )
    class_scans = [
        scan_class(ambient, class_row, determinant_bound, minimum_mw_rank)
        for class_row in residual["residual_group"]["conjugacy_classes"]
        if class_row["action_order"] > 1 and class_row["fixed_rank"] >= 7
    ]
    assert len(class_scans) == 4
    assert sum(
        row["accounting"]["coordinate_subsets_tested"] for row in class_scans
    ) == 13032
    assert [
        (
            row["residual_conjugacy_class_id"],
            row["accounting"]["high_mw_mod2_accepted_seeds"],
        )
        for row in class_scans
    ] == [
        ("3A8-C02", 189),
        ("3A8-C03", 0),
        ("3A8-C04", 0),
        ("3A8-C05", 0),
    ]
    assert class_scans[0]["accounting"][
        "accepted_seed_mw_rank_distribution"
    ] == {"12": 135, "13": 54}
    return {
        "schema": "elkies-k3.3a8-fixed-coordinate-shell-probe.v1",
        "status": "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_3A8_COORDINATE_SHELL_SCAN",
        "proof_scope": {
            "proved": (
                "all rank-seven coordinate summands of pinned LLL bases for "
                "every nonidentity residual matrix conjugacy class of fixed "
                "rank at least seven pass through the declared lattice gates"
            ),
            "not_proved": (
                "residual or full-Weyl embedding orbits, ternary-genus "
                "admissibility, T/NS classes, or determinant-band completeness"
            ),
        },
        "parameters": {
            "ambient_label": "3A8",
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": (
                "7-of-fixed-rank coordinate direct summands of one pinned "
                "LLL basis per eligible residual matrix conjugacy class"
            ),
        },
        "residual_group_order": residual["residual_group"]["order"],
        "component_permutation_image_order": residual["residual_group"][
            "component_permutation_image_order"
        ],
        "class_scans": class_scans,
        "accounting": {
            "conjugacy_classes_scanned": len(class_scans),
            "coordinate_subsets_tested": sum(
                row["accounting"]["coordinate_subsets_tested"]
                for row in class_scans
            ),
            "high_mw_mod2_accepted_seeds_before_residual_dedup": sum(
                row["accounting"]["high_mw_mod2_accepted_seeds"]
                for row in class_scans
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--residual", type=Path, default=RESIDUAL)
    parser.add_argument("--determinant-bound", type=int, default=5000)
    parser.add_argument("--minimum-mw-rank", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    assert arguments.determinant_bound == 5000
    assert arguments.minimum_mw_rank == 12
    payload = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.residual.read_text()),
        arguments.determinant_bound,
        arguments.minimum_mw_rank,
    )
    payload["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(arguments.residual.resolve().relative_to(ROOT)): digest(arguments.residual),
        str(SCAN_COMMON_SOURCE.resolve().relative_to(ROOT)): digest(SCAN_COMMON_SOURCE),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_3a8_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("3A8 coordinate-shell probe is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "3A8SHELL|classes={}|subsets={}|seeds={}|status=PASS_EXACT".format(
            payload["accounting"]["conjugacy_classes_scanned"],
            payload["accounting"]["coordinate_subsets_tested"],
            payload["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
