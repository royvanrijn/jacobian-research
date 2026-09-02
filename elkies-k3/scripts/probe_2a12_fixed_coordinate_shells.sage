#!/usr/bin/env sage-python
"""Scan every eligible residual-class fixed-coordinate shell in N(2A12).

The exact cyclic order-four residual group has one nonidentity class of fixed
rank at least seven: its order-two simultaneous diagram reversal. Test every
rank-seven coordinate direct summand of a pinned LLL basis of that primitive
rank-twelve fixed lattice. Apply the determinant-5000,
discriminant-length-three, MW12--17, and nontrivial mod-two complement-action
gates.

This is exhaustive before residual quotienting for the declared coordinate
language. Residual canonicalization and ternary/T-NS gates are separate.
"""

from __future__ import annotations

import argparse
import json
import runpy
from pathlib import Path

from sage.all import ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = ROOT / "artifacts/generated-results/elkies-k3-2a12-residual-group-v1.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a12-fixed-coordinate-shell-probe-v1.json"
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
    assert residual["schema"] == "elkies-k3.2a12-residual-group.v1"
    assert residual["status"] == "PASS_EXACT_2A12_GLUE_AND_RESIDUAL_GROUP"
    assert residual["residual_group"]["order"] == 4
    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == "2A12"
        ),
    )
    class_scans = [
        scan_class(ambient, class_row, determinant_bound, minimum_mw_rank)
        for class_row in residual["residual_group"]["conjugacy_classes"]
        if class_row["action_order"] > 1 and class_row["fixed_rank"] >= 7
    ]
    assert len(class_scans) == 1
    assert class_scans[0]["accounting"]["coordinate_subsets_tested"] == 792
    assert class_scans[0]["residual_conjugacy_class_id"] == "2A12-C02"
    assert class_scans[0]["accounting"][
        "high_mw_mod2_accepted_seeds"
    ] == 0
    return {
        "schema": "elkies-k3.2a12-fixed-coordinate-shell-probe.v1",
        "status": "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_2A12_COORDINATE_SHELL_SCAN",
        "proof_scope": {
            "proved": (
                "all rank-seven coordinate summands of a pinned LLL basis for "
                "every nonidentity residual matrix conjugacy class of fixed "
                "rank at least seven pass through the declared lattice gates"
            ),
            "not_proved": (
                "residual or full-Weyl embedding orbits, ternary-genus "
                "admissibility, T/NS classes, or determinant-band completeness"
            ),
        },
        "parameters": {
            "ambient_label": "2A12",
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": (
                "7-of-12 coordinate direct summands of one pinned LLL basis "
                "for the order-two fixed lattice"
            ),
        },
        "residual_group_order": residual["residual_group"]["order"],
        "component_permutation_image_order": residual["residual_group"][
            "component_permutation_image_order"
        ],
        "component_diagram_kernel_order": residual["residual_group"][
            "component_diagram_kernel_order"
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
        "elkies-k3/scripts/probe_2a12_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("2A12 coordinate-shell probe is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "2A12SHELL|classes={}|subsets={}|seeds={}|status=PASS_EXACT".format(
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
