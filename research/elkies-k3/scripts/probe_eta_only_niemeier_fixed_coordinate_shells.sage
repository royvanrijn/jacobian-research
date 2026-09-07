#!/usr/bin/env sage-python
"""Scan every eligible eta-only residual fixed-coordinate shell.

Consume the exact residual groups for D24, D16+E8, A24, A17+E7, A15+D9,
and A11+D7+E6. For every nonidentity class of fixed rank at least seven,
test all rank-seven coordinate direct summands of a pinned LLL basis of the
primitive fixed lattice. Apply the determinant-5000,
discriminant-length-three, MW12--17, and nontrivial mod-two gates.

This exhausts the declared coordinate languages before residual quotienting;
it is not an enumeration of every primitive invariant rank-seven sublattice.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-eta-only-niemeier-residual-groups-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-eta-only-niemeier-fixed-coordinate-shell-probe-v1.json"
)
SCAN_COMMON_SOURCE = (
    Path(__file__).resolve().parent / "probe_8a3_fixed_coordinate_shells.sage"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
SCAN_COMMON = runpy.run_path(
    str(SCAN_COMMON_SOURCE), run_name="_eta_only_fixed_coordinate_scan_common"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")

scan_class = SCAN_COMMON["scan_class"]
digest = COMMON["digest"]


def build(catalog, residual, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == (
        "elkies-k3.eta-only-niemeier-residual-groups.v1"
    )
    assert residual["status"] == (
        "PASS_EXACT_ETA_ONLY_SIX_NIEMEIER_RESIDUAL_GROUPS"
    )
    ambient_by_label = {
        row["label"]: matrix(ZZ, row["gram"])
        for row in catalog["rooted_niemeier_lattices"]
    }
    backends = []
    accepted_mw_distribution = Counter()
    for residual_backend in residual["backends"]:
        label = residual_backend["ambient_label"]
        class_scans = [
            scan_class(
                ambient_by_label[label],
                class_row,
                determinant_bound,
                minimum_mw_rank,
            )
            for class_row in residual_backend["residual_group"][
                "conjugacy_classes"
            ]
            if class_row["action_order"] > 1 and class_row["fixed_rank"] >= 7
        ]
        for class_scan in class_scans:
            accepted_mw_distribution.update(
                {
                    int(rank): count
                    for rank, count in class_scan["accounting"][
                        "accepted_seed_mw_rank_distribution"
                    ].items()
                }
            )
        backends.append(
            {
                "ambient_label": label,
                "residual_group_order": residual_backend["residual_group"][
                    "order"
                ],
                "component_permutation_image_order": residual_backend[
                    "residual_group"
                ]["component_permutation_image_order"],
                "component_diagram_kernel_order": residual_backend[
                    "residual_group"
                ]["component_diagram_kernel_order"],
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
        )
    assert [
        (
            row["ambient_label"],
            row["residual_group_order"],
            row["accounting"]["conjugacy_classes_scanned"],
            row["accounting"]["coordinate_subsets_tested"],
        )
        for row in backends
    ] == [
        ("D24", 1, 0, 0),
        ("D16_E8", 1, 0, 0),
        ("A24", 2, 1, 792),
        ("A17_E7", 2, 1, 11440),
        ("A15_D9", 2, 1, 11440),
        ("A11_D7_E6", 2, 1, 11440),
    ]
    assert sum(
        row["accounting"]["coordinate_subsets_tested"] for row in backends
    ) == 35112
    assert all(
        row["accounting"][
            "high_mw_mod2_accepted_seeds_before_residual_dedup"
        ]
        == 0
        for row in backends
    )
    scan_by_label = {
        row["ambient_label"]: row["class_scans"][0]
        for row in backends
        if row["class_scans"]
    }
    assert {
        label: (
            row["accounting"]["determinant_rejected"],
            row["accounting"]["discriminant_length_rejected"],
            row["accounting"]["mw_rank_below_factory_floor_rejected"],
            row["accounting"]["mod2_trivial_rejected"],
        )
        for label, row in scan_by_label.items()
    } == {
        "A24": (0, 792, 0, 0),
        "A17_E7": (0, 6056, 5384, 0),
        "A15_D9": (0, 1332, 10108, 0),
        "A11_D7_E6": (0, 1625, 9815, 0),
    }
    assert accepted_mw_distribution == Counter()
    return {
        "schema": "elkies-k3.eta-only-niemeier-fixed-coordinate-shell-probe.v1",
        "status": "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_ETA_ONLY_COORDINATE_SHELL_SCAN",
        "proof_scope": {
            "proved": (
                "all rank-seven coordinate summands of pinned LLL bases for "
                "every nonidentity eta-only residual class of fixed rank at "
                "least seven pass through the declared lattice gates"
            ),
            "not_proved": (
                "residual or full-Weyl embedding orbits, ternary-genus "
                "admissibility, T/NS classes, or determinant-band completeness"
            ),
        },
        "parameters": {
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": (
                "7-of-fixed-rank coordinate direct summands of one pinned "
                "LLL basis per eligible residual matrix conjugacy class"
            ),
        },
        "backends": backends,
        "accounting": {
            "backends": len(backends),
            "conjugacy_classes_scanned": sum(
                row["accounting"]["conjugacy_classes_scanned"]
                for row in backends
            ),
            "coordinate_subsets_tested": sum(
                row["accounting"]["coordinate_subsets_tested"]
                for row in backends
            ),
            "high_mw_mod2_accepted_seeds_before_residual_dedup": sum(
                row["accounting"][
                    "high_mw_mod2_accepted_seeds_before_residual_dedup"
                ]
                for row in backends
            ),
            "accepted_seed_mw_rank_distribution": {
                str(rank): count
                for rank, count in sorted(accepted_mw_distribution.items())
            },
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
        "elkies-k3/scripts/probe_eta_only_niemeier_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("eta-only coordinate-shell probe is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "ETAONLYSHELL|backends={}|classes={}|subsets={}|seeds={}|status=PASS_EXACT".format(
            payload["accounting"]["backends"],
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
