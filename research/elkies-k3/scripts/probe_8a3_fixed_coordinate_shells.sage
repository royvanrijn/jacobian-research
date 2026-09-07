#!/usr/bin/env sage-python
"""Scan every eligible residual-class fixed-coordinate shell in N(8A3).

For each nonidentity matrix conjugacy class of the exact glue-code residual
group with fixed rank at least seven, test every rank-seven coordinate direct
summand of a pinned LLL basis of the primitive fixed lattice.  Apply the
determinant-5000, discriminant-length-three, MW12--17, and nontrivial mod-two
complement-action gates.

This is an exhaustive pre-residual-quotient scan of the declared coordinate
languages.  Residual canonicalization and ternary/T-NS gates are separate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import GF, ZZ, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-glue-code-residual-group-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shell-probe-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")

rows = COMMON["rows"]
digest = COMMON["digest"]
row_module_basis = COMMON["row_module_basis"]
primitive_closure_index = COMMON["primitive_closure_index"]
discriminant_invariants = COMMON["discriminant_invariants"]
induced_action = COMMON["induced_action"]
root_type = COMMON["root_type"]


def scan_class(ambient, class_row, determinant_bound, minimum_mw_rank):
    action = matrix(ZZ, class_row["representative_matrix"])
    identity24 = identity_matrix(ZZ, 24)
    fixed = row_module_basis(
        (action - identity24)
        .transpose()
        .right_kernel_matrix()
        .change_ring(ZZ)
    )
    assert fixed.nrows() == class_row["fixed_rank"]
    assert primitive_closure_index(fixed) == 1
    fixed_gram = fixed * ambient * fixed.transpose()
    assert fixed_gram.det() == class_row["fixed_determinant"]
    lll_change = fixed_gram.LLL_gram().transpose()
    assert abs(lll_change.det()) == 1
    reduced_basis = lll_change * fixed
    reduced_gram = reduced_basis * ambient * reduced_basis.transpose()

    counters = Counter()
    all_mw = Counter()
    accepted_mw = Counter()
    moved = Counter()
    seeds = []
    for combination in itertools.combinations(range(fixed.nrows()), 7):
        counters["coordinate_subsets_tested"] += 1
        auxiliary_basis = matrix(
            ZZ, [reduced_basis.row(index) for index in combination]
        )
        assert primitive_closure_index(auxiliary_basis) == 1
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        determinant = int(auxiliary.det())
        if determinant > determinant_bound:
            counters["determinant_rejected"] += 1
            continue
        invariants = discriminant_invariants(auxiliary)
        if len(invariants) > 3:
            counters["discriminant_length_rejected"] += 1
            continue
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        frame = complement_basis * ambient * complement_basis.transpose()
        assert frame.det() == determinant
        roots = root_type(frame)
        mw_rank = roots["mw_rank_for_rho_19"]
        all_mw[mw_rank] += 1
        if mw_rank < minimum_mw_rank:
            counters["mw_rank_below_factory_floor_rejected"] += 1
            continue
        complement_action = induced_action(complement_basis, action)
        moved_dimension = int(
            matrix(
                GF(2),
                complement_action - identity_matrix(ZZ, 17),
            ).rank()
        )
        if moved_dimension == 0:
            counters["mod2_trivial_rejected"] += 1
            continue
        accepted_mw[mw_rank] += 1
        moved[moved_dimension] += 1
        seeds.append(
            {
                "coordinate_subset_zero_based": list(combination),
                "determinant": determinant,
                "discriminant_invariants_greater_than_one": invariants,
                "mw_rank_for_rho_19": mw_rank,
                "selected_action_moved_dimension_mod_2": moved_dimension,
                "auxiliary_basis_in_ambient": rows(auxiliary_basis),
            }
        )
    return {
        "residual_conjugacy_class_id": class_row["class_id"],
        "class_size": class_row["class_size"],
        "action_order": class_row["action_order"],
        "fixed_lattice": {
            "rank": fixed.nrows(),
            "determinant": int(fixed_gram.det()),
            "basis_in_ambient": rows(fixed),
            "gram": rows(fixed_gram),
            "lll_change": rows(lll_change),
            "lll_basis_in_ambient": rows(reduced_basis),
            "lll_gram": rows(reduced_gram),
            "primitive_in_ambient": True,
        },
        "accounting": {
            "coordinate_subsets_tested": counters[
                "coordinate_subsets_tested"
            ],
            "determinant_rejected": counters["determinant_rejected"],
            "discriminant_length_rejected": counters[
                "discriminant_length_rejected"
            ],
            "all_length_admissible_mw_rank_distribution": {
                str(key): value for key, value in sorted(all_mw.items())
            },
            "mw_rank_below_factory_floor_rejected": counters[
                "mw_rank_below_factory_floor_rejected"
            ],
            "mod2_trivial_rejected": counters["mod2_trivial_rejected"],
            "high_mw_mod2_accepted_seeds": len(seeds),
            "accepted_seed_mw_rank_distribution": {
                str(key): value for key, value in sorted(accepted_mw.items())
            },
            "selected_action_moved_dimension_mod_2_distribution": {
                str(key): value for key, value in sorted(moved.items())
            },
        },
        "accepted_seeds": seeds,
    }


def build(catalog, residual, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == (
        "elkies-k3.8a3-glue-code-residual-group.v1"
    )
    assert residual["status"] == (
        "PASS_EXACT_8A3_GLUE_CODE_AND_RESIDUAL_GROUP"
    )
    assert residual["residual_group"]["order"] == 2688
    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == "8A3"
        ),
    )
    class_scans = [
        scan_class(
            ambient,
            class_row,
            determinant_bound,
            minimum_mw_rank,
        )
        for class_row in residual["residual_group"]["conjugacy_classes"]
        if class_row["action_order"] > 1 and class_row["fixed_rank"] >= 7
    ]
    assert len(class_scans) == 7
    assert sum(
        row["accounting"]["coordinate_subsets_tested"]
        for row in class_scans
    ) == 24600
    assert [
        (
            row["residual_conjugacy_class_id"],
            row["accounting"]["high_mw_mod2_accepted_seeds"],
        )
        for row in class_scans
    ] == [
        ("8A3-C02", 880),
        ("8A3-C03", 128),
        ("8A3-C04", 0),
        ("8A3-C05", 96),
        ("8A3-C06", 58),
        ("8A3-C07", 0),
        ("8A3-C09", 4),
    ]
    return {
        "schema": "elkies-k3.8a3-fixed-coordinate-shell-probe.v1",
        "status": "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_8A3_COORDINATE_SHELL_SCAN",
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
            "ambient_label": "8A3",
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
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(
            arguments.catalog
        ),
        str(arguments.residual.resolve().relative_to(ROOT)): digest(
            arguments.residual
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_8a3_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("8A3 coordinate-shell probe is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "8A3SHELL|classes={}|subsets={}|seeds={}|status=PASS_EXACT".format(
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
