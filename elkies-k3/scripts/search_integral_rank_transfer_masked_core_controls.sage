#!/usr/bin/env sage-python
"""Run the fixed Golay mask-aware core search on prospective controls.

The default corridors are H3, NS0024, and Q80.  No search parameter is tuned
per corridor: each run resets the same Sage/Python seed and uses the same
good-prime beam rule as the certified Golay-720 construction.  A run either
stores an exact neighbour path to a zero-mask core or an explicit bounded
failure through the declared generation limit.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random as python_random
import runpy

from sage.all import Genus, QQ, ZZ, block_diagonal_matrix, matrix, pari, set_random_seed, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
BASE_SCRIPT = ROOT / "elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-controls-v1.json"

DEFAULT_CORRIDORS = ("H3", "NS0024", "Q80")


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configure_order(base, order):
    """The certified helper uses one module-level cyclic bridge order."""

    for name in ("primary_generator", "mask_profile"):
        base[name].__globals__["BRIDGE_ORDER"] = order


def prepare_corridor(corridor, bridge_artifact, theta_artifact, base, core, reverse):
    edge = next(
        row
        for row in bridge_artifact["edges"]
        if row["corridor"] == corridor and int(row["target_root_rank"]) == 0
    )
    theta_record = next(
        row for row in theta_artifact["corridors"] if row["corridor"] == corridor
    )
    order = int(theta_record["cyclic_glue_order"])
    configure_order(base, order)
    bridges = base["bridge_data"](theta_record, reverse)
    viable_bridges = [
        row for row in bridges if int(pari(row["gram"]).qfminim(2)[0]) == 0
    ]
    assert viable_bridges

    historical_core = matrix(ZZ, edge["core"]["gram"])
    historical_bridge = matrix(ZZ, edge["new_frame"]["bridge_gram"])
    historical_glue = vector(
        QQ,
        [
            QQ(value)
            for value in edge["new_frame"]["glue_generators"][0][
                "K_plus_C_dual_coordinates"
            ]
        ],
    )
    target_frame = core["glued_frame"](
        historical_core, historical_bridge, historical_glue
    )
    core_key = base["finite_form_key"](Genus(historical_core).discriminant_form())
    forced_keys = {
        json.dumps(
            base["finite_form_key"](
                Genus(block_diagonal_matrix(target_frame, -row["gram"]))
                .discriminant_form()
            ),
            sort_keys=True,
        )
        for row in bridges
    }
    assert forced_keys == {json.dumps(core_key, sort_keys=True)}
    all_genera, forced_genus = base["forced_genus"](
        target_frame, bridges[0]["gram"], core_key
    )
    seed = matrix(ZZ, forced_genus.representative())
    assert seed.nrows() == 15 and seed.det() == historical_core.det()
    return {
        "order": order,
        "bridges": bridges,
        "viable_bridges": viable_bridges,
        "historical_core": historical_core,
        "target_frame": target_frame,
        "all_genera": all_genera,
        "forced_genus": forced_genus,
        "seed": seed,
    }


def capped_mask_profile(gram, bridges, reverse, base, order, cap):
    """Return exact zero decisions and occupancy counts truncated at ``cap``."""

    core_generator = base["primary_generator"](gram, order)
    oracle = reverse["CoreCellOracle"](gram)
    rows = []
    for bridge in bridges:
        bridge_gram = bridge["gram"]
        bridge_generator = bridge["generator"]
        multipliers = [
            value
            for value in range(1, order)
            if base["modulo_two"](
                (value * core_generator)
                * gram
                * (value * core_generator)
                + bridge_generator * bridge_gram * bridge_generator
            )
            == 0
        ]
        assert len(multipliers) == 2 and sum(multipliers) == order
        multiplier = multipliers[0]
        violations = 0
        occupied_support = []
        for label in range(order):
            core_class = reverse["discriminant_class"](
                label * multiplier * core_generator
            )
            bridge_class = reverse["discriminant_class"](
                label * bridge_generator
            )
            for bridge_norm in bridge["theta_profile"].get(bridge_class, {}):
                core_norm = QQ(2) - bridge_norm
                if not 0 <= core_norm <= 2:
                    continue
                if oracle.occupied(core_class, core_norm):
                    violations += 1
                    occupied_support.append(
                        [min(label, order - label), str(core_norm)]
                    )
                    if violations >= cap:
                        break
            if violations >= cap:
                break
        rows.append(
            {
                "bridge_class_index": bridge["bridge_class_index"],
                "isotropic_multipliers": multipliers,
                "occupied_forbidden_cells_capped": violations,
                "occupied_support_capped": occupied_support,
                "zero_mask_accepts": violations == 0,
            }
        )
    return rows, len(oracle.cache), core_generator


def search_corridor(
    corridor,
    prepared,
    base,
    core,
    reverse,
    root_descent,
    mask_cap,
    generation_limit,
    persistent_elite,
    support_diversity,
):
    configure_order(base, prepared["order"])
    seed = prepared["seed"]
    bridges = prepared["viable_bridges"]
    set_random_seed(base["SEARCH_SEED"])
    python_random.seed(base["SEARCH_SEED"])
    frontier = [(base["quadratic_form"](seed), [])]
    seen = set()
    generations = []
    hit = None
    archive = []

    for generation in range(1, generation_limit + 1):
        candidates = []
        generated = 0
        tested_rootless = 0
        root_histogram = Counter()
        for parent, path in frontier:
            for _ in range(base["SEARCH_SAMPLES_PER_PARENT"]):
                prime = python_random.choice(base["SEARCH_PRIMES"])
                if seed.det() % prime == 0:
                    continue
                witness = parent.find_primitive_p_divisible_vector__random(prime)
                neighbor = parent.find_p_neighbor_from_vec(prime, witness)
                gram = neighbor.Hessian_matrix()
                key = hashlib.sha256(repr(tuple(gram.list())).encode()).hexdigest()
                new_path = path + [(int(prime), tuple(map(int, witness)))]
                if key in seen:
                    continue
                seen.add(key)
                generated += 1
                roots = int(pari(gram).qfminim(2)[0])
                root_histogram[roots] += 1
                masks = None
                support_signature = None
                if roots == 0:
                    tested_rootless += 1
                    if mask_cap:
                        mask_rows, _, _ = capped_mask_profile(
                            gram,
                            bridges,
                            reverse,
                            base,
                            prepared["order"],
                            mask_cap,
                        )
                        masks = tuple(
                            row["occupied_forbidden_cells_capped"]
                            for row in mask_rows
                        )
                        support_signature = tuple(
                            (
                                row["bridge_class_index"],
                                tuple(tuple(cell) for cell in row["occupied_support_capped"]),
                            )
                            for row in mask_rows
                        )
                    else:
                        mask_rows, _, _ = base["mask_profile"](
                            gram, bridges, reverse
                        )
                        masks = tuple(
                            row["occupied_forbidden_cells"] for row in mask_rows
                        )
                    if min(masks) == 0:
                        hit = {
                            "form": neighbor,
                            "gram": gram,
                            "path": new_path,
                            "mask_rows": mask_rows,
                            "mask_profile": masks,
                        }
                        break
                if root_descent or roots <= 4:
                    candidates.append(
                        (
                            0 if roots == 0 else 100 + 10 * roots,
                            min(masks) if masks else 99,
                            sum(masks) if masks else 999,
                            python_random.random(),
                            neighbor,
                            new_path,
                            roots,
                            masks,
                            support_signature,
                        )
                    )
            if hit is not None:
                break
        if hit is not None:
            record = {
                "generation": generation,
                "unique_raw_neighbors": generated,
                "tested_rootless_neighbors": tested_rootless,
                "root_count_histogram": {
                    str(key): root_histogram[key] for key in sorted(root_histogram)
                },
                "hit": True,
                "hit_mask_profile": list(hit["mask_profile"]),
            }
            generations.append(record)
            print(f"{corridor} generation {generation}: HIT {hit['mask_profile']}", flush=True)
            break

        candidates.sort(key=lambda row: row[:4])
        if support_diversity and any(row[6] == 0 for row in candidates):
            elite = []
            support_keys = set()
            for row in candidates:
                if row[6] != 0:
                    continue
                if row[8] in support_keys:
                    continue
                support_keys.add(row[8])
                elite.append(row)
                if len(elite) == base["SEARCH_ELITE"]:
                    break
            elite_paths = {tuple(row[5]) for row in elite}
            pool = [row for row in candidates if tuple(row[5]) not in elite_paths]
        elif persistent_elite:
            archive = sorted(archive + candidates, key=lambda row: row[:4])[
                : base["SEARCH_ELITE"]
            ]
            elite = archive
            elite_paths = {tuple(row[5]) for row in elite}
            pool = [row for row in candidates if tuple(row[5]) not in elite_paths]
        else:
            elite = candidates[: base["SEARCH_ELITE"]]
            pool = candidates[base["SEARCH_ELITE"] :]
        python_random.shuffle(pool)
        selected = elite + pool[: base["SEARCH_DIVERSITY"]]
        frontier = [(row[4], row[5]) for row in selected]
        record = {
            "generation": generation,
            "unique_raw_neighbors": generated,
            "tested_rootless_neighbors": tested_rootless,
            "root_count_histogram": {
                str(key): root_histogram[key] for key in sorted(root_histogram)
            },
            "eligible_neighbor_count": len(candidates),
            "frontier_size": len(frontier),
            "best_root_and_mask_profile": (
                [
                    selected[0][6],
                    None if selected[0][7] is None else list(selected[0][7]),
                ]
                if selected
                else None
            ),
            "hit": False,
        }
        generations.append(record)
        print(
            f"{corridor} generation {generation}: "
            f"generated={generated} rootless={tested_rootless} "
            f"root-range="
            f"{(min(root_histogram), max(root_histogram)) if root_histogram else None} "
            f"best={record['best_root_and_mask_profile']}",
            flush=True,
        )
        if not frontier:
            break

    result = {
        "corridor": corridor,
        "cyclic_bridge_order": prepared["order"],
        "target_frame_determinant": int(prepared["target_frame"].det()),
        "forced_core_determinant": int(seed.det()),
        "even_rank15_genera_at_determinant": len(prepared["all_genera"]),
        "matching_finite_form_genera": 1,
        "canonical_representative_signed_roots": int(pari(seed).qfminim(2)[0]),
        "all_bridge_class_indices": [
            row["bridge_class_index"] for row in prepared["bridges"]
        ],
        "rootless_bridge_class_indices_used_by_search": [
            row["bridge_class_index"] for row in bridges
        ],
        "generations": generations,
        "unique_raw_neighbors": len(seen),
        "zero_mask_core_found": hit is not None,
        "root_descent_enabled": root_descent,
        "search_mask_violation_cap": mask_cap or None,
        "persistent_elite_enabled": persistent_elite,
        "support_signature_diversity_enabled": support_diversity,
    }
    if hit is None:
        best = archive[0] if persistent_elite and archive else (selected[0] if selected else None)
        if best is not None:
            best_gram = base["lll_reduce"](best[4].Hessian_matrix())
            exact_masks, best_lazy_queries, _ = base["mask_profile"](
                best_gram, bridges, reverse
            )
            result["best_bounded_candidate"] = {
                "path": [
                    {"prime": prime, "witness": list(witness)}
                    for prime, witness in best[5]
                ],
                "reduced_core_gram": [
                    list(map(int, row)) for row in best_gram.rows()
                ],
                "signed_root_count": best[6],
                "exact_mask_results": exact_masks,
                "lazy_core_cells_queried": best_lazy_queries,
            }
        result["proof_boundary"] = (
            "No zero-mask core occurred in this fixed finite beam. This is not "
            "a genus-wide nonexistence result."
        )
        return result

    reduced = base["lll_reduce"](hit["gram"])
    configure_order(base, prepared["order"])
    masks, lazy_queries, core_generator = base["mask_profile"](
        reduced, bridges, reverse
    )
    accepted = next(row for row in masks if row["zero_mask_accepts"])
    bridge = next(
        row
        for row in bridges
        if row["bridge_class_index"] == accepted["bridge_class_index"]
    )
    multiplier = accepted["isotropic_multipliers"][0]
    glue = vector(QQ, list(multiplier * core_generator) + list(bridge["generator"]))
    child = core["glued_frame"](reduced, bridge["gram"], glue)
    reduced_child = base["lll_reduce"](child)
    reduced_target = base["lll_reduce"](prepared["target_frame"])
    reduced_historical_core = base["lll_reduce"](prepared["historical_core"])
    assert int(pari(reduced).qfminim(2)[0]) == 0
    assert core["minimum_norm"](reduced) == 4
    assert int(pari(child).qfminim(2)[0]) == 0
    assert core["minimum_norm"](child) == 4
    assert child.det() == prepared["target_frame"].det()
    child_isometric_to_target = bool(
        pari(reduced_child).qfisom(pari(reduced_target))
    )
    child_discriminant_matches_target = (
        core["discriminant_form_key"](child)
        == core["discriminant_form_key"](prepared["target_frame"])
    )
    assert child_discriminant_matches_target
    result["hit"] = {
        "generation": len(hit["path"]),
        "path": [
            {"prime": prime, "witness": list(witness)}
            for prime, witness in hit["path"]
        ],
        "reduced_core_gram": [list(map(int, row)) for row in reduced.rows()],
        "core_minimum": 4,
        "core_signed_roots": 0,
        "core_automorphism_group_order": int(pari(reduced).qfauto()[0]),
        "historical_core_isometric": bool(
            pari(reduced).qfisom(pari(reduced_historical_core))
        ),
        "mask_results": masks,
        "lazy_core_cells_queried": lazy_queries,
        "accepted_bridge_class_index": accepted["bridge_class_index"],
        "child_rank": child.nrows(),
        "child_determinant": int(child.det()),
        "child_minimum": 4,
        "child_signed_roots": 0,
        "child_discriminant_form_matches_target": child_discriminant_matches_target,
        "child_isometric_to_target_frame": child_isometric_to_target,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corridor",
        action="append",
        choices=DEFAULT_CORRIDORS,
        dest="corridors",
    )
    parser.add_argument(
        "--support-diversity",
        action="store_true",
        help="reserve elite slots for distinct occupied mask-support signatures",
    )
    parser.add_argument(
        "--persistent-elite",
        action="store_true",
        help="retain the globally best elite cores as parents across generations",
    )
    parser.add_argument(
        "--mask-cap",
        type=int,
        default=0,
        help="cap nonzero mask counts during ranking; zero tests remain exact",
    )
    parser.add_argument(
        "--generations",
        type=int,
        default=None,
        help="override the fixed eight-generation bound",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--root-descent",
        action="store_true",
        help="retain the least-rooted neighbors instead of imposing roots <= 4",
    )
    arguments = parser.parse_args()
    corridors = tuple(arguments.corridors or DEFAULT_CORRIDORS)
    if arguments.mask_cap < 0:
        parser.error("--mask-cap must be nonnegative")
    if arguments.support_diversity and not arguments.mask_cap:
        parser.error("--support-diversity requires a positive --mask-cap")

    bridge_artifact = json.loads(BRIDGES.read_text())
    theta_artifact = json.loads(THETA.read_text())
    base = runpy.run_path(str(BASE_SCRIPT))
    core = runpy.run_path(str(CORE_SCRIPT))
    reverse = runpy.run_path(str(REVERSE_SCRIPT))
    results = []
    for corridor in corridors:
        prepared = prepare_corridor(
            corridor, bridge_artifact, theta_artifact, base, core, reverse
        )
        results.append(
            search_corridor(
                corridor,
                prepared,
                base,
                core,
                reverse,
                arguments.root_descent,
                arguments.mask_cap,
                arguments.generations or base["SEARCH_GENERATIONS"],
                arguments.persistent_elite,
                arguments.support_diversity,
            )
        )

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-masked-core-controls.v1",
        "status": "PASS_FIXED_RULE_PROSPECTIVE_CONTROLS",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(BASE_SCRIPT): digest(BASE_SCRIPT),
            relative(CORE_SCRIPT): digest(CORE_SCRIPT),
            relative(REVERSE_SCRIPT): digest(REVERSE_SCRIPT),
        },
        "fixed_search_rule": {
            "random_seed_reset_per_corridor": base["SEARCH_SEED"],
            "primes": list(base["SEARCH_PRIMES"]),
            "generations": arguments.generations or base["SEARCH_GENERATIONS"],
            "samples_per_parent": base["SEARCH_SAMPLES_PER_PARENT"],
            "elite_slots": base["SEARCH_ELITE"],
            "diversity_slots": base["SEARCH_DIVERSITY"],
            "score": "root count first, then minimum and total occupied reverse-mask cells",
            "root_descent_enabled": arguments.root_descent,
            "mask_violation_cap": arguments.mask_cap or None,
            "zero_mask_acceptance_remains_exact_under_cap": True,
            "persistent_elite_enabled": arguments.persistent_elite,
            "support_signature_diversity_enabled": arguments.support_diversity,
        },
        "corridors": results,
        "proof_boundary": (
            "Each hit is checked by exact root enumeration, reverse-mask support, "
            "graph glue, and target-frame isometry. A miss closes only the fixed "
            "finite beam, not the forced genus."
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(relative(output))


if __name__ == "__main__":
    main()
