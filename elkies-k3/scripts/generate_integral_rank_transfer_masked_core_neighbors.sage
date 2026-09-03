#!/usr/bin/env sage-python
"""Generate a zero-mask rank-15 core from a forced genus representative.

The Golay-720 control is deliberately started from Sage's canonical
representative of the finite-form-forced genus, not from the historical core
class.  A deterministic mask-aware Kneser beam finds a compatible rootless
core after seven neighbor steps.  Normal mode replays the short path;
``--search`` reruns the bounded beam that discovered it.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random as python_random
import runpy

from sage.all import (
    Genus,
    QQ,
    ZZ,
    QuadraticForm,
    block_diagonal_matrix,
    lcm,
    matrix,
    pari,
    set_random_seed,
    vector,
)
from sage.quadratic_forms.genera.genus import genera


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-neighbors-v1.json"

CORRIDOR = "Golay720"
BRIDGE_ORDER = 23
SEARCH_SEED = 314159
SEARCH_PRIMES = (7, 11, 13, 17, 19)
SEARCH_GENERATIONS = 8
SEARCH_SAMPLES_PER_PARENT = 300
SEARCH_ELITE = 12
SEARCH_DIVERSITY = 8

# Exact short certificate extracted from the deterministic beam.  Each vector
# is expressed in the basis of the preceding quadratic form.
NEIGHBOR_PATH = (
    (17, (13, 10, 16, 13, 12, 9, 11, 2, 7, 11, 6, 7, 5, 3, 10)),
    (19, (12, 14, 18, 13, 4, 16, 0, 1, 18, 11, 12, 6, 1, 9, 4)),
    (11, (9, 2, 1, 1, 2, 8, 4, 6, 6, 9, 0, 3, 10, 2, 0)),
    (11, (1, 4, 8, 2, 0, 7, 8, 3, 9, 10, 2, 3, 4, 7, 6)),
    (17, (0, 0, 0, 13, 10, 7, 2, 2, 0, 15, 12, 0, 0, 14, 3)),
    (19, (0, 13, 14, 9, 18, 9, 6, 16, 9, 15, 15, 8, 10, 4, 13)),
    (19, (2, 16, 0, 17, 17, 14, 6, 11, 14, 15, 17, 1, 9, 10, 12)),
)

EXPECTED_REDUCED_CORE = matrix(
    ZZ,
    [
        [4, -2, 2, 2, 2, -1, -2, 0, -1, 1, 2, 2, 1, 0, -2],
        [-2, 4, 0, 0, -2, 2, 0, -1, 0, 0, -1, -1, -2, 1, 0],
        [2, 0, 4, 2, 1, 1, -2, 1, -2, 2, 0, 0, 0, 1, -2],
        [2, 0, 2, 4, 2, 0, -1, 1, -2, 0, 0, 0, 1, 0, -1],
        [2, -2, 1, 2, 4, -2, -1, 2, 0, 0, 1, 0, 1, -1, -1],
        [-1, 2, 1, 0, -2, 4, -1, -1, -1, 1, -1, -1, -1, 0, 0],
        [-2, 0, -2, -1, -1, -1, 4, 1, 0, -2, -1, -1, 0, 0, 1],
        [0, -1, 1, 1, 2, -1, 1, 4, 0, -1, 0, -1, 0, -1, 0],
        [-1, 0, -2, -2, 0, -1, 0, 0, 4, 0, 0, 1, 0, -1, 1],
        [1, 0, 2, 0, 0, 1, -2, -1, 0, 4, -1, 0, 0, 0, -1],
        [2, -1, 0, 0, 1, -1, -1, 0, 0, -1, 4, 2, -1, 0, -1],
        [2, -1, 0, 0, 0, -1, -1, -1, 1, 0, 2, 4, 0, 1, -1],
        [1, -2, 0, 1, 1, -1, 0, 0, 0, 0, -1, 0, 4, 0, 1],
        [0, 1, 1, 0, -1, 0, 0, -1, -1, 0, 0, 1, 0, 4, -1],
        [-2, 0, -2, -1, -1, 0, 1, 0, 1, -1, -1, -1, 1, -1, 4],
    ],
)


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def quadratic_form(gram):
    coefficients = [
        gram[left, left] // 2 if left == right else gram[left, right]
        for left in range(gram.nrows())
        for right in range(left, gram.ncols())
    ]
    value = QuadraticForm(ZZ, gram.nrows(), coefficients)
    assert value.Hessian_matrix() == gram
    return value


def finite_form_key(discriminant_form):
    normal = discriminant_form.normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rational_rows(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def fractional_part(value):
    value = QQ(value)
    return value - value.floor()


def modulo_two(value):
    value = QQ(value)
    return value - 2 * (value / 2).floor()


def primary_generator(gram, prime):
    """Find a generator of the unique cyclic prime-primary subgroup."""

    for row in gram.inverse().rows():
        order = lcm(value.denominator() for value in row)
        if order % prime:
            continue
        candidate = vector(
            QQ,
            [fractional_part(value) for value in (order // prime) * row],
        )
        if candidate not in ZZ**gram.nrows() and prime * candidate in ZZ**gram.nrows():
            return candidate
    raise ArithmeticError(f"no order-{prime} discriminant generator")


def lll_reduce(gram):
    transform = matrix(ZZ, pari(gram).qflllgram()).transpose()
    assert abs(transform.det()) == 1
    return transform * gram * transform.transpose()


def bridge_data(theta_record, reverse):
    answer = []
    for row in theta_record["classes"]:
        bridge = matrix(ZZ, row["bridge_gram"])
        generator = reverse["bridge_generator"](bridge)
        answer.append(
            {
                "bridge_class_index": int(row["bridge_class_index"]),
                "gram": bridge,
                "generator": generator,
                "theta_profile": reverse["theta_profile"](bridge),
            }
        )
    return answer


def mask_profile(gram, bridges, reverse, stop_at_first=False):
    """Return occupied forbidden cells for every admissible bridge class."""

    core_generator = primary_generator(gram, BRIDGE_ORDER)
    oracle = reverse["CoreCellOracle"](gram)
    rows = []
    for bridge in bridges:
        bridge_gram = bridge["gram"]
        bridge_generator = bridge["generator"]
        multipliers = [
            value
            for value in range(1, BRIDGE_ORDER)
            if modulo_two(
                (value * core_generator)
                * gram
                * (value * core_generator)
                + bridge_generator * bridge_gram * bridge_generator
            )
            == 0
        ]
        assert len(multipliers) == 2
        assert sum(multipliers) == BRIDGE_ORDER
        multiplier = multipliers[0]
        violations = 0
        mask_cells = set()
        for label in range(BRIDGE_ORDER):
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
                cell = (
                    reverse["canonical_signed_class"](core_class),
                    core_norm,
                )
                mask_cells.add(cell)
                if oracle.occupied(core_class, core_norm):
                    violations += 1
                    if stop_at_first:
                        break
            if stop_at_first and violations:
                break
        rows.append(
            {
                "bridge_class_index": bridge["bridge_class_index"],
                "isotropic_multipliers": multipliers,
                "sign_reduced_mask_cells": len(mask_cells),
                "occupied_forbidden_cells": violations,
                "zero_mask_accepts": violations == 0,
            }
        )
    return rows, len(oracle.cache), core_generator


def forced_genus(target_frame, selected_bridge, core_module_key):
    generated_key = finite_form_key(
        Genus(block_diagonal_matrix(target_frame, -selected_bridge)).discriminant_form()
    )
    determinant = abs(int(target_frame.det() * selected_bridge.det()))
    candidates = genera((15, 0), determinant, even=True)
    matches = [
        candidate
        for candidate in candidates
        if finite_form_key(candidate.discriminant_form()) == generated_key
    ]
    assert len(matches) == 1
    assert generated_key == core_module_key
    return candidates, matches[0]


def replay_path(seed_gram, bridges, reverse):
    current = quadratic_form(seed_gram)
    states = []
    for step, (prime, raw_vector) in enumerate(NEIGHBOR_PATH, start=1):
        assert seed_gram.det() % prime
        witness = vector(ZZ, raw_vector)
        assert any(value % prime for value in witness)
        assert current(witness) % prime == 0
        current = current.find_p_neighbor_from_vec(prime, witness)
        gram = current.Hessian_matrix()
        assert gram.det() == seed_gram.det()
        assert all(value % 2 == 0 for value in gram.diagonal())
        roots = int(pari(gram).qfminim(2)[0])
        masks = None
        if roots == 0:
            masks, _, _ = mask_profile(gram, bridges, reverse)
        states.append(
            {
                "step": step,
                "prime": prime,
                "witness": list(raw_vector),
                "signed_root_count": roots,
                "mask_violation_profile": (
                    None
                    if masks is None
                    else [row["occupied_forbidden_cells"] for row in masks]
                ),
            }
        )
    reduced = lll_reduce(current.Hessian_matrix())
    assert reduced == EXPECTED_REDUCED_CORE
    return reduced, states


def bounded_search(seed_gram, bridges, reverse):
    """Rerun the deterministic beam that produced NEIGHBOR_PATH."""

    set_random_seed(SEARCH_SEED)
    python_random.seed(SEARCH_SEED)
    frontier = [(quadratic_form(seed_gram), [])]
    seen = set()
    generation_records = []
    for generation in range(1, SEARCH_GENERATIONS + 1):
        candidates = []
        generated = 0
        for parent, path in frontier:
            for _ in range(SEARCH_SAMPLES_PER_PARENT):
                prime = python_random.choice(SEARCH_PRIMES)
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
                masks = None
                if roots == 0:
                    mask_rows, _, _ = mask_profile(gram, bridges, reverse)
                    masks = tuple(
                        row["occupied_forbidden_cells"] for row in mask_rows
                    )
                    if min(masks) == 0:
                        reduced = lll_reduce(gram)
                        assert reduced == EXPECTED_REDUCED_CORE
                        assert tuple(new_path) == NEIGHBOR_PATH
                        generation_records.append(
                            {
                                "generation": generation,
                                "unique_raw_neighbors": generated,
                                "hit": True,
                            }
                        )
                        return reduced, new_path, generation_records, len(seen)
                if roots <= 4:
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
                        )
                    )
        candidates.sort(key=lambda row: row[:4])
        elite = candidates[:SEARCH_ELITE]
        pool = candidates[SEARCH_ELITE:]
        python_random.shuffle(pool)
        selected = elite + pool[:SEARCH_DIVERSITY]
        frontier = [(row[4], row[5]) for row in selected]
        generation_records.append(
            {
                "generation": generation,
                "unique_raw_neighbors": generated,
                "eligible_root_count_at_most_four": len(candidates),
                "frontier_size": len(frontier),
                "best_root_and_mask_profile": [
                    selected[0][6],
                    None if selected[0][7] is None else list(selected[0][7]),
                ]
                if selected
                else None,
                "hit": False,
            }
        )
    raise ArithmeticError("bounded search did not reproduce the stored hit")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    bridge_artifact = json.loads(BRIDGES.read_text())
    theta_artifact = json.loads(THETA.read_text())
    reverse = runpy.run_path(str(REVERSE_SCRIPT))
    core_helpers = runpy.run_path(str(CORE_SCRIPT))
    edge = next(
        row
        for row in bridge_artifact["edges"]
        if row["corridor"] == CORRIDOR and int(row["target_root_rank"]) == 0
    )
    theta_record = next(
        row for row in theta_artifact["corridors"] if row["corridor"] == CORRIDOR
    )
    bridges = bridge_data(theta_record, reverse)
    selected_bridge = next(
        row["gram"] for row in bridges if row["bridge_class_index"] == 2
    )

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
    target_frame = core_helpers["glued_frame"](
        historical_core, historical_bridge, historical_glue
    )
    all_genera, generated_genus = forced_genus(
        target_frame,
        selected_bridge,
        finite_form_key(Genus(historical_core).discriminant_form()),
    )
    seed = matrix(ZZ, generated_genus.representative())
    assert seed.nrows() == 15 and seed.det() == EXPECTED_REDUCED_CORE.det()
    seed_roots = int(pari(seed).qfminim(2)[0])
    assert seed_roots == 96

    if arguments.search:
        generated_core, search_path, generation_records, unique_neighbors = bounded_search(
            seed, bridges, reverse
        )
        assert tuple(search_path) == NEIGHBOR_PATH
    else:
        generated_core, path_states = replay_path(seed, bridges, reverse)
        generation_records = None
        unique_neighbors = None
    if arguments.search:
        # Replaying the extracted certificate independently also checks every
        # intermediate arithmetic assertion and supplies stable path metadata.
        replayed_core, path_states = replay_path(seed, bridges, reverse)
        assert replayed_core == generated_core

    final_masks, lazy_queries, core_generator = mask_profile(
        generated_core, bridges, reverse
    )
    assert [row["occupied_forbidden_cells"] for row in final_masks] == [3, 0]
    assert lazy_queries == 13
    assert int(pari(generated_core).qfminim(2)[0]) == 0
    assert core_helpers["minimum_norm"](generated_core) == 4
    assert not pari(generated_core).qfisom(pari(historical_core))

    completing_bridge = next(
        row for row in bridges if row["bridge_class_index"] == 2
    )
    multiplier = final_masks[1]["isotropic_multipliers"][0]
    glue = vector(
        QQ,
        list(multiplier * core_generator) + list(completing_bridge["generator"]),
    )
    child = core_helpers["glued_frame"](
        generated_core, completing_bridge["gram"], glue
    )
    assert child.det() == target_frame.det() == 720
    assert int(pari(child).qfminim(2)[0]) == 0
    assert core_helpers["minimum_norm"](child) == 4
    assert pari(child).qfisom(pari(target_frame))

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-masked-core-neighbors.v1",
        "status": "PASS_CONSTRUCTIVE_MASK_AWARE_CORE_GENERATION",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(REVERSE_SCRIPT): digest(REVERSE_SCRIPT),
            relative(CORE_SCRIPT): digest(CORE_SCRIPT),
        },
        "forced_genus": {
            "corridor": CORRIDOR,
            "target_frame_determinant": int(target_frame.det()),
            "selected_bridge_class_index": 2,
            "bridge_determinant": int(selected_bridge.det()),
            "forced_core_determinant": int(seed.det()),
            "even_rank15_genera_at_determinant": len(all_genera),
            "matching_finite_form_genera": 1,
            "local_symbols": [str(value) for value in generated_genus.local_symbols()],
            "canonical_representative_signed_roots": seed_roots,
            "canonical_representative_is_rootless": False,
        },
        "neighbor_generation": {
            "method": "mask-aware good-prime Kneser beam",
            "random_seed": SEARCH_SEED,
            "primes": list(SEARCH_PRIMES),
            "samples_per_parent": SEARCH_SAMPLES_PER_PARENT,
            "elite_slots": SEARCH_ELITE,
            "diversity_slots": SEARCH_DIVERSITY,
            "hit_generation": len(NEIGHBOR_PATH),
            "search_reexecuted_for_this_artifact": bool(arguments.search),
            "search_generation_records": generation_records,
            "unique_raw_neighbors_before_hit": unique_neighbors,
            "certified_neighbor_path": path_states,
        },
        "generated_core": {
            "rank": generated_core.nrows(),
            "determinant": int(generated_core.det()),
            "minimum": 4,
            "signed_root_count": 0,
            "automorphism_group_order": int(pari(generated_core).qfauto()[0]),
            "gram": [list(map(int, row)) for row in generated_core.rows()],
            "historical_core_isometric": False,
            "distinct_from_historical_core_class": True,
            "lazy_core_cells_queried": lazy_queries,
            "bridge_mask_results": final_masks,
        },
        "completion_truth_check": {
            "bridge_class_index": 2,
            "core_glue_multiplier": multiplier,
            "child_rank": child.nrows(),
            "child_determinant": int(child.det()),
            "child_minimum": 4,
            "child_signed_root_count": 0,
            "child_isometric_to_target_frame": True,
            "child_discriminant_form_matches_target": (
                core_helpers["discriminant_form_key"](child)
                == core_helpers["discriminant_form_key"](target_frame)
            ),
        },
        "proof_boundary": {
            "proved": (
                "The finite-form-forced genus produces a canonical rootful "
                "representative; the seven certified good-prime neighbors "
                "produce a nonhistorical rootless rank-15 core; its class-2 "
                "reverse mask is empty on theta support; and the independently "
                "constructed rank-17 completion is rootless and isometric to "
                "the declared target frame."
            ),
            "bounded": (
                "The beam parameters are fixed and the discovery is proved "
                "only for the Golay-720 control. Completeness, expected running "
                "time, success probability, and a uniform algorithm for every "
                "forced genus are not proved."
            ),
        },
        "reproduce": {
            "short_certificate": (
                "sage -python elkies-k3/scripts/"
                "generate_integral_rank_transfer_masked_core_neighbors.sage --check"
            ),
            "full_bounded_search": (
                "sage -python elkies-k3/scripts/"
                "generate_integral_rank_transfer_masked_core_neighbors.sage "
                "--search --check"
            ),
        },
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        stored = json.loads(output.read_text())
        # Search telemetry is intentionally stored from the expensive full
        # run.  Short replay proves the path without rewriting those fields.
        if not arguments.search:
            payload["neighbor_generation"]["search_reexecuted_for_this_artifact"] = stored[
                "neighbor_generation"
            ]["search_reexecuted_for_this_artifact"]
            payload["neighbor_generation"]["search_generation_records"] = stored[
                "neighbor_generation"
            ]["search_generation_records"]
            payload["neighbor_generation"]["unique_raw_neighbors_before_hit"] = stored[
                "neighbor_generation"
            ]["unique_raw_neighbors_before_hit"]
        encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS constructive mask-aware core generation")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(relative(output))


if __name__ == "__main__":
    main()
