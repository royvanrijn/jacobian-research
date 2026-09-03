#!/usr/bin/env sage-python
"""Plan prescribed ADE mutations without sampling Kneser neighbours.

The planner works on the rank-15 core plus a fixed rank-two bridge.  A blind
benchmark fixture discloses the parent state, the good prime, the parent root
lines which must survive, and the desired child root metric, but not the
withheld isotropic line.  Candidate lines are generated deterministically on
the survivor-constrained projective quadric.  Death incidences are checked
before a candidate core is materialized; the completion is then classified
against the prescribed birth/root metric.  An optional diagnostic can also
require the withheld target core isometry class.  The separate H0l checker
remains the exact pre-materialization affine-CVP oracle.

This is a bounded target planner for core Kneser moves.  It does not turn a
core move into an elliptic-neighbour equation or a marked rational map.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import runpy
import time

from sage.all import GF, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CONTROL = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-controls-v1.json"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
Q80_COMPLETION = GENERATED / "elkies-k3-integral-rank-transfer-q80-defect-beam-v1.json"
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
INVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_ns0024_inverse_ade_mutation.sage"
SIGNATURE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_root_system_signature.sage"
OUTPUT = GENERATED / "elkies-k3-inverse-ade-target-planner-benchmark-v1.json"


def relative(path):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_projective(values, prime):
    values = tuple(int(value) % prime for value in values)
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return None
    inverse = pow(pivot, -1, prime)
    return tuple((value * inverse) % prime for value in values)


def signature_key(signature):
    return (
        signature["ade_type"],
        int(signature["root_rank"]),
        int(signature["root_line_count"]),
        tuple(
            sorted(
                (row["type"], int(row["rank"]), int(row["signed_root_count"]))
                for row in signature["components"]
            )
        ),
    )


def line_pairing(witness, gram, line, prime):
    integral_pairing = witness["core"] * gram
    assert integral_pairing in ZZ ** gram.nrows()
    return int(integral_pairing.dot_product(line) % prime)


def graph_state(gram, bridge, order, base, reverse, inverse):
    base["BRIDGE_ORDER"] = order
    search["configure_order"](base, order)
    masks, _, generator = base["mask_profile"](gram, [bridge], reverse)
    multiplier = int(masks[0]["isotropic_multipliers"][0])
    bridge_vectors = inverse["dual_vectors_through_two"](
        bridge["gram"], reverse["discriminant_class"]
    )
    core_vectors = inverse["dual_vectors_through_two"](
        gram, reverse["discriminant_class"]
    )
    witnesses = inverse["graph_root_witnesses"](
        core_vectors,
        bridge_vectors,
        generator,
        bridge["generator"],
        multiplier,
        order,
        gram.nrows(),
    )
    return {
        "generator": generator,
        "multiplier": multiplier,
        "bridge_vectors": bridge_vectors,
        "witnesses": witnesses,
    }


def completion_signature(gram, bridge, order, base, core, reverse, signatures):
    """Classify the canonical graph completion by direct rank-17 root enumeration."""

    search["configure_order"](base, order)
    masks, _, generator = base["mask_profile"](gram, [bridge], reverse)
    multiplier = int(masks[0]["isotropic_multipliers"][0])
    glue = vector(
        QQ,
        list(multiplier * generator) + list(bridge["generator"]),
    )
    frame = core["glued_frame"](gram, bridge["gram"], glue)
    return signatures["metric_root_signature"](frame)


def predict_child(gram, bridge, order, state, prime, line, inverse, birth_death, reverse):
    """Predict exactly the completion roots, querying only bridge-required cells."""

    adjusted = birth_death["adjusted_isotropic_lift"](gram, prime, line)
    kernel, pivot, pivot_inverse = birth_death["congruence_kernel_basis"](
        gram, prime, adjusted
    )
    kernel_gram = kernel * gram * kernel.transpose()
    kernel_lll = matrix(ZZ, pari(kernel_gram).qflllgram()).transpose()
    reduced_kernel = kernel_lll * kernel
    assert abs(reduced_kernel.det()) == prime
    graph_by_bridge_class = {}
    graph_labels = {}
    for label in range(order):
        core_class = reverse["discriminant_class"](
            label * state["multiplier"] * state["generator"]
        )
        bridge_class = reverse["discriminant_class"](
            label * bridge["generator"]
        )
        graph_by_bridge_class[bridge_class] = core_class
        graph_labels[bridge_class] = label

    shell_cache = {}
    signed = {}
    for bridge_value, bridge_norm, bridge_class in state["bridge_vectors"]:
        bridge_class = tuple(bridge_class)
        if bridge_class not in graph_by_bridge_class:
            continue
        core_class = vector(QQ, graph_by_bridge_class[bridge_class])
        target_norm = QQ(2) - bridge_norm
        if not 0 <= target_norm <= 2:
            continue
        cache_key = (tuple(core_class), target_norm)
        if cache_key not in shell_cache:
            dual_pairing = core_class * gram
            assert dual_pairing in ZZ ** gram.nrows()
            correction = vector(ZZ, [0] * gram.nrows())
            correction[pivot] = (
                -ZZ(dual_pairing.dot_product(adjusted)) * pivot_inverse
            ) % prime
            layer_vectors = []
            for layer in range(prime):
                shift = core_class + correction + QQ(layer) * adjusted / prime
                for core_value in birth_death["affine_vectors_of_norm"](
                    gram, reduced_kernel, shift, target_norm
                ):
                    layer_vectors.append((core_value, layer))
            shell_cache[cache_key] = layer_vectors
        for core_value, layer in shell_cache[cache_key]:
            canonical = inverse["canonical_line"](core_value, bridge_value)
            positive = tuple(core_value) + tuple(bridge_value)
            canonical_layer = layer if positive == canonical else (-layer) % prime
            canonical_bridge_class = reverse["discriminant_class"](
                canonical[gram.nrows():]
            )
            signed[canonical] = {
                "core": vector(QQ, canonical[: gram.nrows()]),
                "bridge": vector(QQ, canonical[gram.nrows():]),
                "graph_label": graph_labels[canonical_bridge_class],
                "affine_layer": canonical_layer,
            }
    witnesses = [signed[key] for key in sorted(signed)]
    signature = inverse["root_metric_signature"](
        witnesses, gram, bridge["gram"]
    )
    return adjusted, witnesses, signature, len(shell_cache)


def make_withheld_template(gram, bridge, order, prime, hidden_line, base, core, reverse, inverse, birth_death, signatures):
    """Compile a target template, deliberately omitting the source line."""

    state = graph_state(gram, bridge, order, base, reverse, inverse)
    line = vector(ZZ, hidden_line)
    survivor_indices = [
        index
        for index, witness in enumerate(state["witnesses"])
        if line_pairing(witness, gram, line, prime) == 0
    ]
    hidden_child = base["quadratic_form"](gram).find_p_neighbor_from_vec(
        prime, line
    ).Hessian_matrix()
    hidden_reduced = base["lll_reduce"](hidden_child)
    child_signature = completion_signature(
        hidden_reduced, bridge, order, base, core, reverse, signatures
    )
    return state, {
        "prime": int(prime),
        "surviving_parent_root_line_indices": survivor_indices,
        "parent_root_lines": len(state["witnesses"]),
        "desired_child_signature": {
            "ade_type": child_signature["ade_type"],
            "root_rank": child_signature["root_rank"],
            "root_line_count": child_signature["root_line_count"],
            "components": [
                {
                    "type": row["type"],
                    "rank": row["rank"],
                    "signed_root_count": row["signed_root_count"],
                }
                for row in child_signature["components"]
            ],
        },
        "desired_child_signature_key": signature_key(child_signature),
        "withheld_fields": [
            "selected_isotropic_line",
            "adjusted_isotropic_lift",
            "truth_child_gram",
        ],
    }, hidden_child


def sparse_projective_vectors(field, dimension, maximum_support):
    one = field.one()
    nonzero = list(field)[1:]
    for support_size in range(1, min(maximum_support, dimension) + 1):
        for support in combinations(range(dimension), support_size):
            for tail in product(nonzero, repeat=support_size - 1):
                coefficients = [field.zero()] * dimension
                coefficients[support[0]] = one
                for index, value in zip(support[1:], tail):
                    coefficients[index] = value
                yield vector(field, coefficients)


def constrained_quadric_lines(
    gram, prime, equality_rows, maximum_support, dense_probes, seed_material
):
    """Deterministically enumerate a bounded rational parametrization."""

    field = GF(prime)
    ambient_rank = gram.nrows()
    equations = matrix(field, equality_rows) if equality_rows else matrix(field, 0, ambient_rank)
    subspace = equations.right_kernel()
    basis = subspace.basis_matrix()
    restricted = basis * gram.change_ring(field) * basis.transpose()
    dimension = basis.nrows()

    seed = None
    seed_probes = 0
    for coefficients in sparse_projective_vectors(field, dimension, min(3, dimension)):
        seed_probes += 1
        candidate = coefficients * basis
        if candidate * gram.change_ring(field) * candidate == 0:
            seed = candidate
            break
    if seed is None:
        raise RuntimeError("failed to find an isotropic seed in the survivor subspace")

    seen = set()
    seed_key = canonical_projective(seed, prime)
    seen.add(seed_key)
    yield vector(ZZ, seed_key), {
        "survivor_subspace_dimension": dimension,
        "seed_probes": seed_probes,
        "parameter_probes": 0,
    }

    parameter_probes = 0
    for coefficients in sparse_projective_vectors(field, dimension, maximum_support):
        parameter_probes += 1
        parameter = coefficients * basis
        parameter_norm = parameter * gram.change_ring(field) * parameter
        pairing = seed * gram.change_ring(field) * parameter
        # For the bilinear Gram convention, this is the second intersection
        # of the line through the isotropic seed and ``parameter`` with Q=0.
        candidate = parameter_norm * seed - 2 * pairing * parameter
        key = canonical_projective(candidate, prime)
        if key is None or key in seen:
            continue
        seen.add(key)
        yield vector(ZZ, key), {
            "survivor_subspace_dimension": dimension,
            "seed_probes": seed_probes,
            "parameter_probes": parameter_probes,
        }

    for counter in range(dense_probes):
        encoded = hashlib.shake_256(
            seed_material + counter.to_bytes(8, "big")
        ).digest(2 * dimension)
        coefficients = vector(
            field,
            [int.from_bytes(encoded[2 * i : 2 * i + 2], "big") % prime for i in range(dimension)],
        )
        if not coefficients:
            continue
        parameter_probes += 1
        parameter = coefficients * basis
        parameter_norm = parameter * gram.change_ring(field) * parameter
        pairing = seed * gram.change_ring(field) * parameter
        candidate = parameter_norm * seed - 2 * pairing * parameter
        key = canonical_projective(candidate, prime)
        if key is None or key in seen:
            continue
        seen.add(key)
        yield vector(ZZ, key), {
            "survivor_subspace_dimension": dimension,
            "seed_probes": seed_probes,
            "parameter_probes": parameter_probes,
            "dense_parameter_probes": counter + 1,
        }


def plan_template(gram, bridge, order, state, template, target_core_gram, base, core, reverse, inverse, birth_death, signatures, maximum_support, dense_probes, require_core_isometry, max_materialized):
    prime = template["prime"]
    survivors = set(template["surviving_parent_root_line_indices"])
    equality_rows = []
    incidence_rows = []
    for index, witness in enumerate(state["witnesses"]):
        row = witness["core"] * gram
        finite_row = [int(value % prime) for value in row]
        incidence_rows.append(finite_row)
        if index in survivors:
            equality_rows.append(finite_row)
    incidence = matrix(GF(prime), incidence_rows)

    statistics = Counter()
    start = time.monotonic()
    last_meta = {}
    seed_material = json.dumps(
        {
            "prime": prime,
            "equalities": equality_rows,
            "target": template["desired_child_signature"],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    for line, meta in constrained_quadric_lines(
        gram, prime, equality_rows, maximum_support, dense_probes, seed_material
    ):
        last_meta = meta
        statistics["isotropic_lines_proposed"] += 1
        values = incidence * vector(GF(prime), line)
        actual_survivors = {index for index, value in enumerate(values) if value == 0}
        if actual_survivors != survivors:
            statistics["rejected_by_survival_death_incidence"] += 1
            continue
        # Incidence is the cheap target generator.  Only its survivors are
        # materialized, and then the bridge completion supplies the exact
        # born/additional-root metric used for acceptance.
        statistics["materialized_neighbors"] += 1
        if max_materialized and statistics["materialized_neighbors"] > max_materialized:
            statistics["materialized_neighbors"] -= 1
            statistics["stopped_at_materialized_budget"] = 1
            break
        child = base["quadratic_form"](gram).find_p_neighbor_from_vec(
            prime, line
        ).Hessian_matrix()
        reduced_child = base["lll_reduce"](child)
        if template["desired_child_signature"]["ade_type"] == "rootless":
            if int(pari(reduced_child).qfminim(2)[0]):
                statistics["rejected_by_birth_or_metric_signature"] += 1
                continue
            masks, _, _ = base["mask_profile"](
                reduced_child, [bridge], reverse, stop_at_first=True
            )
            if not masks[0]["zero_mask_accepts"]:
                statistics["rejected_by_birth_or_metric_signature"] += 1
                continue
            predicted_signature = {
                "ade_type": "rootless",
                "root_rank": 0,
                "root_line_count": 0,
                "components": [],
            }
        else:
            predicted_signature = completion_signature(
                reduced_child, bridge, order, base, core, reverse, signatures
            )
            if signature_key(predicted_signature) != tuple(template["desired_child_signature_key"]):
                statistics["rejected_by_birth_or_metric_signature"] += 1
                continue
        if require_core_isometry:
            assert target_core_gram is not None
            if not pari(reduced_child).qfisom(
                pari(base["lll_reduce"](target_core_gram))
            ):
                statistics["rejected_by_target_core_isometry"] += 1
                continue
        return {
            "status": "HIT",
            "selected_isotropic_line": list(map(int, line)),
            "predicted_child_ade_type": predicted_signature["ade_type"],
            "predicted_child_root_lines": predicted_signature["root_line_count"],
            "target_core_isometry_required": require_core_isometry,
            "statistics": dict(statistics),
            "quadric_enumerator": last_meta,
            "elapsed_seconds": time.monotonic() - start,
        }
    return {
        "status": "BOUNDED_MISS",
        "statistics": dict(statistics),
        "quadric_enumerator": last_meta,
        "elapsed_seconds": time.monotonic() - start,
    }


def replay_core_path(seed, path, base, continuation=()):
    current = base["quadratic_form"](seed)
    parents = []
    for row in path:
        gram = current.Hessian_matrix()
        line = vector(ZZ, row["witness"])
        parents.append((gram, int(row["prime"]), line))
        current = current.find_p_neighbor_from_vec(int(row["prime"]), line)
    if continuation:
        current = base["quadratic_form"](
            base["lll_reduce"](current.Hessian_matrix())
        )
        for row in continuation:
            gram = current.Hessian_matrix()
            line = vector(ZZ, row["witness"])
            parents.append((gram, int(row["prime"]), line))
            current = current.find_p_neighbor_from_vec(int(row["prime"]), line)
    return parents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corridor", action="append", choices=("H3", "NS0024", "Q80"))
    parser.add_argument("--maximum-parameter-support", type=int, default=3)
    parser.add_argument(
        "--maximum-edges",
        type=int,
        default=0,
        help="benchmark only the first N withheld edges (zero means all)",
    )
    parser.add_argument(
        "--dense-probes",
        type=int,
        default=150000,
        help="deterministic dense parameters after the sparse quadric prefix",
    )
    parser.add_argument(
        "--require-core-isometry",
        action="store_true",
        help="also require the withheld child core class (can be much more expensive)",
    )
    parser.add_argument(
        "--max-materialized",
        type=int,
        default=3000,
        help="stop an edge after this many completion classifications (zero is unlimited)",
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.maximum_parameter_support < 1:
        parser.error("--maximum-parameter-support must be positive")

    controls = json.loads(CONTROL.read_text())
    q80_completion = json.loads(Q80_COMPLETION.read_text())
    bridge_artifact = json.loads(BRIDGES.read_text())
    theta_artifact = json.loads(THETA.read_text())
    global search
    search = runpy.run_path(str(SEARCH_SCRIPT))
    base = runpy.run_path(str(search["BASE_SCRIPT"]))
    core = runpy.run_path(str(search["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(search["REVERSE_SCRIPT"]))
    inverse = runpy.run_path(str(INVERSE_SCRIPT))
    birth_death = runpy.run_path(str(inverse["BIRTH_DEATH_SCRIPT"]))
    signatures = runpy.run_path(str(SIGNATURE_SCRIPT))

    selected = arguments.corridor or ["H3", "NS0024", "Q80"]
    historical_lower_bounds = {
        "H3": 42300,
        "NS0024": 7477,
        "Q80": 42300 + 30228,
    }
    results = []
    for corridor in selected:
        prepared = search["prepare_corridor"](
            corridor, bridge_artifact, theta_artifact, base, core, reverse
        )
        control = next(row for row in controls["corridors"] if row["corridor"] == corridor)
        bridge_index = 2 if corridor == "Q80" else int(control["completion"]["bridge_class_index"])
        bridge = next(
            row for row in prepared["viable_bridges"]
            if int(row["bridge_class_index"]) == bridge_index
        )
        continuation = (
            q80_completion["hit"]["directed_path_after_near_miss"]
            if corridor == "Q80"
            else ()
        )
        parent_rows = replay_core_path(
            prepared["seed"], control["path"], base, continuation
        )
        edge_results = []
        for edge_index, (gram, prime, hidden_line) in enumerate(parent_rows, start=1):
            if arguments.maximum_edges and edge_index > arguments.maximum_edges:
                break
            state, template, truth_child_gram = make_withheld_template(
                gram,
                bridge,
                prepared["order"],
                prime,
                hidden_line,
                base,
                core,
                reverse,
                inverse,
                birth_death,
                signatures,
            )
            planned = plan_template(
                gram,
                bridge,
                prepared["order"],
                state,
                template,
                truth_child_gram if arguments.require_core_isometry else None,
                base,
                core,
                reverse,
                inverse,
                birth_death,
                signatures,
                arguments.maximum_parameter_support,
                arguments.dense_probes,
                arguments.require_core_isometry,
                arguments.max_materialized,
            )
            edge_results.append(
                {
                    "edge": edge_index,
                    "prime": prime,
                    "parent_ade_type": inverse["root_metric_signature"](
                        state["witnesses"], gram, bridge["gram"]
                    )["ade_type"],
                    "target_ade_type": template["desired_child_signature"]["ade_type"],
                    "surviving_parent_root_lines": len(
                        template["surviving_parent_root_line_indices"]
                    ),
                    "parent_root_lines": template["parent_root_lines"],
                    "withheld_line_not_in_planner_input": True,
                    "planner": planned,
                }
            )
            print(
                f"{corridor} edge {edge_index}: {planned['status']} "
                f"proposed={planned['statistics'].get('isotropic_lines_proposed', 0)} "
                f"cvp={planned['statistics'].get('affine_cvp_predictions', 0)} "
                f"materialized={planned['statistics'].get('materialized_neighbors', 0)}",
                flush=True,
            )
        recovered = all(row["planner"]["status"] == "HIT" for row in edge_results)
        materialized = sum(
            row["planner"]["statistics"].get("materialized_neighbors", 0)
            for row in edge_results
        )
        baseline = historical_lower_bounds[corridor]
        results.append(
            {
                "corridor": corridor,
                "edges": edge_results,
                "summary": {
                    "all_withheld_edges_recovered": recovered,
                    "materialized_completion_candidates": materialized,
                    "historical_raw_neighbor_candidates_lower_bound": baseline,
                    "candidate_reduction_factor_if_recovered": (
                        str(QQ(baseline) / materialized)
                        if recovered and materialized
                        else None
                    ),
                    "orders_of_magnitude_gate_at_least_10x": bool(
                        recovered and materialized and baseline >= 10 * materialized
                    ),
                },
            }
        )

    payload = {
        "schema": "elkies-k3.inverse-ade-target-planner-benchmark.v1",
        "status": (
            "PASS_BOUNDED_BLIND_TARGET_RECOVERY"
            if all(edge["planner"]["status"] == "HIT" for row in results for edge in row["edges"])
            else "BOUNDED_TARGET_PLANNER_MISSES_REMAIN"
        ),
        "inputs": {
            relative(CONTROL): digest(CONTROL),
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(Q80_COMPLETION): digest(Q80_COMPLETION),
            relative(SEARCH_SCRIPT): digest(SEARCH_SCRIPT),
            relative(INVERSE_SCRIPT): digest(INVERSE_SCRIPT),
            relative(SIGNATURE_SCRIPT): digest(SIGNATURE_SCRIPT),
        },
        "planner_rule": {
            "random_sampling": False,
            "maximum_parameter_support": arguments.maximum_parameter_support,
            "dense_probes_per_edge": arguments.dense_probes,
            "require_withheld_child_core_isometry": arguments.require_core_isometry,
            "maximum_materialized_per_edge": arguments.max_materialized,
            "selection_order": [
                "survivor linear subspace",
                "deterministic projective-quadric parametrization",
                "exact survivor/death incidence",
                "materialize only incidence survivors",
                "bridge-completion birth and full root-metric classification",
            ]
            + (
                ["target-core integral isometry"]
                if arguments.require_core_isometry
                else []
            ),
        },
        "corridors": results,
        "proof_boundary": (
            "Each edge is a retrospective blind-to-line benchmark: its prime, exact "
            "surviving parent root lines, and child ADE metric are prescribed; the "
            "historical isotropic line and child core Gram matrix are withheld. The "
            "bounded sparse parametrization is not a complete quadric enumeration. "
            "The search uses exact survival/death incidence before construction, but "
            "uses completion classification rather than the slower affine-CVP oracle "
            "for birth acceptance. These are rank-15 core Kneser moves, not "
            "elliptic-neighbour equations."
        ),
        "reproduce": (
            "sage -python elkies-k3/scripts/plan_inverse_ade_targets.sage "
            "--maximum-parameter-support {} --dense-probes {} "
            "--max-materialized {}{}"
        ).format(
            arguments.maximum_parameter_support,
            arguments.dense_probes,
            arguments.max_materialized,
            " --require-core-isometry" if arguments.require_core_isometry else "",
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        old = json.loads(output.read_text())
        # Timings are informative and deliberately excluded from stale checks.
        def strip_timings(value):
            if isinstance(value, dict):
                return {key: strip_timings(item) for key, item in value.items() if key != "elapsed_seconds"}
            if isinstance(value, list):
                return [strip_timings(item) for item in value]
            return value
        if strip_timings(old) != strip_timings(payload):
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
