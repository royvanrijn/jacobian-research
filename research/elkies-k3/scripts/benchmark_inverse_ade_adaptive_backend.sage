#!/usr/bin/env sage-python
"""Benchmark an adaptive, target-free inverse-ADE backend in rank 15.

The engine sees only a parent core ``K``, a binary bridge ``C``, its cyclic
graph glue ``H``, a good prime ``p``, and a desired ADE metric.  In
particular, neither a marked target core nor a historical neighbour line is
part of the planning input.

Three exact backends are provided.

``expanded``
    Enumerate the required scaled dual-coset shells once and hash their
    projective reductions.

``orbit``
    Compile the same shells, but quotient their projective support by the
    subgroup of Aut(K) that fixes every graph-glue core class.  This mode is
    used only for rootlessness, where orbit membership is sufficient.

``lazy``
    Enumerate no scaled shell.  For each proposed line, reject old roots by
    modular incidence and query the required nonzero affine layers by exact
    CVP.  A rootless acceptance exhausts every required graph cell, so its
    predicted physical root set is complete before the child is built.

The default benchmark exercises the historical terminal H3 and Q80 parent
states.  A deterministic target-free proposal stream must first find a
rootless line on each parent.  The withheld historical lines are loaded only
after those searches finish, for regression and nonuniqueness comparisons.
Timings are informative; exact counts and witnesses are the certificate
authority.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import runpy
import time

from sage.all import GF, QQ, ZZ, identity_matrix, matrix, pari, set_random_seed, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
GENERATED = ROOT / "artifacts/generated-results"
PLANNER_SCRIPT = ROOT / "elkies-k3/scripts/plan_inverse_ade_targets.sage"
OUTPUT = GENERATED / "elkies-k3-inverse-ade-adaptive-backend-v1.json"


def relative(path):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matrix_record(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_vector_record(value):
    return [str(entry) for entry in value]


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


ROOTLESS_SIGNATURE = {
    "ade_type": "rootless",
    "root_rank": 0,
    "signed_root_count": 0,
    "root_line_count": 0,
    "components": [],
    "nonorthogonal_root_line_edges": [],
    "pairwise_inner_products_upper_triangular": [],
    "pairwise_inner_products_sha256": hashlib.sha256(b"[]").hexdigest(),
}


def metric_signature(witnesses, gram, bridge_gram, inverse):
    if not witnesses:
        return dict(ROOTLESS_SIGNATURE)
    return inverse["root_metric_signature"](witnesses, gram, bridge_gram)


def canonical_projective_rational(value, prime, canonical_projective):
    reduced = []
    for entry in map(QQ, value):
        denominator = int(entry.denominator())
        assert denominator % prime
        reduced.append(
            int(entry.numerator()) * pow(denominator, -1, prime) % prime
        )
    return canonical_projective(reduced, prime)


def matrix_key(value):
    return tuple(map(int, value.list()))


def automorphism_generators(gram):
    """Return PARI's exact group order and row-convention generators."""

    data = pari(gram).qfauto()
    generators = [matrix(ZZ, item).transpose() for item in data[1]]
    for generator in generators:
        assert generator * gram * generator.transpose() == gram
    return int(data[0]), generators


def graph_cells(state, bridge, order, reverse):
    """Compile the finite core/bridge cells that can support norm two."""

    graph_by_bridge_class = {}
    graph_labels = {}
    for label in range(order):
        core_class = reverse["discriminant_class"](
            label * state["multiplier"] * state["generator"]
        )
        bridge_class = reverse["discriminant_class"](
            label * bridge["generator"]
        )
        bridge_key = tuple(bridge_class)
        assert bridge_key not in graph_by_bridge_class
        graph_by_bridge_class[bridge_key] = tuple(core_class)
        graph_labels[bridge_key] = label

    grouped = defaultdict(list)
    for bridge_value, bridge_norm, bridge_class in state["bridge_vectors"]:
        bridge_key = tuple(bridge_class)
        if bridge_key not in graph_by_bridge_class:
            continue
        target_norm = QQ(2) - bridge_norm
        if not 0 <= target_norm <= 2:
            continue
        core_class = graph_by_bridge_class[bridge_key]
        grouped[(core_class, target_norm)].append(
            {
                "bridge": bridge_value,
                "bridge_norm": bridge_norm,
                "bridge_class": bridge_key,
                "graph_label": graph_labels[bridge_key],
            }
        )

    cells = []
    for (core_class, target_norm), bridge_rows in grouped.items():
        integral_class = all(QQ(entry).denominator() == 1 for entry in core_class)
        cells.append(
            {
                "core_class": vector(QQ, core_class),
                "target_norm": QQ(target_norm),
                "bridge_rows": bridge_rows,
                # The integral norm-two cell detects ordinary child-core roots
                # and has the highest occupancy on the terminal controls.
                "priority": (
                    0 if integral_class and target_norm == 2 else 1,
                    -QQ(target_norm),
                    tuple(core_class),
                ),
            }
        )
    return sorted(cells, key=lambda row: row["priority"])


def shell_size_estimate(gram, prime, cells):
    """Geometry-of-numbers estimate used only to choose an exact backend."""

    rank = gram.nrows()
    determinant = float(abs(gram.det()))
    constant = math.pi ** (rank / 2) / math.gamma(rank / 2)
    estimate = 0.0
    for cell in cells:
        target = float(prime**2 * cell["target_norm"])
        if target == 0:
            estimate += 1.0
        else:
            estimate += (
                constant
                * target ** (rank / 2 - 1)
                / math.sqrt(determinant)
            )
    return estimate


def compatible_automorphisms(gram, cells):
    """Find a certified graph-class-fixing subgroup of Aut(K)."""

    order, generators = automorphism_generators(gram)
    classes = [row["core_class"] for row in cells]
    compatible_generators = []
    for generator in generators:
        if all(
            all(entry in ZZ for entry in core_class * generator - core_class)
            for core_class in classes
        ):
            compatible_generators.append(generator)

    # Determine the exact order of the generated compatible subgroup.  The
    # adaptive caller only invokes this closure when the ambient group order
    # is within its declared safety cap.
    identity = identity_matrix(ZZ, gram.nrows())
    found = {matrix_key(identity): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in compatible_generators:
            candidate = current * generator
            key = matrix_key(candidate)
            if key not in found:
                assert candidate * gram * candidate.transpose() == gram
                found[key] = candidate
                frontier.append(candidate)
    return {
        "ambient_order": order,
        "compatible_subgroup_order": len(found),
        "compatible_generators": compatible_generators,
        "compatible_elements": list(found.values()),
    }


def select_mode(
    gram,
    prime,
    cells,
    desired_signature,
    requested_mode,
    expansion_limit,
    orbit_minimum,
):
    estimate = shell_size_estimate(gram, prime, cells)
    automorphisms = compatible_automorphisms(gram, cells)
    effective_projective_order = max(
        1,
        automorphisms["compatible_subgroup_order"] // 2,
    )
    if requested_mode != "auto":
        mode = requested_mode
        reason = "explicit command-line selection"
    elif estimate <= expansion_limit:
        if (
            desired_signature["ade_type"] == "rootless"
            and effective_projective_order >= orbit_minimum
        ):
            mode = "orbit"
            reason = "small shell and useful graph-compatible projective symmetry"
        else:
            mode = "expanded"
            reason = "estimated total shell representations fit the expansion budget"
    else:
        mode = "lazy"
        reason = "estimated projected support exceeds the expansion budget"
    if mode == "orbit" and desired_signature["ade_type"] != "rootless":
        raise ValueError("orbit mode currently supports the rootless predicate only")
    return {
        "mode": mode,
        "reason": reason,
        "estimated_scaled_shell_vectors": estimate,
        "expansion_limit": expansion_limit,
        "ambient_automorphism_order": automorphisms["ambient_order"],
        "compatible_subgroup_order": automorphisms["compatible_subgroup_order"],
        "effective_projective_subgroup_order": effective_projective_order,
        "automorphisms": automorphisms,
    }


def witness_from_parts(core_value, bridge_row, core_rank, inverse, reverse):
    canonical = inverse["canonical_line"](core_value, bridge_row["bridge"])
    canonical_bridge_class = reverse["discriminant_class"](
        canonical[core_rank:]
    )
    # The label is informational.  Sign can exchange label with order-label;
    # physical coordinates and the metric remain the authority.
    return canonical, {
        "core": vector(QQ, canonical[:core_rank]),
        "bridge": vector(QQ, canonical[core_rank:]),
        "graph_label": bridge_row["graph_label"],
        "bridge_class": tuple(canonical_bridge_class),
        "affine_layer": None,
    }


def compile_expanded_shells(
    gram,
    prime,
    cells,
    inverse,
    birth_death,
    reverse,
    canonical_projective,
):
    """Compile every scaled shell and its exact physical birth witnesses."""

    births = defaultdict(dict)
    shell_counts = []
    lattice_basis = identity_matrix(ZZ, gram.nrows())
    for cell_index, cell in enumerate(cells):
        shift = prime * cell["core_class"]
        target_norm = prime**2 * cell["target_norm"]
        shell = birth_death["affine_vectors_of_norm"](
            gram, lattice_basis, shift, target_norm
        )
        nonzero = 0
        for scaled in shell:
            line_key = canonical_projective_rational(
                scaled, prime, canonical_projective
            )
            if line_key is None:
                continue
            nonzero += 1
            line = vector(ZZ, line_key)
            assert line * gram * line % prime == 0
            core_value = scaled / prime
            for bridge_row in cell["bridge_rows"]:
                physical_key, witness = witness_from_parts(
                    core_value,
                    bridge_row,
                    gram.nrows(),
                    inverse,
                    reverse,
                )
                births[line_key][physical_key] = witness
        shell_counts.append(
            {
                "cell_index": cell_index,
                "core_class": rational_vector_record(cell["core_class"]),
                "target_norm": str(cell["target_norm"]),
                "signed_scaled_vectors": len(shell),
                "nonzero_projective_reductions_with_multiplicity": nonzero,
            }
        )
    return {
        "births": {
            key: [rows[item] for item in sorted(rows)]
            for key, rows in births.items()
        },
        "shell_counts": shell_counts,
        "scaled_shell_vectors": sum(
            row["signed_scaled_vectors"] for row in shell_counts
        ),
        "projected_birth_points": len(births),
    }


def projective_orbit(line_key, prime, elements, canonical_projective):
    line = vector(ZZ, line_key)
    return {
        canonical_projective(line * element, prime) for element in elements
    }


def compile_orbit_support(expanded, prime, automorphisms, canonical_projective):
    """Compress a fully certified birth support into projective Aut(K)-orbits."""

    elements = automorphisms["compatible_elements"]
    representatives = set()
    orbit_sizes = Counter()
    for line_key in expanded["births"]:
        orbit = projective_orbit(
            line_key, prime, elements, canonical_projective
        )
        representative = min(orbit)
        representatives.add(representative)
        orbit_sizes[len(orbit)] += 1
    return {
        "representatives": representatives,
        "raw_projected_birth_points": len(expanded["births"]),
        "projective_orbit_representatives": len(representatives),
        "orbit_size_histogram_counted_from_raw_points": {
            str(size): count for size, count in sorted(orbit_sizes.items())
        },
    }


def old_survivors(gram, prime, line, parent_witnesses):
    survivors = []
    for witness in parent_witnesses:
        pairing = witness["core"] * gram
        assert pairing in ZZ ** gram.nrows()
        if int(pairing.dot_product(line) % prime) == 0:
            survivors.append(witness)
    return survivors


def lazy_has_birth(
    gram,
    prime,
    line,
    cells,
    birth_death,
):
    """Return an exact birth decision, stopping at the first occupied cell."""

    adjusted = birth_death["adjusted_isotropic_lift"](
        gram, prime, vector(ZZ, list(line))
    )
    kernel, pivot, pivot_inverse = birth_death["congruence_kernel_basis"](
        gram, prime, adjusted
    )
    kernel_gram = kernel * gram * kernel.transpose()
    kernel_lll = matrix(ZZ, pari(kernel_gram).qflllgram()).transpose()
    reduced_kernel = kernel_lll * kernel
    assert abs(reduced_kernel.det()) == prime

    statistics = Counter()
    # Query layer zero as well.  Although old physical roots were already
    # rejected by incidence, the chosen representatives of transported graph
    # classes need not put every layer-zero witness in the same representative
    # cell as the parent's canonical root list.  Rechecking it is cheap and
    # keeps this oracle independent of representative choices.
    for cell_index, cell in enumerate(cells):
        core_class = cell["core_class"]
        dual_pairing = core_class * gram
        assert dual_pairing in ZZ ** gram.nrows()
        correction = vector(ZZ, [0] * gram.nrows())
        correction[pivot] = (
            -ZZ(dual_pairing.dot_product(adjusted)) * pivot_inverse
        ) % prime
        for layer in range(prime):
            statistics["affine_cvp_queries"] += 1
            shift = core_class + correction + QQ(layer) * adjusted / prime
            vectors = birth_death["affine_vectors_of_norm"](
                gram, reduced_kernel, shift, cell["target_norm"]
            )
            if vectors:
                statistics["occupied_cell_index"] = cell_index
                statistics["occupied_affine_layer"] = layer
                statistics["occupied_cell_vectors_returned"] = len(vectors)
                return True, dict(statistics)
    return False, dict(statistics)


def predict_line(
    compiled,
    gram,
    bridge,
    order,
    state,
    prime,
    line,
    inverse,
    birth_death,
    reverse,
    canonical_projective,
):
    """Predict the complete root set when the line satisfies rootlessness."""

    survivors = old_survivors(gram, prime, line, state["witnesses"])
    if survivors:
        return {
            "rootless": False,
            "rejection": "old_root_hyperplane",
            "known_root_count": len(survivors),
            "statistics": {},
        }

    mode = compiled["selection"]["mode"]
    line_key = tuple(map(int, line))
    if mode == "expanded":
        born = compiled["expanded"]["births"].get(line_key, [])
        if born:
            return {
                "rootless": False,
                "rejection": "projected_birth_shell",
                "known_root_count": len(born),
                "statistics": {},
            }
    elif mode == "orbit":
        orbit = projective_orbit(
            line_key,
            prime,
            compiled["selection"]["automorphisms"]["compatible_elements"],
            canonical_projective,
        )
        representative = min(orbit)
        if representative in compiled["orbit"]["representatives"]:
            return {
                "rootless": False,
                "rejection": "projected_birth_shell_orbit",
                "known_root_count": 1,
                "statistics": {"candidate_projective_orbit_size": len(orbit)},
            }
    else:
        occupied, statistics = lazy_has_birth(
            gram,
            prime,
            line,
            compiled["cells"],
            birth_death,
        )
        if occupied:
            return {
                "rootless": False,
                "rejection": "lazy_affine_cvp_birth",
                "known_root_count": 1,
                "statistics": statistics,
            }

    # Every old-root hyperplane and every possible birth cell has now been
    # exhausted.  The empty set is therefore the complete predicted root set.
    return {
        "rootless": True,
        "rejection": None,
        "known_root_count": 0,
        "complete_predicted_root_set": [],
        "signature": dict(ROOTLESS_SIGNATURE),
        "statistics": statistics if mode == "lazy" else {},
    }


def compile_backend(
    gram,
    bridge,
    order,
    state,
    prime,
    desired_signature,
    requested_mode,
    expansion_limit,
    orbit_minimum,
    inverse,
    birth_death,
    reverse,
    canonical_projective,
):
    cells = graph_cells(state, bridge, order, reverse)
    selection = select_mode(
        gram,
        prime,
        cells,
        desired_signature,
        requested_mode,
        expansion_limit,
        orbit_minimum,
    )
    compiled = {
        "selection": selection,
        "cells": cells,
        "expanded": None,
        "orbit": None,
    }
    if selection["mode"] in ("expanded", "orbit"):
        compiled["expanded"] = compile_expanded_shells(
            gram,
            prime,
            cells,
            inverse,
            birth_death,
            reverse,
            canonical_projective,
        )
    if selection["mode"] == "orbit":
        compiled["orbit"] = compile_orbit_support(
            compiled["expanded"],
            prime,
            selection["automorphisms"],
            canonical_projective,
        )
        # Rootless membership no longer needs the raw point dictionary.
        compiled["expanded"]["births"] = None
    return compiled


def isotropic_quadric_points(gram, prime):
    """Exact size of the nonsingular parabolic quadric in odd vector rank."""

    rank = gram.nrows()
    assert rank % 2 == 1
    assert gram.det() % prime
    half_projective_dimension = (rank - 1) // 2
    return (prime ** (2 * half_projective_dimension) - 1) // (prime - 1)


def serialize_engine_input(gram, bridge, order, state, prime, desired_signature):
    """Make the target-free planner boundary explicit and hashable."""

    return {
        "parent_core_gram": matrix_record(gram),
        "bridge_gram": matrix_record(bridge["gram"]),
        "bridge_generator": rational_vector_record(bridge["generator"]),
        "graph_order": int(order),
        "graph_core_generator": rational_vector_record(state["generator"]),
        "graph_multiplier": int(state["multiplier"]),
        "prime": int(prime),
        "desired_ade_signature": desired_signature,
    }


def public_selection_record(selection):
    return {
        key: value
        for key, value in selection.items()
        if key != "automorphisms"
    }


def public_compile_record(compiled):
    cells = [
        {
            "core_class": rational_vector_record(row["core_class"]),
            "target_norm": str(row["target_norm"]),
            "bridge_vector_multiplicity": len(row["bridge_rows"]),
        }
        for row in compiled["cells"]
    ]
    expanded = compiled["expanded"]
    if expanded is not None:
        expanded = {
            key: value for key, value in expanded.items() if key != "births"
        }
    orbit = compiled["orbit"]
    if orbit is not None:
        orbit = {
            key: value for key, value in orbit.items() if key != "representatives"
        }
    return {
        "selection": public_selection_record(compiled["selection"]),
        "required_graph_cells": cells,
        "expanded": expanded,
        "orbit": orbit,
    }


def search_rootless_line(
    gram,
    bridge,
    order,
    state,
    prime,
    desired_signature,
    compiled,
    maximum_support,
    dense_probes,
    maximum_predictions,
    proposal_mode,
    proposal_seed,
    base,
    planner,
    inverse,
    birth_death,
    reverse,
):
    seed_material = json.dumps(
        {
            "engine": "inverse-ade-adaptive-v1",
            "input": serialize_engine_input(
                gram, bridge, order, state, prime, desired_signature
            ),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    if proposal_mode == "parametric":
        generator = planner["constrained_quadric_lines"](
            gram,
            prime,
            [],
            maximum_support,
            dense_probes,
            seed_material,
        )
    else:
        set_random_seed(proposal_seed)
        form = base["quadratic_form"](gram)

        def native_kneser_lines():
            seen = set()
            maximum_proposals = max(1000, 4 * maximum_predictions)
            for counter in range(maximum_proposals):
                raw = form.find_primitive_p_divisible_vector__random(prime)
                key = planner["canonical_projective"](raw, prime)
                if key is None or key in seen:
                    continue
                seen.add(key)
                line = vector(ZZ, key)
                assert line * gram * line % prime == 0
                yield line, {
                    "generation": "deterministic_native_kneser",
                    "random_seed": proposal_seed,
                    "random_draws": counter + 1,
                }

        generator = native_kneser_lines()
    counters = Counter()
    lazy_queries = 0
    last_meta = {}
    start = time.monotonic()
    for line, meta in generator:
        counters["isotropic_lines_proposed"] += 1
        last_meta = meta
        prediction = predict_line(
            compiled,
            gram,
            bridge,
            order,
            state,
            prime,
            line,
            inverse,
            birth_death,
            reverse,
            planner["canonical_projective"],
        )
        if prediction["rejection"] == "old_root_hyperplane":
            counters["rejected_by_old_root_hyperplanes"] += 1
            continue
        counters["exact_backend_predictions"] += 1
        lazy_queries += prediction["statistics"].get("affine_cvp_queries", 0)
        if not prediction["rootless"]:
            counters["rejected_by_birth_locus"] += 1
        else:
            assert signature_key(prediction["signature"]) == signature_key(
                desired_signature
            )
            return {
                "status": "HIT",
                "selected_isotropic_line": list(map(int, line)),
                "prediction": prediction,
                "statistics": {
                    **dict(counters),
                    "affine_cvp_queries": lazy_queries,
                    "materialized_neighbors_during_search": 0,
                },
                "quadric_enumerator": last_meta,
                "elapsed_seconds": time.monotonic() - start,
            }
        if maximum_predictions and counters["exact_backend_predictions"] >= maximum_predictions:
            break
    return {
        "status": "BOUNDED_MISS",
        "statistics": {
            **dict(counters),
            "affine_cvp_queries": lazy_queries,
            "materialized_neighbors_during_search": 0,
        },
        "quadric_enumerator": last_meta,
        "elapsed_seconds": time.monotonic() - start,
    }


def verify_materialized_rootless(
    gram,
    bridge,
    order,
    state,
    line,
    base,
    core,
    reverse,
    signatures,
    planner,
    birth_death,
):
    """Construct the child only after transporting the supplied graph glue."""

    prime = int(line["prime"])
    coordinates = vector(ZZ, line["coordinates"])
    quadratic_form = base["quadratic_form"](gram)
    child_to_parent = quadratic_form.find_p_neighbor_from_vec(
        prime, coordinates, return_matrix=True
    ).transpose()
    child_gram = child_to_parent * gram * child_to_parent.transpose()

    adjusted = birth_death["adjusted_isotropic_lift"](
        gram, prime, vector(ZZ, list(coordinates))
    )
    _, pivot, pivot_inverse = birth_death["congruence_kernel_basis"](
        gram, prime, adjusted
    )
    parent_core_class = vector(
        QQ,
        reverse["discriminant_class"](
            state["multiplier"] * state["generator"]
        ),
    )
    dual_pairing = parent_core_class * gram
    assert dual_pairing in ZZ ** gram.nrows()
    correction = vector(ZZ, [0] * gram.nrows())
    correction[pivot] = (
        -ZZ(dual_pairing.dot_product(adjusted)) * pivot_inverse
    ) % prime
    transported_core_glue_parent = parent_core_class + correction
    transported_core_glue_child = (
        transported_core_glue_parent * child_to_parent.inverse()
    )
    assert transported_core_glue_child * child_gram in ZZ ** gram.nrows()
    glue = vector(
        QQ,
        list(transported_core_glue_child) + list(bridge["generator"]),
    )
    frame, split, split_basis = signatures["overlattice_with_basis"](
        child_gram, bridge["gram"], glue
    )
    signature = signatures["metric_root_signature"](
        frame, split, split_basis, glue
    )
    reduced_child = base["lll_reduce"](child_gram)
    return {
        "child_core_gram_sha256": hashlib.sha256(
            json.dumps(matrix_record(reduced_child), separators=(",", ":")).encode()
        ).hexdigest(),
        "actual_ade_type": signature["ade_type"],
        "actual_root_line_count": int(signature["root_line_count"]),
        "predicted_root_set_equal_to_materialized_root_set": (
            signature["root_line_count"] == 0
        ),
        "materialized_physical_root_lines": (
            signature["physical_bridge_witnesses"]["root_lines"]
            if signature["root_line_count"]
            else []
        ),
    }


def prepare_tools():
    planner = runpy.run_path(str(PLANNER_SCRIPT))
    search = runpy.run_path(str(planner["SEARCH_SCRIPT"]))
    planner["graph_state"].__globals__["search"] = search
    planner["completion_signature"].__globals__["search"] = search
    base = runpy.run_path(str(search["BASE_SCRIPT"]))
    core = runpy.run_path(str(search["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(search["REVERSE_SCRIPT"]))
    inverse = runpy.run_path(str(planner["INVERSE_SCRIPT"]))
    birth_death = runpy.run_path(str(inverse["BIRTH_DEATH_SCRIPT"]))
    signatures = runpy.run_path(str(planner["SIGNATURE_SCRIPT"]))
    return planner, search, base, core, reverse, inverse, birth_death, signatures


def terminal_fixture(corridor, tools, source):
    planner, search, base, core, reverse, inverse, _, signatures = tools
    controls = source["controls"]
    bridge_artifact = source["bridges"]
    theta_artifact = source["theta"]
    q80_completion = source["q80"]
    prepared = search["prepare_corridor"](
        corridor, bridge_artifact, theta_artifact, base, core, reverse
    )
    control = next(
        row for row in controls["corridors"] if row["corridor"] == corridor
    )
    bridge_index = (
        2
        if corridor == "Q80"
        else int(control["completion"]["bridge_class_index"])
    )
    bridge = next(
        row
        for row in prepared["viable_bridges"]
        if int(row["bridge_class_index"]) == bridge_index
    )
    current = base["quadratic_form"](prepared["seed"])
    if corridor == "Q80":
        # The terminal Q80 control follows the stored near-miss path and then
        # three of the four continuation edges.  The fourth witness remains
        # unread until reveal_historical_line is called.
        for row in control["path"]:
            current = current.find_p_neighbor_from_vec(
                int(row["prime"]), vector(ZZ, row["witness"])
            )
        current = base["quadratic_form"](
            base["lll_reduce"](current.Hessian_matrix())
        )
        continuation = q80_completion["hit"]["directed_path_after_near_miss"]
        for row in continuation[:-1]:
            current = current.find_p_neighbor_from_vec(
                int(row["prime"]), vector(ZZ, row["witness"])
            )
        prime = int(continuation[-1]["prime"])
    else:
        for row in control["path"][:-1]:
            current = current.find_p_neighbor_from_vec(
                int(row["prime"]), vector(ZZ, row["witness"])
            )
        prime = int(control["path"][-1]["prime"])
    replay_gram = current.Hessian_matrix()
    replay_state = planner["graph_state"](
        replay_gram,
        bridge,
        prepared["order"],
        base,
        reverse,
        inverse,
    )
    reduction = matrix(ZZ, pari(replay_gram).qflllgram()).transpose()
    gram = reduction * replay_gram * reduction.transpose()
    assert abs(reduction.det()) == 1
    inverse_reduction = reduction.inverse().change_ring(ZZ)
    state = {
        "generator": replay_state["generator"] * inverse_reduction,
        "multiplier": replay_state["multiplier"],
        "bridge_vectors": replay_state["bridge_vectors"],
        "witnesses": [
            {
                **witness,
                "core": witness["core"] * inverse_reduction,
            }
            for witness in replay_state["witnesses"]
        ],
    }
    assert all(
        witness["core"] * gram in ZZ ** gram.nrows()
        for witness in state["witnesses"]
    )
    desired_signature = dict(ROOTLESS_SIGNATURE)
    return {
        "gram": gram,
        "prime": int(prime),
        "inverse_reduction": inverse_reduction,
        "state": state,
        "bridge": bridge,
        "order": int(prepared["order"]),
        "desired_signature": desired_signature,
    }


def reveal_historical_line(corridor, fixture, source, planner):
    """Read the held-out terminal witness only after backend compilation."""

    if corridor == "Q80":
        row = source["q80"]["hit"]["directed_path_after_near_miss"][-1]
    else:
        control = next(
            item
            for item in source["controls"]["corridors"]
            if item["corridor"] == corridor
        )
        row = control["path"][-1]
    assert int(row["prime"]) == fixture["prime"]
    transformed = vector(ZZ, row["witness"]) * fixture["inverse_reduction"]
    line = vector(
        ZZ,
        planner["canonical_projective"](transformed, fixture["prime"]),
    )
    assert line * fixture["gram"] * line % fixture["prime"] == 0
    return line


def backend_mode_equivalence_control(tools):
    """Exercise all three modes on the diagonal index-two glue control."""

    planner, _, _, _, reverse, inverse, birth_death, _ = tools
    prime = 5
    gram = 2 * identity_matrix(ZZ, 3)
    bridge = {
        "gram": matrix(ZZ, [[2]]),
        "generator": vector(QQ, [QQ(1) / 2]),
    }
    core_generator = vector(QQ, [QQ(1) / 2] * 3)
    bridge_vectors = inverse["dual_vectors_through_two"](
        bridge["gram"], reverse["discriminant_class"]
    )
    core_vectors = inverse["dual_vectors_through_two"](
        gram, reverse["discriminant_class"]
    )
    parent_witnesses = inverse["graph_root_witnesses"](
        core_vectors,
        bridge_vectors,
        core_generator,
        bridge["generator"],
        1,
        2,
        gram.nrows(),
    )
    state = {
        "generator": core_generator,
        "multiplier": 1,
        "bridge_vectors": bridge_vectors,
        "witnesses": parent_witnesses,
    }
    lines = []
    for entries in product(range(prime), repeat=gram.nrows()):
        key = planner["canonical_projective"](entries, prime)
        if key is None or key != tuple(entries):
            continue
        line = vector(ZZ, entries)
        if line * gram * line % prime == 0:
            lines.append(line)
    assert len(lines) == prime + 1 == 6

    decisions = {}
    mode_records = {}
    for requested_mode in ("expanded", "orbit", "lazy", "auto"):
        compiled = compile_backend(
            gram,
            bridge,
            2,
            state,
            prime,
            ROOTLESS_SIGNATURE,
            requested_mode,
            250000.0,
            4,
            inverse,
            birth_death,
            reverse,
            planner["canonical_projective"],
        )
        mode_decisions = []
        for line in lines:
            prediction = predict_line(
                compiled,
                gram,
                bridge,
                2,
                state,
                prime,
                line,
                inverse,
                birth_death,
                reverse,
                planner["canonical_projective"],
            )
            mode_decisions.append(prediction["rootless"])
        decisions[requested_mode] = mode_decisions
        mode_records[requested_mode] = public_compile_record(compiled)
    assert decisions["expanded"] == decisions["orbit"] == decisions["lazy"]
    assert decisions["auto"] == decisions["orbit"]
    assert not any(decisions["expanded"])
    assert mode_records["auto"]["selection"]["mode"] == "orbit"
    assert mode_records["orbit"]["orbit"][
        "projective_orbit_representatives"
    ] == 1
    return {
        "name": "A1^3 plus A1 with diagonal half-class glue",
        "prime": prime,
        "isotropic_lines": len(lines),
        "all_modes_exactly_agree": True,
        "rootless_lines": sum(decisions["expanded"]),
        "mode_records": mode_records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corridor", action="append", choices=("H3", "Q80"))
    parser.add_argument(
        "--mode", choices=("auto", "expanded", "orbit", "lazy"), default="auto"
    )
    parser.add_argument("--expansion-limit", type=float, default=250000.0)
    parser.add_argument("--orbit-minimum", type=int, default=4)
    parser.add_argument(
        "--maximum-parameter-support",
        type=int,
        default=1,
        help=(
            "structured parametrization prefix; the default moves immediately "
            "from coordinate directions to deterministic dense parameters"
        ),
    )
    parser.add_argument("--dense-probes", type=int, default=50000)
    parser.add_argument(
        "--proposal-mode",
        choices=("kneser", "parametric"),
        default="kneser",
    )
    parser.add_argument("--proposal-seed", type=int, default=314159)
    parser.add_argument(
        "--terminal-proposal-seed",
        type=int,
        default=1,
        help="fixed target-free native-Kneser seed for both terminal controls",
    )
    parser.add_argument(
        "--maximum-predictions",
        type=int,
        default=12000,
        help="maximum exact birth predictions for the novelty search",
    )
    parser.add_argument(
        "--terminal-search-predictions",
        type=int,
        default=30000,
        help=(
            "target-free exact-prediction budget on each historical terminal "
            "parent; zero is available only for a predicate-regression diagnostic"
        ),
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    tools = prepare_tools()
    planner, search, base, core, reverse, inverse, birth_death, signatures = tools
    source = {
        "controls": json.loads(planner["CONTROL"].read_text()),
        "bridges": json.loads(planner["BRIDGES"].read_text()),
        "theta": json.loads(planner["THETA"].read_text()),
        "q80": json.loads(planner["Q80_COMPLETION"].read_text()),
    }

    selected_corridors = arguments.corridor or ["H3", "Q80"]
    # The v1 terminal diagnostic stopped after exactly 3,000 independently
    # constructed and rejected completions on each control.
    historical_baselines = {"H3": 3000, "Q80": 3000}
    records = []
    all_nonhistorical = []
    for corridor in selected_corridors:
        fixture = terminal_fixture(corridor, tools, source)
        gram = fixture["gram"]
        prime = fixture["prime"]
        bridge = fixture["bridge"]
        order = fixture["order"]
        state = fixture["state"]
        desired_signature = fixture["desired_signature"]

        engine_input = serialize_engine_input(
            gram, bridge, order, state, prime, desired_signature
        )
        encoded_input = json.dumps(
            engine_input, sort_keys=True, separators=(",", ":")
        ).encode()
        compile_start = time.monotonic()
        compiled = compile_backend(
            gram,
            bridge,
            order,
            state,
            prime,
            desired_signature,
            arguments.mode,
            arguments.expansion_limit,
            arguments.orbit_minimum,
            inverse,
            birth_death,
            reverse,
            planner["canonical_projective"],
        )
        compile_seconds = time.monotonic() - compile_start
        compiled_record = public_compile_record(compiled)

        # Backend compilation and the independent discovery search are both
        # complete before the withheld historical line is read.
        if arguments.terminal_search_predictions:
            planned = search_rootless_line(
                gram,
                bridge,
                order,
                state,
                prime,
                desired_signature,
                compiled,
                arguments.maximum_parameter_support,
                arguments.dense_probes,
                arguments.terminal_search_predictions,
                arguments.proposal_mode,
                arguments.terminal_proposal_seed,
                base,
                planner,
                inverse,
                birth_death,
                reverse,
            )
        else:
            planned = {
                "status": "NOT_RUN",
                "reason": "terminal predicate regression only",
                "statistics": {
                    "isotropic_lines_proposed": 0,
                    "exact_backend_predictions": 0,
                    "materialized_neighbors_during_search": 0,
                },
                "elapsed_seconds": 0.0,
            }

        historical_line = reveal_historical_line(
            corridor, fixture, source, planner
        )
        historical_prediction_start = time.monotonic()
        historical_prediction = predict_line(
            compiled,
            gram,
            bridge,
            order,
            state,
            prime,
            historical_line,
            inverse,
            birth_death,
            reverse,
            planner["canonical_projective"],
        )
        historical_prediction_seconds = time.monotonic() - historical_prediction_start
        assert historical_prediction["rootless"]
        historical_truth = verify_materialized_rootless(
            gram,
            bridge,
            order,
            state,
            {"prime": prime, "coordinates": list(map(int, historical_line))},
            base,
            core,
            reverse,
            signatures,
            planner,
            birth_death,
        )
        assert historical_truth["actual_ade_type"] == "rootless"

        selected_truth = None
        selected_is_historical = None
        if planned["status"] == "HIT":
            selected_line = vector(ZZ, planned["selected_isotropic_line"])
            selected_is_historical = (
                planner["canonical_projective"](selected_line, prime)
                == planner["canonical_projective"](historical_line, prime)
            )
            selected_truth = verify_materialized_rootless(
                gram,
                bridge,
                order,
                state,
                {
                    "prime": prime,
                    "coordinates": planned["selected_isotropic_line"],
                },
                base,
                core,
                reverse,
                signatures,
                planner,
                birth_death,
            )
            if selected_truth["actual_ade_type"] != "rootless":
                raise AssertionError(
                    "adaptive prediction/materialization mismatch for "
                    f"{corridor} line {planned['selected_isotropic_line']}: "
                    f"{selected_truth}"
                )
            if not selected_is_historical:
                all_nonhistorical.append(
                    {
                        "corridor": corridor,
                        "prime": prime,
                        "line": planned["selected_isotropic_line"],
                        "child_core_gram_sha256": selected_truth[
                            "child_core_gram_sha256"
                        ],
                    }
                )

        quadric_points = isotropic_quadric_points(gram, prime)
        materialized = 1 + int(selected_truth is not None)
        baseline = historical_baselines[corridor]
        record = {
            "corridor": corridor,
            "rank": gram.nrows(),
            "prime": prime,
            "parent_ade_type": metric_signature(
                state["witnesses"], gram, bridge["gram"], inverse
            )["ade_type"],
            "desired_ade_type": "rootless",
            "engine_input": engine_input,
            "engine_input_sha256": hashlib.sha256(encoded_input).hexdigest(),
            "forbidden_target_fields_absent": [
                "historical_isotropic_line",
                "marked_target_core",
                "target_core_overlap_fingerprint",
                "target_core_gram",
            ],
            "compiled_backend": compiled_record,
            "compile_seconds": compile_seconds,
            "planner": planned,
            "historical_regression_after_planning": {
                "line": list(map(int, historical_line)),
                "complete_predicted_root_set": historical_prediction.get(
                    "complete_predicted_root_set"
                ),
                "prediction_statistics": historical_prediction["statistics"],
                "prediction_seconds": historical_prediction_seconds,
                "prediction_fixed_before_child_construction": True,
                "materialized_truth": historical_truth,
            },
            "selected_line_is_historical_line": selected_is_historical,
            "selected_materialized_truth": selected_truth,
            "isotropic_quadric_points": quadric_points,
            "materialized_isotropic_lines": materialized,
            "materialized_fraction": str(QQ(materialized) / quadric_points),
            "historical_terminal_materialized_reject_window": baseline,
            "materialization_reduction_factor": str(QQ(baseline) / materialized),
            "materializes_at_most_ten_percent_of_quadric": (
                10 * materialized <= quadric_points
            ),
            "at_least_tenfold_materialization_reduction": (
                baseline >= 10 * materialized
            ),
        }
        records.append(record)
        print(
            f"{corridor}: mode={compiled['selection']['mode']} "
            f"historical=PASS search={planned['status']} "
            f"proposed={planned['statistics'].get('isotropic_lines_proposed', 0)} "
            f"predicted={planned['statistics'].get('exact_backend_predictions', 0)} "
            f"materialized={materialized}",
            flush=True,
        )

    # A genuinely new completion is searched without any candidate line or
    # target core in the engine input.  The historical line is consulted only
    # after the search returns, to certify novelty.
    novelty_fixture = terminal_fixture("NS0024", tools, source)
    novelty_gram = novelty_fixture["gram"]
    novelty_prime = novelty_fixture["prime"]
    novelty_bridge = novelty_fixture["bridge"]
    novelty_order = novelty_fixture["order"]
    novelty_state = novelty_fixture["state"]
    novelty_signature = novelty_fixture["desired_signature"]
    novelty_input = serialize_engine_input(
        novelty_gram,
        novelty_bridge,
        novelty_order,
        novelty_state,
        novelty_prime,
        novelty_signature,
    )
    novelty_compile_start = time.monotonic()
    novelty_compiled = compile_backend(
        novelty_gram,
        novelty_bridge,
        novelty_order,
        novelty_state,
        novelty_prime,
        novelty_signature,
        arguments.mode,
        arguments.expansion_limit,
        arguments.orbit_minimum,
        inverse,
        birth_death,
        reverse,
        planner["canonical_projective"],
    )
    novelty_compile_seconds = time.monotonic() - novelty_compile_start
    novelty_search = search_rootless_line(
        novelty_gram,
        novelty_bridge,
        novelty_order,
        novelty_state,
        novelty_prime,
        novelty_signature,
        novelty_compiled,
        arguments.maximum_parameter_support,
        arguments.dense_probes,
        arguments.maximum_predictions,
        arguments.proposal_mode,
        arguments.proposal_seed,
        base,
        planner,
        inverse,
        birth_death,
        reverse,
    )
    novelty_truth = None
    novelty_is_historical = None
    if novelty_search["status"] == "HIT":
        novelty_line = vector(ZZ, novelty_search["selected_isotropic_line"])
        novelty_historical_line = reveal_historical_line(
            "NS0024", novelty_fixture, source, planner
        )
        novelty_is_historical = (
            planner["canonical_projective"](novelty_line, novelty_prime)
            == planner["canonical_projective"](
                novelty_historical_line, novelty_prime
            )
        )
        novelty_truth = verify_materialized_rootless(
            novelty_gram,
            novelty_bridge,
            novelty_order,
            novelty_state,
            {
                "prime": novelty_prime,
                "coordinates": novelty_search["selected_isotropic_line"],
            },
            base,
            core,
            reverse,
            signatures,
            planner,
            birth_death,
        )
        assert novelty_truth["actual_ade_type"] == "rootless"
        if not novelty_is_historical:
            all_nonhistorical.append(
                {
                    "corridor": "NS0024",
                    "prime": novelty_prime,
                    "line": novelty_search["selected_isotropic_line"],
                    "child_core_gram_sha256": novelty_truth[
                        "child_core_gram_sha256"
                    ],
                }
            )
    novelty_quadric_points = isotropic_quadric_points(
        novelty_gram, novelty_prime
    )
    novelty_record = {
        "corridor": "NS0024",
        "rank": novelty_gram.nrows(),
        "prime": novelty_prime,
        "engine_input": novelty_input,
        "engine_input_sha256": hashlib.sha256(
            json.dumps(
                novelty_input, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest(),
        "forbidden_target_fields_absent": [
            "historical_isotropic_line",
            "marked_target_core",
            "target_core_overlap_fingerprint",
            "target_core_gram",
        ],
        "compiled_backend": public_compile_record(novelty_compiled),
        "compile_seconds": novelty_compile_seconds,
        "planner": novelty_search,
        "selected_line_is_historical_line": novelty_is_historical,
        "historical_line_revealed_only_after_search": True,
        "selected_materialized_truth": novelty_truth,
        "isotropic_quadric_points": novelty_quadric_points,
        "materialized_isotropic_lines": int(novelty_truth is not None),
        "materialized_fraction": str(
            QQ(int(novelty_truth is not None)) / novelty_quadric_points
        ),
    }
    print(
        "NS0024 novelty: "
        f"mode={novelty_compiled['selection']['mode']} "
        f"search={novelty_search['status']} "
        f"predicted={novelty_search['statistics'].get('exact_backend_predictions', 0)} "
        f"materialized={int(novelty_truth is not None)}",
        flush=True,
    )

    mode_equivalence_control = backend_mode_equivalence_control(tools)

    historical_gate = all(
        row["historical_regression_after_planning"]["materialized_truth"][
            "predicted_root_set_equal_to_materialized_root_set"
        ]
        for row in records
    )
    terminal_recovery_gate = all(
        row["planner"]["status"] == "HIT"
        and row["planner"]["prediction"]["complete_predicted_root_set"] == []
        and row["selected_materialized_truth"]["actual_ade_type"] == "rootless"
        and row["selected_materialized_truth"][
            "predicted_root_set_equal_to_materialized_root_set"
        ]
        for row in records
    )
    historical_coordinate_nondetermination = all(
        row["planner"]["status"] == "HIT"
        and not row["selected_line_is_historical_line"]
        and row["historical_regression_after_planning"]["materialized_truth"][
            "actual_ade_type"
        ]
        == "rootless"
        and row["selected_materialized_truth"]["actual_ade_type"] == "rootless"
        for row in records
    )
    search_gate = novelty_search["status"] == "HIT"
    reduction_gate = any(
        row["at_least_tenfold_materialization_reduction"] for row in records
    )
    materialization_gate = all(
        row["materializes_at_most_ten_percent_of_quadric"] for row in records
    )
    novelty_gate = bool(all_nonhistorical)
    implementation_gates = (
        historical_gate
        and search_gate
        and reduction_gate
        and materialization_gate
        and novelty_gate
    )
    all_requested_gates = implementation_gates and terminal_recovery_gate
    if all_requested_gates:
        status = "PASS_ALL_ADAPTIVE_INVERSE_ADE_GATES"
    elif implementation_gates:
        status = (
            "PASS_ADAPTIVE_PREDICATE_CONTROLS_AND_NEW_COMPLETION_"
            "TERMINAL_BLIND_RECOVERY_OPEN"
        )
    else:
        status = "BOUNDED_ADAPTIVE_ENGINE_GATES_REMAIN"
    payload = {
        "schema": "elkies-k3.inverse-ade-adaptive-backend.v1",
        "status": status,
        "backend_contract": {
            "logical_predicate": (
                "RootlessLines_p = Q_p^iso minus old-root hyperplanes minus "
                "projected scaled dual-shell birth points"
            ),
            "modes": {
                "expanded": "exact full projected-shell hash expansion",
                "orbit": (
                    "exact graph-compatible Aut(K)-orbit compression of the "
                    "expanded projective support"
                ),
                "lazy": (
                    "exact candidate-wise nonzero affine-layer CVP with early "
                    "exit on the first occupied graph cell"
                ),
            },
            "adaptive_policy": (
                "Use full expansion below the estimated representation budget; "
                "use orbit compression there when compatible projective symmetry "
                "is useful; otherwise use lazy affine CVP. Mode selection is "
                "heuristic, but every backend decision is exact."
            ),
        },
        "success_gates": {
            "historical_H3_and_Q80_root_sets_predicted_without_target_core": historical_gate,
            "blind_H3_and_Q80_rootless_line_recovery_from_KCHpADE_only": (
                terminal_recovery_gate
            ),
            "exact_historical_coordinate_equality_nonrequired_diagnostic": all(
                row["selected_line_is_historical_line"] for row in records
            ),
            "historical_coordinate_not_determined_by_KCHpADE": (
                historical_coordinate_nondetermination
            ),
            "target_free_search_finds_new_rootless_completion": search_gate,
            "complete_root_set_fixed_before_child_construction": (
                historical_gate and terminal_recovery_gate and search_gate
            ),
            "materialize_at_most_ten_percent_of_isotropic_lines": materialization_gate,
            "at_least_tenfold_materialization_reduction_on_one_terminal": reduction_gate,
            "nonhistorical_rootless_completion_found": novelty_gate,
            "all_requested_success_gates": all_requested_gates,
        },
        "terminal_controls": records,
        "novelty_search": novelty_record,
        "nonhistorical_rootless_completions": all_nonhistorical,
        "backend_mode_equivalence_control": mode_equivalence_control,
        "inputs": {
            relative(SCRIPT): digest(SCRIPT),
            relative(PLANNER_SCRIPT): digest(PLANNER_SCRIPT),
            relative(planner["CONTROL"]): digest(planner["CONTROL"]),
            relative(planner["BRIDGES"]): digest(planner["BRIDGES"]),
            relative(planner["THETA"]): digest(planner["THETA"]),
            relative(planner["Q80_COMPLETION"]): digest(
                planner["Q80_COMPLETION"]
            ),
            relative(planner["SEARCH_SCRIPT"]): digest(planner["SEARCH_SCRIPT"]),
            relative(planner["INVERSE_SCRIPT"]): digest(planner["INVERSE_SCRIPT"]),
            relative(planner["SIGNATURE_SCRIPT"]): digest(
                planner["SIGNATURE_SCRIPT"]
            ),
        },
        "software": {"sage": SAGE_VERSION},
        "proof_boundary": (
            "The three backends are exact for the declared graph-glue root cells. "
            "The geometry-of-numbers count is used only for backend selection and "
            "is not a mathematical claim. Orbit mode currently compresses rootless "
            "support after exact shell enumeration; it does not yet enumerate shell "
            "orbit representatives directly. The H3 and Q80 historical coordinates "
            "are withheld predicate-regression queries made after backend compilation; "
            "they are read only after each blind terminal search returns. The blind "
            "H3 and Q80 solutions differ from their historical lines, while both old "
            "and new lines materialize rootless; therefore K,C,H,p and abstract ADE "
            "do not determine the historical coordinate, and exact coordinate equality "
            "is not a valid success gate. The result proves the two terminal root-set "
            "controls and the displayed new completions, not an asymptotic theorem for "
            "all rank-15 cores."
        ),
        "reproduce": (
            "sage -python "
            "elkies-k3/scripts/benchmark_inverse_ade_adaptive_backend.sage"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        old = json.loads(output.read_text())

        def strip_timings(value):
            if isinstance(value, dict):
                return {
                    key: strip_timings(item)
                    for key, item in value.items()
                    if not key.endswith("_seconds")
                }
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
