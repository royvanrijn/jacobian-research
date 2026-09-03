#!/usr/bin/env sage-python
"""Mass-close small genera and analyze exact defect-directed p-neighbor graphs.

The defect in this calibration is the complete signed norm-two set.  Thus a
zero-support state is a rootless lattice, and a directed line is one whose
Kneser neighbor kills every current physical root.  New roots may be born in
the nonzero affine layers of the neighbor dual, exactly as in Theorem H0i.1.

The examples are intentionally ternary: every isotropic line can be listed,
every child can be classified by exact integral isometry, and completeness of
the vertex set is certified by the Minkowski--Siegel mass.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
from itertools import product
import json
from pathlib import Path

from sage.all import QQ, ZZ, Genus, QuadraticForm, matrix, pari, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v1.json"
)

FIXTURES = (
    {
        "name": "ternary_det112",
        "seed": [[2, -1, 0], [-1, 4, 0], [0, 0, 16]],
        "discovery_prime": 3,
        "analysis_primes": (3, 5),
        "expected_classes": 4,
        "expected_mass": QQ(3) / 4,
        "expected_zero_states": 1,
    },
    {
        "name": "ternary_det126",
        "seed": [[2, 0, 0], [0, 4, -1], [0, -1, 16]],
        "discovery_prime": 5,
        "analysis_primes": (5,),
        "expected_classes": 3,
        "expected_mass": QQ(3) / 4,
        "expected_zero_states": 1,
    },
    {
        "name": "ternary_det316",
        "seed": [[2, -1, -1], [-1, 12, -4], [-1, -4, 16]],
        "discovery_prime": 3,
        "analysis_primes": (3, 5),
        "expected_classes": 9,
        "expected_mass": QQ(39) / 16,
        "expected_zero_states": 6,
    },
)


def rational(value):
    value = QQ(value)
    return {
        "numerator": int(value.numerator()),
        "denominator": int(value.denominator()),
        "text": str(value),
    }


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def gram_digest(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quadratic_form(gram):
    coefficients = [
        gram[left, left] // 2 if left == right else gram[left, right]
        for left in range(gram.nrows())
        for right in range(left, gram.ncols())
    ]
    form = QuadraticForm(ZZ, gram.nrows(), coefficients)
    assert form.Hessian_matrix() == gram
    return form


def lll_reduce(gram):
    transform = matrix(ZZ, pari(gram).qflllgram())
    reduced = transform.transpose() * gram * transform
    assert abs(transform.det()) == 1
    assert reduced.det() == gram.det()
    return reduced


def projective_isotropic_lines(gram, prime):
    """List every isotropic line once, normalized at its first nonzero entry."""

    rank = gram.nrows()
    answer = []
    for entries in product(range(prime), repeat=rank):
        if not any(entries):
            continue
        pivot = next(index for index, entry in enumerate(entries) if entry)
        if entries[pivot] != 1:
            continue
        line = vector(ZZ, entries)
        if int((line * gram * line) // 2) % prime == 0:
            answer.append(line)
    # In odd dimension three the nonsingular projective quadric has p+1 points.
    assert rank != 3 or len(answer) == prime + 1
    return tuple(answer)


def isometry_index(gram, representatives):
    for index, representative in enumerate(representatives):
        if pari(gram).qfisom(pari(representative)) != 0:
            return index
    return None


def enumerate_genus(seed, prime):
    """Close one p-neighbor component; the caller certifies it by exact mass."""

    representatives = [lll_reduce(seed)]
    pending = deque([0])
    while pending:
        source = representatives[pending.popleft()]
        form = quadratic_form(source)
        for line in projective_isotropic_lines(source, prime):
            child = form.find_p_neighbor_from_vec(prime, line).Hessian_matrix()
            child = lll_reduce(child)
            if isometry_index(child, representatives) is None:
                representatives.append(child)
                pending.append(len(representatives) - 1)
    return representatives


def physical_roots(gram):
    enumeration = pari(gram).qfminim(2)
    assert int(enumeration[0]) % 2 == 0
    positive = matrix(ZZ, enumeration[2].sage()).columns()
    signed = set()
    for root in positive:
        signed.add(tuple(map(int, root)))
        signed.add(tuple(map(int, -root)))
    return tuple(sorted(signed))


def theta_shells(gram, maximum=20):
    cumulative = 0
    result = []
    for norm in range(2, maximum + 1, 2):
        new_cumulative = int(pari(gram).qfminim(norm)[0])
        result.append(
            {"norm": norm, "signed_vectors": new_cumulative - cumulative}
        )
        cumulative = new_cumulative
    return result


def root_complement(gram, roots):
    if not roots:
        return None
    root_matrix = matrix(ZZ, [list(root) for root in roots])
    equations = root_matrix * gram
    complement_basis = equations.right_kernel_matrix()
    complement = complement_basis * gram * complement_basis.transpose()
    return {
        "rank": complement.nrows(),
        "determinant": int(complement.det()) if complement.nrows() else 1,
        "gram": rows(complement),
    }


def state_record(label, gram):
    roots = physical_roots(gram)
    root_matrix = matrix(ZZ, [list(root) for root in roots]) if roots else None
    root_gram = root_matrix * gram * root_matrix.transpose() if roots else None
    minimum = int(pari(gram).qfminim()[1])
    assert (minimum == 4) == (not roots)
    return {
        "state": label,
        "rank": gram.nrows(),
        "determinant": int(gram.det()),
        "minimum": minimum,
        "zero_support": not roots,
        "automorphism_group_order": int(pari(gram).qfauto()[0]),
        "gram": rows(gram),
        "gram_sha256": gram_digest(gram),
        "physical_defect_signature": {
            "signed_root_count": len(roots),
            "signed_root_vectors": [list(root) for root in roots],
            "signed_root_gram": rows(root_gram) if root_gram is not None else [],
            "root_span_rank": int(root_matrix.rank()) if root_matrix is not None else 0,
            "root_complement": root_complement(gram, roots),
        },
        "theta_shells_through_norm_20": theta_shells(gram),
    }


def tarjan_scc(adjacency):
    index = 0
    stack = []
    on_stack = set()
    indices = {}
    low = {}
    components = []

    def visit(vertex):
        nonlocal index
        indices[vertex] = index
        low[vertex] = index
        index += 1
        stack.append(vertex)
        on_stack.add(vertex)
        for child in sorted(adjacency.get(vertex, ())):
            if child not in indices:
                visit(child)
                low[vertex] = min(low[vertex], low[child])
            elif child in on_stack:
                low[vertex] = min(low[vertex], indices[child])
        if low[vertex] == indices[vertex]:
            component = []
            while True:
                child = stack.pop()
                on_stack.remove(child)
                component.append(child)
                if child == vertex:
                    break
            components.append(tuple(sorted(component)))

    for vertex in sorted(adjacency):
        if vertex not in indices:
            visit(vertex)
    return tuple(sorted(components, key=lambda item: (item[0], len(item), item)))


def shortest_zero_paths(adjacency, zero_states):
    reverse = {vertex: set() for vertex in adjacency}
    for source, destinations in adjacency.items():
        for destination in destinations:
            reverse[destination].add(source)
    distance = {vertex: 0 for vertex in zero_states}
    next_step = {}
    pending = deque(sorted(zero_states))
    while pending:
        destination = pending.popleft()
        for source in sorted(reverse[destination]):
            if source in distance:
                continue
            distance[source] = distance[destination] + 1
            next_step[source] = destination
            pending.append(source)

    paths = {}
    for source in sorted(adjacency):
        if source not in distance:
            paths[source] = None
            continue
        path = [source]
        while path[-1] not in zero_states:
            path.append(next_step[path[-1]])
        paths[source] = path
    return distance, paths


def graph_summary(state_records, edge_records, primes):
    state_names = [row["state"] for row in state_records]
    zero_states = {row["state"] for row in state_records if row["zero_support"]}
    selected = [row for row in edge_records if row["prime"] in primes]
    full_weights = Counter((row["source"], row["destination"]) for row in selected)
    directed = [row for row in selected if row["defect_directed"]]
    directed_weights = Counter((row["source"], row["destination"]) for row in directed)

    full_adjacency = {state: set() for state in state_names}
    directed_adjacency = {state: set() for state in state_names}
    for source, destination in full_weights:
        full_adjacency[source].add(destination)
    for source, destination in directed_weights:
        directed_adjacency[source].add(destination)

    distance, paths = shortest_zero_paths(directed_adjacency, zero_states)
    defective_states = set(state_names) - zero_states
    defective_adjacency = {
        state: directed_adjacency[state] & defective_states
        for state in sorted(defective_states)
    }
    defective_components = tarjan_scc(defective_adjacency)
    closed_traps = []
    for component in defective_components:
        component_set = set(component)
        exits = set().union(*(directed_adjacency[state] for state in component))
        if exits <= component_set:
            closed_traps.append(component)

    state_by_name = {row["state"]: row for row in state_records}
    reachable_defective = defective_states & set(distance)
    unreachable_defective = defective_states - set(distance)
    separator = None
    if reachable_defective and unreachable_defective:
        reachable_orders = {
            state_by_name[state]["automorphism_group_order"]
            for state in reachable_defective
        }
        unreachable_orders = {
            state_by_name[state]["automorphism_group_order"]
            for state in unreachable_defective
        }
        if reachable_orders.isdisjoint(unreachable_orders):
            separator = {
                "field": "automorphism_group_order",
                "reachable_values": sorted(reachable_orders),
                "unreachable_values": sorted(unreachable_orders),
                "boundary": "exact separator in this finite graph, not a universal invariant",
            }

    def weighted_rows(counter):
        return [
            {"source": source, "destination": destination, "line_count": count}
            for (source, destination), count in sorted(counter.items())
        ]

    one_step_zero_lines = {
        state: sum(
            count
            for (source, destination), count in directed_weights.items()
            if source == state and destination in zero_states
        )
        for state in state_names
    }
    destination_defects = {}
    for state in state_names:
        histogram = Counter()
        for row in directed:
            if row["source"] == state:
                histogram[row["child_signed_root_count"]] += 1
        destination_defects[state] = {
            str(defect): count for defect, count in sorted(histogram.items())
        }

    return {
        "primes": list(primes),
        "full_graph": {
            "weighted_adjacency": weighted_rows(full_weights),
            "sccs": [list(component) for component in tarjan_scc(full_adjacency)],
            "strongly_connected": len(tarjan_scc(full_adjacency)) == 1,
        },
        "defect_directed_graph": {
            "definition": (
                "An edge is directed when the defining isotropic line is "
                "nonorthogonal modulo p to every signed physical root of the parent."
            ),
            "weighted_adjacency": weighted_rows(directed_weights),
            "defective_sccs": [list(component) for component in defective_components],
            "closed_defect_traps": [list(component) for component in closed_traps],
            "zero_states": sorted(zero_states),
            "reachable_defective_states": sorted(reachable_defective),
            "unreachable_defective_states": sorted(unreachable_defective),
            "shortest_distance_to_zero": {
                state: distance.get(state) for state in state_names
            },
            "shortest_paths_to_zero": {
                state: paths[state] for state in state_names
            },
            "one_step_zero_line_count": one_step_zero_lines,
            "directed_destination_defect_histogram": destination_defects,
            "observed_separator": separator,
        },
    }


def analyze_fixture(fixture):
    seed = matrix(ZZ, fixture["seed"])
    assert seed.is_positive_definite()
    assert all(seed[index, index] % 2 == 0 for index in range(seed.nrows()))
    assert seed.det() % fixture["discovery_prime"]
    representatives = enumerate_genus(seed, fixture["discovery_prime"])

    # Relabel only after closure, keeping defective states before zero states.
    representatives.sort(
        key=lambda gram: (
            not bool(physical_roots(gram)),
            len(physical_roots(gram)),
            tuple(gram.list()),
        )
    )
    labels = []
    defect_index = 0
    zero_index = 0
    for gram in representatives:
        if physical_roots(gram):
            labels.append(f"D{defect_index}")
            defect_index += 1
        else:
            labels.append(f"Z{zero_index}")
            zero_index += 1

    state_records = [
        state_record(label, gram) for label, gram in zip(labels, representatives)
    ]
    automorphism_orders = [
        row["automorphism_group_order"] for row in state_records
    ]
    observed_mass = sum(QQ(1) / order for order in automorphism_orders)
    target_genus = Genus(seed)
    target_mass = target_genus.mass()
    assert len(representatives) == fixture["expected_classes"]
    assert target_mass == fixture["expected_mass"]
    assert observed_mass == target_mass
    assert zero_index == fixture["expected_zero_states"]

    edge_records = []
    for prime in fixture["analysis_primes"]:
        assert seed.det() % prime
        for source_index, (source_label, source) in enumerate(
            zip(labels, representatives)
        ):
            source_roots = physical_roots(source)
            form = quadratic_form(source)
            for line in projective_isotropic_lines(source, prime):
                child = form.find_p_neighbor_from_vec(prime, line).Hessian_matrix()
                assert child.det() == source.det()
                assert Genus(child) == target_genus
                reduced_child = lll_reduce(child)
                destination_index = isometry_index(reduced_child, representatives)
                assert destination_index is not None
                survivor_roots = [
                    root
                    for root in source_roots
                    if int(vector(ZZ, root) * source * line) % prime == 0
                ]
                child_roots = physical_roots(child)
                assert len(child_roots) >= len(survivor_roots)
                edge_records.append(
                    {
                        "prime": prime,
                        "source": source_label,
                        "destination": labels[destination_index],
                        "isotropic_line": list(map(int, line)),
                        "parent_signed_root_count": len(source_roots),
                        "old_survivor_signed_root_count": len(survivor_roots),
                        "old_death_signed_root_count": (
                            len(source_roots) - len(survivor_roots)
                        ),
                        "new_birth_signed_root_count": (
                            len(child_roots) - len(survivor_roots)
                        ),
                        "child_signed_root_count": len(child_roots),
                        "defect_directed": not survivor_roots,
                    }
                )

    # Exact self-adjointness check for the unquotiented line multiplicities.
    state_by_name = {row["state"]: row for row in state_records}
    for prime in fixture["analysis_primes"]:
        weights = Counter(
            (row["source"], row["destination"])
            for row in edge_records
            if row["prime"] == prime
        )
        for (source, destination), count in weights.items():
            reverse_count = weights[destination, source]
            source_order = state_by_name[source]["automorphism_group_order"]
            destination_order = state_by_name[destination]["automorphism_group_order"]
            assert QQ(count) / source_order == QQ(reverse_count) / destination_order

    prime_sets = [(prime,) for prime in fixture["analysis_primes"]]
    if len(fixture["analysis_primes"]) > 1:
        prime_sets.append(tuple(fixture["analysis_primes"]))
    graph_records = [
        graph_summary(state_records, edge_records, primes) for primes in prime_sets
    ]
    return {
        "name": fixture["name"],
        "genus": {
            "rank": seed.nrows(),
            "determinant": int(seed.det()),
            "even": True,
            "discovery_prime": fixture["discovery_prime"],
            "class_count": len(representatives),
            "zero_support_class_count": zero_index,
            "target_mass": rational(target_mass),
            "enumerated_mass": rational(observed_mass),
            "mass_closed": True,
            "local_symbols": [str(symbol) for symbol in target_genus.local_symbols()],
        },
        "states": state_records,
        "line_edges": edge_records,
        "graphs": graph_records,
    }


def graph_for(record, primes):
    return next(row for row in record["graphs"] if row["primes"] == list(primes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    genera = [analyze_fixture(fixture) for fixture in FIXTURES]
    by_name = {row["name"]: row for row in genera}

    det112_p3 = graph_for(by_name["ternary_det112"], (3,))[
        "defect_directed_graph"
    ]
    assert len(det112_p3["unreachable_defective_states"]) == 2
    assert sorted(map(len, det112_p3["closed_defect_traps"])) == [2]
    assert graph_for(by_name["ternary_det112"], (3,))["full_graph"][
        "strongly_connected"
    ]
    det112_p5 = graph_for(by_name["ternary_det112"], (5,))[
        "defect_directed_graph"
    ]
    assert not det112_p5["unreachable_defective_states"]
    assert max(
        distance
        for state, distance in det112_p5["shortest_distance_to_zero"].items()
        if state.startswith("D")
    ) == 1

    det126_p5 = graph_for(by_name["ternary_det126"], (5,))[
        "defect_directed_graph"
    ]
    det126_distances = sorted(
        distance
        for state, distance in det126_p5["shortest_distance_to_zero"].items()
        if state.startswith("D")
    )
    assert det126_distances == [1, 2]
    distance_two_state = next(
        state
        for state, distance in det126_p5["shortest_distance_to_zero"].items()
        if distance == 2
    )
    distance_two_path = det126_p5["shortest_paths_to_zero"][distance_two_state]
    det126_states = {
        row["state"]: row for row in by_name["ternary_det126"]["states"]
    }
    assert [
        det126_states[state]["physical_defect_signature"]["signed_root_count"]
        for state in distance_two_path
    ] == [2, 2, 0]

    det316_p3 = graph_for(by_name["ternary_det316"], (3,))[
        "defect_directed_graph"
    ]
    assert len(det316_p3["unreachable_defective_states"]) == 2
    assert sorted(map(len, det316_p3["closed_defect_traps"])) == [1, 1]
    assert graph_for(by_name["ternary_det316"], (3,))["full_graph"][
        "strongly_connected"
    ]
    det316_p5 = graph_for(by_name["ternary_det316"], (5,))[
        "defect_directed_graph"
    ]
    assert not det316_p5["unreachable_defective_states"]

    payload = {
        "schema": "elkies-k3.small-genus-defect-graphs.v1",
        "status": "PASS_EXACT_MASS_COMPLETE_SMALL_GENUS_DEFECT_GRAPHS",
        "definition": {
            "state": "integral isometry class in one positive even ternary genus",
            "physical_defect": "the complete signed set of norm-two vectors",
            "zero_support": "the lattice is rootless",
            "directed_edge": (
                "a good-p Kneser line nonorthogonal modulo p to every current "
                "physical root; all old defects die, while child roots are births"
            ),
        },
        "inputs": {
            str(SCRIPT.relative_to(ROOT)): file_digest(SCRIPT),
            "fixture_grams": "embedded exactly in the checker",
        },
        "genera": genera,
        "findings": {
            "closed_fixed_prime_traps": (
                "The mass-complete determinant-112 directed 3-graph has a closed "
                "two-state defective SCC.  The determinant-316 directed 3-graph "
                "has two closed singleton traps, one a self-trap and one with no "
                "directed outgoing line."
            ),
            "finite_prime_escape": (
                "Prime 5 alone sends every defective state to zero in at most one "
                "step in determinants 112 and 316.  Thus the observed 3-traps are "
                "prime-set traps, not all-good-prime traps."
            ),
            "nontrivial_distance": (
                "The determinant-126 directed 5-graph has defective distances one "
                "and two.  The distance-two path has signed defect counts 2->2->0, "
                "so defect cardinality does not determine distance."
            ),
            "stronger_observed_state_data": (
                "Automorphism-group order separates the reachable and unreachable "
                "defective regions in both directed 3-controls.  Full theta shells, "
                "root complements, and the per-prime directed destination profile "
                "distinguish equal-defect states, but no universal monotone is claimed."
            ),
        },
        "proof_boundary": {
            "proved": (
                "Every genus list closes the exact Minkowski--Siegel mass; every "
                "projective isotropic line at the declared good primes is enumerated; "
                "every child is classified by exact integral isometry; SCCs and "
                "shortest paths are computed from the resulting complete graphs."
            ),
            "not_proved": (
                "These root-defect ternary controls do not establish an all-good-prime "
                "trap, a rank-15 completion-mask theorem, a universal finite prime "
                "bound, or a scalar Lyapunov function."
            ),
        },
        "software_assumptions": {
            "sage": SAGE_VERSION,
            "pari": ".".join(map(str, pari.version())),
            "exact_arithmetic": True,
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage --check"
        ),
    }

    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS exact mass-complete small-genus defect graphs")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
