#!/usr/bin/env python3
"""Independent full orbit marking and geometric/pointwise relator replay.

Does not import the search permutation model, tuple augmentation, group
enumerator or propagation. The old PD parser supplies only legacy arc names;
geometric under-to-under paths independently check that transport.
"""

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from itertools import permutations
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram

HERE = Path(__file__).resolve().parent


def involution_type(p):
    return sorted(p) == list(range(len(p))) and all(p[p[i]] == i for i in range(len(p))) and sum(p[i] != i for i in range(len(p))) == 4


def connected(generators):
    n = len(generators[0])
    # Connected components of the literal transposition-edge graph.
    components = [{i} for i in range(n)]
    for p in generators:
        for i, j in enumerate(p):
            if i >= j:
                continue
            left = next(c for c in components if i in c)
            right = next(c for c in components if j in c)
            if left is not right:
                left.update(right)
                components.remove(right)
    return len(components) == 1


def orbit_data(row, actions, elements, size):
    a, b, c, d = row["tuple"]
    assert a == 0 and all(0 <= i < size for i in [b, c, d])
    orbit = {(g[b] * size + g[c]) * size + g[d] for g in actions}
    assert len(orbit) * row["stabilizer"] == len(actions)
    assert len(orbit) == row["raw_weight"]
    assert row["transitive"] == connected([elements[i] for i in [a, b, c, d]])
    return orbit


def verify_catalogue(catalogue, degree):
    elements = [tuple(p) for p in catalogue["elements"]]
    # Independent enumeration by filtering all permutations, rather than
    # choosing four-point supports and one of three pairings.
    all_perms = list(permutations(range(degree)))
    expected = {p for p in all_perms if involution_type(p)}
    assert len(elements) == len(set(elements)) == len(expected)
    assert set(elements) == expected
    index = {p: i for i, p in enumerate(elements)}
    size = len(elements)
    tau = elements[0]
    assert tau == (1, 0, 3, 2) + tuple(range(4, degree))
    centralizer = [p for p in all_perms if all(p[tau[i]] == tau[p[i]] for i in range(degree))]
    assert [list(p) for p in centralizer] == catalogue["centralizer"]
    actions = []
    for p in centralizer:
        action = []
        for q in elements:
            conjugate = [None] * degree
            for i in range(degree):
                conjugate[p[i]] = p[q[i]]
            action.append(index[tuple(conjugate)])
        assert action[0] == 0
        actions.append(action)
    # Mark every actual orbit on the entire normalized raw tuple cube.
    # Disjointness plus no unmarked tuple proves coverage without relying on
    # the producer's successive-stabilizer recursion or just a count sum.
    marked = bytearray(size ** 3)
    transitive = []
    raw, connected_raw = 0, 0
    for row in catalogue["representatives"]:
        orbit = orbit_data(row, actions, elements, size)
        for encoded in orbit:
            assert marked[encoded] == 0, "overlapping tuple orbits"
            marked[encoded] = 1
        raw += len(orbit)
        if row["transitive"]:
            transitive.append(row)
            connected_raw += len(orbit)
    assert all(marked) and raw == size ** 3
    counts = catalogue["counts"]
    assert counts["class_size"] == size and counts["centralizer_size"] == len(actions)
    assert counts["four_tuple_orbits"] == len(catalogue["representatives"])
    assert counts["raw_normalized_tuples"] == raw
    assert counts["transitive_orbits"] == len(transitive)
    assert counts["transitive_raw_normalized_tuples"] == connected_raw
    assert {int(k): v for k, v in counts["stabilizer_histogram"].items()} == dict(Counter(r["stabilizer"] for r in catalogue["representatives"]))
    # Reverse pointwise evaluation of conjugation; all meridians involutive.
    table = [[index[tuple(a[b[a[i]]] for i in range(degree))] for b in elements] for a in elements]
    return elements, actions, table, transitive


def geometric_crossings(pd):
    ends = defaultdict(list)
    for c, row in enumerate(pd):
        for h, edge in enumerate(row):
            ends[edge].append((c, h))
    assert all(len(pair) == 2 for pair in ends.values())
    visited, groups = set(), []
    for c in range(len(pd)):
        for h in [0, 2]:
            if (c, h) in visited:
                continue
            current, edges = (c, h), []
            visited.add(current)
            for _ in range(2 * len(pd) + 1):
                edge = pd[current[0]][current[1]]
                edges.append(edge)
                pair = ends[edge]
                d, k = pair[1] if pair[0] == current else pair[0]
                if k % 2 == 0:
                    visited.add((d, k))
                    break
                current = d, (k + 2) % 4
            else:
                raise AssertionError("nonterminating geometric arc")
            groups.append(sorted(edges))
    assert len(groups) == len(pd) and len(visited) == 2 * len(pd)
    physical = {e: i for i, g in enumerate(groups) for e in g}
    legacy = oriented_diagram(pd)
    names = {}
    for row, (o, a, b, _) in zip(pd, legacy):
        for edge, name in [(row[1], o), (row[0], a), (row[2], b)]:
            p = physical[edge]
            if p in names:
                assert names[p] == name
            names[p] = name
    assert set(names) == set(range(len(pd))) == set(names.values())
    return [(names[physical[row[1]]], names[physical[row[0]]], names[physical[row[2]]]) for row in pd]


def check_trace(crossings, seeds, trace):
    assert len(seeds) == len(set(seeds)) == 4
    seen = set(seeds)
    valid = {(o, a, b) for o, a, b in crossings} | {(o, b, a) for o, a, b in crossings}
    for o, source, target in trace:
        assert (o, source, target) in valid
        assert o in seen and source in seen and target not in seen
        seen.add(target)
    assert seen == set(range(len(crossings)))


def recheck_targets(report, panel, elements, table, transitive):
    assert report["status"] == "COMPLETE" and report["sharp_names"] == []
    assert [r["name"] for r in report["results"]] == [c["name"] for c in panel]
    checks = []
    for result, case in zip(report["results"], panel):
        crossings = geometric_crossings(case["pd"])
        seeds, trace = result["seeds"], result["seed_trace"]
        check_trace(crossings, seeds, trace)
        successes = 0
        for row in transitive:
            labels = [-1] * len(crossings)
            for seed, value in zip(seeds, row["tuple"]):
                labels[seed] = value
            # Exactly one pass through the geometric upper trace, then test
            # all relators. Producer propagates iteratively with early failure.
            for o, source, target in trace:
                labels[target] = table[labels[o]][labels[source]]
            assert -1 not in labels
            if all(table[labels[o]][labels[a]] == labels[b] for o, a, b in crossings):
                successes += 1
        assert successes == len(result["certificates"]) == 0
        assert result["status"] == "COMPLETE" and result["transitive_seed_orbits_tested"] == len(transitive)
        checks.append({"name": case["name"], "orbits_replayed": len(transitive), "hits": successes})
    return checks


def group_order(generators):
    n = len(generators[0])
    identity = tuple(range(n))
    reached, frontier = {identity}, [identity]
    while frontier:
        following = []
        for p in frontier:
            for q in generators:
                # Left multiplication, separate from producer right BFS.
                value = tuple(q[p[i]] for i in range(n))
                if value not in reached:
                    assert len(reached) < 20160
                    reached.add(value)
                    following.append(value)
        frontier = following
    return len(reached)


def check_control(item, q, exponent):
    x, y, m = [tuple(item[k]) for k in ["x", "y", "m"]]
    n = len(x)
    assert involution_type(m) and connected([x, y])
    def apply_power(p, point, count):
        for _ in range(count):
            point = p[point]
        return point
    for i in range(n):
        assert apply_power(x, i, 4) == apply_power(y, i, q) == i
        # x*m = y^exponent is the meridian equation without inversion code.
        assert x[m[i]] == apply_power(y, i, exponent)
    assert group_order([x, y]) == item["image_order"]


def main():
    started = time.monotonic()
    panel = json.loads((HERE / "parabolic_panel.json").read_text())["cases"]
    results, rejects = [], []
    for degree, prefix in [(8, "double_transposition"), (9, "double_transposition9")]:
        catalogue = json.loads((HERE / (prefix + "_catalogue.json")).read_text())
        report = json.loads((HERE / (prefix + "_results.json")).read_text())
        assert hashlib.sha256((HERE / (prefix + "_catalogue.json")).read_bytes()).hexdigest() == report["catalogue_sha256"]
        for p, digest in report["input_sha256"].items():
            assert hashlib.sha256((HERE / p).read_bytes()).hexdigest() == digest
        elements, actions, table, transitive = verify_catalogue(catalogue, degree)
        rows = recheck_targets(report, panel, elements, table, transitive)
        results.append({"degree": degree, "catalogue": catalogue["counts"], "targets": rows})
        print(json.dumps({"degree": degree, "status": "PASS full raw tuple marking and independent relators", "transitive_orbits": len(transitive)}), flush=True)
        for kind in ["false_stabilizer", "lost_transitivity"]:
            bad = deepcopy(transitive[0])
            if kind == "false_stabilizer":
                bad["stabilizer"] *= 2
            else:
                bad["transitive"] = False
            try:
                orbit_data(bad, actions, elements, len(elements))
            except AssertionError:
                rejects.append(f"degree{degree}:{kind}")
            else:
                raise AssertionError("accepted altered orbit input")
    controls = json.loads((HERE / "double_transposition_results.json").read_text())["controls"]
    checked_controls = []
    for control, q, exponent in zip(controls, [7, 15], [2, 4]):
        for item in control["transitive_pairs"]:
            check_control(item, q, exponent)
        checked_controls.append({"name": control["name"], "pairs": len(control["transitive_pairs"]), "image_orders": control["image_orders"]})
    bad = deepcopy(controls[1]["transitive_pairs"][0])
    bad["m"] = list(range(8))
    try:
        check_control(bad, 15, 4)
    except AssertionError:
        rejects.append("false_torus_meridian")
    else:
        raise AssertionError("accepted false torus meridian")
    bad_trace = deepcopy(json.loads((HERE / "double_transposition_results.json").read_text())["results"][0])
    bad_trace["seed_trace"][0][2] = bad_trace["seeds"][0]
    try:
        check_trace(geometric_crossings(panel[0]["pd"]), bad_trace["seeds"], bad_trace["seed_trace"])
    except AssertionError:
        rejects.append("false_upper_trace")
    else:
        raise AssertionError("accepted false upper trace")
    # Elementary class counts: all nonzero v and nonzero annihilating f
    # produce distinct transvections I+v*f over F2.
    transvections = set()
    for v in range(1, 16):
        for f in range(1, 16):
            if (v & f).bit_count() % 2 == 0:
                transvections.add(tuple(j ^ (v if (f & j).bit_count() % 2 else 0) for j in range(16)))
    assert len(transvections) == 105
    files = ["double_transposition_plan.json", "double_transposition_degree9_plan.json",
        "double_transposition_results.json", "double_transposition9_results.json",
        "double_transposition_catalogue.json", "double_transposition9_catalogue.json",
        "verify_double_transpositions.py", "parabolic_lifts.py", "parabolic_panel.json"]
    receipt = {"status": "PASS independent degree8/9 double-transposition gate", "degrees": results,
        "controls": checked_controls, "rejected_corruptions": rejects,
        "class_boundary": {"A8_double_transpositions": 210, "GL4_F2_transvections": len(transvections), "A8_order": 20160, "AGL3_F2_order": 1344},
        "seconds": round(time.monotonic() - started, 6),
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in files},
        "scope": "Full orbit coverage and independent pointwise/geometric negative checks on exactly16 targets. Torus topology/meridian convention and uniform graph bound are written arguments; parser reused only for legacy arc index transport. Full MRC and all16 targets remain open."}
    (HERE / "double_transposition_verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "controls": checked_controls, "rejected_corruptions": rejects, "seconds": receipt["seconds"]}, indent=2))


if __name__ == "__main__":
    main()
