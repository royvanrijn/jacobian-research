#!/usr/bin/env python3
"""Independent dense/permutation replay of lift dimensions and frozen cases."""

from collections import Counter, deque
from itertools import combinations, permutations, product
import hashlib
import json
from pathlib import Path

from radical_towers import apply_generator, build
from signed_lifts import solve

HERE = Path(__file__).resolve().parent


def dense_rank(matrix):
    if not matrix:
        return 0
    a = [list(row) for row in matrix]
    pivot = 0
    for column in range(len(a[0])):
        hit = next((i for i in range(pivot, len(a)) if a[i][column]), None)
        if hit is None:
            continue
        a[pivot], a[hit] = a[hit], a[pivot]
        for i in range(pivot + 1, len(a)):
            if a[i][column]:
                a[i] = [x ^ y for x, y in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def swap_bits(x, p, q):
    bp, bq = (x >> p) & 1, (x >> q) & 1
    return x ^ ((bp ^ bq) << p) ^ ((bp ^ bq) << q)


def direct_dimension(data):
    r, edges = data["degree"], data["transpositions"]
    s = len(edges)
    matrix = []
    for word in data["relators"]:
        block = [[0] * s for _ in range(r)]
        # Independently evaluate one variable at a time, instead of propagating
        # symbolic linear forms or using the production row reducer.
        for j in range(s):
            x = 0
            for letter in reversed(word):
                i = abs(letter) - 1
                p, q = edges[i]
                x = swap_bits(x, p, q)
                if i == j:
                    x ^= (1 << p) | (1 << q)
            for coordinate in range(r):
                block[coordinate][j] = (x >> coordinate) & 1
        matrix.extend(block)
    return s - dense_rank(matrix)


def pretzel_s4(r):
    edges = [(2, 3), (0, 1), (1, 3)] + [(0, 3) if i % 2 else (1, 3) for i in range(3, r)]
    parameters = [2] + [3] * (r - 1)
    return {"label": f"rank-two-form pretzel P(2,3 repeated {r-1} times)",
            "degree": 4, "transpositions": edges,
            "relators": [[i + 1, i + 1] for i in range(r)]
                + [[i + 1, (i + 1) % r + 1] * p for i, p in enumerate(parameters)]}


def permutation_group_size(tower):
    size = tower["sheets"]
    generators = [tuple(apply_generator(x, u, a) for x in range(size))
                  for u, a in zip(tower["directions_bits"], tower["offsets"])]
    identity = tuple(range(size))
    seen, queue = {identity}, deque([identity])
    while queue:
        g = queue.popleft()
        for s in generators:
            h = tuple(s[i] for i in g)
            if h not in seen:
                seen.add(h)
                queue.append(h)
    return len(seen)


def independent_crossings(code):
    starts = [i for i, x in enumerate(code) if x < 0]
    segments = []
    for i, start in enumerate(starts):
        finish = starts[(i + 1) % len(starts)]
        positions = []
        j = start
        while True:
            positions.append(j)
            if j == finish and len(positions) > 1:
                break
            j = (j + 1) % len(code)
        segments.append([code[j] for j in positions])
    crossings = []
    for c in range(1, len(starts) + 1):
        o = next(i for i, segment in enumerate(segments) if c in segment)
        a = next(i for i, segment in enumerate(segments) if segment[-1] == -c)
        b = next(i for i, segment in enumerate(segments) if segment[0] == -c)
        crossings.append((o, a, b))
    return crossings


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def replay_panel_case(case, saved):
    crossings = independent_crossings(case["gauss_code"])
    n = len(crossings)
    seeds = saved["seed_strands"]
    seen = set(seeds)
    for o, old, new in saved["seed_propagation"]:
        assert o in seen and old in seen and new not in seen
        assert (o, old, new) in crossings or (o, new, old) in crossings
        seen.add(new)
    assert seen == set(range(n)) and len(seeds) == 5
    summaries = []
    for r, published in zip([3, 4], saved["by_degree"]):
        edges = list(combinations(range(r), 2))
        transpositions = []
        for p, q in edges:
            permutation = list(range(r))
            permutation[p], permutation[q] = permutation[q], permutation[p]
            transpositions.append(tuple(permutation))
        counters = Counter()
        for tail in product(range(len(edges)), repeat=4):
            assignment = [0] + list(tail)
            mapping = dict(zip(seeds, (transpositions[i] for i in assignment)))
            for o, old, new in saved["seed_propagation"]:
                mapping[new] = compose(mapping[o], compose(mapping[old], mapping[o]))
            if any(mapping[b] != compose(mapping[o], compose(mapping[a], mapping[o]))
                   for o, a, b in crossings):
                continue
            labels = []
            reached = {0}
            for i in range(n):
                labels.append(tuple(j for j in range(r) if mapping[i][j] != j))
            changed = True
            while changed:
                before = set(reached)
                for p, q in labels:
                    if p in reached or q in reached:
                        reached.update([p, q])
                changed = reached != before
            if len(reached) != r:
                continue
            data = {"degree": r, "transpositions": labels,
                    "relators": [[o + 1, a + 1, o + 1, b + 1] for o, a, b in crossings]}
            counters[direct_dimension(data)] += 1
        assert dict(counters) == {int(k): v for k, v in published["surjective_colouring_bound_counts"].items()}
        best = published["best_certificate"]
        if best:
            assert direct_dimension(best["presentation"]) == best["bound"]
            vectors = [[(a >> i) & 1 for i in range(n)] for a in best["lift_basis_bits"]]
            assert dense_rank(vectors) == best["bound"]
        summaries.append({"degree": r, "bound_counts": dict(counters)})
    tower = build(saved["by_degree"][1]["best_certificate"]["presentation"])
    assert tower["dimension"] == 3 and tower["independent_orbit_size"] == 8
    for o, a, b in crossings:
        for x in range(8):
            y = x
            for j in [b, o, a, o]:
                y = apply_generator(y, tower["directions_bits"][j], tower["offsets"][j])
            assert y == x
    return {"name": case["name"], "independent_replay": "PASS", "degrees": summaries,
            "eight_sheet_affine_certificate": tower}


def main():
    controls, orders = [], []
    for r in range(3, 10):
        data = pretzel_s4(r)
        assert direct_dimension(data) == r
        assert solve(data)["meridional_rank_lower_bound"] == r
        tower = build(data)
        assert tower["form_rank"] == 2 and tower["dimension"] == r - 1
        assert tower["independent_orbit_size"] == 1 << (r - 1)
        for word in data["relators"]:
            for x in range(tower["sheets"]):
                y = x
                for letter in reversed(word):
                    j = abs(letter) - 1
                    y = apply_generator(y, tower["directions_bits"][j], tower["offsets"][j])
                assert y == x
        if r <= 5:
            actual = permutation_group_size(tower)
            assert actual == tower["image_order"]
            orders.append({"bound": r, "sheets": tower["sheets"], "group_order": actual})
        controls.append({"meridional_bound": r, "affine_dimension": r - 1,
                         "form_rank": 2, "image_order": tower["image_order"]})
    # Connected sums use the standard trefoil meridian presentation; the count
    # grows while the base group S3 and the target W(D3)=S4 stay fixed.
    sums = []
    for t in range(1, 9):
        data = {"degree": 3, "transpositions": [(0, 1)] + [(0, 2)] * t,
                "relators": [[1, j + 2, 1, -(j + 2), -1, -(j + 2)] for j in range(t)]}
        assert direct_dimension(data) == t + 1
        assert solve(data)["meridional_rank_lower_bound"] == t + 1
        sums.append({"trefoil_summands": t, "bound": t + 1})
    base = pretzel_s4(5)
    for permutation in permutations(range(4)):
        relabelled = {**base, "transpositions": [(permutation[p], permutation[q]) for p, q in base["transpositions"]]}
        assert direct_dimension(relabelled) == solve(relabelled)["meridional_rank_lower_bound"] == 5
    # Add a redundant conjugate meridian. Intrinsic lift dimension must not grow.
    augmented = {**base, "transpositions": base["transpositions"] + [(1, 2)],
                 "relators": base["relators"] + [[6, 3, -1, -3]]}
    assert direct_dimension(augmented) == solve(augmented)["meridional_rank_lower_bound"] == 5

    panel = json.loads((HERE / "lift_panel.json").read_text())
    saved = json.loads((HERE / "lift_panel_results.json").read_text())
    source_receipt = json.loads((HERE / "panel_source_verification.json").read_text())
    assert source_receipt["status"] == "PASS" and source_receipt["workbook_unchanged"]
    assert source_receipt["panel_sha256"] == saved["panel_sha256"]
    assert source_receipt["reader_sha256"] == hashlib.sha256((HERE / "verify_panel_source.py").read_bytes()).hexdigest()
    assert source_receipt["source_sha256"] == hashlib.sha256((HERE / panel["source"]).read_bytes()).hexdigest()
    assert hashlib.sha256((HERE / "lift_panel.json").read_bytes()).hexdigest() == saved["panel_sha256"]
    assert len(panel["cases"]) == len(saved["results"]) == 16
    replays = []
    for case, result in zip(panel["cases"], saved["results"]):
        assert case["name"] == result["name"] and result["status"] == "COMPLETE"
        replays.append(replay_panel_case(case, result))
    inputs = ["signed_lifts.py", "radical_towers.py", "search_lift_panel.py", "verify_lift_dimension.py",
              "lift_panel.json", "lift_panel_results.json", "LIFT_DIMENSION_AND_RADICAL_TOWERS.tex",
              "verify_panel_source.py", "panel_source_verification.json", "sources/wirt-hm/manifest.json"]
    report = {"status": "PASS", "written_proof_is_not_formal_verification": True,
              "radical_tower_controls": controls, "independent_group_orders": orders,
              "connected_sum_controls": sums, "base_relabellings": 24,
              "redundant_meridian_dimension_unchanged": True,
              "frozen_panel_independent_replay": replays,
              "new_sharp_panel_certificates": sum(r["sharp_counting_certificate"] for r in saved["results"]),
              "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}}
    (HERE / "lift_dimension_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    (HERE / "example_radical_pretzel_5.json").write_text(json.dumps(pretzel_s4(5), indent=2) + "\n")
    print(json.dumps({"status": "PASS", "panel_cases_independently_replayed": len(replays),
                      "new_sharp_panel_certificates": report["new_sharp_panel_certificates"],
                      "group_orders": orders}, indent=2))


if __name__ == "__main__":
    main()
