#!/usr/bin/env python3
"""Seven-sheet branched-cover certificates, with an independent full edge chain.

Default replay uses only the standard library and retained literal inputs.
--census additionally re-reads the pinned independently packaged named census.
No quotient search, historical numerical topology test, or Smith producer is used.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def diagram(pd):
    n = len(pd)
    ends, adjacency = {}, {e: set() for row in pd for e in row}
    for i, row in enumerate(pd):
        assert len(row) == 4 and len(set(row)) == 4
        for h, e in enumerate(row):
            ends.setdefault(e, []).append(4 * i + h)
        a, b = row[1], row[3]
        adjacency[a].add(b)
        adjacency[b].add(a)
    assert len(ends) == 2 * n and all(len(v) == 2 for v in ends.values())
    pairs = {a: b for pair in ends.values() for a, b in (pair, pair[::-1])}
    current, visits = 0, set()
    for _ in range(2 * n):
        assert current not in visits
        visits.add(current)
        current = pairs[4 * (current // 4) + (current + 2) % 4]
    assert current == 0 and len({v // 4 for v in visits}) == n
    unseen, faces = set(pairs), 0
    while unseen:
        start = current = min(unseen)
        while True:
            unseen.remove(current)
            v = pairs[current]
            current = 4 * (v // 4) + (v + 1) % 4
            if current == start:
                break
        faces += 1
    assert faces == n + 2, "PD is not a sphere diagram"
    unseen, components = set(adjacency), []
    while unseen:
        queue, comp = [min(unseen)], set()
        for e in queue:
            if e in comp:
                continue
            comp.add(e)
            queue.extend(adjacency[e] - comp)
        unseen -= comp
        components.append(sorted(comp))
    components.sort()
    assert len(components) == n
    arc = {e: i for i, comp in enumerate(components) for e in comp}
    crossings = [(arc[b], arc[a], arc[c]) for a, b, c, d in pd]
    assert all(arc[b] == arc[d] for a, b, c, d in pd)
    return crossings, pairs


def census_bijection(pd, census_pd, cert):
    _, a = diagram(pd)
    _, b = diagram(census_pd)
    f, orientation, swap = cert["dart_images"], cert["orientation"], cert["over_under_swap"]
    assert orientation in (-1, 1) and swap in (0, 1)
    assert sorted(f) == list(range(4 * len(pd))) and len(pd) == len(census_pd)
    for x, y in enumerate(f):
        assert f[a[x]] == b[y]
        assert (x % 2) ^ (y % 2) == swap
        assert f[4 * (x // 4) + (x + 1) % 4] == 4 * (y // 4) + (y + orientation) % 4


def matvec(rows, v):
    assert len(rows) == 3 and all(isinstance(r, int) and 0 <= r < 8 for r in rows)
    return sum(((r & v).bit_count() % 2) << i for i, r in enumerate(rows))


def permutations(matrices):
    out = [tuple(matvec(a, v) - 1 for v in range(1, 8)) for a in matrices]
    for p in out:
        assert sorted(p) == list(range(7))
        assert sum(p[i] != i for i in range(7)) == 4
        assert all(p[p[i]] == i for i in range(7))
    return out


def full_chain(perms, crossings):
    """Keep all directed edges, lift all relators, fill every meridian cycle."""
    n, count = 7, len(perms)
    edges = [(g, v) for g in range(count) for v in range(n)]
    parent, tree = list(range(n)), set()
    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v
    # A different, reverse spanning tree from the folded discovery computation.
    for j in reversed(range(len(edges))):
        g, a = edges[j]
        b = perms[g][a]
        if root(a) != root(b):
            parent[root(a)] = root(b)
            tree.add(j)
    assert len(tree) == 6, "monodromy is not transitive"
    rows = []
    for o, a, b in crossings:
        # All images are involutions; square fillings allow this word in pi1.
        word = (o, a, o, b)
        for start in range(n):
            row, v = [0] * len(edges), start
            for g in word:
                row[g * n + v] += 1
                v = perms[g][v]
            assert v == start, "crossing relation fails"
            rows.append(row)
    for g, p in enumerate(perms):
        unused = set(range(n))
        while unused:
            row, start = [0] * len(edges), min(unused)
            v = start
            while True:
                unused.remove(v)
                row[g * n + v] += 1
                v = p[v]
                if v == start:
                    break
            rows.append(row)
    chords = [j for j in range(len(edges)) if j not in tree]
    return [[r[j] for j in chords] for r in rows], len(chords)


def rank_mod3(rows, columns):
    a = [[v % 3 for v in row] for row in rows]
    pivot = 0
    for j in range(columns):
        found = next((i for i in range(pivot, len(a)) if a[i][j]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        factor = a[pivot][j]
        a[pivot] = [(v * factor) % 3 for v in a[pivot]]
        for i in range(pivot + 1, len(a)):
            if a[i][j]:
                factor = a[i][j]
                a[i] = [(v - factor * w) % 3 for v, w in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def check_case(case):
    pd = case["pd"]
    crossings, _ = diagram(pd)
    census_bijection(pd, case["census_pd"], case["census_map_witness"])
    seen = set(case["seed_arcs"])
    assert len(seen) == len(case["seed_arcs"]) == 4
    for o, old, new in case["seed_trace"]:
        assert o in seen and old in seen and new not in seen
        assert (o, old, new) in crossings or (o, new, old) in crossings
        seen.add(new)
    assert seen == set(range(len(pd)))
    dimensions = []
    for action in case["actions"]:
        assert len(action["meridian_matrices_F2"]) == len(pd)
        perms = permutations(action["meridian_matrices_F2"])
        rows, columns = full_chain(perms, crossings)
        dimension = columns - rank_mod3(rows, columns)
        assert dimension == action["expected_H1_F3_dimension"]
        dimensions.append(dimension)
    assert len(dimensions) == 2 and dimensions[0] == dimensions[1]
    return {"name": case["name"], "cover_degree": 7, "meridian_index": 2,
            "H1_F3_dimensions": dimensions, "bridge_upper": 4,
            "meridional_lower": (dimensions[0] + 6 + 1) // 2,
            "sharp_four": dimensions[0] == 1}


def check_external_census(data):
    import importlib.metadata as metadata
    import snappy_15_knots
    from spherogram import Link
    source = data["census_provenance"]
    assert metadata.version("snappy_15_knots") == source["version"]
    assert metadata.version("spherogram") == source["spherogram_version"]
    root = Path(snappy_15_knots.__file__).parent
    for pin in source["source_files"]:
        content = (root / pin["path"]).read_bytes()
        assert len(content) == pin["bytes"]
        assert hashlib.sha256(content).hexdigest() == pin["sha256"]
    for case in data["cases"]:
        link = Link("K" + case["name"].replace("_", ""))
        assert [list(row) for row in link.PD_code()] == case["census_pd"]
        assert json.loads(json.dumps(link.DT_code())) == case["census_DT"]
    return len(data["cases"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", action="store_true")
    args = parser.parse_args()
    path = HERE / "seven_sheet_cover_inputs.json"
    data = json.loads(path.read_text())
    source = data["prior_uncertainty_source"]
    snapshot = HERE / source["snapshot"]
    raw = snapshot.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == source["sha256"]
    assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == source["git_blob"]
    names = json.loads(raw)
    assert len(names) == 293 and all(c["name"] in names for c in data["cases"])
    results = [check_case(case) for case in data["cases"]]
    assert len(results) == 29
    assert len({r["name"] for r in results}) == 29
    assert sum(r["sharp_four"] for r in results) == 24
    # Relevant corruptions must fail; neither a class typo nor a missing filling is accepted.
    corruptions = []
    for kind in ("identity", "matrix_entry", "seed_trace", "census_map"):
        bad = copy.deepcopy(data["cases"][0])
        if kind == "identity":
            bad["actions"][0]["meridian_matrices_F2"][0] = [1, 2, 4]
        elif kind == "matrix_entry":
            bad["actions"][0]["meridian_matrices_F2"][0][0] ^= 1
        elif kind == "seed_trace":
            bad["seed_trace"][0][2] = bad["seed_arcs"][0]
        else:
            bad["census_map_witness"]["dart_images"][0] = bad["census_map_witness"]["dart_images"][1]
        try:
            check_case(bad)
        except (AssertionError, KeyError, ValueError):
            corruptions.append(kind)
        else:
            raise AssertionError("corrupted certificate accepted: " + kind)
    census = check_external_census(data) if args.census else None
    files = [path, Path(__file__), snapshot]
    receipt = {"status": "PASS", "scope": data["scope"],
               "checks": "Full directed edge-cover matrices over F3, every lifted crossing and all fixed/moved meridian-cycle fillings; reverse tree, independent from folded Smith discovery. Exact retained named-census map bijections and four-seed traces.",
               "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
               "results": results, "sharp_four_count": 24,
               "corruptions_rejected": corruptions,
               "external_census_reread_count": census,
               "assurance": "Internal exact computation. No formal verification, external review, exhaustive priority or other-family exclusion is inferred."}
    dest = HERE / ("seven_sheet_census_verification.json" if args.census else "seven_sheet_cover_verification.json")
    dest.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "sharp_four": 24, "zero_H1_controls": 5,
                      "actions": 58, "census_rereads": census,
                      "corruptions_rejected": len(corruptions)}))


if __name__ == "__main__":
    main()
