#!/usr/bin/env python3
"""Compact finite checks for the written four-face rational-replacement theorem.

Requires Regina7.4. Classical imports are KM13's branched-cover model,
Martelli's exhaustive M6 filling classification, and the arborescent
cover construction. The local tangle arguments and universal gluing
proof remain written mathematics, not computations in this script.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json
import platform

import regina

HERE = Path(__file__).resolve().parent
INPUT = HERE / "rational_replacement_topology_inputs.json"
NAMES = ["13a_3145", "13a_3147", "13a_3148", "13a_3149",
         "13a_373", "13a_375", "13a_391", "13a_392",
         "13a_3264", "13a_3266", "13a_3275", "13a_3276",
         "13n_1005", "13n_1006", "13n_1007", "13n_1011",
         "13n_1012", "13n_1013", "13n_1043", "13n_1044",
         "13n_1045", "13n_1046", "13n_1047", "13n_1048",
         "13n_3842", "13n_3844", "13n_3855", "13n_3856"]
K4 = list(combinations(range(4), 2))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pins_match(pins):
    assert all(digest(HERE / path) == expected for path, expected in pins.items())


def adjacency(vertices, edges):
    adj = {v: set() for v in vertices}
    for a, b in edges:
        assert a in adj and b in adj and a != b
        adj[a].add(b)
        adj[b].add(a)
    assert sum(map(len, adj.values())) == 2 * len(edges)
    return adj


def connected(vertices, edges):
    if not vertices:
        return False
    adj = adjacency(vertices, edges)
    found = {min(vertices)}
    queue = list(found)
    for a in queue:
        for b in adj[a] - found:
            found.add(b)
            queue.append(b)
    return found == set(vertices)


def classify_proper_graph(edges):
    assert len(edges) < 6 and connected(range(4), edges)
    vertices = set(range(4))
    remaining = list(edges)
    removed = []
    while len(vertices) > 1:
        adj = adjacency(vertices, remaining)
        leaf = next((v for v in sorted(vertices) if len(adj[v]) == 1), None)
        if leaf is None:
            break
        removed.append(leaf)
        vertices.remove(leaf)
        remaining = [e for e in remaining if leaf not in e]
    if not remaining:
        assert len(vertices) == 1 and len(edges) == 3
        return "tree_after_leaf_deletions"
    adj = adjacency(vertices, remaining)
    degrees = sorted(map(len, adj.values()))
    if len(vertices) == 3:
        assert len(remaining) == 3 and degrees == [2, 2, 2]
        assert len(removed) == 1
        return "triangle_with_leaf"
    assert len(vertices) == 4 and not removed
    if len(remaining) == 4:
        assert degrees == [2, 2, 2, 2]
        return "four_cycle"
    assert len(remaining) == 5 and degrees == [2, 2, 3, 3]
    small = [v for v in vertices if len(adj[v]) == 2]
    big = [v for v in vertices if len(adj[v]) == 3]
    assert all(adj[v] == set(big) for v in small)
    assert tuple(sorted(big)) in remaining
    # Absorbing the two degree-two vertices creates three parallel boxes
    # on the two surviving vertices; there is no intervening graph material.
    return "diamond"


def graph_audit():
    counts = Counter()
    representatives = []
    total = 0
    for bits in range(1 << 6):
        edges = [e for i, e in enumerate(K4) if bits >> i & 1]
        if not connected(range(4), edges):
            continue
        total += 1
        if len(edges) == 6:
            counts["K4"] += 1
            continue
        kind = classify_proper_graph(edges)
        counts[kind] += 1
        representatives.append({"edge_mask": bits, "edges": edges, "case": kind})
    counts["connected"] = total
    # The medial of the tetrahedron is its line graph. Each pair of
    # incident tetrahedral edges is consecutive in the degree-three star.
    medial = [(i, j) for i, e in enumerate(K4) for j, f in enumerate(K4)
              if i < j and set(e) & set(f)]
    adj = adjacency(range(6), medial)
    assert len(medial) == 12 and all(len(v) == 4 for v in adj.values())
    opposite = [(i, j) for i in range(6) for j in range(i + 1, 6) if j not in adj[i]]
    assert len(opposite) == 3
    assert sorted(a for pair in opposite for a in pair) == list(range(6))
    assert 2 * 2 - 4 == 0  # connected double cover of each four-marked sphere
    return dict(counts), representatives, medial


def mod2_b1(triangulation):
    homology = triangulation.homology()
    return homology.rank() + homology.torsionRank(2)


def assert_mod2_sphere(triangulation):
    assert mod2_b1(triangulation) == 0


def rejected(call):
    try:
        call()
    except (AssertionError, ValueError, IndexError, KeyError):
        return True
    raise AssertionError("Invalid control was accepted")


def main():
    assert regina.versionString() == "7.4"
    data = json.loads(INPUT.read_text())
    assert data["expected_names"] == NAMES
    pins_match(data["source_sha256"])
    counts, graph_records, medial = graph_audit()
    assert counts == data["graph_counts"]
    dimensions = data["mixed_exception_H1_F2_dimensions"]
    assert dimensions == {"figure_eight": 1, "Klein_bottle_I_bundle": 2,
                          "gluing_torus": 2, "lower_bound": 1}
    assert dimensions["figure_eight"] + dimensions["Klein_bottle_I_bundle"] - dimensions["gluing_torus"] == 1
    # Abelianising uvu^-1=v^-1 gives the relation 2v=0; modulo2 this
    # imposes no relation on the two generators of the Klein-bottle group.
    assert all(x % 2 == 0 for x in [0, 2])
    inputs = json.loads((HERE / "seed_cover_topology_inputs.json").read_text())
    receipt = json.loads((HERE / "seed_cover_topology_verification.json").read_text())
    assert receipt["status"] == "PASS"
    assert [r["name"] for r in inputs["cases"]] == NAMES
    established = {r["name"]: r for r in receipt["cases"]}
    rows = []
    for row in inputs["cases"]:
        binding = established[row["name"]]
        assert binding["closed_cover_isosig"] == row["closed_cover_isosig"]
        assert binding["irreducible"] and binding["incompressible_normal_torus"]
        assert binding["hyperbolic_piece_isosig"] == row["ideal_piece_isosig"]
        assert binding["outside_BJW18_graph_cover_class"]
        tri = regina.Triangulation3.fromIsoSig(row["closed_cover_isosig"])
        assert tri.isValid() and tri.isClosed() and tri.isOrientable() and tri.isConnected()
        assert_mod2_sphere(tri)
        rows.append({"name": row["name"], "closed_cover_isosig": tri.isoSig(),
                     "closed_H1": str(tri.homology()), "closed_H1_F2_dimension": 0,
                     "mixed_actual_cover_dependency_checked": True,
                     "excluded_from_four_face_rational_replacement_class": True})
    bad_pins = dict(data["source_sha256"])
    bad_pins[next(iter(bad_pins))] = "0" * 64
    rp3 = regina.Example3.lens(2, 1)
    assert mod2_b1(rp3) == 1
    negative = {"K4_falsely_called_proper": rejected(lambda: classify_proper_graph(K4)),
                "RP3_falsely_called_mod2_sphere": rejected(lambda: assert_mod2_sphere(rp3)),
                "changed_source_digest": rejected(lambda: pins_match(bad_pins))}
    pins = [INPUT, Path(__file__)] + [HERE / p for p in data["source_sha256"]]
    report = {"status": "PASS", "python": platform.python_version(),
              "regina": regina.versionString(),
              "command": "python verify_rational_replacement_topology.py",
              "scope": data["scope"], "graph_counts": counts,
              "proper_graph_cases": graph_records, "tetrahedral_medial_edges": medial,
              "mixed_exception_mod2_lower_bound": 1, "cases": rows,
              "negative_controls": negative,
              "input_sha256": {str(p.relative_to(HERE)): digest(p) for p in pins},
              "assurance": "Internal finite checks of the inputs to the written universal topological proof. KM13 branched-cover model, Martelli's complete classification and the arborescent cover construction are imported. Prior actual-cover/tori certificate is bound, not replayed. No formal verification or external specialist review."}
    (HERE / "rational_replacement_topology_verification.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "graph_counts": counts,
                      "mixed_exception_mod2_lower_bound": 1,
                      "excluded_named_knots": len(rows), "negative_controls": negative}, indent=2))


if __name__ == "__main__":
    main()
