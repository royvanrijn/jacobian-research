#!/usr/bin/env python3
"""Exact finite Jones/crosscap inputs for the written BBK21 exclusion.

Standard-library replay. Imported mathematical inputs are KL16 Theorem1.1,
its stated Menasco diagram criterion, and BBK21 Theorem1; this checker does
not reprove them or replay the earlier bridge-number computations.
"""
from collections import defaultdict
from copy import deepcopy
from math import comb
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
INPUT = HERE / "twisted_link_crosscap_inputs.json"
EXPECTED = ["13a_373", "13a_375", "13a_391", "13a_392",
            "13a_3264", "13a_3266", "13a_3275", "13a_3276",
            "13a_3145", "13a_3147", "13a_3148", "13a_3149",
            "13a_304", "13a_362", "13a_623", "13a_629",
            "13a_3232", "13a_3259", "13a_3374", "13a_3375"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins(pins):
    assert all(digest(HERE / name) == expected for name, expected in pins.items())


def endpoints(pd):
    ends = defaultdict(list)
    for c, row in enumerate(pd):
        assert len(row) == 4
        for h, edge in enumerate(row):
            ends[edge].append((c, h))
    assert len(ends) == 2 * len(pd) and all(len(pair) == 2 for pair in ends.values())
    return ends


def alternating_knot(pd):
    ends = endpoints(pd)
    current = start = (0, 0)
    seen = []
    for _ in range(2 * len(pd)):
        assert current not in seen
        seen.append(current)
        c, h = current
        outgoing = (c, (h + 2) % 4)
        pair = ends[pd[c][outgoing[1]]]
        current = pair[1] if pair[0] == outgoing else pair[0]
    assert current == start
    assert len({c for c, h in seen}) == len(pd)
    assert all(seen[i][1] % 2 != seen[(i + 1) % len(seen)][1] % 2
               for i in range(len(seen)))


def projection_edge_connectivity(pd):
    edges = [[c for c, h in pair] for pair in endpoints(pd).values()]
    n = len(pd)
    # Fix the last vertex outside the cut, removing complementary duplication.
    return min(sum(((bits >> u) ^ (bits >> v)) & 1 for u, v in edges)
               for bits in range(1, 1 << (n - 1)))


def assert_prime_diagram(pd):
    assert projection_edge_connectivity(pd) >= 4


def bracket(pd):
    n = len(pd)
    labels = sorted(endpoints(pd))
    poly = defaultdict(int)
    for bits in range(1 << n):
        parent = {e: e for e in labels}

        def find(e):
            while e != parent[e]:
                e = parent[e]
            return e

        def union(a, b):
            parent[find(a)] = find(b)

        for i, (a, b, c, d) in enumerate(pd):
            if bits >> i & 1:
                union(b, c)
                union(d, a)
            else:
                union(a, b)
                union(c, d)
        circles = len({find(e) for e in labels})
        base = n - 2 * bits.bit_count()
        for k in range(circles):
            exponent = base + 2 * (circles - 1) - 4 * k
            poly[exponent] += (-1) ** (circles - 1) * comb(circles - 1, k)
    return {e: c for e, c in sorted(poly.items()) if c}


def checkerboard_coefficients(pd):
    """Rebuild both checkerboard region graphs by corner identifications."""
    n = len(pd)
    parent = {(c, h): (c, h) for c in range(n) for h in range(4)}

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a

    def union(a, b):
        parent[find(a)] = find(b)

    for (c, h), (d, k) in endpoints(pd).values():
        union((c, (h - 1) % 4), (d, k))
        union((c, h), (d, (k - 1) % 4))
    faces = {find(a) for a in parent}
    assert len(faces) == n + 2
    adj = {f: set() for f in faces}
    for c in range(n):
        for h in range(4):
            a, b = find((c, h)), find((c, (h + 1) % 4))
            adj[a].add(b)
            adj[b].add(a)
    colour = {min(faces): 0}
    queue = list(colour)
    for a in queue:
        for b in adj[a]:
            if b not in colour:
                colour[b] = 1 - colour[a]
                queue.append(b)
            assert colour[b] != colour[a]
    assert len(colour) == len(faces)
    result = []
    for shade in range(2):
        vertices = {a for a in faces if colour[a] == shade}
        edges = set()
        for c in range(n):
            h = 0 if colour[find((c, 0))] == shade else 1
            a, b = find((c, h)), find((c, (h + 2) % 4))
            assert a != b  # no nugatory crossing in this colour
            edges.add(tuple(sorted([a, b])))
        result.append({"vertices": len(vertices), "reduced_edges": len(edges),
                       "second_coefficient_abs": len(edges) - len(vertices) + 1})
    return result


def check_case(row, named_pds):
    pd = row["pd"]
    assert len(pd) == 13 and pd == named_pds[row["name"]]
    alternating_knot(pd)
    assert_prime_diagram(pd)
    assert projection_edge_connectivity(pd) == row["projection_edge_connectivity"]
    poly = bracket(pd)
    assert [[e, c] for e, c in poly.items()] == row["bracket_coefficients_A"]
    lo, hi = min(poly), max(poly)
    assert abs(poly[lo]) == abs(poly[hi]) == 1
    beta, beta_prime = abs(poly.get(hi - 4, 0)), abs(poly.get(lo + 4, 0))
    assert (beta, beta_prime) == (row["beta_abs"], row["beta_prime_abs"])
    total = beta + beta_prime
    assert total == row["T"] and (hi - lo) // 4 == row["Jones_span"] == 13
    colours = checkerboard_coefficients(pd)
    assert sorted(r["second_coefficient_abs"] for r in colours) == sorted([beta, beta_prime])
    assert sorted((r["vertices"], r["reduced_edges"]) for r in colours) == sorted(
        (r["vertices"], r["reduced_edges"]) for r in row["checkerboard_replay"])
    assert sorted(r["vertices"] for r in colours) != [2, len(pd)]
    assert total > 1  # unlike the standard odd two-braid torus diagram
    lower = (total + 2) // 3 + 1
    assert lower == row["crosscap_lower_KL16"] and lower > 3
    return {"name": row["name"], "prime_alternating_diagram": True,
            "complete_bracket_states": 1 << len(pd),
            "Jones_second_coefficient_absolute_values": [beta, beta_prime],
            "T": total, "crosscap_lower_KL16": lower,
            "crosscap_upper_if_twisted_and_bridge_four": 3,
            "outside_standard_BBK21_twisted_link_class": True}


def rejection(call):
    try:
        call()
    except (AssertionError, ValueError, IndexError, KeyError):
        return True
    raise AssertionError("Invalid control was accepted")


def main():
    data = json.loads(INPUT.read_text())
    assert [r["name"] for r in data["cases"]] == EXPECTED
    verify_pins(data["source_sha256"])
    named_pds = {r["name"]: r["pd"] for r in json.loads(
        (HERE / "seed_cover_topology_inputs.json").read_text())["cases"]}
    seven = json.loads((HERE / "seven_sheet_cover_verification.json").read_text())
    periodic = json.loads((HERE / "periodic_twist_family_verification.json").read_text())
    assert seven["status"] == "PASS" and periodic["status"].startswith("PASS")
    a6 = json.loads((HERE / "a6_affine_inputs.json").read_text())
    a6_proof = json.loads((HERE / "a6_affine_verification.json").read_text())
    assert a6_proof["status"] == "PASS" and a6_proof["named_cases"] == 24
    assert a6_proof["input_sha256"]["a6_affine_inputs.json"] == digest(HERE / "a6_affine_inputs.json")
    named_pds.update({r["name"]: r["pd"] for r in a6["cases"]})
    bridge_rows = seven["results"] + periodic["families"] + a6_proof["results"]
    bridge_four = {r["name"] for r in bridge_rows
                   if r["bridge_upper"] == r["meridional_lower"] == 4}
    assert set(EXPECTED) <= bridge_four
    rows = [check_case(row, named_pds) for row in data["cases"]]
    trefoil = [[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]]
    figure_eight = [[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]]
    assert bracket(trefoil) == {-5: -1, 3: -1, 7: 1}
    assert bracket(figure_eight) == {-8: 1, -4: -1, 0: 1, 4: -1, 8: 1}
    composite = [[1, 4, 2, 5], [5, 2, 6, 3], [3, 6, 4, 7],
                 [7, 10, 8, 11], [11, 8, 0, 9], [9, 0, 10, 1]]
    assert projection_edge_connectivity(composite) == 2
    bad = deepcopy(data["cases"][0])
    bad["bracket_coefficients_A"][0][1] += 1
    bad_pins = dict(data["source_sha256"])
    bad_pins[next(iter(bad_pins))] = "0" * 64
    negative = {"changed_bracket_polynomial": rejection(lambda: check_case(bad, named_pds)),
                "composite_diagram_as_prime": rejection(lambda: assert_prime_diagram(composite)),
                "changed_source_digest": rejection(lambda: verify_pins(bad_pins))}
    pins = [INPUT, Path(__file__)] + [HERE / p for p in data["source_sha256"]]
    report = {"status": "PASS", "python": platform.python_version(),
              "command": "python3 verify_twisted_link_crosscap.py",
              "scope": data["scope"], "cases": rows,
              "positive_controls": ["trefoil exact bracket", "figure-eight exact bracket"],
              "negative_controls": negative,
              "input_sha256": {str(p.relative_to(HERE)): digest(p) for p in pins},
              "assurance": "Internal exact finite replay; KL16, BBK21 and the stated diagram criterion imported. Earlier knot identity and bridge proofs are bound by retained inputs/receipts, not recomputed here. No formal or external review."}
    (HERE / "twisted_link_crosscap_verification.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "globally_excluded_standard_twisted_links": EXPECTED,
                      "negative_controls": negative}, indent=2))


if __name__ == "__main__":
    main()
