#!/usr/bin/env python3
"""Finite controls for the written all-rank proofs; no knot-table search."""

from collections import deque
from itertools import combinations, product
from math import factorial, gcd
import hashlib
import json
from pathlib import Path

from signed_lifts import nullspace, orbit, relator_rows, sign_potentials, signed_swap, solve

HERE = Path(__file__).resolve().parent


def pretzel_quotient(n, first=2):
    edges = [(0, n - 1), (0, n - 1)] + [(i, n - 1) for i in range(1, n - 1)]
    parameters = [first] + [3] * (n - 1)
    relators = [[i + 1, i + 1] for i in range(n)]
    relators += [[i + 1, (i + 1) % n + 1] * abs(p) for i, p in enumerate(parameters)]
    return {"label": f"pretzel Coxeter quotient P({','.join(map(str, parameters))})",
            "degree": n, "transpositions": edges, "relators": relators,
            "provenance": "Baader Blair Kjuchukova 2021 section 2; cycle quotient, not a full knot presentation"}


def generated_group_order(n, edges, signs):
    points = [x for x in range(1 << n) if x.bit_count() % 2 == 0]
    index = {x: i for i, x in enumerate(points)}
    generators = [tuple(index[signed_swap(x, p, q, (signs >> j) & 1)] for x in points)
                  for j, (p, q) in enumerate(edges)]
    identity = tuple(range(len(points)))
    seen, queue = {identity}, deque([identity])
    while queue:
        g = queue.popleft()
        for s in generators:
            h = tuple(s[i] for i in g)
            if h not in seen:
                seen.add(h)
                queue.append(h)
    return len(seen)


def main():
    report = {"status": "PASS", "formal_verification": False, "external_review": False,
              "scope": "finite controls; uniform arguments and published knot quotient are in PROOF.tex"}
    # Exercise nontrivial elimination beyond the rank-zero pretzel lift spaces.
    matrix_cases = 0
    for rows in product(range(8), repeat=3):
        kernel = nullspace(list(rows), 3)
        generated = {0}
        for v in kernel:
            generated |= {x ^ v for x in list(generated)}
        expected = {x for x in range(8) if all((row & x).bit_count() % 2 == 0 for row in rows)}
        assert generated == expected
        matrix_cases += 1
    report["exhaustive_three_by_three_binary_matrices"] = matrix_cases

    # Independently evaluate affine words on every sheet against symbolic rows.
    edges = pretzel_quotient(4)["transpositions"]
    word_cases = 0
    for length in range(7):
        for word in product(range(1, 5), repeat=length):
            try:
                rows = relator_rows(4, edges, [list(word)])
            except ValueError:
                continue  # its permutation part is not a relator
            for signs in range(16):
                predicted = all((row & signs).bit_count() % 2 == 0 for row in rows)
                actual = True
                for x in range(16):
                    if x.bit_count() % 2:
                        continue
                    y = x
                    for letter in reversed(word):
                        p, q = edges[letter - 1]
                        y = signed_swap(y, p, q, (signs >> (letter - 1)) & 1)
                    actual = actual and y == x
                assert actual == predicted
                word_cases += 1
    report["relator_word_sign_cases_through_length_six"] = word_cases
    # Exhaust every simple graph and every signing on three and four vertices.
    cases = 0
    for n in (3, 4):
        complete = list(combinations(range(n), 2))
        for mask in range(1 << len(complete)):
            edges = [e for i, e in enumerate(complete) if (mask >> i) & 1]
            if sign_potentials(n, edges, 0)[0] != 1:
                continue
            for signs in range(1 << len(edges)):
                balanced = sign_potentials(n, edges, signs)[1] is not None
                assert (len(orbit(n, edges, signs)) == 1 << (n - 1)) == (not balanced)
                cases += 1
    report["exhaustive_connected_simple_signed_graphs_n3_n4"] = cases

    examples, negative = [], []
    for n in range(3, 15):
        data = pretzel_quotient(n)
        result = solve(data)
        assert result["lift_space_dimension"] == n
        assert result["unbalanced_lift_dimension"] == 1
        assert result["transitive_signed_lift_exists"]
        # Independent even-parity bit-vector orbit, without the linear solver.
        assert len(orbit(n, data["transpositions"], 1)) == 1 << (n - 1)
        result["form_radical_dimension"] = int(n % 2 == 0)
        result["nondegenerate_original_theorem_applies"] = n % 2 == 1
        examples.append(result)
        altered = solve(pretzel_quotient(n, first=3))
        assert altered["lift_space_dimension"] == n - 1
        assert not altered["transitive_signed_lift_exists"]
        negative.append({"degree": n, "fixed_symmetric_quotient_has_no_unbalanced_lift": True})
    report["pretzel_quotients"] = examples
    report["changed_even_twist_to_odd_negative_controls"] = negative

    # All signings of graphs with a repeated edge, including zero translations.
    parallel_cases = 0
    for n in range(3, 9):
        edges = pretzel_quotient(n)["transpositions"]
        for signs in range(1 << n):
            balanced = sign_potentials(n, edges, signs)[1] is not None
            assert balanced == (((signs >> 0) ^ (signs >> 1)) & 1 == 0)
            assert (len(orbit(n, edges, signs)) == 1 << (n - 1)) == (not balanced)
            parallel_cases += 1
    report["parallel_edge_signings"] = parallel_cases

    orders = []
    for n in (3, 4, 5):
        edges = pretzel_quotient(n)["transpositions"]
        size = generated_group_order(n, edges, 1)
        assert size == (1 << (n - 1)) * factorial(n)
        orders.append({"degree": n, "permutation_group_order": size})
    report["independent_group_enumerations"] = orders

    # For even n, fixed-point failure alone can give a false positive.
    n = 4
    edges = [(0, 3), (1, 3), (2, 3)]
    signs = 1  # ambient potential b=e0 has odd parity
    fixed = [x for x in range(1 << n) if x.bit_count() % 2 == 0
             and all(signed_swap(x, p, q, (signs >> i) & 1) == x
                     for i, (p, q) in enumerate(edges))]
    assert not fixed
    assert sign_potentials(n, edges, signs)[1] is not None
    assert len(orbit(n, edges, signs)) == 4 < 8
    report["balanced_even_degree_without_common_fixed_point"] = {
        "degree": n, "sheets": 8, "zero_orbit": 4, "common_fixed_points": 0,
        "ambient_potential_exists": True}

    # Reject a false symmetric quotient instead of constructing an affine lift.
    bad = pretzel_quotient(4)
    bad["relators"].append([1, 3])
    try:
        solve(bad)
    except ValueError as error:
        assert "fails relator" in str(error)
    else:
        raise AssertionError("invalid symmetric quotient accepted")
    report["invalid_symmetric_relator_rejected"] = True

    # Exact group-order inputs to the uniform torus obstruction.
    torus = []
    for p, q in [(3, 5), (4, 11), (5, 17), (6, 37), (7, 67), (8, 131)]:
        assert all(q % d for d in range(2, q))
        assert q > (1 << (p - 1)) - 1 and gcd(p, q) == 1
        orders = []
        for k in range(2, p):
            agl = (1 << k)
            for i in range(k):
                agl *= (1 << k) - (1 << i)
            assert gcd(q, agl) == 1
            orders.append({"dimension": k, "AGL_order": agl, "gcd_with_q": 1})
        torus.append({"p": p, "q": q, "excluded_admissible_dimensions": orders})
    report["torus_affine_obstruction_controls"] = torus

    inputs = ["signed_lifts.py", "verify_signed_lifts.py", "PROOF.tex", "LITERATURE.tex", "sources.json"]
    report["input_sha256"] = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}
    (HERE / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    (HERE / "example_pretzel_5.json").write_text(json.dumps(pretzel_quotient(5), indent=2) + "\n")
    print(json.dumps({"status": "PASS", "simple_signed_graph_cases": cases,
                      "parallel_sign_cases": parallel_cases, "pretzel_degrees": [3, 14],
                      "group_orders": report["independent_group_enumerations"],
                      "report": str(HERE / "verification.json")}, indent=2))


if __name__ == "__main__":
    main()
