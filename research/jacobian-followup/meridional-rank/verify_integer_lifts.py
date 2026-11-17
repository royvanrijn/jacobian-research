#!/usr/bin/env python3
"""Independent affine evaluation, Smith-step replay, and modular elimination.

Panel verification does not call the producer's matrix or Smith algorithms.
Knot names are imported from the pinned primary table, not identified here.
"""

from collections import Counter
from itertools import combinations, product
from math import gcd, factorial
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def crossings_from_walk(code):
    n = len(code) // 2
    under, over, arc = {}, {}, -1
    for c in code:
        if c < 0:
            under[-c] = (arc % n, (arc + 1) % n)
            arc += 1
        else:
            over[c] = arc % n
    assert arc + 1 == n and set(under) == set(over) == set(range(1, n + 1))
    return [(over[c], *under[c]) for c in range(1, n + 1)]


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def independent_colourings(code, seeds, degree):
    cr = crossings_from_walk(code)
    permutations = []
    for p, q in combinations(range(degree), 2):
        s = list(range(degree))
        s[p], s[q] = q, p
        permutations.append(tuple(s))
    identity = tuple(range(degree))
    index = {s: i for i, s in enumerate(permutations)}
    out = {}
    for tail in product(range(len(permutations)), repeat=len(seeds) - 1):
        assignment = (0,) + tail
        group, frontier = {identity}, [identity]
        gens = [permutations[i] for i in assignment]
        while frontier:
            a = frontier.pop()
            for b in gens:
                c = compose(a, b)
                if c not in group:
                    group.add(c)
                    frontier.append(c)
        if len(group) != factorial(degree):
            continue
        labels = {s: permutations[i] for s, i in zip(seeds, assignment)}
        changed, consistent = True, True
        while changed and consistent:
            changed = False
            for o, a, b in cr:
                if o in labels and (a in labels or b in labels):
                    old, new = (a, b) if a in labels else (b, a)
                    value = compose(compose(labels[o], labels[old]), labels[o])
                    if new in labels and labels[new] != value:
                        consistent = False
                        break
                    if new not in labels:
                        labels[new] = value
                        changed = True
        if consistent:
            assert len(labels) == len(code) // 2
            out[assignment] = [index[labels[i]] for i in range(len(labels))]
    return out


def independent_matrix(code, degree, labels):
    edges = list(combinations(range(degree), 2))
    count = len(labels)
    output = []
    for o, a, b in crossings_from_walk(code):
        word = (o, a, o, b)
        block = [[0] * count for _ in range(degree)]
        for variable in range(count):
            point = [0] * degree
            for generator in reversed(word):
                p, q = edges[labels[generator]]
                point[p], point[q] = point[q], point[p]
                if generator == variable:
                    point[p] += 1
                    point[q] -= 1
            for i in range(degree):
                block[i][variable] = point[i]
        output.extend(row for row in block if any(row))
    return output


def replay(matrix, cert):
    a = [row[:] for row in matrix]
    m, n = cert["rows"], cert["columns"]
    assert len(a) == m and all(len(row) == n for row in a)
    for step in cert["steps"]:
        axis, operation, *args = step
        assert axis in ("r", "c")
        limit = m if axis == "r" else n
        assert all(0 <= v < limit for v in args[:2] if operation != "neg")
        if operation == "swap":
            i, j = args
            if axis == "r":
                a[i], a[j] = a[j], a[i]
            else:
                a = [[row[j] if k == i else row[i] if k == j else value
                      for k, value in enumerate(row)] for row in a]
        elif operation == "add":
            target, source, multiple = args
            assert target != source and isinstance(multiple, int)
            if axis == "r":
                a[target] = [a[target][j] + multiple * a[source][j] for j in range(n)]
            else:
                a = [[value + multiple * row[source] if j == target else value
                      for j, value in enumerate(row)] for row in a]
        elif operation == "neg":
            assert axis == "r" and len(args) == 1 and 0 <= args[0] < m
            a[args[0]] = [-v for v in a[args[0]]]
        else:
            raise AssertionError("non-unimodular or unknown operation")
    d = cert["diagonal"]
    assert all(v > 0 for v in d) and all(y % x == 0 for x, y in zip(d, d[1:]))
    assert all(a[i][j] == (d[i] if i == j and i < len(d) else 0)
               for i in range(m) for j in range(n))
    assert cert["characteristic_zero_dimension"] == n - len(d)
    assert cert["maximum_dimension_all_characteristics"] == n - sum(v == 1 for v in d)
    return d


def nullspace(matrix, columns, prime):
    a = [[v % prime for v in row] for row in matrix]
    pivot_columns, row = [], 0
    for c in range(columns):
        pivot = next((i for i in range(row, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        factor = pow(a[row][c], -1, prime)
        a[row] = [x * factor % prime for x in a[row]]
        for i in range(len(a)):
            if i != row:
                factor = a[i][c]
                a[i] = [(x - factor * y) % prime for x, y in zip(a[i], a[row])]
        pivot_columns.append(c)
        row += 1
    basis = []
    for free in set(range(columns)) - set(pivot_columns):
        v = [0] * columns
        v[free] = 1
        for i, c in enumerate(pivot_columns):
            v[c] = -a[i][free] % prime
        assert all(sum(x * y for x, y in zip(r, v)) % prime == 0 for r in matrix)
        basis.append(v)
    return basis


def check_seeds(code, result):
    cr = crossings_from_walk(code)
    seen = set(result["seeds"])
    assert len(seen) == result["bridge_upper_bound"]
    for o, old, new in result["seed_trace"]:
        assert o in seen and old in seen and new not in seen
        assert (o, old, new) in cr or (o, new, old) in cr
        seen.add(new)
    assert seen == set(range(len(code) // 2))


def matrix_controls():
    # Independent gcds of minors determine Smith factors for 2x2 matrices.
    from integer_lifts import smith_certificate
    count = 0
    for values in product(range(-2, 3), repeat=4):
        matrix = [list(values[:2]), list(values[2:])]
        cert = smith_certificate(matrix)
        diagonal = replay(matrix, cert)
        common = 0
        for value in values:
            common = gcd(common, value)
        determinant = abs(values[0] * values[3] - values[1] * values[2])
        expected = [] if common == 0 else [common]
        if determinant:
            expected.append(determinant // common)
        assert diagonal == expected
        count += 1
    # A diagonal with distinct prime factors checks the all-prime maximum.
    cert = smith_certificate([[2, 0], [0, 3]])
    assert replay([[2, 0], [0, 3]], cert) == [1, 6]
    assert cert["maximum_dimension_all_characteristics"] == 1
    assert smith_certificate([], 4)["maximum_dimension_all_characteristics"] == 4
    return count


def main():
    controls = matrix_controls()
    summaries, sharp = [], []
    inputs = {"verify_integer_lifts.py", "integer_lifts.py", "search_integer_lifts.py",
              "search_old_panel_integer.py", "sources/current13/manifest.json",
              "sources/wirt-hm/manifest.json", "panel_source_verification.json"}
    all_certificates = 0
    for panel_name, result_name in [("current13_panel.json", "current13_integer_results.json"),
                                    ("lift_panel.json", "old_panel_integer_results.json")]:
        panel = json.loads((HERE / panel_name).read_text())
        report = json.loads((HERE / result_name).read_text())
        inputs.update([panel_name, result_name])
        for path, digest in report["input_sha256"].items():
            assert hashlib.sha256((HERE / path).read_bytes()).hexdigest() == digest
        assert len(panel["cases"]) == len(report["results"]) == 16
        for case, result in zip(panel["cases"], report["results"]):
            assert case["name"] == result["name"] and result["status"] == "COMPLETE"
            check_seeds(case["gauss_code"], result)
            summary = {"name": case["name"], "panel": panel_name, "degrees": []}
            chosen = False
            for degree_result in result["by_degree"]:
                r = degree_result["degree"]
                colourings = independent_colourings(case["gauss_code"], result["seeds"], r)
                certificates = degree_result["certificates"]
                assert set(colourings) == {tuple(c["seed_assignment"]) for c in certificates}
                assert len(colourings) == len(certificates)
                distributions = Counter()
                for c in certificates:
                    assert colourings[tuple(c["seed_assignment"])] == c["arc_labels"]
                    matrix = independent_matrix(case["gauss_code"], r, c["arc_labels"])
                    d = replay(matrix, c["smith"])
                    columns = len(case["gauss_code"]) // 2
                    for p in [2, 3, 5, 7, 11]:
                        basis = nullspace(matrix, columns, p)
                        assert len(basis) == columns - sum(v % p != 0 for v in d)
                    bound = columns - d.count(1)
                    distributions[bound] += 1
                    all_certificates += 1
                    if not chosen and bound == result["bridge_upper_bound"]:
                        assert r == 4 and d == [1] * 10 + [5] and columns == 15
                        basis = nullspace(matrix, columns, 5)
                        assert len(basis) == 5
                        # Evaluation on the certified five seeds is itself injective.
                        evaluated = [[v[s] for s in result["seeds"]] for v in basis]
                        assert not nullspace(evaluated, 5, 5)
                        sharp.append({"name": case["name"], "gauss_code": case["gauss_code"],
                                      "source_row": case["source_row"], "degree": r, "prime": 5,
                                      "arc_transpositions": [list(combinations(range(r), 2))[i] for i in c["arc_labels"]],
                                      "lift_basis": basis, "smith_diagonal": d,
                                      "seeds": result["seeds"], "seed_trace": result["seed_trace"],
                                      "certified_bounds": {"meridional_lower": 5, "bridge_upper": 5},
                                      "novel_resolution": "NOT CLAIMED; prior class and independent knot-label audits incomplete"})
                        chosen = True
                assert dict(sorted(distributions.items())) == {
                    int(k): v for k, v in degree_result["maximum_over_all_characteristics_counts"].items()}
                summary["degrees"].append({"degree": r, "surjective_colourings": len(colourings),
                                            "maximum_dimensions": dict(sorted(distributions.items()))})
            summaries.append(summary)
    assert len(sharp) == 16 and all_certificates == 512
    for directory in ["sources/current13", "sources/wirt-hm"]:
        manifest = json.loads((HERE / directory / "manifest.json").read_text())
        for item in manifest["files"]:
            assert hashlib.sha256((HERE / directory / item["path"]).read_bytes()).hexdigest() == item["sha256"]
    sealed = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in sorted(inputs)}
    (HERE / "sharp_five_certificates.json").write_text(json.dumps({"claim_boundary":
        "Exact equalities for the pinned table diagrams, using the standard Wirtinger-number bridge theorem. Names follow the primary workbook. No new-priority claim or full MRC proof.",
        "input_sha256": sealed, "certificates": sharp}, indent=2) + "\n")
    verification = {"status": "PASS", "input_sha256": sealed,
                    "sharp_certificates_sha256": hashlib.sha256((HERE / "sharp_five_certificates.json").read_bytes()).hexdigest(),
                    "smith_two_by_two_gcd_of_minors_controls": controls,
                    "independently_replayed_lift_matrices": all_certificates,
                    "modular_cross_checks": [2, 3, 5, 7, 11],
                    "sharp_five_certificates": len(sharp), "panel_results": summaries,
                    "full_conjecture_closed": False,
                    "independent_of_producer": "Separate diagram walk, full permutation-group enumeration, reverse affine point evaluation, elementary-step replay, modular elimination. Uniform proofs remain unformalized and externally unreviewed."}
    (HERE / "integer_lift_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps({k: verification[k] for k in ["status", "smith_two_by_two_gcd_of_minors_controls",
          "independently_replayed_lift_matrices", "sharp_five_certificates", "full_conjecture_closed"]}, indent=2))


if __name__ == "__main__":
    main()
