"""Verify the FullOC certificate strengthened by grouped-root collision bounds."""

from fractions import Fraction
from itertools import product
import gzip
import json
from math import gcd
from pathlib import Path

from verify_abc_pair_deviation_certificate import (
    load_repaired_certificate,
    require,
    verify_repaired_certificate,
)


DELTA_FILE = Path(__file__).with_name("abc_grouped_fulloc_delta.json.gz")


def sparse_branches(checker, node):
    """Expand the diagonal bound and all terms of the maximum for the off-diagonal bound."""
    require(set(node) == {"type", "A", "B", "G", "H", "I", "O", "cap", "children"},
            "unexpected grouped-root node keys")
    a, b = node["A"], node["B"]
    letters = ("a", "b", "c")
    require(a in letters and b in letters and a != b, "invalid letter pair")
    c = next(letter for letter in letters if letter not in (a, b))
    first, second, free_b, fixed_c = (
        checker.subset(node[field], field == "O") for field in ("G", "H", "I", "O")
    )
    require(not set(first).intersection(second), "overlapping root groups")

    add, coords, var = checker.add, checker.coords, checker.var

    def weight(letter, indices, divisor=1):
        return add(*[(Fraction(index, divisor), var(f"{letter}{index}")) for index in indices])

    frozen_a = add((1, var(a)), (-1, coords(a, checker.R)))
    mass_b = coords(b, free_b)
    frozen_b = add((1, var(b)), (-1, mass_b))
    frozen_c = coords(c, fixed_c)
    modulus_c = weight(c, fixed_c)
    modulus_b = add((7, var(b)), *[
        (-7 if index in free_b else index - 7, var(f"{b}{index}")) for index in checker.R
    ])
    energy = add((2, coords(a, checker.R)), (-1, coords(a, first)), (-1, coords(a, second)))
    positive_parts = [
        (Fraction(1, 4), add((2, weight(a, group, gcd(*group))), (-1, modulus_b)))
        for group in (first, second)
    ]

    cap = node["cap"]
    require(isinstance(cap, list) and cap, "invalid residue or difference bound")
    if cap == ["baseline"]:
        difference = mass_b
    elif len(cap) == 2 and cap[0] == "diff" and type(cap[1]) is int and 2 <= cap[1] <= 6:
        difference = add((2, coords(b, [index for index in free_b if index % cap[1]])))
    elif cap == ["value"]:
        difference = add()
        positive_parts.append((Fraction(1, 2), add((1, weight(b, free_b)), (-1, modulus_c))))
    elif len(cap) == 2 and cap[0] == "group":
        group = checker.subset(cap[1], False)
        require(set(group).issubset(free_b), "residue group outside free coordinates")
        difference = coords(b, [index for index in free_b if index not in group])
        positive_parts.append((Fraction(1, 2), add((1, weight(b, group, gcd(*group))), (-1, modulus_c))))
    else:
        raise ValueError("unknown residue or difference bound")

    diag = add((1, var("D")), (-1, frozen_a), (-1, frozen_b), (-1, var(c)))
    off = add((1, var("D")), (-1, frozen_a), (-1, frozen_b), (-Fraction(1, 2), var(c)),
              (-Fraction(1, 2), frozen_c), (-Fraction(1, 2), energy), (-Fraction(1, 2), difference))
    branches = [[checker.row(diag)]]
    for switches in product((0, 1), repeat=len(positive_parts)):
        row = add((1, off), *[
            (-coefficient, term)
            for selected, (coefficient, term) in zip(switches, positive_parts, strict=True)
            if selected
        ])
        branches.append([checker.row(row)])
    return branches


def check_symmetric_profile(checker):
    node = {"type": "sparse", "A": "a", "B": "b", "G": [2, 4], "H": [3],
            "I": [3], "O": [], "cap": ["diff", 3], "children": []}
    coordinates = (19, 18, 54, 36, 0, 0, 132)
    point = [Fraction(value, 396) for value in coordinates * 3] + [Fraction(0)] * 4
    rows = sparse_branches(checker, node)
    require(max(-checker.dot(row, point) for [(row, _)] in rows) == Fraction(37, 66),
            "unexpected grouped-root bound on the symmetric profile")


def iter_leaves(node, path=()):
    if node["type"] == "leaf":
        yield path, node
        return
    for index, child in enumerate(node["children"]):
        yield from iter_leaves(child, path + (index,))


def ancestor_rows(checker, tree, path):
    rows = checker.base(Fraction(tree["lambda"]))
    node = tree["tree"]
    for index in path:
        rows.extend(checker.branches(node, tree["mode"])[index])
        node = node["children"][index]
    return rows


def check_subtree(checker, mode, node, rows, limit, statistics):
    kind = node["type"]
    require(kind in statistics, "unknown replacement node type")
    statistics[kind] += 1
    if kind == "leaf":
        require(set(node) == {"type", "y", "bound"}, "unexpected replacement leaf keys")
        require(isinstance(node["y"], dict) and node["y"], "empty leaf multipliers")
        coefficients = [Fraction(0)] * checker.N
        rhs = Fraction(0)
        for index_text, multiplier in node["y"].items():
            require(index_text.isdigit(), "invalid leaf row index")
            index = int(index_text)
            value = Fraction(multiplier)
            require(str(index) == index_text and index < len(rows) and value > 0,
                    "invalid leaf row or multiplier")
            row, bound = rows[index]
            coefficients = [a + value * b for a, b in zip(coefficients, row, strict=True)]
            rhs += value * bound
        require(coefficients[checker.D] == 1 and min(coefficients) >= 0,
                "invalid leaf dual coefficients")
        require(rhs == Fraction(node["bound"]) <= limit, "replacement leaf bound mismatch")
        return rhs

    branches = sparse_branches(checker, node) if kind == "sparse" else checker.branches(node, mode)
    require(isinstance(node["children"], list) and len(node["children"]) == len(branches),
            "replacement branch coverage mismatch")
    return max(
        check_subtree(checker, mode, child, rows + extra, limit, statistics)
        for child, extra in zip(node["children"], branches, strict=True)
    )


def check_budget_envelope(originals, replacements):
    def line(leaf):
        slope = Fraction(leaf["y"].get("0", "0"))
        return Fraction(leaf["bound"]) - slope, slope

    effective = {line(leaf) for path, leaf in originals.items() if path not in replacements}
    for subtree in replacements.values():
        effective.update(line(leaf) for _, leaf in iter_leaves(subtree))

    flat = Fraction(33678, 57833), Fraction(0)
    rising = Fraction(65607, 413150), Fraction(87492, 206575)
    require({flat, rising}.issubset(effective), "budget envelope lines are absent")
    crossing = (flat[0] - rising[0]) / rising[1]
    for budget in (Fraction(0), crossing, Fraction(1)):
        maximum = max(intercept + slope * budget for intercept, slope in (flat, rising))
        require(all(intercept + slope * budget <= maximum for intercept, slope in effective),
                "budget envelope fails on an endpoint")
    return {"flat": str(flat[0]), "intercept": str(rising[0]), "slope": str(rising[1]),
            "switch": str(crossing)}


def main():
    checker, tree, _ = load_repaired_certificate()
    original = verify_repaired_certificate(checker, tree)
    print(json.dumps({"repaired_fulloc_tree": original}, sort_keys=True), flush=True)
    check_symmetric_profile(checker)

    delta = json.loads(gzip.decompress(DELTA_FILE.read_bytes()))
    require(delta["schema_version"] == 1, "unknown grouped-root delta schema")
    limit = Fraction(delta["upper_bound"])
    require(0 < limit < Fraction(original["upper"]), "expected a strict improvement")
    originals = dict(iter_leaves(tree["tree"]))
    old = {path: Fraction(node["bound"]) for path, node in originals.items()}
    require(len(old) == original["leaf"], "unexpected baseline leaf count")
    replacements = {}
    for entry in delta["replacements"]:
        require(all(type(i) is int and i >= 0 for i in entry["path"]), "invalid replacement path")
        path = tuple(entry["path"])
        require(path not in replacements and old.get(path) == Fraction(entry["original_bound"]),
                "duplicate path or unexpected original leaf bound")
        replacements[path] = entry["subtree"]
    require(replacements.keys() == {path for path, bound in old.items() if bound > limit},
            "replacements must cover precisely all baseline leaves above the new threshold")

    statistics = {"leaf": 0, "geo": 0, "sparse": 0}
    upper = max(bound for path, bound in old.items() if path not in replacements)
    for path, subtree in replacements.items():
        rows = ancestor_rows(checker, tree, path)
        upper = max(upper, check_subtree(checker, tree["mode"], subtree, rows, limit, statistics))
    print(json.dumps({"grouped_fulloc_tree": {"upper": str(upper), "threshold": str(limit),
          "replaced_leaves": len(replacements), "checked_total_leaves": len(old) - len(replacements)
                                                        + statistics["leaf"], **statistics}},
          sort_keys=True))
    print(json.dumps({"budget_envelope": check_budget_envelope(originals, replacements)}, sort_keys=True))


if __name__ == "__main__":
    main()
