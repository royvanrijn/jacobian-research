"""Verify the branch correction and the 7/12 FullOC exponent certificate."""

from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from types import ModuleType


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "sources" / "zeng-abc-modulus-locking-24407593"
TREE_FILE = SOURCE_DIR / "cert_fulloc_5840_improved_73_125.json.gz"
CHECKER_FILE = SOURCE_DIR / "check_extended_pairdev.py"
DELTA_FILE = ROOT / "abc_pair_deviation_7_12_delta.json"
TREE_SHA256 = "a94bbf1c0a8bd1267f9be28cd862d980e2fe1eab7b74f2b6760d76d342ae8069"
CHECKER_SHA256 = "a1ac4cd83a296c55e2c9a245e845d66dc7935f8624febdd25a61a81104b31437"
REVERSE_PAIRS = (("a", "b"), ("a", "c"), ("b", "c"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def source_bytes(path, sha256):
    contents = path.read_bytes()
    require(hashlib.sha256(contents).hexdigest() == sha256, f"SHA-256 mismatch for {path}")
    return contents


def follow(node, path):
    for index in path:
        node = node["children"][index]
    return node


def leaf_count(node):
    if node["type"] == "leaf":
        return 1
    return sum(leaf_count(child) for child in node["children"])


def check_gap(checker, published_tree, pairdev_path):
    numerators = (
        2365, 5523, 11266, 5500, 3388, 0, 28042,
        1573, 5754, 11266, 5720, 3278, 0, 27591,
        5971, 3388, 11277, 9829, 51, 0, 30516,
        11266, 11474, 13217, 50309,
    )
    point = tuple(Fraction(numerator, 86149) for numerator in numerators)
    mode = published_tree["mode"]
    rows = checker.base(Fraction(published_tree["lambda"]))
    node = published_tree["tree"]
    for index in pairdev_path:
        rows.extend(checker.branches(node, mode)[index])
        node = node["children"][index]

    def satisfies(inequalities):
        return all(checker.dot(row, point) <= bound for row, bound in inequalities)

    require(node["type"] == "pairdev" and Fraction(node["eps"]) == Fraction(1, 100),
            "unexpected branch")
    require(len(rows) == 115 and min(point) >= 0 and satisfies(rows),
            "point must obey the root and ancestor rows")
    branches = checker.branches(node, mode)
    require(len(branches) == 4 and not any(map(satisfies, branches)),
            "point must miss all four original children")
    excess = point[checker.D] - Fraction(published_tree["max_certified_D"])
    require(excess == Fraction(268, 123180763), "unexpected excess")
    return node, {"ancestor_rows": len(rows), "original_children_satisfied": 0,
                  "D_excess": str(excess)}


def main():
    checker = ModuleType("zeng_24407593_checker")
    exec(compile(source_bytes(CHECKER_FILE, CHECKER_SHA256), str(CHECKER_FILE), "exec"),
         checker.__dict__)
    certificate = json.loads(gzip.decompress(source_bytes(TREE_FILE, TREE_SHA256)))
    delta = json.loads(DELTA_FILE.read_text(encoding="utf-8"))
    require(delta["schema_version"] == 1, "unknown certificate delta version")
    require(certificate["lambda"] == "1" and certificate["mode"] == "all_deg",
            "unexpected source exponent program")

    pairdev, gap = check_gap(checker, certificate, delta["pairdev_path"])
    print(json.dumps({"published_tree_gap": gap}, sort_keys=True), flush=True)
    reverse_children = delta["reverse_children"]
    require(len(pairdev["children"]) == 4 and len(reverse_children) == 3,
            "wrong number of deviation children")
    pairdev["children"].extend(reverse_children)
    added_leaves = sum(map(leaf_count, reverse_children))
    require(added_leaves == 19, "unexpected reverse subtree leaf count")

    replacements = delta["leaf_replacements"]
    require([item["original_bound"] for item in replacements] == ["1689/2893", "5845/10009"],
            "unexpected replacement list")
    for item in replacements:
        parent = follow(certificate["tree"], item["path"][:-1])
        index = item["path"][-1]
        old_leaf = parent["children"][index]
        require(old_leaf["type"] == "leaf" and old_leaf["bound"] == item["original_bound"],
                "replacement path does not identify the expected original leaf")
        parent["children"][index] = item["subtree"]
        added_leaves += leaf_count(item["subtree"]) - 1

    original_branches = checker.branches

    def repaired_branches(node, mode):
        branches = original_branches(node, mode)
        if node["type"] == "pairdev":
            epsilon = Fraction(node["eps"])
            require(epsilon >= 0, "negative pair deviation threshold")
            for low, high in REVERSE_PAIRS:
                coeffs = checker.add((1, checker.coords(low, [3, 4])),
                                     (-1, checker.coords(high, [3, 4])))
                branches.append([checker.row(coeffs, -epsilon)])
        return branches

    checker.branches = repaired_branches
    certificate["leaves"] += added_leaves
    certificate["max_certified_D"] = delta["replacement_upper_bound"]
    result = checker.verify(certificate)
    require(result["pairdev"] == 1 and result["leaf"] == 84177 and result["upper"] == "7/12",
            "unexpected repaired exponent certificate result")
    print(json.dumps({"repaired_tree": result}, sort_keys=True))


if __name__ == "__main__":
    main()
