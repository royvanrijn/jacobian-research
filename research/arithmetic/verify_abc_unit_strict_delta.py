"""Verify the sparse exact improvement of the pinned 107/200 abc certificate."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import lzma
from pathlib import Path

from verify_abc_height_fulloc_certificate import (
    IMPROVED_UNIT_CERTIFICATE_PARTS,
    branch_rows,
    certificate_root_rows,
    check_leaf,
    load_checker,
    require,
    sparse_rows,
    verify,
)


BASELINE_SHA256 = "09038226bc3dbcd5c5339f88fb4a73bebd0fd62764599733ad266111af08538b"
BASELINE_BOUND = Fraction(107, 200)
BASELINE_LEAVES = 1822745
DIRECTORY = Path(__file__).resolve().parent
DELTA = DIRECTORY / "abc_unit_strict_delta.json.xz"


def pinned_baseline():
    payload = b"".join(lzma.decompress((DIRECTORY / name).read_bytes())
                       for name in IMPROVED_UNIT_CERTIFICATE_PARTS)
    require(sha256(payload).hexdigest() == BASELINE_SHA256, "unit baseline does not match the checked bytes")
    certificate = json.loads(payload)
    require(certificate["schema_version"] == 7 and Fraction(certificate["budget"]) == 1
            and Fraction(certificate["threshold"]) == BASELINE_BOUND, "unexpected unit baseline metadata")
    return certificate


def verify_delta(certificate, delta):
    require(set(delta) == {"schema_version", "baseline_sha256", "source_bound", "target", "replacements"}
            and delta["schema_version"] == 1 and delta["baseline_sha256"] == BASELINE_SHA256
            and Fraction(delta["source_bound"]) == BASELINE_BOUND, "invalid sparse certificate header")
    limit = Fraction(delta["target"])
    require(0 < limit < BASELINE_BOUND, "the sparse certificate must improve the baseline bound")
    checker = load_checker()
    root_rows = sparse_rows(certificate_root_rows(checker, Fraction(1), certificate["schema_version"]))

    @lru_cache(maxsize=8192)
    def expand(metadata):
        return tuple(sparse_rows(rows) for rows in branch_rows(checker, json.loads(metadata), True))

    def regions(node):
        metadata = json.dumps({key: value for key, value in node.items() if key != "children"},
                              sort_keys=True, separators=(",", ":"))
        rows = expand(metadata)
        require(isinstance(node.get("children"), list) and len(node["children"]) == len(rows),
                "missing affine branch")
        return rows

    patched_counts = Counter()

    def verify_replacement(node, rows):
        kind = node["type"]
        patched_counts[kind] += 1
        if kind == "leaf":
            return check_leaf(checker, node, rows, limit)
        return max(verify_replacement(child, rows + appended)
                   for child, appended in zip(node["children"], regions(node), strict=True))

    replacements = {}
    upper_patched = Fraction(0)
    for entry in delta["replacements"]:
        require(set(entry) == {"path", "tree"} and isinstance(entry["path"], list)
                and all(type(step) is int and step >= 0 for step in entry["path"]),
                "invalid replacement path")
        path = tuple(entry["path"])
        require(path not in replacements, "duplicate replacement path")
        node, rows = certificate["tree"], root_rows
        for step in path:
            require(node["type"] != "leaf", "replacement path descends beyond its leaf")
            branches = regions(node)
            require(step < len(branches), "replacement child index is out of range")
            rows += branches[step]
            node = node["children"][step]
        require(node["type"] == "leaf" and limit < Fraction(node["bound"]) <= BASELINE_BOUND,
                "replacement does not target an over-limit baseline leaf")
        upper_patched = max(upper_patched, verify_replacement(entry["tree"], rows))
        replacements[path] = (node, entry["tree"])

    baseline_count = 0
    upper_untouched = Fraction(0)
    stored_bounds = {}
    remaining = set(replacements)
    stack = [((), certificate["tree"])]
    while stack:
        path, node = stack.pop()
        if node["type"] != "leaf":
            stack.extend((path + (i,), child) for i, child in enumerate(node["children"]))
            continue
        baseline_count += 1
        if path in remaining:
            remaining.remove(path)
            continue
        bound_text = node["bound"]
        if bound_text not in stored_bounds:
            stored_bounds[bound_text] = Fraction(bound_text)
        upper_untouched = max(upper_untouched, stored_bounds[bound_text])
    require(baseline_count == BASELINE_LEAVES and not remaining and upper_untouched <= limit,
            "the sparse replacements do not cover every over-limit leaf")
    upper = max(upper_untouched, upper_patched)
    for node, tree in replacements.values():
        node.clear()
        node.update(tree)
    certificate["threshold"] = delta["target"]
    result = {"budget": "1", "upper": str(upper), "upper_patched": str(upper_patched),
              "upper_untouched": str(upper_untouched), "replacements": len(replacements),
              "leaf": baseline_count - len(replacements) + patched_counts["leaf"]}
    return result


def reweight_verified_tree(certificate, at_budget):
    budget_change = at_budget - Fraction(certificate["budget"])
    pending = [certificate["tree"]]
    upper = None
    while pending:
        node = pending.pop()
        if node["type"] != "leaf":
            pending.extend(node["children"])
            continue
        bound = Fraction(node["bound"]) + budget_change * Fraction(node["y"].get("0", 0))
        if upper is None or bound > upper:
            upper = bound
    require(upper is not None, "the reconstructed tree has no leaves")
    return upper


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="also replay every leaf of the reconstructed tree")
    parser.add_argument("--at-budget", type=Fraction,
                        help="reweight the verified strict unit tree at this radical budget")
    args = parser.parse_args()
    if args.at_budget is not None and not 0 < args.at_budget <= 1:
        parser.error("--at-budget must be in (0, 1]")
    certificate = pinned_baseline()
    delta = json.loads(lzma.decompress(DELTA.read_bytes()))
    result = verify_delta(certificate, delta)
    if args.full:
        entire = verify(certificate)
        require(result["upper"] == entire["upper"] and result["leaf"] == entire["leaf"],
                "the full replay disagrees with the sparse replay")
        result["full_replay"] = True
    if args.at_budget is not None:
        result["evaluated_budget"] = str(args.at_budget)
        result["reweighted_upper"] = str(reweight_verified_tree(certificate, args.at_budget))
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
