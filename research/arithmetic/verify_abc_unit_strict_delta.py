"""Verify sparse exact improvements of the pinned unit abc certificate."""

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
    return replay_replacements(certificate, delta["replacements"], limit,
                               Fraction(1), BASELINE_LEAVES)


def replay_replacements(certificate, entries, limit, budget, source_leaves):
    """Check replacement branches and the exact bounds of all unchanged leaves."""
    checker = load_checker()
    source_budget = Fraction(certificate["budget"])
    source_limit = Fraction(certificate["threshold"])
    budget_change = budget - source_budget
    schema = certificate["schema_version"]
    rows = certificate_root_rows(checker, budget, schema)
    if budget_change:
        source_rows = certificate_root_rows(checker, source_budget, schema)
        require(rows[1:] == source_rows[1:] and rows[0][0] == source_rows[0][0]
                and rows[0][1] == budget and source_rows[0][1] == source_budget,
                "the budget must affect only the first row's right-hand side")
    root_rows = sparse_rows(rows)

    @lru_cache(maxsize=131072)
    def leaf_bound(bound, slope):
        return Fraction(bound) + budget_change * Fraction(slope)

    def source_bound(node):
        return leaf_bound(node["bound"], node["y"].get("0", 0) if budget_change else 0)

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
    for entry in entries:
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
        require(node["type"] == "leaf" and limit < source_bound(node) <= source_limit,
                "replacement does not target an over-limit source leaf")
        upper_patched = max(upper_patched, verify_replacement(entry["tree"], rows))
        replacements[path] = (node, entry["tree"])

    baseline_count = 0
    upper_untouched = Fraction(0)
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
        bound = source_bound(node)
        upper_untouched = max(upper_untouched, bound)
        if budget_change:
            node["bound"] = str(bound)
    require(baseline_count == source_leaves and not remaining and upper_untouched <= limit,
            "the sparse replacements do not cover every over-limit leaf")
    upper = max(upper_untouched, upper_patched)
    for node, tree in replacements.values():
        node.clear()
        node.update(tree)
    certificate["budget"] = str(budget)
    certificate["threshold"] = str(limit)
    result = {"budget": str(budget), "upper": str(upper), "upper_patched": str(upper_patched),
              "upper_untouched": str(upper_untouched), "replacements": len(replacements),
              "leaf": baseline_count - len(replacements) + patched_counts["leaf"]}
    return result


def verify_overlay(certificate, overlay, strict_sha256, strict_leaf_count):
    """Replay a downward retarget of the verified strict unit tree."""
    fields = {"schema_version", "unit_baseline_sha256", "unit_strict_delta_sha256",
              "source_budget", "target_budget", "target", "replacements"}
    require(set(overlay) == fields and overlay["schema_version"] == 1,
            "invalid unit overlay header")
    require(overlay["unit_baseline_sha256"] == BASELINE_SHA256
            and overlay["unit_strict_delta_sha256"] == strict_sha256,
            "unit overlay source hash mismatch")
    old, new = Fraction(overlay["source_budget"]), Fraction(overlay["target_budget"])
    limit = Fraction(overlay["target"])
    require(old == Fraction(certificate["budget"]) == 1 and 0 < new < old
            and 0 < limit < Fraction(certificate["threshold"]),
            "invalid unit overlay budgets or exponent")
    require(isinstance(overlay["replacements"], list), "invalid unit overlay replacements")
    result = replay_replacements(certificate, overlay["replacements"], limit, new, strict_leaf_count)
    require(Fraction(result["upper"]) == limit, "unit overlay replay disagrees with its target")
    return result


def reweight_verified_tree(certificate, at_budget):
    budget_change = at_budget - Fraction(certificate["budget"])

    @lru_cache(maxsize=131072)
    def reweighted_bound(bound, slope):
        return Fraction(bound) + budget_change * Fraction(slope)

    pending = [certificate["tree"]]
    upper = None
    while pending:
        node = pending.pop()
        if node["type"] != "leaf":
            pending.extend(node["children"])
            continue
        bound = reweighted_bound(node["bound"], node["y"].get("0", 0))
        if upper is None or bound > upper:
            upper = bound
    require(upper is not None, "the reconstructed tree has no leaves")
    return upper


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="also replay every leaf of the reconstructed tree")
    parser.add_argument("--overlay", type=Path, help="replay a pinned downward unit-budget overlay")
    parser.add_argument("--at-budget", type=Fraction,
                        help="reweight the verified tree at this radical budget")
    args = parser.parse_args()
    if args.at_budget is not None and not 0 < args.at_budget <= 1:
        parser.error("--at-budget must be in (0, 1]")
    certificate = pinned_baseline()
    delta_bytes = DELTA.read_bytes()
    delta = json.loads(lzma.decompress(delta_bytes))
    result = verify_delta(certificate, delta)
    if args.overlay is not None:
        overlay = json.loads(lzma.decompress(args.overlay.read_bytes()))
        result = verify_overlay(certificate, overlay, sha256(delta_bytes).hexdigest(), result["leaf"])
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
