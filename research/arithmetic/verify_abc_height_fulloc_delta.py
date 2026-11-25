"""Replay the exact height refinements of the grouped FullOC certificate."""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import json
from pathlib import Path
import sys

from verify_abc_grouped_fulloc_certificate import (
    DELTA_FILE as GROUPED_DELTA_FILE,
    load_repaired_certificate,
    sparse_branches,
)
from verify_abc_height_fulloc_certificate import LETTERS, branch_rows, check_leaf, sparse_rows
from verify_abc_pair_deviation_certificate import require


DELTA_FILE = Path(__file__).with_name("abc_height_fulloc_delta.jsonl.gz")
BUDGET = Fraction(1)
THRESHOLD = Fraction(29, 50)
SCHEMA = {"schema_version": 1, "budget": str(BUDGET), "threshold": str(THRESHOLD)}


def replacement_index(entries):
    replacements = {}
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {"path", "original_bound", "subtree"},
                "unexpected replacement metadata")
        indices = entry["path"]
        require(isinstance(indices, list) and all(type(i) is int and i >= 0 for i in indices),
                "invalid replacement path")
        path = tuple(indices)
        require(path not in replacements, f"duplicate replacement path {path}")
        replacements[path] = (Fraction(entry["original_bound"]), entry["subtree"])
    return replacements


def read_delta(path):
    with gzip.open(path, "rt", encoding="utf-8") as source:
        require(json.loads(next(source)) == SCHEMA, "unexpected height-delta header")
        return replacement_index(map(json.loads, source))


def verify(height_replacements):
    checker, repaired, _ = load_repaired_certificate()
    require(repaired["mode"] == "all_deg" and Fraction(repaired["lambda"]) == BUDGET
            and Fraction(repaired["max_certified_D"]) == Fraction(7, 12)
            and Fraction(repaired["D_target"]) > Fraction(7, 12),
            "unexpected repaired FullOC program")
    with gzip.open(GROUPED_DELTA_FILE, "rt", encoding="utf-8") as source:
        grouped = json.load(source)
    require(set(grouped) == {"schema_version", "upper_bound", "replacements"}
            and grouped["schema_version"] == 1, "unknown grouped-root delta schema")
    grouped_limit = Fraction(grouped["upper_bound"])
    require(THRESHOLD < grouped_limit < Fraction(7, 12), "unexpected grouped-root threshold")
    grouped_replacements = replacement_index(grouped["replacements"])
    height_replacements = dict(height_replacements)

    counts = {layer: Counter() for layer in ("repaired", "grouped", "height")}
    uppers = {}
    effective = Counter()

    @lru_cache(maxsize=16384)
    def expand(metadata, short):
        data = json.loads(metadata)
        kind = data["type"]
        if kind in ("tailgeo", "tail", "residue"):
            parts = branch_rows(checker, data, short=short)
        else:
            data["children"] = []
            parts = sparse_branches(checker, data) if kind == "sparse" else checker.branches(
                data, repaired["mode"],
            )
        return tuple(sparse_rows(added) for added in parts)

    def replacement(replacements, path, bound, label):
        information = replacements.pop(path, None)
        require(information is not None and information[0] == bound,
                f"missing {label} replacement or mismatched leaf bound at {path}")
        effective[label] += 1
        return information[1]

    def walk(node, rows, path, layer, short=None):
        kind = node["type"]
        counts[layer][kind] += 1
        if kind == "leaf":
            limit = {"repaired": Fraction(7, 12), "grouped": grouped_limit,
                     "height": THRESHOLD}[layer]
            bound = check_leaf(checker, node, rows, limit)
            if layer not in uppers or bound > uppers[layer]:
                uppers[layer] = bound
            if layer == "repaired" and bound > grouped_limit:
                child = replacement(grouped_replacements, path, bound, "grouped")
                walk(child, rows, path, "grouped")
            elif layer != "height" and bound > THRESHOLD:
                child = replacement(height_replacements, path, bound, "height")
                require(child.get("type") in ("height", "leaf"), "missing height disjunction")
                walk(child, rows, path, "height")
            else:
                effective["leaf"] += 1
                if "effective" not in uppers or bound > uppers["effective"]:
                    uppers["effective"] = bound
            return

        children = node.get("children")
        require(isinstance(children, list), f"missing children at {path}")
        if kind == "height":
            require(layer == "height" and short is None and set(node) == {"type", "children"}
                    and len(children) == len(LETTERS), "invalid height disjunction")
            # At least two normalized terms have exponent one. Each choice of the
            # possible short letter is an exhaustive case, with no added affine row.
            for letter, child in zip(LETTERS, children, strict=True):
                walk(child, rows, path, layer, letter)
            return

        allowed = {"repaired": {"geo", "ss", "pairdev"}, "grouped": {"geo", "sparse"},
                   "height": {"tailgeo", "tail", "residue"}}
        require(kind in allowed[layer] and (layer != "height" or short in LETTERS),
                f"invalid {layer} node type {kind}")
        metadata = json.dumps({key: value for key, value in node.items() if key != "children"},
                              sort_keys=True, separators=(",", ":"))
        branches = expand(metadata, short)
        require(len(children) == len(branches), f"incomplete {layer} branch")
        for i, (child, added) in enumerate(zip(children, branches, strict=True)):
            walk(child, rows + added, path + (i,) if layer != "height" else path, layer, short)

    walk(repaired["tree"], sparse_rows(checker.base(BUDGET)), (), "repaired")
    require(counts["repaired"]["leaf"] == repaired["leaves"] == 84177
            and counts["repaired"]["pairdev"] == 1 and uppers["repaired"] == Fraction(7, 12),
            "unexpected repaired FullOC replay result")
    require(not grouped_replacements and not height_replacements, "unused replacement paths")
    require(uppers["effective"] <= THRESHOLD, "height-delta upper bound exceeds its threshold")
    return {"budget": str(BUDGET), "threshold": str(THRESHOLD), "upper": str(uppers["effective"]),
            "replaced_grouped_leaves": effective["grouped"],
            "replaced_height_leaves": effective["height"], "checked_total_leaves": effective["leaf"],
            **{layer: dict(totals) for layer, totals in counts.items()}}


def main():
    require(len(sys.argv) <= 2, "expected at most one height-delta archive")
    path = Path(sys.argv[1]) if len(sys.argv) == 2 else DELTA_FILE
    print(json.dumps(verify(read_delta(path)), sort_keys=True))


if __name__ == "__main__":
    main()
