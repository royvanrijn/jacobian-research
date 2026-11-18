#!/usr/bin/env python3
"""Exact finite A4 action controls; not a replay of Domrina's geometric proof."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import platform


HERE = Path(__file__).resolve().parent
PDF = HERE / "sources/domrina-four-sheeted-general-2000.pdf"
PDF_SHA256 = "0be24c5c6a9cc9423c3049c3402f2c98e44bac5f1c6018b3999989e0f32e5018"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    result = [0] * len(a)
    for i, v in enumerate(a):
        result[v] = i
    return tuple(result)


def alternating_group(n):
    return frozenset(
        p for p in permutations(range(n))
        if sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n)) % 2 == 0
    )


def closure(generators, degree):
    identity = tuple(range(degree))
    generators = set(generators)
    generators |= {inverse(g) for g in generators}
    result = {identity}
    frontier = [identity]
    while frontier:
        a = frontier.pop()
        for b in generators:
            c = compose(a, b)
            if c not in result:
                result.add(c)
                frontier.append(c)
    return frozenset(result)


def conjugate_subgroup(g, subgroup):
    gi = inverse(g)
    return frozenset(compose(compose(g, h), gi) for h in subgroup)


def subgroup_core(group, subgroup):
    return frozenset.intersection(*(conjugate_subgroup(g, subgroup) for g in group))


def normal_closure(group, subgroup):
    generators = set().union(*(conjugate_subgroup(g, subgroup) for g in group))
    return closure(generators, len(next(iter(group))))


def all_subgroups(group):
    n = len(next(iter(group)))
    identity = tuple(range(n))
    result = {frozenset({identity})}
    frontier = list(result)
    while frontier:
        subgroup = frontier.pop()
        for g in group - subgroup:
            enlarged = closure(set(subgroup) | {g}, n)
            if enlarged not in result:
                result.add(enlarged)
                frontier.append(enlarged)
    return result


def conjugacy_classes_of_subgroups(group, subgroups):
    remaining = set(subgroups)
    result = []
    while remaining:
        h = min(remaining, key=lambda s: (len(s), sorted(s)))
        orbit = {conjugate_subgroup(g, h) for g in group}
        require(orbit <= subgroups, "Subgroup list not closed under conjugation")
        result.append((h, orbit))
        remaining -= orbit
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "plane_a4_exclusion_verification.json")
    args = parser.parse_args()

    pdf_bytes = PDF.read_bytes()
    require(pdf_bytes.startswith(b"%PDF"), "Primary source is not a PDF")
    require(sha256(pdf_bytes).hexdigest() == PDF_SHA256, "Primary PDF pin changed")

    a4 = alternating_group(4)
    identity4 = tuple(range(4))
    require(len(a4) == 12, "Wrong literal A4 order")
    subgroups = all_subgroups(a4)
    counts = Counter(len(h) for h in subgroups)
    require(dict(counts) == {1: 1, 2: 3, 3: 4, 4: 1, 12: 1}, "A4 subgroup classification changed")

    v4 = next(h for h in subgroups if len(h) == 4)
    require(subgroup_core(a4, v4) == v4, "V4 is not normal")
    require(all(len(h) != 6 for h in subgroups), "Unexpected order-six subgroup")

    core_free = {h for h in subgroups if subgroup_core(a4, h) == frozenset({identity4})}
    require(len(core_free) == 8, "Unexpected number of core-free subgroups")
    actions = []
    for h, orbit in conjugacy_classes_of_subgroups(a4, core_free):
        normal = normal_closure(a4, h)
        degree = len(a4) // len(h)
        expected_normal_order = {1: 1, 2: 4, 3: 12}[len(h)]
        require(len(normal) == expected_normal_order, "Wrong stabilizer normal closure")
        if len(h) == 2:
            require(normal == v4, "C2 normal closure is not V4")
        actions.append({
            "stabilizer_order": len(h),
            "conjugate_stabilizers": len(orbit),
            "action_degree": degree,
            "core_order": len(subgroup_core(a4, h)),
            "normal_closure_order": len(normal),
            "passes_KMFI1": normal == a4,
        })
    require(sorted(a["action_degree"] for a in actions) == [4, 6, 12], "Wrong faithful transitive action degrees")
    require([a["action_degree"] for a in actions if a["passes_KMFI1"]] == [4], "Wrong KMFI1 survivor")

    a4_stabilizer = frozenset(g for g in a4 if g[0] == 0)
    require(len(a4_stabilizer) == 3 and normal_closure(a4, a4_stabilizer) == a4,
            "Natural A4 action must pass elementary normal generation")

    a5 = alternating_group(5)
    a5_stabilizer = frozenset(g for g in a5 if g[0] == 0)
    dt = (1, 0, 3, 2, 4)
    dt_class = {compose(compose(g, dt), inverse(g)) for g in a5}
    require(len(a5) == 60 and len(a5_stabilizer) == 12, "Wrong natural A5 group/stabilizer orders")
    require(normal_closure(a5, a5_stabilizer) == a5, "Natural A5 stabilizer must normally generate")
    require(len(dt_class) == 15 and closure(dt_class, 5) == a5, "A5 double transpositions do not generate")
    require(all(sum(g[i] == i for i in range(5)) == 1 for g in dt_class), "Wrong A5 fixed-sheet count")

    result = {
        "schema_version": 1,
        "result": "PASS",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "primary_pdf_sha256": PDF_SHA256,
        "A4_subgroup_counts_by_order": dict(sorted(counts.items())),
        "A4_core_free_subgroup_count": len(core_free),
        "A4_faithful_transitive_actions": actions,
        "only_KMFI1_surviving_degree": 4,
        "natural_A5_controls": {
            "order": len(a5), "stabilizer_order": len(a5_stabilizer),
            "normal_closure_order": len(normal_closure(a5, a5_stabilizer)),
            "double_transposition_class_size": len(dt_class),
            "double_transposition_fixed_sheets": 1,
            "passes_necessary_group_conditions": True,
            "plane_realization_or_exclusion_certified_here": False,
        },
        "proof_boundary": (
            "Exact enumeration certifies finite subgroup/core/normal-closure controls and source identity. "
            "The universal plane exclusion combines the written elementary proof, retained KMFI1, "
            "and Domrina's published degree-four theorem. This computation does not replay that theorem. "
            "No degree-five source theorem is imported."
        ),
        "formal_verification": False,
        "independent_replay": False,
        "external_review": False,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"result": "PASS", "subgroups": len(subgroups), "action_degrees": [a["action_degree"] for a in actions],
                      "KMFI1_survivor": 4, "receipt": str(args.output)}))


if __name__ == "__main__":
    main()
