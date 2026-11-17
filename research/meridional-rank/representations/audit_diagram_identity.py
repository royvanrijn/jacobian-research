#!/usr/bin/env python3
"""Exact census transfer and bounded, frozen Montesinos diagram matching.

Requires the versions in diagram_audit_requirements.txt; acceptance uses
integer diagram data only. No hyperbolic numerical decision is used.
"""

import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", message="Plink failed to import tkinter.*")
import spherogram
from spherogram import Link, RationalTangle
import snappy_15_knots

HERE = Path(__file__).resolve().parent


def pd_walk(pd):
    ends = {}
    for c, crossing in enumerate(pd):
        assert len(crossing) == 4
        for h, edge in enumerate(crossing):
            ends.setdefault(edge, []).append((c, h))
    assert all(len(v) == 2 for v in ends.values())
    current = start = (0, 0)
    visits, incoming = [], {}
    for _ in range(2 * len(pd)):
        c, h = current
        visits.append((c, h % 2))  # odd local slots form the overstrand
        incoming.setdefault(c, []).append(h)
        outgoing = (c, (h + 2) % 4)
        pair = ends[pd[outgoing[0]][outgoing[1]]]
        current = pair[1] if pair[0] == outgoing else pair[0]
    assert current == start and len(incoming) == len(pd)
    signs = {}
    for c, slots in incoming.items():
        assert len(slots) == 2 and {h % 2 for h in slots} == {0, 1}
        under = next(h for h in slots if h % 2 == 0)
        over = next(h for h in slots if h % 2 == 1)
        signs[c] = 1 if (under - over) % 4 == 1 else -1
    return [(c, over, signs[c]) for c, over in visits]


def normalize(sequence, full):
    best = None
    best_change = None
    for reverse in [False, True]:
        oriented = list(reversed(sequence)) if reverse else sequence
        for shift in range(len(sequence)):
            shifted = oriented[shift:] + oriented[:shift]
            for over_flip in [False, True]:
                for sign_flip in ([False, True] if full else [False]):
                    ids, result = {}, []
                    for c, over, sign in shifted:
                        if c not in ids:
                            ids[c] = len(ids)
                        value = [ids[c], int(over) ^ int(over_flip)]
                        if full:
                            value.append(-sign if sign_flip else sign)
                        result.append(tuple(value))
                    key = tuple(result)
                    if best is None or key < best:
                        best, best_change = key, {"reverse": reverse, "shift": shift,
                                                  "over_flip": over_flip, "sign_flip": sign_flip}
    return best, best_change


def main():
    versions = {name: metadata.version(name) for name in
                ["snappy", "spherogram", "snappy_15_knots", "snappy_manifolds"]}
    assert versions == {"snappy": "3.3.2", "spherogram": "2.4.1",
                        "snappy_15_knots": "1.2.1", "snappy_manifolds": "1.4"}
    plan_path = HERE / "montesinos_audit_plan.json"
    plan = json.loads(plan_path.read_text())
    sharp_path = HERE / plan["target_source"]
    assert hashlib.sha256(sharp_path.read_bytes()).hexdigest() == plan["target_sha256"]
    sharp = json.loads(sharp_path.read_text())["certificates"]
    assert [c["name"] for c in sharp] == plan["targets"]
    records, target_keys = [], {}
    for certificate in sharp:
        name = certificate["name"]
        link = Link("K" + name)
        pd = [list(c) for c in link.PD_code()]
        walk = pd_walk(pd)
        source_walk = [(abs(c) - 1, int(c > 0), 1) for c in certificate["gauss_code"]]
        source_plain, source_change = normalize(source_walk, False)
        census_plain, census_change = normalize(walk, False)
        assert source_plain == census_plain, f"over/under diagrams differ for {name}"
        key, change = normalize(walk, True)
        record = {"name": name, "lookup": "K" + name, "census_pd": pd,
                  "census_DT": link.DT_code(), "source_plain_transform": source_change,
                  "census_plain_transform": census_change,
                  "census_full_transform": change, "exact_over_under_transfer": True,
                  "family_matches": []}
        records.append(record)
        target_keys.setdefault(key, []).append(record)
    checkpoint = HERE / "diagram_identity_audit.json"
    report = {"status": "IN_PROGRESS", "versions": versions,
              "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [plan_path, sharp_path, Path(__file__)]},
              "census_data": [], "records": records, "tested_candidates": 0,
              "acceptance": "Exact over/under transfer from BKM to independently packaged HT named diagrams. Family acceptance additionally requires exact full signed-Gauss planar-map match up to traversal and global spatial reflections. No numerical or invariant-only acceptance.",
              "claim_boundary": plan["claim_boundary"]}
    package = Path(snappy_15_knots.__file__).parent
    for path in sorted(package.rglob("*.sqlite")):
        report["census_data"].append({"file": str(path.relative_to(package)),
                                     "bytes": path.stat().st_size,
                                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    checkpoint.write_text(json.dumps(report, indent=2) + "\n")
    for index, candidate in enumerate(plan["candidates"]):
        tangle = RationalTangle(candidate["integer_twists"])
        for numerator, denominator in candidate["fractions"]:
            tangle = tangle + RationalTangle(numerator, denominator)
        link = tangle.numerator_closure()
        if len(link.link_components) == 1:
            pd = [list(c) for c in link.PD_code()]
            key, change = normalize(pd_walk(pd), True)
            for record in target_keys.get(key, []):
                record["family_matches"].append({"candidate_index": index, **candidate,
                                                  "pd": pd, "full_transform": change})
        report["tested_candidates"] = index + 1
        checkpoint.write_text(json.dumps(report, indent=2) + "\n")
    report["status"] = "COMPLETE"
    report["exact_census_transfers"] = len(records)
    report["exact_Montesinos_five_tangle_matches"] = sum(bool(r["family_matches"]) for r in records)
    checkpoint.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "census_transfers": len(records),
                      "candidates": report["tested_candidates"],
                      "matches": [{"name": r["name"], "fractions": r["family_matches"][0]["fractions"]}
                                  for r in records if r["family_matches"]]}, indent=2))


if __name__ == "__main__":
    main()
