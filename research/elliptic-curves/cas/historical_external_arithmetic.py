#!/usr/bin/env python3
"""Separate, retrospective equation-only controls; never mutate the frozen census.

prepare replays existing rank certificates (no point search), freezes the entire
inventory selection roster, and snapshots every prospective file hash. run calls
the unmodified census BASE/LOCAL controller with its frozen runtime policies.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as F
from pathlib import Path
import subprocess

import compact_atlas_specialization as atlas
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
import run_wide_arithmetic_profile as worker
import wide_arithmetic_profile_core as core

ROOT = worker.ROOT
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/wide_arithmetic_historical_external_v1"
INVENTORY = ROOT / "elliptic-curves/data/research_curves/database.json"
INDEX = ROOT / "artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v22.json"
HISTORICAL = ("historical_external", "retrospective_known_high_rank")
PROSPECTIVE = ("prospective_broad_2080", "frozen_prospective")
REQUIRED = {"-2448/11": 27, "2012/211": 27, "2828/2015": 27,
            "4286/1881": 27, "110314/102227": 28}
read, sha, save, require = worker.read, worker.sha, worker.save, worker.require


def labels(kind):
    require(kind in (HISTORICAL, PROSPECTIVE), "unrecognized cohort")
    return dict(zip(("cohort", "selection_mode"), kind))


def tree_hashes(path):
    return {str(p.relative_to(path)): sha(p) for p in sorted(path.rglob("*")) if p.is_file()}


def select(roster):
    selected = [r for r in roster if str(r["family"]) == "11952" and int(r["rank_lower_bound"]) >= 27]
    require(len({r["id"] for r in selected}) == len(selected), "duplicate inventory IDs")
    require(len({str(F(r["parameter"])) for r in selected}) == len(selected), "duplicate historical parameters")
    found = {str(F(r["parameter"])): int(r["rank_lower_bound"]) for r in selected}
    require(all(found.get(t, 0) >= rank for t, rank in REQUIRED.items()), "mandatory historical controls missing")
    return sorted(selected, key=lambda r: F(r["parameter"]))


def verify_admission(current, original, family):
    require(current["id"] == original["id"], "inventory ID mismatch")
    require(current["rank_lower_bound"] == original["rank_lower_bound"], "rank metadata mismatch")
    require(str(current["family"]) == str(original["family"]) == "11952", "family mismatch")
    require(F(current["parameter"]) == F(original["parameter"]), "parameter mismatch")
    model = tuple(map(F, original["curve"]))
    minimal = tuple(map(F, current["ainvs"]))
    specialized, _ = atlas.specialize(family, current["parameter"])
    require(cert.isomorphic(specialized, model) and cert.isomorphic(model, minimal), "specialization/model mismatch")
    points = tuple(tuple(map(F, p)) for p in original["points"])
    witness = original["rank_certificate"]
    replay = checked_rank(model, points, [s["prime"] for s in witness["signatures"]],
                          witness["no_rational_2_torsion_prime"])
    require(core.stable_json(replay) == core.stable_json(witness), "finite-reduction certificate mismatch")
    require(replay["rank_lower_bound"] == current["rank_lower_bound"], "rank replay mismatch")
    u, r, s, t = map(F, current["original_to_minimal_isomorphism"])
    transported = tuple(((x-r)/u**2, (y-s*(x-r)-t)/u**3) for x, y in points)
    require(transported == tuple(tuple(map(F, p)) for p in current["points"]), "point transport mismatch")
    require(all(cert.is_on_weierstrass_curve(minimal, p) for p in transported), "minimal-model point failure")
    return {"id": current["id"], "status": "PASS_EXACT_EXISTING_CERTIFICATE_REPLAY",
            "rank_lower_bound": replay["rank_lower_bound"], "parameter": str(F(current["parameter"])),
            "rank_provenance": current["rank_provenance"],
            "local_search_rank_lower_bound": current["local_search_rank_lower_bound"],
            "finite_reduction_primes": [s["prime"] for s in replay["signatures"]],
            "no_rational_2_torsion_prime": replay["no_rational_2_torsion_prime"],
            "specialization_isomorphism": "PASS", "point_transport": "PASS"}


def prepare(out, census):
    require(not out.exists(), "external output already exists; use run/check to resume")
    require(out != census and census not in out.parents, "external output must be outside census")
    plan, population = worker.check_plan(census)
    report = read(census / "REPORT.json")
    require(len(population) == 2080 and report["full_population_checkpointed"], "census not complete")
    before = tree_hashes(census)
    inventory = read(INVENTORY)
    roster = [{k: r[k] for k in ("id", "family", "parameter", "rank_lower_bound", "rank_provenance")}
              for r in inventory["curves"]]
    selected_ids = {r["id"] for r in select(roster)}
    current = {r["id"]: r for r in inventory["curves"] if r["id"] in selected_ids}
    original = {r["id"]: r for r in read(INDEX)["curves"] if r["id"] in selected_ids}
    family = next(r for r in read(atlas.ATLAS)["families"] if r["family"] == "11952")
    admitted, rows, evidence = [], [], []
    for r in select(roster):
        c, o = current[r["id"]], original[r["id"]]
        admitted.append(verify_admission(c, o, family))
        row = {"id": c["id"], "family": "11952", "t": str(F(c["parameter"])),
               "ainvs": [int(F(a)) for a in c["ainvs"]],
               "final_rank_lower_bound": c["rank_lower_bound"],
               "initial_rank_lower_bound": None, "improved_since_initial": None,
               "rank_provenance": c["rank_provenance"],
               "local_search_rank_lower_bound": c["local_search_rank_lower_bound"],
               **labels(HISTORICAL)}
        require(all(F(a).denominator == 1 for a in c["ainvs"]), "nonintegral inventory model")
        row["curve_key"] = core.curve_key(row)
        rows.append(row)
        evidence.append({"inventory": c, "original": o})
        print(f"HISTORICAL_ADMISSION|id={c['id']}|t={row['t']}|rank_LB={row['final_rank_lower_bound']}|PASS", flush=True)
    require(len({r["curve_key"] for r in rows}) == len(rows), "duplicate historical equations")
    policies = {mode: read(census / f"{mode}-runtime.json") for mode in ("base", "local")}
    for mode, policy in policies.items():
        require(sha(Path(policy["sage"])) == policy["sage_sha256"], "Sage executable changed")
        require(subprocess.check_output([policy["sage"], "--version"], text=True).strip() == policy["sage_version"], "Sage version changed")
    sources = {worker.rel(p): sha(p) for p in (INVENTORY, INDEX, atlas.ATLAS, Path(__file__),
               Path(atlas.__file__), Path(cert.__file__), Path(__import__("memory_rank_certificate").__file__))}
    sources.update(plan["sources"])
    require(tree_hashes(census) == before, "prospective files changed during admission")
    for name in ("inputs", "base", "local", "logs"):
        (out / name).mkdir(parents=True)
    save(out / "inventory_roster.json", {"rows": roster})
    save(out / "rank_evidence.json", {"rows": evidence})
    save(out / "admission.json", {"rows": admitted, "point_searches_launched": 0})
    for row in rows:
        save(out / "inputs" / f"{row['curve_key']}.json", row)
    for mode, policy in policies.items():
        save(out / f"{mode}-runtime.json", policy)
    save(out / "plan.json", {"schema": "elliptic-curves.historical-external-plan.v1",
         "census": str(census), "prospective_files": before, "sources": sources,
         "rows": rows, "policies": policies, "max_jobs": 2, "class_workers": 0,
         "selection": "All inventory family 11952 rows with certified rank lower bound >=27; frozen before arithmetic",
         "evidence_hashes": {name: sha(out / name) for name in ("inventory_roster.json", "rank_evidence.json", "admission.json")}})
    print(f"HISTORICAL_PREPARE|PASS|rows={len(rows)}|prospective_files_preserved={len(before)}", flush=True)


def check(out, complete=False):
    plan = read(out / "plan.json")
    require(tree_hashes(Path(plan["census"])) == plan["prospective_files"], "prospective file set or bytes changed")
    for path, digest in plan["sources"].items():
        require(sha(ROOT / path) == digest, f"bound source changed: {path}")
    for name, digest in plan["evidence_hashes"].items():
        require(sha(out / name) == digest, f"evidence changed: {name}")
    require({r["id"] for r in select(read(out / "inventory_roster.json")["rows"])} == {r["id"] for r in plan["rows"]}, "selection changed")
    for mode, policy in plan["policies"].items():
        require(read(out / f"{mode}-runtime.json") == policy, "runtime policy changed")
        require(sha(Path(policy["sage"])) == policy["sage_sha256"], "Sage executable changed")
        expected = {f"{r['curve_key']}.json" for r in plan["rows"]}
        actual = {p.name for p in (out / mode).glob("*.json")}
        require(actual <= expected and (not complete or actual == expected), f"unexpected or missing {mode} checkpoints")
        for row in plan["rows"]:
            inp = out / "inputs" / f"{row['curve_key']}.json"
            require(inp.read_text() == core.stable_json(row), "historical input changed")
            p = out / mode / inp.name
            if not p.exists():
                continue
            result = read(p)
            require(result["input_sha256"] == sha(inp) and result["curve_key"] == row["curve_key"], "checkpoint binding failed")
            require(result["status"] == "PASS" or result["status"].startswith("UNKNOWN"), "invalid arithmetic status")
            for key in ("timeout_seconds", "memory_gb"):
                require(result["runtime"][key] == policy[key], "worker budget mismatch")
            require(result["runtime"]["mode"] == mode, "worker mode mismatch")
            command = [policy["sage"], "-python", str(worker.WORKER), "--input", str(inp),
                       "--output", str(p.with_suffix(".tmp.json")), "--mode", mode,
                       "--memory-gb", str(policy["memory_gb"])]
            require(result["runtime"]["command"] == command, "worker invocation mismatch")
    require(not (out / "class").exists(), "CLASS lane must remain off")
    return plan


def run(out, jobs):
    plan = check(out)
    require(1 <= jobs <= plan["max_jobs"], "jobs outside frozen cap")
    for mode, policy in plan["policies"].items():
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(worker.run_one, out, row, mode, policy["timeout_seconds"],
                                  policy["memory_gb"], policy["sage"]): row for row in plan["rows"]}
            for future in as_completed(futures):
                print(f"HISTORICAL_{mode.upper()}|t={futures[future]['t']}|status={future.result()}", flush=True)
    check(out, complete=True)
    print("HISTORICAL_RUN|PASS_CHECKPOINTED|prospective_bytes_unchanged=PASS|point_searches=0", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("prepare", "run", "check"))
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--census", type=Path, default=worker.DEFAULT_OUTPUT)
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    out = args.output.resolve()
    if args.command == "prepare":
        prepare(out, args.census.resolve())
    elif args.command == "run":
        run(out, args.jobs)
    else:
        check(out, complete=True)
        print("HISTORICAL_CHECK|PASS|complete=PASS|prospective_bytes_unchanged=PASS")


if __name__ == "__main__":
    main()
