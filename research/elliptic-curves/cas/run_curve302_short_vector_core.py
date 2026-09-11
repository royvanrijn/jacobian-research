#!/usr/bin/env python3
"""Fail-closed static geometry experiments for the Curve302 exceptional quotient.

Three stages, no point searches:
  1. completely enumerate primitive quotient directions through the largest
     static quotient norm actually acquired in the fourteen sealed V3 runs;
  2. build the intrinsic short-vector saturation filtration and compare it with
     the exact common integral cores reconstructed from all trajectory prefixes;
  3. rank all 180 acquisitions in the complete vocabulary and measure one-step
     short-vector basin multiplicity for every observed common-core landmark.

The source is the already sealed ``curve302-closure-structure`` result folder.
All source files and policy limits are hash-frozen at prepare time.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from curve302_short_vector_core import (
    enumerate_primitive_directions, fraction, in_rational_span,
    intersection_saturated, lattice_contains, primitive, qnorm, rank_interval,
    rational_rank, saturation_basis, saturation_index,
    RationalBasis, IntegerVocabulary, target_rank_intervals,
)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT / "artifacts/local/elliptic-curves"
SELF = Path(__file__).resolve()
CORE = CAS / "curve302_short_vector_core.py"
DEFAULT = LOCAL / "curve302-short-vector-core-v1"
STAGES = ("enumeration", "filtration", "ranks-basins")
NAMES = tuple([f"recovered-local-{i:02d}" for i in range(1, 5)] +
              [f"recovered-strict-{i:02d}" for i in range(1, 4)] +
              [f"residual-strict-{i:02d}" for i in range(1, 8)])


def require(test, message):
    if not test:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def atomic(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(obj, stream, indent=2, sort_keys=True, allow_nan=False); stream.write("\n")
            stream.flush(); os.fsync(stream.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def sha(path):
    h = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def ftext(value: F) -> str:
    value = fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def parse_form(payload):
    q = tuple(tuple(fraction(v) for v in row) for row in payload["schur_quotient"])
    require(len(q) == 14 and all(len(row) == 14 for row in q), "wrong quotient dimension")
    # enumeration performs the exact positive-definite LDL check
    return q


def find_source(explicit=None):
    if explicit:
        folder = Path(explicit).resolve()
        require((folder / "REPORT.json").is_file(), "explicit source has no REPORT.json")
        return folder
    exact = LOCAL / "curve302-closure-structure-v1"
    if (exact / "REPORT.json").is_file() and read(exact / "REPORT.json").get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS":
        return exact.resolve()
    candidates = []
    for folder in LOCAL.glob("curve302-closure-structure-v1*"):
        if (folder / "REPORT.json").is_file() and read(folder / "REPORT.json").get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS":
            candidates.append(folder.resolve())
    require(len(candidates) == 1, f"expected one passed closure-structure source, found {len(candidates)}; pass --source")
    return candidates[0]


def source_data(folder):
    folder = Path(folder)
    report = read(folder / "REPORT.json")
    require(report.get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS", "source experiment is not passed")
    for stage in ("quotient-relations", "trajectories"):
        path = folder / f"{stage}.json"
        require(path.is_file(), f"source missing {path.name}")
        require(report["outputs"][stage] == sha(path), f"source hash mismatch: {stage}")
    rel = read(folder / "quotient-relations.json")
    traj = read(folder / "trajectories.json")
    require(rel.get("status") == "PASS_QUOTIENT_RELATION_ANALYSIS", "quotient source not passed")
    require(traj.get("status") == "PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED", "trajectory source not passed")
    require(tuple(rel["direction_ids"]) == NAMES == tuple(traj["direction_ids"]), "direction roster changed")
    form = parse_form(rel)
    runs = []
    total = 0
    finals = []
    for raw in traj["runs"]:
        seed = raw["seed"]
        require(seed in NAMES, "unknown seed")
        events = []
        for stage in raw["stages"]:
            for item in stage["new"]:
                require(item["integral"] is True and int(item["denominator"]) == 1, "nonintegral acquisition")
                word = tuple(int(v) for v in item["quotient_word"])
                p = primitive(word)
                require(tuple(item["primitive_quotient_word"]) == p, "primitive acquisition changed")
                events.append({"word": word, "primitive": p, "epoch": int(stage["epoch"])})
        qfinal = int(raw["final_rank"]) - 17
        require(qfinal == 1 + len(events), "trajectory final dimension mismatch")
        runs.append({"seed": seed, "seed_index": NAMES.index(seed), "events": events, "final_dimension": qfinal})
        total += len(events); finals.append(qfinal)
    require(len(runs) == 14 and total == 180, f"expected fourteen runs / 180 acquisitions, got {len(runs)} / {total}")
    require(sorted(finals) == [12] + [14] * 13, "terminal dimensions no longer match 13x31 plus rank29 exception")
    return {"form": form, "runs": sorted(runs, key=lambda r: r["seed_index"]), "source": folder}


def reusable_enumeration(donor, source_hashes, data, max_directions, max_nodes):
    """Import only an already sealed complete enumeration of identical inputs.

    A later failed stage in the donor is preserved and does not invalidate its
    completed enumeration. No donor file or receipt is edited or relabelled.
    """
    donor = Path(donor).resolve()
    old = read(donor / "plan.json")
    require(sha(donor / "plan.json") == read(donor / "manifest.json")["plan_sha256"], "donor plan changed")
    require(old["schema"] == "curve302-short-vector-core-plan.v1", "wrong donor schema")
    require(old["source_hashes"] == source_hashes, "donor source differs")
    seal = read(donor / "phases/enumeration/seal.json")
    meta = read(donor / "enumeration.json")
    require(meta["status"] == "PASS_COMPLETE_EXACT_ENUMERATION", "donor enumeration is incomplete")
    require(seal["output_sha256"] == sha(donor / "enumeration.json"), "donor enumeration metadata changed")
    require(seal["enumeration_sha256"] == meta["enumeration_sha256"] == sha(enumeration_path(donor)), "donor enumeration TSV changed")
    observed = [e["primitive"] for run in data["runs"] for e in run["events"]]
    bound = max(qnorm(data["form"], v) for v in observed)
    require(fraction(meta["bound"]) == bound and meta["observed_acquisitions"] == len(observed), "donor bound or exposure differs")
    require(meta["distinct_observed_directions"] == len(set(observed)), "donor observed roster differs")
    require(meta["direction_count"] <= max_directions and meta["nodes"] <= max_nodes, "donor exceeds current limits")
    return {"folder": str(donor), "producer_plan_sha256": sha(donor / "plan.json"),
            "producer_code_hashes": old["code_hashes"],
            "files_sha256": {name: sha(donor/name) for name in ("enumeration.json", "primitive-directions.tsv")}}


def prepare(folder, source, max_directions, max_nodes, stage_seconds, reuse_enumeration=None):
    folder = Path(folder).resolve(); source = find_source(source)
    require(not folder.exists(), "output folder exists; use resume/check or a new --folder")
    data = source_data(source)
    source_files = [source / "REPORT.json", source / "quotient-relations.json", source / "trajectories.json"]
    source_hashes = {p.name: sha(p) for p in source_files}
    reused = None if reuse_enumeration is None else reusable_enumeration(
        reuse_enumeration, source_hashes, data, max_directions, max_nodes)
    folder.mkdir(parents=True)
    plan = {
        "schema": "curve302-short-vector-core-plan.v1", "status": "SEALED_NOT_RUN",
        "source": str(source), "source_hashes": source_hashes,
        "code_hashes": {SELF.name: sha(SELF), CORE.name: sha(CORE)},
        "policy": {"max_directions": int(max_directions), "max_nodes": int(max_nodes), "stage_seconds": int(stage_seconds)},
        "expected": {"directions": 14, "runs": 14, "acquisitions": 180},
        "reused_enumeration": reused,
        "boundary": "Retrospective static quotient-lattice analysis only; no point search, rank claim, prospective selector or propagation theorem.",
    }
    atomic(folder / "plan.json", plan)
    atomic(folder / "manifest.json", {"plan_sha256": sha(folder / "plan.json"), "created_unix": time.time()})
    if reused is not None:
        for name, digest in reused["files_sha256"].items():
            original = Path(reused["folder"])/name
            shutil.copyfile(original, folder/name)
            require(sha(original) == sha(folder/name) == digest, "donor changed during copy")
        phase = folder/"phases/enumeration"
        atomic(phase/"import.json", reused)
        atomic(phase/"seal.json", {"output_sha256": sha(folder/"enumeration.json"),
                                  "enumeration_sha256": sha(enumeration_path(folder)),
                                  "import_sha256": sha(phase/"import.json")})
        print("SHORT_CORE_REUSED_ENUMERATION|complete=PASS|donor="+reused["folder"], flush=True)
    print(f"SHORT_CORE_PREPARED|source={source}|folder={folder}|events=180", flush=True)


def guard(folder):
    folder = Path(folder).resolve()
    require((folder / "plan.json").is_file() and (folder / "manifest.json").is_file(), "missing sealed plan")
    plan = read(folder / "plan.json"); manifest = read(folder / "manifest.json")
    require(sha(folder / "plan.json") == manifest["plan_sha256"], "plan changed")
    require(plan.get("schema") == "curve302-short-vector-core-plan.v1", "wrong plan schema")
    source = Path(plan["source"])
    for name, digest in plan["source_hashes"].items():
        require((source / name).is_file() and sha(source / name) == digest, f"source changed: {name}")
    for path, digest in ((SELF, plan["code_hashes"][SELF.name]), (CORE, plan["code_hashes"][CORE.name])):
        require(sha(path) == digest, f"code changed after prepare: {path.name}")
    return plan, source_data(source)


def enumeration_path(folder): return Path(folder) / "primitive-directions.tsv"


def read_enumeration(folder):
    rows = []
    with enumeration_path(folder).open() as stream:
        header = next(stream).rstrip("\n")
        require(header == "norm\tvector", "wrong enumeration header")
        for line in stream:
            norm, vector = line.rstrip("\n").split("\t")
            rows.append((fraction(norm), tuple(map(int, vector.split(",")))))
    require(all(rows[i-1] < rows[i] for i in range(1, len(rows))), "enumeration order changed or duplicate row")
    return rows


def stage_enumeration(folder, data, policy):
    observed = [e["primitive"] for run in data["runs"] for e in run["events"]]
    observed_norms = [qnorm(data["form"], v) for v in observed]
    bound = max(observed_norms)
    rows, stats = enumerate_primitive_directions(data["form"], bound,
        max_directions=policy["max_directions"], max_nodes=policy["max_nodes"])
    present = {v for _, v in rows}
    require(all(v in present for v in observed), "complete bound failed to contain an observed primitive acquisition")
    out = enumeration_path(folder)
    tmp = out.with_suffix(".tmp")
    with tmp.open("w") as stream:
        stream.write("norm\tvector\n")
        for norm, vector in rows:
            stream.write(f"{ftext(norm)}\t{','.join(map(str, vector))}\n")
    os.replace(tmp, out)
    result = {
        "schema": "curve302-short-vector-enumeration.v1", "status": "PASS_COMPLETE_EXACT_ENUMERATION",
        "bound": ftext(bound), "direction_count": len(rows), "distinct_observed_directions": len(set(observed)),
        "observed_acquisitions": len(observed), "nodes": stats.nodes, "leaves": stats.leaves,
        "enumeration_sha256": sha(out),
        "definition": "All primitive integer directions modulo sign with v^T Q v <= max static Q-norm among the 180 acquired primitive directions.",
    }
    atomic(Path(folder) / "enumeration.json", result)
    print(f"SHORT_CORE_ENUM|directions={len(rows)}|bound={ftext(bound)}|nodes={stats.nodes}", flush=True)


def prefix_lattices(data):
    prefixes = defaultdict(list)
    per_run = {}
    count = 0
    n = 14
    for run in data["runs"]:
        rows = [tuple(int(i == run["seed_index"]) for i in range(n))]
        states = {1: saturation_basis(rows, n)}
        require(saturation_index(rows, n) == 1, "seed axis not saturated")
        count += 1; prefixes[1].append(states[1])
        for event in run["events"]:
            rows.append(event["word"])
            require(rational_rank(rows, n) == len(rows), "prefix lost rational independence")
            require(saturation_index(rows, n) == 1, "trajectory prefix saturation index is not one")
            sat = saturation_basis(rows, n)
            d = len(rows); states[d] = sat; prefixes[d].append(sat); count += 1
        per_run[run["seed"]] = states
    require(count == 194, f"expected 194 prefix lattices, got {count}")
    return prefixes, per_run


def stage_filtration(folder, data, policy):
    enum = read_enumeration(folder); n = 14
    echelon = RationalBasis(n)
    filtration = []
    i = 0
    while i < len(enum):
        norm = enum[i][0]; j = i
        before = len(echelon.independent)
        while j < len(enum) and enum[j][0] == norm:
            v = enum[j][1]
            echelon.add(v)
            j += 1
        after = len(echelon.independent)
        if after > before:
            sat = saturation_basis(echelon.independent, n)
            filtration.append({"norm": ftext(norm), "rank_before": before, "rank_after": after,
                               "basis": [list(v) for v in sat], "shell_first_rank": i + 1, "shell_last_rank": j})
        i = j
        # The final shell is consumed in full (including ties). Once rank=n,
        # saturation is Z^n at every later shell, so no further change exists.
        if after == n:
            break
    require(filtration and filtration[-1]["rank_after"] == 14, "enumerated short vectors do not span the full quotient")

    prefixes, per_run = prefix_lattices(data)
    common = []
    for d in sorted(prefixes):
        core = intersection_saturated(prefixes[d], n)
        common.append({"quotient_dimension": d, "run_count": len(prefixes[d]), "rank": len(core),
                       "basis": [list(v) for v in core]})
    landmarks = []
    prior = -1
    for row in common:
        if row["rank"] > prior:
            core = tuple(tuple(v) for v in row["basis"])
            first_containing = next((f for f in filtration if lattice_contains(tuple(tuple(v) for v in f["basis"]), core, n)), None)
            equality = next((f for f in filtration if tuple(tuple(v) for v in f["basis"]) == core), None)
            landmarks.append({"quotient_dimension": row["quotient_dimension"], "core_rank": row["rank"],
                              "core_basis": row["basis"], "intrinsic_first_containing_norm": None if first_containing is None else first_containing["norm"],
                              "intrinsic_exact_equality_norm": None if equality is None else equality["norm"]})
            prior = row["rank"]
    result = {
        "schema": "curve302-short-vector-filtration.v1", "status": "PASS_INTRINSIC_AND_OBSERVED_FILTRATIONS",
        "prefix_lattices_checked": 194, "intrinsic_filtration": filtration,
        "observed_common_integral_cores": common, "common_core_landmarks": landmarks,
        "boundary": "Intrinsic filtration depends only on the rounded Schur quotient form and the complete primitive-vector bound; observed cores are recomputed independently from exact trajectory words.",
    }
    atomic(Path(folder) / "filtration.json", result)
    print("SHORT_CORE_FILTRATION|intrinsic_steps={}|common_landmarks={}".format(len(filtration), len(landmarks)), flush=True)


def dot_zero(functionals, vector):
    return all(sum(a * b for a, b in zip(row, vector)) == 0 for row in functionals)


def rational_annihilator(rows, n):
    from sympy import Matrix
    A = Matrix(rows) if rows else Matrix.zeros(0, n)
    import math
    answer = []
    for v in A.nullspace():
        fractions = [fraction(x) for x in v]
        denominator = math.lcm(*(x.denominator for x in fractions))
        answer.append(primitive(tuple(int(x*denominator) for x in fractions)))
    return tuple(answer)


def stage_ranks_basins(folder, data, policy):
    enum = read_enumeration(folder); filtration = read(Path(folder) / "filtration.json"); n = 14
    rank_rows = []
    frequency = Counter(e["primitive"] for run in data["runs"] for e in run["events"])
    intervals = target_rank_intervals(enum, frequency)
    for run in data["runs"]:
        for step, event in enumerate(run["events"]):
            lo, hi, norm = intervals[event["primitive"]]
            rank_rows.append({"seed": run["seed"], "step": step, "primitive": list(event["primitive"]),
                              "norm": ftext(norm), "rank_best": lo, "rank_worst": hi, "multiplicity_across_180": frequency[event["primitive"]]})
    require(len(rank_rows) == 180, "acquisition rank census incomplete")
    worsts = sorted(r["rank_worst"] for r in rank_rows)
    rank_summary = {
        "events": 180, "coverage": 180, "vocabulary": len(enum),
        "top10": sum(r["rank_worst"] <= 10 for r in rank_rows),
        "top100": sum(r["rank_worst"] <= 100 for r in rank_rows),
        "top1000": sum(r["rank_worst"] <= 1000 for r in rank_rows),
        "median_rank_worst": worsts[len(worsts)//2], "maximum_rank_worst": max(worsts),
    }

    common_by_dim = {int(row["quotient_dimension"]): tuple(tuple(v) for v in row["basis"])
                     for row in filtration["observed_common_integral_cores"]}
    landmarks = [int(row["quotient_dimension"]) for row in filtration["common_core_landmarks"] if int(row["quotient_dimension"]) >= 2]
    prefixes, per_run = prefix_lattices(data)
    vocabulary = IntegerVocabulary(enum)
    basins = []
    for d in landmarks:
        core = common_by_dim[d]
        for run in data["runs"]:
            if d > run["final_dimension"]: continue
            prefix = per_run[run["seed"]][d - 1]
            actual = run["events"][d - 2]["primitive"]
            actual_norm = qnorm(data["form"], actual)
            rank_prefix = len(prefix)
            deficit = rational_rank((*prefix, *core), n) - rank_prefix
            require(deficit in (0, 1), "one-step common-core deficit exceeds one")
            ann_prefix = rational_annihilator(prefix, n)
            extended_space = saturation_basis((*prefix, *core), n)
            ann_extended = rational_annihilator(extended_space, n)
            eligible_full, hit_full, eligible_actual, hit_actual, hit_min = vocabulary.basin_counts(
                ann_prefix, ann_extended, actual_norm, already_contained=deficit == 0)
            require(eligible_full > 0 and eligible_actual > 0, "empty basin denominator")
            require(deficit == 0 or dot_zero(ann_extended, actual), "actual acquisition does not enter its observed common core")
            basins.append({
                "landmark_dimension": d, "seed": run["seed"], "core_rank": len(core), "deficit_before": deficit,
                "actual_next_norm": ftext(actual_norm), "eligible_at_actual_norm": eligible_actual,
                "hits_at_actual_norm": hit_actual, "hit_fraction_at_actual_norm": hit_actual / eligible_actual,
                "eligible_at_complete_bound": eligible_full, "hits_at_complete_bound": hit_full,
                "hit_fraction_at_complete_bound": hit_full / eligible_full,
                "minimum_hit_norm": None if hit_min is None else ftext(hit_min),
            })
    by_landmark = []
    for d in landmarks:
        rows = [r for r in basins if r["landmark_dimension"] == d]
        by_landmark.append({"landmark_dimension": d, "runs": len(rows),
                            "mean_hit_fraction_at_actual_norm": sum(r["hit_fraction_at_actual_norm"] for r in rows)/len(rows),
                            "mean_hit_fraction_at_complete_bound": sum(r["hit_fraction_at_complete_bound"] for r in rows)/len(rows),
                            "already_contained_before_count": sum(r["deficit_before"] == 0 for r in rows)})

    repeated = [{"primitive": list(v), "count": c, "rank": intervals[v][:2], "norm": ftext(qnorm(data["form"], v))}
                for v, c in frequency.most_common() if c > 1]
    result = {
        "schema": "curve302-short-vector-ranks-basins.v1", "status": "PASS_COMPLETE_RANK_AND_BASIN_CENSUS",
        "rank_summary": rank_summary, "acquisition_ranks": rank_rows, "repeated_directions": repeated,
        "basin_landmark_summary": by_landmark, "basins": basins,
        "boundary": "Ranks are exact within the complete primitive vocabulary through the maximum observed static quotient norm. Basin fractions are retrospective one-step saturated-span counts, not V3 probabilities.",
    }
    atomic(Path(folder) / "ranks-basins.json", result)
    print("SHORT_CORE_RANKS|coverage=180|vocabulary={}|median_worst={}|basins={}".format(len(enum), rank_summary["median_rank_worst"], len(basins)), flush=True)


def output_for(stage, folder):
    return Path(folder) / {"enumeration":"enumeration.json", "filtration":"filtration.json", "ranks-basins":"ranks-basins.json"}[stage]


def run_stage(folder, stage):
    plan, _ = guard(folder); phase = Path(folder) / "phases" / stage
    seal = phase / "seal.json"; started = phase / "started.json"; log = phase / "worker.log"
    if seal.exists():
        record = read(seal); out = output_for(stage, folder)
        require(out.is_file() and record["output_sha256"] == sha(out), "sealed stage output changed")
        if stage == "enumeration": require(record["enumeration_sha256"] == sha(enumeration_path(folder)), "sealed enumeration TSV changed")
        if "import_sha256" in record:
            require(record["import_sha256"] == sha(phase/"import.json"), "import receipt changed")
            require(read(phase/"import.json") == plan["reused_enumeration"], "import provenance differs from plan")
        return
    require(not started.exists(), "interrupted unsealed stage; preserve folder and use a fresh --folder")
    atomic(started, {"stage":stage, "time":time.time()})
    command = [sys.executable, str(SELF), "_stage", "--folder", str(Path(folder).resolve()), "--stage", stage]
    with log.open("w") as stream:
        try:
            proc = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                  timeout=plan["policy"]["stage_seconds"], env={**os.environ, "PYTHONUNBUFFERED":"1"})
        except subprocess.TimeoutExpired:
            atomic(phase / "failure.json", {"status":"UNKNOWN_TIMEOUT"}); raise RuntimeError(stage + " timed out; evidence retained")
    require(proc.returncode == 0, stage + " failed; inspect " + str(log))
    out = output_for(stage, folder); require(out.is_file(), "stage returned success without output")
    record = {"output_sha256":sha(out)}
    if stage == "enumeration": record["enumeration_sha256"] = sha(enumeration_path(folder))
    atomic(seal, record)


def execute(folder):
    guard(folder)
    for stage in STAGES: run_stage(folder, stage)
    report = {"status":"PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS",
              "outputs":{stage:sha(output_for(stage, folder)) for stage in STAGES},
              "enumeration_sha256":sha(enumeration_path(folder))}
    atomic(Path(folder) / "REPORT.json", report); print(json.dumps(report, indent=2, sort_keys=True))


def check(folder):
    plan, data = guard(folder); report = read(Path(folder) / "REPORT.json")
    require(report.get("status") == "PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS", "run is not complete")
    for stage in STAGES: run_stage(folder, stage); require(report["outputs"][stage] == sha(output_for(stage, folder)), "report hash mismatch")
    require(report["enumeration_sha256"] == sha(enumeration_path(folder)), "enumeration report hash mismatch")
    with tempfile.TemporaryDirectory(prefix="curve302-short-core-check-") as tmp:
        tmp = Path(tmp)
        # Recompute mathematical stages directly; source/code bindings were already checked.
        stage_enumeration(tmp, data, plan["policy"])
        stage_filtration(tmp, data, plan["policy"])
        stage_ranks_basins(tmp, data, plan["policy"])
        for name in ("enumeration.json","primitive-directions.tsv","filtration.json","ranks-basins.json"):
            require(sha(tmp / name) == sha(Path(folder) / name), "deterministic recomputation differs: " + name)
    print("SHORT_CORE_CHECK|status=PASS|deterministic=PASS", flush=True)


def status(folder):
    folder = Path(folder)
    if not (folder / "plan.json").is_file(): print(json.dumps({"status":"NOT_PREPARED"})); return
    rows = []
    for stage in STAGES:
        phase = folder / "phases" / stage
        rows.append({"stage":stage, "status":"SEALED" if (phase/"seal.json").exists() else "INTERRUPTED" if (phase/"started.json").exists() else "PENDING"})
    overall = read(folder / "REPORT.json")["status"] if (folder / "REPORT.json").is_file() else "INCOMPLETE"
    print(json.dumps({"status":overall, "stages":rows}, indent=2))


def smoke():
    from curve302_short_vector_core import enumerate_primitive_directions, saturation_basis, intersection_saturated
    q = ((F(2),F(1),F(0)),(F(1),F(2),F(0)),(F(0),F(0),F(5)))
    rows, stats = enumerate_primitive_directions(q, F(5), max_directions=1000, max_nodes=10000)
    require(rows and all(qnorm(q,v) <= 5 for _,v in rows), "synthetic enumeration failed")
    sat = saturation_basis(((2,0,0),(0,2,0)),3); require(sat == ((1,0,0),(0,1,0)), "synthetic saturation failed")
    inter = intersection_saturated((((1,0,0),(0,1,0)), ((1,0,0),(0,0,1))),3)
    require(inter == ((1,0,0),), "synthetic intersection failed")
    print(f"SHORT_CORE_SMOKE|status=PASS|directions={len(rows)}|nodes={stats.nodes}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("prepare","run","resume","status","check","smoke","_stage"))
    p.add_argument("--source", type=Path)
    p.add_argument("--reuse-enumeration", type=Path, help="new folder only: import a hash-verified complete enumeration of identical source inputs")
    p.add_argument("--folder", type=Path, default=DEFAULT)
    p.add_argument("--stage", choices=STAGES)
    p.add_argument("--max-directions", type=int, default=2_000_000)
    p.add_argument("--max-nodes", type=int, default=100_000_000)
    p.add_argument("--stage-seconds", type=int, default=3600)
    a = p.parse_args()
    if a.command == "smoke": smoke(); return
    if a.command == "status": status(a.folder); return
    if a.command == "prepare": prepare(a.folder, a.source, a.max_directions, a.max_nodes, a.stage_seconds, a.reuse_enumeration); return
    if a.command == "check": check(a.folder); return
    if a.command == "_stage":
        plan, data = guard(a.folder); require(a.stage is not None, "_stage requires --stage")
        if a.stage == "enumeration": stage_enumeration(a.folder, data, plan["policy"])
        elif a.stage == "filtration": stage_filtration(a.folder, data, plan["policy"])
        else: stage_ranks_basins(a.folder, data, plan["policy"])
        return
    if a.command == "run":
        require(not a.folder.exists(), "folder exists; use resume")
        prepare(a.folder, a.source, a.max_directions, a.max_nodes, a.stage_seconds, a.reuse_enumeration)
    execute(a.folder)

if __name__ == "__main__": main()
