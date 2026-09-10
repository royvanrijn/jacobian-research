#!/usr/bin/env python3
"""Run three fail-closed Curve302 closure-structure experiments.

Stages:
  1. exhaustive 2^14 closure/hypergraph analysis;
  2. rounded-height quotient near-relation analysis;
  3. exact reconciliation of all fourteen completed seeded V3 trajectories.

No point search is launched.  Existing campaigns are never modified.  ``prepare``
freezes committed inputs plus the exact local terminal/audit files used by stage 3.
A failed/censored/unreceipted stage is preserved as UNKNOWN and is not silently
rerun in place.
"""
from __future__ import annotations

from collections import Counter
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from curve302_closure_structure_core import (
    DIM, FULL, GENERIC_RANK, SCALE, EXPECTED_DIRECTIONS, arrival_thresholds,
    atomic, closure_law_audit, closure_mask, closure_waves, inclusion_minimal,
    mask_names, read, require, sha, trigger_report, validate_landscape,
)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL = ROOT / "artifacts/local/elliptic-curves"
SELF = Path(__file__).resolve()
DEFAULT = LOCAL / "curve302-closure-structure-v1"

LANDSCAPE = ART / "curve302_exceptional_subgroup_landscape_v1.json"
UNLOCK = ART / "curve302_unlock_seed_closure_v2.json"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
M24 = ART / "curve302_recovered_followup_wave_03_mod2_v1.json"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
FILTRATION = ART / "curve302_recovered_quotient_local_filtration_v1.json"
PUBLIC_SPAN = ART / "curve302_recovered_public_span_v1/result.json"
PANEL_EXPORT = ART / "curve302_seed_universality_panel_v1.json"
TWO_SEED_EXPORT = ART / "curve302_seeded_v3_amplifier_v1.json"
LANDSCAPE_SOURCE = CAS / "audit_curve302_exceptional_subgroup_landscape.sage"
GROUP_SOURCE = CAS / "half_lattice_pointed_sieve.py"
GEOMETRY_SOURCE = CAS / "prospective_half_lattice_v3.sage"
CORE_SOURCE = CAS / "curve302_closure_structure_core.py"
SAGE_SOURCE = CAS / "curve302_closure_structure_sage.py"
TWO_SEED_LOCAL = LOCAL / "curve302-seeded-v3-amplifier-v1"
PANEL_LOCAL = LOCAL / "curve302-seed-universality-panel-v1"
STAGES = ("closure-laws", "quotient-relations", "trajectories")


def relative(path):
    return str(Path(path).relative_to(ROOT))


def git_head():
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def sage_launcher():
    for value in (os.environ.get("V3_SAGE"), shutil.which("sage"), str(Path.home()/".local/bin/sage"), "/usr/bin/sage"):
        if value and Path(value).is_file() and os.access(value, os.X_OK):
            return str(Path(value).resolve())
    raise RuntimeError("Sage launcher not found; set V3_SAGE=/absolute/path/to/sage")


def seed_folder(seed):
    return (TWO_SEED_LOCAL if seed in {"recovered-strict-02", "recovered-strict-03"} else PANEL_LOCAL) / seed


def stage_audit_path(folder, stage):
    return folder / "replay-M17" / f"epoch-{int(stage['epoch']):02d}" / stage["audit"]


def raw_seed_bindings(seed):
    folder = seed_folder(seed)
    required = [folder/"seed-input.json", folder/"seeded-verified.json", folder/"replay-M17/terminal.json"]
    missing = [str(path) for path in required if not path.is_file()]
    require(not missing, f"{seed}: missing completed raw evidence: {missing}")
    verified = read(folder/"seeded-verified.json")
    terminal = read(folder/"replay-M17/terminal.json")
    require(verified.get("status") == "PASS_INDEPENDENT_SEEDED_V3_REPLAY", f"{seed}: replay is not passed")
    require(int(verified["initial_rank"]) == 18 and int(terminal["initial_rank"]) == 18, f"{seed}: wrong initial rank")
    paths = list(required)
    for stage in terminal.get("stages", []):
        audit = stage_audit_path(folder, stage)
        require(audit.is_file(), f"{seed}: missing stage-final audit {audit}")
        paths.append(audit)
    return {relative(path): sha(path) for path in paths}


def committed_bindings():
    paths = [LANDSCAPE, UNLOCK, VISIBILITY, M24, PARENT, FILTRATION, PUBLIC_SPAN,
             PANEL_EXPORT, TWO_SEED_EXPORT, LANDSCAPE_SOURCE, GROUP_SOURCE,
             GEOMETRY_SOURCE, CORE_SOURCE, SAGE_SOURCE, SELF]
    missing = [str(path) for path in paths if not path.is_file()]
    require(not missing, f"missing committed inputs: {missing}")
    return {relative(path): sha(path) for path in paths}


def prepare(folder):
    folder = Path(folder).resolve()
    require(not folder.exists(), "experiment folder exists; use resume/status/check or choose a new --folder")
    source = read(LANDSCAPE)
    names, _ = validate_landscape(source)
    require(tuple(names) == EXPECTED_DIRECTIONS, "direction naming/order changed")
    unlock = read(UNLOCK)
    require(unlock.get("status") == "PASS_RETROSPECTIVE_SINGLE_SEED_CLOSURE_CENSUS", "unlock census is not passed")
    require(tuple(unlock["direction_ids"]) == EXPECTED_DIRECTIONS, "unlock direction roster changed")
    raw = {seed: raw_seed_bindings(seed) for seed in EXPECTED_DIRECTIONS}
    sage = sage_launcher()
    software = json.loads(subprocess.check_output([
        sage, "-python", "-c",
        "import json,sage.version; from sage.all import pari; print(json.dumps({'sage':sage.version.version,'pari':str(pari.version())}))",
    ], cwd=ROOT, text=True, timeout=60).strip().splitlines()[-1])
    folder.mkdir(parents=True)
    plan = {
        "schema": "curve302-closure-structure.v1",
        "status": "SEALED_NOT_RUN",
        "source_commit": git_head(),
        "sage": sage,
        "software": software,
        "direction_ids": list(EXPECTED_DIRECTIONS),
        "stages": list(STAGES),
        "limits": {"stage_seconds": 3600, "rss_bytes": 6*1024**3},
        "committed_inputs": committed_bindings(),
        "raw_seed_inputs": raw,
        "boundary": (
            "Retrospective structure analysis only. All fourteen exceptional directions and completed seeded runs are known inputs. "
            "No point search is launched and no rank upper bound or prospective success probability is inferred."
        ),
    }
    atomic(folder/"plan.json", plan)
    atomic(folder/"manifest.json", {"plan_sha256": sha(folder/"plan.json"), "created_unix": time.time()})
    print("CLOSURE_LAB_PREPARED|folder={}|raw_runs=14".format(folder), flush=True)


def guard(folder):
    folder = Path(folder).resolve()
    require((folder/"plan.json").is_file() and (folder/"manifest.json").is_file(), "missing sealed plan")
    plan, manifest = read(folder/"plan.json"), read(folder/"manifest.json")
    require(plan.get("schema") == "curve302-closure-structure.v1", "wrong plan schema")
    require(sha(folder/"plan.json") == manifest["plan_sha256"], "plan changed")
    for name, digest in plan["committed_inputs"].items():
        path = ROOT/name
        require(path.is_file() and sha(path) == digest, "committed input changed: " + name)
    for seed, bindings in plan["raw_seed_inputs"].items():
        for name, digest in bindings.items():
            path = ROOT/name
            require(path.is_file() and sha(path) == digest, f"raw seed evidence changed: {seed}: {name}")
    return plan


def experiment_closure(folder):
    source = read(LANDSCAPE)
    names, states = validate_landscape(source)
    global_threshold = int(source["optimization"]["minimax"]["objective_numerator"])
    historical_threshold = max(int(row["retained_numerator"]) for row in source["historical_comparison"]["attested_tail_rows"])
    existing = read(UNLOCK)
    expected_minimax = {row["seed_direction"]: int(row["minimax_bottleneck_numerator"])
                        for row in existing["single_seed_ranking"]}
    thresholds = {}
    for label, threshold in (("global_minimax", global_threshold), ("historical_tail", historical_threshold)):
        closures = [closure_mask(states, mask, threshold) for mask in range(1 << DIM)]
        laws = closure_law_audit(closures)
        full_generators = inclusion_minimal([mask for mask, closed in enumerate(closures) if closed == FULL])
        singleton = {}
        for i, name in enumerate(names):
            final, waves = closure_waves(states, 1 << i, threshold)
            singleton[name] = {
                "final_rank": GENERIC_RANK + final.bit_count(), "full": final == FULL,
                "waves": [[names[j] for j in wave] for wave in waves],
            }
        thresholds[label] = {
            "threshold_numerator": threshold,
            "threshold_scaled_height": threshold/(4*SCALE),
            "laws": laws,
            "singleton_closures": singleton,
            "full_singleton_count": sum(row["full"] for row in singleton.values()),
            "minimal_full_generator_count": len(full_generators),
            "minimal_full_generator_size_histogram": dict(sorted(Counter(mask.bit_count() for mask in full_generators).items())),
            "minimal_full_generators": [mask_names(mask, names) for mask in full_generators],
            "minimal_triggers": trigger_report(states, names, threshold),
        }

    critical = []
    for i, name in enumerate(names):
        arrivals, full_threshold = arrival_thresholds(states, 1 << i)
        require(full_threshold == expected_minimax[name], f"{name}: minimax differs from existing closure census")
        critical.append({
            "seed": name, "full_closure_threshold": full_threshold,
            "arrival_thresholds": {names[j]: int(arrivals[j]) for j in range(DIM)},
        })

    base = states[0]["retained_numerators"]
    pairwise = []
    for i, name in enumerate(names):
        child = states[1 << i]["retained_numerators"]
        targets = []
        for j, target in enumerate(names):
            if i == j:
                targets.append(None)
                continue
            before, after = int(base[j]), int(child[j])
            targets.append({"target": target, "before": before, "after": after, "drop": before-after,
                            "ratio": None if before == 0 else after/before})
        pairwise.append({"added": name, "targets": targets})

    result = {
        "schema": "curve302-closure-laws.v1",
        "status": "PASS_EXHAUSTIVE_16384_CLOSURE_ANALYSIS",
        "direction_ids": names,
        "thresholds": thresholds,
        "critical_singleton_arrivals": critical,
        "pairwise_immediate_unlock": pairwise,
        "input_sha256": sha(LANDSCAPE),
        "existing_single_seed_census_sha256": sha(UNLOCK),
        "boundary": "Exact combinatorics of retained finite-atlas costs; all directions are retrospective known M31/M17 coordinates.",
    }
    atomic(Path(folder)/"closure-laws.json", result)
    print("CLOSURE_LAB_EXP1|states=16384|full_singletons_global={}|full_singletons_tail={}".format(
        thresholds["global_minimax"]["full_singleton_count"], thresholds["historical_tail"]["full_singleton_count"]), flush=True)


def phase_output(folder, stage):
    return Path(folder)/{
        "closure-laws": "closure-laws.json",
        "quotient-relations": "quotient-relations.json",
        "trajectories": "trajectories.json",
    }[stage]


def run_phase(folder, stage):
    from research_runtime.supervisor import Limits, run as supervise
    plan = guard(folder)
    phase = Path(folder)/"phases"/stage
    receipt, seal = phase/"supervisor.json", phase/"seal.json"
    if receipt.exists():
        row = read(receipt)
        require(row["outcome"] == "completed" and row["returncode"] == 0,
                "previous phase failed/censored; preserve it and use a fresh folder")
        require(seal.is_file(), "completed phase has no seal")
        bound = read(seal)
        require(bound["receipt_sha256"] == sha(receipt), "phase receipt changed")
        output = phase_output(folder, stage)
        require(output.is_file() and bound["output_sha256"] == sha(output), "phase output changed")
        return
    require(not (phase/"started.json").exists(), "unreceipted interrupted phase; preserve it and use a fresh folder")
    atomic(phase/"started.json", {"stage": stage, "time": time.time()})
    command = [plan["sage"], "-python", str(SELF), "_phase", "--folder", str(Path(folder).resolve()), "--stage", stage]
    result = supervise(
        command,
        limits=Limits(wall_seconds=plan["limits"]["stage_seconds"], rss_bytes=plan["limits"]["rss_bytes"]),
        log_path=phase/"worker.log", checkpoint_path=receipt, cwd=ROOT,
        env={**os.environ, "PYTHONUNBUFFERED":"1", "OPENBLAS_NUM_THREADS":"1", "OMP_NUM_THREADS":"1", "MKL_NUM_THREADS":"1"},
    )
    if result["outcome"] == "completed" and result["returncode"] == 0:
        output = phase_output(folder, stage)
        require(output.is_file(), "phase returned success without output")
        atomic(seal, {"receipt_sha256": sha(receipt), "output_sha256": sha(output)})
        return
    atomic(Path(folder)/"REPORT.json", {
        "status": "UNKNOWN_FAILED_OR_CENSORED_STAGE", "stage": stage,
        "outcome": result["outcome"], "returncode": result["returncode"],
    })
    raise RuntimeError(f"{stage} failed/censored; evidence retained")


def execute(folder):
    plan = guard(folder)
    for stage in STAGES:
        run_phase(folder, stage)
    outputs = {stage: sha(phase_output(folder, stage)) for stage in STAGES}
    report = {
        "status": "PASS_THREE_CLOSURE_EXPERIMENTS",
        "outputs": outputs,
        "summary": {
            "closure": read(Path(folder)/"closure-laws.json")["status"],
            "relations": read(Path(folder)/"quotient-relations.json")["status"],
            "trajectories": read(Path(folder)/"trajectories.json")["status"],
        },
        "boundary": plan["boundary"],
    }
    atomic(Path(folder)/"REPORT.json", report)
    print(json.dumps(report, indent=2, sort_keys=True))


def check(folder):
    plan = guard(folder)
    report = read(Path(folder)/"REPORT.json")
    require(report.get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS", "experiment is not complete")
    for stage in STAGES:
        run_phase(folder, stage)
        require(report["outputs"][stage] == sha(phase_output(folder, stage)), "report/output hash differs")
    closure = read(Path(folder)/"closure-laws.json")
    relations = read(Path(folder)/"quotient-relations.json")
    trajectories = read(Path(folder)/"trajectories.json")
    require(len(closure["critical_singleton_arrivals"]) == 14, "singleton closure census incomplete")
    require(len(relations["schur_quotient"]) == 14, "quotient form wrong dimension")
    require(len(trajectories["runs"]) == 14, "trajectory census incomplete")
    require(all(int(row["final_rank"]) >= 29 for row in trajectories["runs"]), "trajectory endpoint below completed panel")
    print("CLOSURE_LAB_CHECK|status=PASS|runs=14|states=16384", flush=True)


def smoke():
    from curve302_closure_structure_core import binary_rank, minimal_triggers, primitive_vector, rref_binary
    n = 3
    states = []
    for mask in range(1 << n):
        costs = []
        for target in range(n):
            if mask >> target & 1:
                costs.append(None)
            elif target == 0:
                costs.append(5)
            elif target == 1:
                costs.append(4 if mask & 1 else 10)
            else:
                costs.append(3 if (mask & 3) == 3 else 20)
        states.append({"state_mask": mask, "retained_numerators": costs})
    closures = [closure_mask(states, mask, 5, n) for mask in range(1 << n)]
    require(closure_mask(states, 1, 5, n) == 7, "synthetic closure did not cascade")
    require(minimal_triggers(states, 2, 5, n) == [3], "synthetic trigger changed")
    require(closure_law_audit(closures, n)["closure_axioms"]["idempotent"], "synthetic closure law failed")
    require(rref_binary([0b011,0b101,0b110], 3) == (5,6) and binary_rank([3,5,6],3) == 2, "binary RREF changed")
    require(primitive_vector([-2,0,2]) == (1,0,-1), "primitive normalization changed")

    source = read(LANDSCAPE)
    names, real_states = validate_landscape(source)
    require(tuple(names) == EXPECTED_DIRECTIONS, "real direction roster changed")
    unlock = read(UNLOCK)
    require(len(unlock["single_seed_ranking"]) == 14, "real singleton census incomplete")
    for i, name in enumerate(names):
        _, full_threshold = arrival_thresholds(real_states, 1 << i)
        expected = next(row["minimax_bottleneck_numerator"] for row in unlock["single_seed_ranking"] if row["seed_direction"] == name)
        require(full_threshold == int(expected), "real minimax smoke mismatch: " + name)
    print("CLOSURE_LAB_SMOKE|status=PASS|synthetic=PASS|real_singletons=14", flush=True)


def status(folder):
    folder = Path(folder)
    if not (folder/"plan.json").exists():
        print(json.dumps({"status":"NOT_PREPARED"}, indent=2)); return
    plan = read(folder/"plan.json")
    stages = []
    for stage in plan["stages"]:
        phase = folder/"phases"/stage
        if (phase/"seal.json").exists(): state = "SEALED"
        elif (phase/"supervisor.json").exists(): state = "FAILED_OR_CENSORED"
        elif (phase/"started.json").exists(): state = "INTERRUPTED_UNRECEIPTED"
        else: state = "PENDING"
        stages.append({"stage":stage,"status":state})
    effective = read(folder/"REPORT.json")["status"] if (folder/"REPORT.json").exists() else "INCOMPLETE"
    print(json.dumps({"status":effective,"stages":stages}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare","run","resume","status","check","smoke","_phase"))
    parser.add_argument("--folder", type=Path, default=DEFAULT)
    parser.add_argument("--stage", choices=STAGES)
    args = parser.parse_args()
    if args.command == "smoke": smoke(); return
    if args.command == "status": status(args.folder); return
    if args.command == "prepare": prepare(args.folder); return
    if args.command == "check": check(args.folder); return
    if args.command == "_phase":
        guard(args.folder); require(args.stage is not None, "_phase requires --stage")
        if args.stage == "closure-laws": experiment_closure(args.folder)
        else:
            import curve302_closure_structure_sage as sage_lab
            if args.stage == "quotient-relations": sage_lab.relation_experiment(args.folder)
            else: sage_lab.trajectory_experiment(args.folder)
        return
    if args.command == "run":
        if args.folder.exists(): raise RuntimeError("folder exists; use resume")
        prepare(args.folder)
    execute(args.folder)


if __name__ == "__main__":
    main()
