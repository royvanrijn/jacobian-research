#!/usr/bin/env python3
"""Historical Curve302 chart exposure / multiplicity / ordering controls.

No V3 point search is launched.  The controller consumes a completed Curve302
short-vector bundle plus the sealed closure-structure evidence and historical
chart audit/replay records.  It recognizes exposure only from explicit quotient
word evidence; unsupported point-only chart schemas fail closed with a probe.

Commands: prepare | run | resume | status | check | probe | _stage
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time

from curve302_short_core_controls import Vocabulary, integer, primitive, qnorm, require
from curve302_chart_exposure import (
    UnsupportedChartSchema, build_candidate_population, candidate_metrics_for_stage,
    counterfactual_stage, multiplicity_comparison, normalize_explicit_ledger,
    normalize_raw_tree, schema_probe, stage_prefixes, summarize_counterfactuals,
    validate_ledger,
)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
SELF = Path(__file__).resolve()
HELPER = CAS / "curve302_chart_exposure.py"
CONTROL_HELPER = CAS / "curve302_short_core_controls.py"
DEFAULT = ROOT / "artifacts/local/elliptic-curves/curve302-chart-exposure-v1"
NAMES = tuple([f"recovered-local-{i:02d}" for i in range(1,5)] +
              [f"recovered-strict-{i:02d}" for i in range(1,4)] +
              [f"residual-strict-{i:02d}" for i in range(1,8)])
STAGES = ("exposure-census", "multiplicity", "counterfactuals")
SHORT_FILES = ("REPORT.json", "enumeration.json", "filtration.json", "ranks-basins.json",
               "primitive-directions.tsv", "plan.json")
STRUCTURE_FILES = ("REPORT.json", "quotient-relations.json", "trajectories.json")


class EvidenceError(ValueError):
    pass


def read(path):
    path = Path(path)
    def unique(pairs):
        out = {}
        for k, v in pairs:
            require(k not in out, f"duplicate JSON key {k} in {path}")
            out[k] = v
        return out
    def reject(v):
        raise EvidenceError(f"nonfinite JSON value {v} in {path}")
    return json.loads(path.read_text(), object_pairs_hook=unique, parse_constant=reject)


def sha(path):
    h = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as out:
            json.dump(payload, out, sort_keys=True, indent=2, allow_nan=False)
            out.write("\n"); out.flush()
            if os.environ.get("CURVE302_TEST_SKIP_FSYNC") != "1":
                os.fsync(out.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def software():
    import sympy
    return {"python": platform.python_version(), "sympy": sympy.__version__,
            "machine": platform.machine(), "byteorder": sys.byteorder}


def discover_short(explicit=None):
    if explicit:
        folder = Path(explicit).resolve()
        require((folder/"REPORT.json").is_file(), "explicit short-vector source has no REPORT.json")
        return folder
    parents = [ROOT/"artifacts/local/elliptic-curves", ROOT/"artifacts/generated-results/elliptic-curves"]
    candidates = []
    for parent in parents:
        if not parent.exists():
            continue
        for report in parent.rglob("REPORT.json"):
            try:
                if read(report).get("status") == "PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS":
                    candidates.append(report.parent.resolve())
            except Exception:
                pass
    candidates = sorted(set(candidates))
    require(len(candidates) == 1, f"found {len(candidates)} completed short-vector bundles; pass --source")
    return candidates[0]


def discover_structure(short, explicit=None):
    if explicit:
        return Path(explicit).resolve()
    plan = read(Path(short)/"plan.json")
    value = plan.get("source")
    if value and Path(value).is_dir():
        return Path(value).resolve()
    # Compatibility with copied/moved experiments: locate byte-identical source
    # by its three committed output hashes.
    hashes = plan.get("source_hashes", {})
    parents = [ROOT/"artifacts/local/elliptic-curves", ROOT/"artifacts/generated-results/elliptic-curves"]
    candidates = []
    for parent in parents:
        if not parent.exists():
            continue
        for report in parent.rglob("REPORT.json"):
            folder = report.parent
            try:
                if read(report).get("status") != "PASS_THREE_CLOSURE_EXPERIMENTS":
                    continue
                if all((folder/name).is_file() and sha(folder/name) == digest for name, digest in hashes.items()):
                    candidates.append(folder.resolve())
            except Exception:
                continue
    candidates = sorted(set(candidates))
    require(len(candidates) == 1, f"closure-structure source moved and {len(candidates)} byte-identical candidates found; pass --structure")
    return candidates[0]


def source_data(short, structure):
    short, structure = Path(short), Path(structure)
    sr = read(short/"REPORT.json")
    require(sr.get("status") == "PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS", "short-vector source not passed")
    for stage in ("enumeration", "filtration", "ranks-basins"):
        require(sr["outputs"][stage] == sha(short/f"{stage}.json"), f"short-vector output hash mismatch: {stage}")
    enumeration = read(short/"enumeration.json")
    require(enumeration.get("status") == "PASS_COMPLETE_EXACT_ENUMERATION", "short-vector enumeration incomplete")
    require(sha(short/"primitive-directions.tsv") == enumeration["enumeration_sha256"] == sr["enumeration_sha256"],
            "primitive vocabulary hash mismatch")

    cr = read(structure/"REPORT.json")
    require(cr.get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS", "closure structure source not passed")
    for stage in ("quotient-relations", "trajectories"):
        require(cr["outputs"][stage] == sha(structure/f"{stage}.json"), f"closure output hash mismatch: {stage}")
    rel, traj = read(structure/"quotient-relations.json"), read(structure/"trajectories.json")
    require(rel.get("status") == "PASS_QUOTIENT_RELATION_ANALYSIS", "quotient relation source invalid")
    require(traj.get("status") == "PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED", "trajectory source invalid")
    require(tuple(rel["direction_ids"]) == NAMES == tuple(traj["direction_ids"]), "direction roster changed")
    q = tuple(tuple(F(str(x)) for x in row) for row in rel["schur_quotient"])
    require(len(q) == 14 and all(len(r) == 14 for r in q), "wrong quotient form shape")

    runs = []
    total = 0; total_charts = 0
    for raw in traj["runs"]:
        seed = raw["seed"]
        require(seed in NAMES, "unknown trajectory seed")
        events = []
        no_gain_epochs = []
        raw_stages = raw["stages"]
        for expected_epoch, stage in enumerate(raw_stages):
            require(integer(stage["epoch"]) == expected_epoch, "trajectory epochs reordered")
            if not stage["new"]:
                require(integer(stage["after"]) == integer(stage["before"]),
                        f"empty acquisition stage changed rank: {seed}/{expected_epoch}")
                require(expected_epoch == len(raw_stages) - 1,
                        f"nonterminal no-gain stage unsupported: {seed}/{expected_epoch}")
                no_gain_epochs.append(expected_epoch)
                continue
            # Historical seeded closure runs in this experiment acquired one new
            # displayed-D direction per accepted stage. Refuse to invent ordering
            # inside a multi-gain stage.
            require(len(stage["new"]) == 1,
                    f"multi-gain stage unsupported for ordering control: {seed}/{expected_epoch}")
            item = stage["new"][0]
            require(item.get("integral") is True and integer(item["denominator"]) == 1, "nonintegral trajectory acquisition")
            word = tuple(integer(v) for v in item["quotient_word"])
            p = primitive(word)
            require(tuple(item["primitive_quotient_word"]) == p, "trajectory primitive word changed")
            events.append({"word": word, "primitive": p, "epoch": expected_epoch})
        qfinal = integer(raw["final_rank"]) - 17
        require(qfinal == 1 + len(events), "trajectory endpoint dimension mismatch")
        charts = integer(raw["charts"])
        require(charts > 0, "trajectory chart count missing")
        runs.append({"seed": seed, "seed_index": NAMES.index(seed), "events": events,
                     "epochs": list(range(len(raw_stages))),
                     "no_gain_epochs": no_gain_epochs,
                     "final_dimension": qfinal, "charts": charts})
        total += len(events); total_charts += charts
    require(len(runs) == 14 and len({r["seed"] for r in runs}) == 14, "expected fourteen seeded runs")
    require(total == 180, f"expected 180 acquisitions, got {total}")
    require(sorted(r["final_dimension"] for r in runs) == [12] + [14]*13, "endpoint census changed")
    runs.sort(key=lambda r: r["seed_index"])

    filtration = read(short/"filtration.json")
    require(filtration.get("status") == "PASS_INTRINSIC_AND_OBSERVED_FILTRATIONS", "short-vector filtration invalid")
    cores = {}
    rows = filtration["observed_common_integral_cores"]
    require([integer(r["quotient_dimension"]) for r in rows] == list(range(1,15)), "observed core dimensions changed")
    for row in rows:
        d = integer(row["quotient_dimension"])
        basis = tuple(tuple(integer(x) for x in v) for v in row["basis"])
        require(len(basis) == integer(row["rank"]), "common core rank/basis mismatch")
        cores[d] = basis
    return {"form": q, "runs": runs, "cores": cores, "enumeration": enumeration,
            "total_charts": total_charts}


def find_raw_root(structure, explicit=None):
    if explicit:
        root = Path(explicit).resolve(); require(root.is_dir(), "--raw-root is not a directory"); return root
    # The closure-structure controller snapshots all seeded V3 evidence into its
    # own experiment folder; recurse from there rather than guessing original
    # campaign paths.
    return Path(structure).resolve()


def copy_with_hash(src, dst, expected=None):
    src, dst = Path(src), Path(dst)
    before = sha(src)
    if expected is not None:
        require(before == expected, f"source hash changed before copy: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    require(sha(dst) == before == sha(src), f"source changed while copying: {src}")
    return before


def prepare(folder, source=None, structure=None, raw_root=None, ledger_path=None, *, static_limit=1000,
            random_orders=256, stage_seconds=3600, memory_bytes=6*1024**3):
    folder = Path(folder).resolve()
    require(not folder.exists(), "output folder exists; use resume/check or a fresh --folder")
    short = discover_short(source); structure = discover_structure(short, structure)
    data = source_data(short, structure)
    require(1 <= static_limit <= integer(data["enumeration"]["direction_count"]), "invalid static candidate limit")
    require(1 <= random_orders <= 10000, "invalid random-order count")
    require(stage_seconds > 0 and memory_bytes >= 512*1024**2, "invalid stage resource budget")
    expected_epochs = {r["seed"]: set(r.get("epochs", range(len(r["events"])))) for r in data["runs"]}

    folder.mkdir(parents=True, exist_ok=False)
    try:
        if ledger_path:
            ledger_file = Path(ledger_path).resolve(); payload = read(ledger_file)
            ledger = normalize_explicit_ledger(payload, NAMES, 14)
            ledger["source_files"] = [str(ledger_file)]
            ledger["normalization"] = "explicit-ledger"
            raw = ledger_file.parent
        else:
            raw = find_raw_root(structure, raw_root)
            ledger = normalize_raw_tree(raw, NAMES, expected_epochs, 14)
        stats = validate_ledger(ledger, NAMES, data["runs"], expected_total_charts=data["total_charts"])
        # Every actual acquisition must be explicitly exposed in at least one
        # historical chart for that stage. Otherwise this is an incomplete
        # chart-result schema and zeros would be meaningless.
        ledger_runs = {r["seed"]: r for r in ledger["runs"]}
        missing_actual = []
        for run in data["runs"]:
            for event in run["events"]:
                stage = ledger_runs[run["seed"]]["stages"][event["epoch"]]
                words = {tuple(e["word"]) for c in stage["charts"] for e in c["exposures"]}
                if event["primitive"] not in words:
                    missing_actual.append((run["seed"], event["epoch"]))
        require(not missing_actual, f"explicit chart evidence does not expose actual acquisition in {len(missing_actual)} stages; first {missing_actual[:5]}")
    except Exception as exc:
        try:
            raw = find_raw_root(structure, raw_root) if not ledger_path else Path(ledger_path).resolve().parent
            probe = schema_probe(raw, NAMES)
        except Exception as probe_exc:
            probe = {"probe_error": repr(probe_exc)}
        atomic(folder/"SCHEMA_PROBE.json", probe)
        atomic(folder/"PREPARE_FAILED.json", {"status": "UNKNOWN_UNSUPPORTED_OR_AMBIGUOUS_CHART_SCHEMA",
                                                "error": repr(exc),
                                                "boundary": "No exposure zero or chart-order result is inferred from an unsupported schema."})
        raise

    # Snapshot bound high-level inputs.
    input_hashes = {}
    for label, origin, files in (("short", short, SHORT_FILES), ("structure", structure, STRUCTURE_FILES)):
        for name in files:
            src = origin/name; require(src.is_file(), f"missing source file {src}")
            dst = folder/"inputs"/label/name
            input_hashes[str(dst.relative_to(folder))] = copy_with_hash(src, dst)
    # Snapshot only raw files actually selected by the normalization adapter.
    raw_hashes = {}
    raw_sources = []
    for i, text in enumerate(sorted(set(ledger.get("source_files", [])))):
        src = Path(text).resolve(); require(src.is_file(), f"selected raw chart source missing: {src}")
        dst = folder/"inputs/raw"/f"{i:04d}-{src.name}"
        raw_hashes[str(dst.relative_to(folder))] = copy_with_hash(src, dst)
        raw_sources.append({"original": str(src), "snapshot": str(dst.relative_to(folder)), "sha256": raw_hashes[str(dst.relative_to(folder))]})
    # Strip absolute source paths from the canonical frozen ledger; provenance is
    # retained separately in raw_sources.
    def scrub(obj):
        if isinstance(obj, dict):
            return {k: (Path(v).name if k == "source_file" and isinstance(v, str) else scrub(v)) for k, v in obj.items()
                    if k != "source_files"}
        if isinstance(obj, list): return [scrub(v) for v in obj]
        if isinstance(obj, tuple): return [scrub(v) for v in obj]
        return obj
    frozen_ledger = scrub(ledger)
    atomic(folder/"inputs/chart-exposure-ledger.json", frozen_ledger)
    input_hashes["inputs/chart-exposure-ledger.json"] = sha(folder/"inputs/chart-exposure-ledger.json")

    code_hashes = {}
    for src in (SELF, HELPER, CONTROL_HELPER):
        dst = folder/"code"/src.name
        code_hashes[str(dst.relative_to(folder))] = copy_with_hash(src, dst)
    policy = {"static_limit": int(static_limit), "random_orders": int(random_orders),
              "master_seed": "curve302-chart-exposure-v1", "stage_seconds": int(stage_seconds),
              "memory_bytes": int(memory_bytes),
              "candidate_rule": "first_N_static_directions_plus_all_actual_outside",
              "norm_match_rule": "frozen_decade_band_with_complete_static_shell_boundaries",
              "exposure_rule": "explicit_14D_quotient_word_historical_chart_evidence_only",
              "absence_rule": "absence_is_comparable_only_when_chart_complete_flag_is_true"}
    plan = {"schema": "curve302-chart-exposure-plan.v1", "status": "SEALED_NOT_RUN",
            "short_source": str(short), "structure_source": str(structure), "raw_root": str(raw),
            "inputs": input_hashes, "raw_inputs": raw_hashes, "raw_sources": raw_sources,
            "code": code_hashes, "software": software(), "policy": policy, "stages": list(STAGES),
            "expected": {"runs": 14, "acquisitions": 180, "charts": data["total_charts"],
                         "candidate_vocabulary_directions": integer(data["enumeration"]["direction_count"])},
            "ledger_stats": stats,
            "boundary": "Retrospective model checking on one known D/M17 and historical chart transcripts. No point search, prospective rank selector, formal significance test, or propagation theorem."}
    atomic(folder/"plan.json", plan)
    atomic(folder/"manifest.json", {"plan_sha256": sha(folder/"plan.json")})
    print(f"CHART_EXPOSURE_PREPARED|charts={stats['charts']}|exposures={stats['exposures']}|raw_files={len(raw_sources)}", flush=True)


def guard(folder):
    folder = Path(folder).resolve()
    require((folder/"plan.json").is_file() and (folder/"manifest.json").is_file(), "missing sealed plan")
    plan = read(folder/"plan.json")
    require(plan.get("schema") == "curve302-chart-exposure-plan.v1", "wrong plan schema")
    require(tuple(plan["stages"]) == STAGES, "stage list changed")
    require(sha(folder/"plan.json") == read(folder/"manifest.json")["plan_sha256"], "plan changed")
    require(software() == plan["software"], "runtime changed since prepare")
    for name, digest in {**plan["inputs"], **plan["raw_inputs"], **plan["code"]}.items():
        file = folder/name
        require(file.is_file() and sha(file) == digest, f"frozen input/code changed: {name}")
    for src in (SELF, HELPER, CONTROL_HELPER):
        key = "code/" + src.name
        require(sha(src) == plan["code"][key], f"controller/helper changed after prepare: {src.name}; use frozen code")
    return plan


def load_frozen(folder):
    folder = Path(folder)
    short, structure = folder/"inputs/short", folder/"inputs/structure"
    data = source_data(short, structure)
    ledger = read(folder/"inputs/chart-exposure-ledger.json")
    stats = validate_ledger(ledger, NAMES, data["runs"], expected_total_charts=data["total_charts"])
    require(stats == read(folder/"plan.json")["ledger_stats"], "ledger stats changed")
    wanted = [e["primitive"] for r in data["runs"] for e in r["events"]]
    enum = data["enumeration"]
    vocab = Vocabulary.load(short/"primitive-directions.tsv", 14, integer(enum["direction_count"]), F(enum["bound"]), wanted, data["form"])
    candidates, positions = build_candidate_population(vocab, data["runs"], static_limit=integer(read(folder/"plan.json")["policy"]["static_limit"]))
    prefixes = stage_prefixes(data["runs"], 14)
    ledger_by_seed = {r["seed"]: r for r in ledger["runs"]}
    return data, ledger, vocab, candidates, positions, prefixes, ledger_by_seed


def stage_rows(folder):
    data, ledger, vocab, candidates, positions, prefixes, ledger_by_seed = load_frozen(folder)
    rows = []
    for run in data["runs"]:
        for event in run["events"]:
            epoch = event["epoch"]
            charts = ledger_by_seed[run["seed"]]["stages"][epoch]["charts"]
            metrics = candidate_metrics_for_stage(charts, candidates, positions, prefixes[(run["seed"], epoch)],
                                                  event["primitive"], vocab, data["form"], 14)
            metrics.update(seed=run["seed"], epoch=epoch, quotient_dimension_before=len(prefixes[(run["seed"], epoch)]),
                           actual_word=list(event["primitive"]))
            rows.append(metrics)
    require(len(rows) == 180, "stage metric count changed")
    return data, positions, prefixes, rows


def do_stage(folder, stage):
    folder = Path(folder); plan = guard(folder)
    data, positions, prefixes, rows = stage_rows(folder)
    if stage == "exposure-census":
        compact = []
        for r in rows:
            special = {}
            for name, idx in (("L2", 1), ("L1", 0), ("L4", 3)):
                axis = tuple(int(i == idx) for i in range(14))
                found = None
                for c in [r["actual"], *r["controls"]]:
                    if tuple(c["word"]) == axis:
                        found = c; break
                # The special axis may be outside the actual rank-matched control
                # band. Recompute from chart rows if necessary.
                if found is None:
                    count = 0; first = None; min_bits = None
                    for ci, chart in enumerate(r["chart_rows"]):
                        if axis in map(tuple, chart["exposed"]):
                            count += 1; first = ci if first is None else first
                            bits = chart.get("quartic_coefficient_bits")
                            if bits is not None: min_bits = bits if min_bits is None else min(min_bits, bits)
                    found = {"word": axis, "count": count, "first_order": first, "min_quartic_bits": min_bits,
                             "static_rank": (positions[axis] + 1 if axis in positions else None)}
                special[name] = found
            compact.append({"seed": r["seed"], "epoch": r["epoch"],
                            "quotient_dimension_before": r["quotient_dimension_before"],
                            "chart_count": len(r["chart_rows"]), "candidate_count": r["candidate_count"],
                            "coverage_complete": r["coverage_complete"],
                            "complete_chart_count": r["complete_chart_count"],
                            "incomplete_or_unknown_chart_count": r["incomplete_or_unknown_chart_count"],
                            "exposed_candidate_count": r["exposed_candidate_count"],
                            "actual": r["actual"], "rank_band": r["rank_band"],
                            "control_count": len(r["controls"]), "special_axes": special})
        result = {"schema": "curve302-chart-exposure-census.v1", "status": "PASS_EXACT_HISTORICAL_EXPOSURE_CENSUS",
                  "stages": compact,
                  "summary": {"stages": len(compact), "charts": sum(x["chart_count"] for x in compact),
                              "actual_exposed": sum(x["actual"]["count"] > 0 for x in compact),
                              "complete_coverage_stages": sum(x["coverage_complete"] for x in compact),
                              "candidate_rule": plan["policy"]["candidate_rule"]},
                  "boundary": "Exposure means explicit quotient-word evidence in the historical chart audit/replay. Unsupported point-only records are not treated as misses. Absence-based comparisons require explicit complete chart coverage and otherwise remain UNKNOWN."}
        atomic(folder/"exposure-census.json", result)
    elif stage == "multiplicity":
        result = multiplicity_comparison(rows)
        result["schema"] = "curve302-chart-multiplicity.v1"
        atomic(folder/"multiplicity.json", result)
    elif stage == "counterfactuals":
        cf = []
        for r in rows:
            cf.append(counterfactual_stage(r, prefixes[(r["seed"], r["epoch"])], data["cores"], NAMES, positions,
                                           random_orders=integer(plan["policy"]["random_orders"]),
                                           master_seed=plan["policy"]["master_seed"]))
        result = summarize_counterfactuals(cf); result["schema"] = "curve302-stage-local-chart-order-controls.v1"
        atomic(folder/"counterfactuals.json", result)
    else:
        raise EvidenceError("unknown stage")


def output_path(folder, stage):
    return Path(folder)/f"{stage}.json"


def seal_path(folder, stage):
    return Path(folder)/"seals"/f"{stage}.json"


def write_seal(folder, stage):
    path = output_path(folder, stage)
    require(path.is_file(), f"stage output missing: {stage}")
    seal = {"stage": stage, "output": path.name, "sha256": sha(path)}
    atomic(seal_path(folder, stage), seal)


def verify_seal(folder, stage):
    seal = read(seal_path(folder, stage)); path = output_path(folder, stage)
    require(seal["stage"] == stage and seal["output"] == path.name and seal["sha256"] == sha(path), f"bad seal: {stage}")


def child_limits(seconds, memory_bytes):
    def apply():
        resource.setrlimit(resource.RLIMIT_CPU, (seconds, seconds + 5))
        resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
    return apply


def launch_stage(folder, stage):
    plan = guard(folder); policy = plan["policy"]
    log = Path(folder)/"logs"/f"{stage}.log"; log.parent.mkdir(exist_ok=True)
    cmd = [sys.executable, str(SELF), "_stage", "--folder", str(Path(folder).resolve()), "--stage", stage]
    with log.open("wb") as out:
        proc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT,
                              timeout=integer(policy["stage_seconds"]) + 30,
                              preexec_fn=child_limits(integer(policy["stage_seconds"]), integer(policy["memory_bytes"])))
    require(proc.returncode == 0, f"stage {stage} failed; see {log}")
    write_seal(folder, stage)


def finalize(folder):
    folder = Path(folder); plan = guard(folder)
    for stage in STAGES:
        verify_seal(folder, stage)
    report = {"status": "PASS_THREE_CHART_EXPOSURE_EXPERIMENTS", "plan_sha256": sha(folder/"plan.json"),
              "outputs": {stage: sha(output_path(folder, stage)) for stage in STAGES},
              "ledger_sha256": plan["inputs"]["inputs/chart-exposure-ledger.json"],
              "boundary": plan["boundary"]}
    atomic(folder/"REPORT.json", report)
    write_summary(folder)


def write_summary(folder):
    folder = Path(folder)
    c, m, f = read(folder/"exposure-census.json"), read(folder/"multiplicity.json"), read(folder/"counterfactuals.json")
    lines = ["# Curve302 historical chart exposure controls", "",
             f"Historical stages: **{c['summary']['stages']}**; charts: **{c['summary']['charts']}**.",
             f"Actual acquisitions explicitly exposed: **{c['summary']['actual_exposed']}/{c['summary']['stages']}**.",
             f"Stages with complete chart coverage: **{c['summary']['complete_coverage_stages']}/{c['summary']['stages']}**.", "",
             "## Norm-matched multiplicity", "",
             f"Usable rank-band stages: **{m['summary']['usable_rank_band_stages']}**.",
             f"Median exposure-count percentile: **{m['summary']['median_exposure_percentile']}**.",
             f"Median parameter-height percentile: **{m['summary']['median_parameter_height_percentile']}**.",
             f"Median quartic-bit percentile: **{m['summary']['median_quartic_bits_percentile']}**.", "",
             "## Stage-local ordering controls", ""]
    lines += ["| Policy | Draws | Gains | Match actual | Next observed core | L2 | L1 | L4 | Unknown stages |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in f["summary"]:
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            row["policy"], row.get("draws",0), row.get("gains",0), row.get("matches_actual",0),
            row.get("contains_observed_next_core_integrally",0), row.get("axis_L2",0), row.get("axis_L1",0),
            row.get("axis_L4",0), row.get("unknown_stages",0)))
    lines += ["", "The reordering controls are stage-local only: a counterfactual gain is never propagated into a later historical stage.",
              "Exposure is transcript/replay evidence, not a quotient-height proxy. These are retrospective controls on one known rank-31 subgroup."]
    (folder/"SUMMARY.md").write_text("\n".join(lines) + "\n")


def resume(folder):
    guard(folder)
    for stage in STAGES:
        out, seal = output_path(folder, stage), seal_path(folder, stage)
        if out.exists() or seal.exists():
            require(out.exists() and seal.exists(), f"interrupted/unsealed stage {stage}; use a fresh folder")
            verify_seal(folder, stage); continue
        launch_stage(folder, stage)
    finalize(folder)


def status(folder):
    folder = Path(folder)
    if (folder/"PREPARE_FAILED.json").exists():
        print(json.dumps(read(folder/"PREPARE_FAILED.json"), indent=2)); return
    plan = guard(folder)
    rows = []
    for stage in STAGES:
        state = "SEALED" if output_path(folder, stage).exists() and seal_path(folder, stage).exists() else "PENDING"
        if state == "SEALED":
            try: verify_seal(folder, stage)
            except Exception: state = "INVALID"
        rows.append({"stage": stage, "status": state})
    print(json.dumps({"plan": plan["status"], "stages": rows, "report": (folder/"REPORT.json").exists()}, indent=2))


def check(folder):
    folder = Path(folder); guard(folder)
    report = read(folder/"REPORT.json")
    require(report.get("status") == "PASS_THREE_CHART_EXPOSURE_EXPERIMENTS", "final report not passed")
    for stage in STAGES:
        verify_seal(folder, stage); require(report["outputs"][stage] == sha(output_path(folder, stage)), f"report hash mismatch: {stage}")
    # Deterministic byte-for-byte recomputation under frozen inputs/code.
    with tempfile.TemporaryDirectory(prefix="curve302-chart-exposure-check-") as tmp:
        tmp = Path(tmp)
        # Reuse original frozen folder through _stage; output files are redirected
        # by copying the sealed control state to a local scratch directory.
        scratch = tmp/"scratch"; shutil.copytree(folder, scratch)
        for stage in STAGES:
            output_path(scratch, stage).unlink(missing_ok=True); seal_path(scratch, stage).unlink(missing_ok=True)
            do_stage(scratch, stage)
            require(sha(output_path(scratch, stage)) == report["outputs"][stage], f"deterministic recomputation differs: {stage}")
    print("CHART_EXPOSURE_CHECK|status=PASS|stages=3|deterministic=PASS", flush=True)


def probe(source=None, structure=None, raw_root=None):
    short = discover_short(source); structure = discover_structure(short, structure)
    root = find_raw_root(structure, raw_root)
    print(json.dumps(schema_probe(root, NAMES), indent=2, sort_keys=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("prepare","run","resume","status","check","probe","_stage"))
    ap.add_argument("--folder", type=Path, default=DEFAULT)
    ap.add_argument("--source", type=Path)
    ap.add_argument("--structure", type=Path)
    ap.add_argument("--raw-root", type=Path)
    ap.add_argument("--ledger", type=Path)
    ap.add_argument("--static-limit", type=int, default=1000)
    ap.add_argument("--random-orders", type=int, default=256)
    ap.add_argument("--stage-seconds", type=int, default=3600)
    ap.add_argument("--memory-gib", type=int, default=6)
    ap.add_argument("--stage", choices=STAGES)
    args = ap.parse_args()
    try:
        if args.command == "probe":
            probe(args.source, args.structure, args.raw_root); return
        if args.command in ("prepare","run"):
            prepare(args.folder, args.source, args.structure, args.raw_root, args.ledger,
                    static_limit=args.static_limit, random_orders=args.random_orders,
                    stage_seconds=args.stage_seconds, memory_bytes=args.memory_gib*1024**3)
            if args.command == "prepare": return
            resume(args.folder); return
        if args.command == "resume": resume(args.folder); return
        if args.command == "status": status(args.folder); return
        if args.command == "check": check(args.folder); return
        if args.command == "_stage":
            require(args.stage is not None, "_stage requires --stage")
            do_stage(args.folder, args.stage); return
    except (EvidenceError, UnsupportedChartSchema, ValueError, OSError, subprocess.TimeoutExpired) as exc:
        print(f"CHART_EXPOSURE_FAIL|{type(exc).__name__}|{exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
