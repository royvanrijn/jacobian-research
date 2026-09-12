#!/usr/bin/env python3
"""Final report layer: labeled views, isolated cohorts, deterministic replay.

No worker, population, control bank, original summary or correlation is modified.
Only the prospective rows enter summarize_profiles; historical rows enter the
explicitly retrospective comparison table, never a population histogram.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from io import StringIO
from pathlib import Path

import historical_external_arithmetic as lane
import wide_arithmetic_profile_core as core

read, sha, require = lane.read, lane.sha, lane.require
METRICS = ("bk_local_term", "phi_m_count", "phi_a_count", "log2_abs_field_discriminant",
           "field_ramified_prime_count", "log2_conductor", "log2_abs_minimal_discriminant",
           "forced_class_2rank_lower_from_known_rank")


def validate_cohorts(prospective, historical):
    for rows, kind in ((prospective, lane.PROSPECTIVE), (historical, lane.HISTORICAL)):
        require(all((r.get("cohort"), r.get("selection_mode")) == kind for r in rows), "cohort contamination")
        require(len({r["curve_key"] for r in rows}) == len(rows), "duplicate cohort keys")
    require(not ({r["curve_key"] for r in prospective} & {r["curve_key"] for r in historical}), "overlapping cohort keys")


def prospective_summary(prospective, historical):
    validate_cohorts(prospective, historical)
    # Deliberately never pass the retrospective panel to this function.
    return core.summarize_profiles(prospective)


def summarize_group(name, rows, kind):
    require(rows and all((r["cohort"], r["selection_mode"]) == kind for r in rows), "mixed/empty comparison group")
    out = {"group": name, **lane.labels(kind), "count": len(rows),
           "base_status_counts": dict(sorted(Counter(r["base_status"] for r in rows).items())),
           "local_status_counts": dict(sorted(Counter(r["local_status"] for r in rows).items()))}
    for metric in METRICS:
        values = [r[metric] for r in rows if r.get(metric) is not None]
        out[f"n_{metric}"] = len(values)
        out[f"median_{metric}"] = core.median(values)
        out[f"min_{metric}"] = min(values) if values else None
        out[f"max_{metric}"] = max(values) if values else None
    roots = [r["root_number"] for r in rows if r.get("root_number") is not None]
    signatures = [r["field_signature"] for r in rows if r.get("field_signature") is not None]
    out.update(root_minus_one_count=roots.count(-1), root_known_count=len(roots),
               totally_real_count=signatures.count([3, 0]), signature_known_count=len(signatures))
    return out


def comparison(prospective, historical):
    validate_cohorts(prospective, historical)
    rows = [summarize_group("Historical 11952, LB >=27", historical, lane.HISTORICAL)]
    broad25 = [r for r in prospective if r["family"] == "11952" and r["t"] == "921/653"]
    require(len(broad25) == 1 and broad25[0]["final_rank_lower_bound"] == 25, "broad 11952 control mismatch")
    groups = [("Broad 11952 @ 921/653, LB 25", broad25),
              ("Prospective LB >=23 tail", [r for r in prospective if r["final_rank_lower_bound"] >= 23])]
    for rank in (23, 24, 25):
        groups.append((f"Prospective LB {rank}", [r for r in prospective if r["final_rank_lower_bound"] == rank]))
    groups.extend([("Full prospective population", prospective),
                   ("Prospective 11952 only (same-family diagnostic)", [r for r in prospective if r["family"] == "11952"])])
    rows.extend(summarize_group(name, group, lane.PROSPECTIVE) for name, group in groups)
    return rows


def labeled_profile(row, directory, kind):
    key = row["curve_key"]
    stages = {mode: read(directory / mode / f"{key}.json") for mode in ("base", "local")}
    profile = core.flatten_profile(row, stages["base"], stages["local"], None)
    profile.update(lane.labels(kind))
    profile.update(family=row["family"], base_status=stages["base"]["status"], local_status=stages["local"]["status"])
    profile["checkpoint_sha256"] = {mode: sha(directory / mode / f"{key}.json") for mode in stages}
    if kind == lane.HISTORICAL:
        for name in ("rank_provenance", "local_search_rank_lower_bound"):
            profile[name] = row[name]
    return profile


def fmt(value):
    if value is None:
        return "UNKNOWN"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def table(headers, rows):
    def cell(value):
        return fmt(value).replace("|", r"\|").replace("\n", " ")
    return ["| " + " | ".join(cell(h) for h in headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"] + [
        "| " + " | ".join(cell(v) for v in row) + " |" for row in rows]


def markdown(summary, groups, historical, plan):
    counts = summary["stage_status_counts"]
    lines = ["# Completed arithmetic census and separate historical 11952 controls", "",
             "All ranks below are **certified lower bounds**, not exact ranks. "
             "This report runs no point searches or class-group probes.", "",
             f"Prospective census: {summary['population_count']} fibres; BASE {counts['base']}; LOCAL {counts['local']}.", "",
             "## Cohort isolation and replay", "",
             "The 2,080 rows are labeled `cohort = prospective_broad_2080`, "
             "`selection_mode = frozen_prospective` in `prospective_profiles.json`. "
             "The historical rows are labeled `cohort = historical_external`, "
             "`selection_mode = retrospective_known_high_rank` in `historical_profiles.json`. "
             "These are final-analysis views: the original worker schemas/checkpoints remain unedited.", "",
             "The full 445-row repository inventory roster was frozen before selecting every "
             "11952 fibre with certified lower bound at least 27. All five qualifying certificates "
             "were replayed exactly: point membership, independent finite-reduction columns, "
             "rational 2-torsion exclusion, specialization isomorphism and point transport. "
             "The rank-28 entry is a public-point reproduction; its original local-search bound is 27.", "",
             "Historical arithmetic calls the **identical unchanged** BASE/LOCAL worker and controller "
             "with the census's Sage version, 60/180-second per-worker budgets and memory policy. "
             "At most two historical workers run together. No cached factorization or known points "
             "are supplied to the workers. Every timeout/error is retained as UNKNOWN; no refills.", "",
             "[Arithmetic replay](ARITHMETIC_REPLAY.json) checks all five BASE results and all three "
             "completed LOCAL results, using integral-monic 2-division fields, certified PARI maximal "
             "orders, prime-factor proofs and exact local reductions. This shares Sage/PARI with the "
             "worker; it is not a separate full descent. Neither timed-out LOCAL result is filled in.", "",
             f"The hashes and file set of all {len(plan['prospective_files'])} existing census files "
             "are checked unchanged. This includes the original plan, source snapshots, worker inputs, "
             "results, summaries, frozen CLASS control selection and launch state. "
             "The original prospective stratification and Spearman diagnostics are reproduced exactly. "
             "Historical rows never enter those computations or any population histogram.", "",
             "## Historical controls: per fibre", ""]
    lines += table(["11952 parameter", "Rank LB", "BASE", "LOCAL", "u+n", "#Phi_m", "#Phi_a",
                    "log2 |D_K|", "Ramified primes", "Signature", "Root", "log2 conductor", "Forced g lower*"],
                   [[r["t"], r["final_rank_lower_bound"], r["base_status"], r["local_status"],
                     r.get("bk_local_term"), r.get("phi_m_count"), r.get("phi_a_count"),
                     r.get("log2_abs_field_discriminant"), r.get("field_ramified_prime_count"),
                     r.get("field_signature"), r.get("root_number"), r.get("log2_conductor"),
                     r.get("forced_class_2rank_lower_from_known_rank")] for r in historical])
    lines += ["", "## Separate retrospective/prospective comparison", "",
              "Medians use available values only, never zero for UNKNOWN. Groups overlap "
              "(the LB-25 fibre also belongs to the tail and full population); they are not independent samples. "
              "The historical panel is selected after knowing high lower bounds and is not a prospective cohort.", ""]
    lines += table(["Group", "Fibres", "LOCAL PASS", "Median u+n", "Median #Phi_m", "Median #Phi_a",
                    "Median log2 |D_K|", "Median ramified primes", "Median log2 conductor", "Median forced g lower*"],
                   [[r["group"], r["count"], r["local_status_counts"].get("PASS", 0),
                     *[r[f"median_{m}"] for m in ("bk_local_term", "phi_m_count", "phi_a_count",
                         "log2_abs_field_discriminant", "field_ramified_prime_count", "log2_conductor",
                         "forced_class_2rank_lower_from_known_rank")]] for r in groups])
    lines += ["", "The JSON/CSV comparison retains the nonmissing count and range for each metric, "
              "BASE/LOCAL status counts, root-number counts and field signatures.", "",
              "*`forced g lower = max(0, rank_LB-(u+n))` uses the known rank lower bound. "
              "It is **not independent equation-derived evidence** and is not a class-group upper bound.", "",
              "## Unchanged prospective strata", ""]
    lines += table(["Rank LB stratum", "Fibres", "Known u+n", "Median u+n", "Median log2 |D_K|", "Median log2 conductor"],
                   [[name, r["count"], r["n_bk_local_term"], r["median_bk_local_term"],
                     r["median_log2_abs_field_discriminant"], r["median_log2_conductor"]]
                    for name, r in summary["bucket_summary"].items()])
    lines += ["", "## Unchanged prospective descriptive diagnostics", ""]
    lines += table(["Metric", "Pairs", "Spearman vs final rank LB", "Spearman vs follow-up gain"],
                   [[m, summary["correlation_pair_counts"][m]["final"],
                     summary["correlations"][f"spearman_{m}_vs_final_lower_bound"],
                     summary["correlations"][f"spearman_{m}_vs_followup_gain"]]
                    for m in ("bk_local_term", "log2_abs_minimal_discriminant", "log2_abs_field_discriminant", "log2_conductor")])
    lines += ["", "The completed-case local-term association is weak. Discriminant/conductor sizes "
              "show a stronger negative association with discovered lower bound. Neither establishes "
              "a rank mechanism: the original search was adaptive, ranks are censored lower bounds, "
              "parent composition can confound pooled diagnostics, and LOCAL completion is heavily censored "
              "by the fixed timeout. No p-values, exact ranks, full Selmer dimensions, or class-group upper "
              "bounds are inferred. A weak local association does not locate a missing signal in g.", "",
              "`check` byte-rebuilds every analysis output and verifies the untouched census and bound sources.", ""]
    return "\n".join(lines)


def csv_text(rows):
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=sorted({k for r in rows for k in r}))
    writer.writeheader()
    writer.writerows({k: core.stable_json(v).strip() if isinstance(v, (dict, list)) else v for k, v in r.items()} for r in rows)
    return output.getvalue()


def build(out):
    plan = lane.check(out, complete=True)
    replay = read(out / "ARITHMETIC_REPLAY.json")
    require(replay["status"] == "PASS_EXISTING_CHECKPOINT_REPLAY", "arithmetic replay missing")
    require(sha(Path(__file__).with_name("verify_historical_external_arithmetic.sage")) == replay["checker_sha256"], "arithmetic checker changed")
    expected_bindings = {f"{mode}/{r['curve_key']}.json": sha(out / mode / f"{r['curve_key']}.json")
                         for r in plan["rows"] for mode in ("inputs", "base", "local")}
    require(expected_bindings == replay["checkpoint_sha256"], "arithmetic replay input mismatch")
    census = Path(plan["census"])
    census_plan, population = lane.worker.check_plan(census)
    sources = [r for r in census_plan["source_files_with_rows"] if Path(r["file"]).name == "queue.json"]
    require(len(sources) == 1 and sha(Path(sources[0]["file"])) == sources[0]["sha256"], "family queue hash mismatch")
    queue = {r["id"]: r for r in read(sources[0]["file"])["rows"]}
    prospective, raw = [], []
    for row in population:
        source = queue[row["id"]]
        require([int(v) for v in source["model"]] == row["ainvs"] and source["parameter"] == row["t"], "family metadata not bound to census equation")
        enriched = {**row, "family": source["family"]}
        prospective.append(labeled_profile(enriched, census, lane.PROSPECTIVE))
        key = row["curve_key"]
        require(not (census / "class" / f"{key}.json").exists(), "CLASS unexpectedly active")
        raw.append(core.flatten_profile(row, read(census / "base" / f"{key}.json"), read(census / "local" / f"{key}.json"), None))
    require(raw == read(census / "profiles.json")["rows"], "original prospective profile mismatch")
    historical = [labeled_profile(r, out, lane.HISTORICAL) for r in plan["rows"]]
    summary = prospective_summary(prospective, historical)
    summary["stage_status_counts"] = {mode: dict(sorted(Counter(r[f"{mode}_status"] for r in prospective).items())) for mode in ("base", "local")}
    summary["stage_status_counts"]["class"] = {}
    require(summary == read(census / "summary.json"), "prospective diagnostics changed")
    groups = comparison(prospective, historical)
    outputs = {
        "prospective_profiles.json": core.stable_json({"schema": "elliptic-curves.labeled-arithmetic-profiles.v1", "rows": prospective}),
        "historical_profiles.json": core.stable_json({"schema": "elliptic-curves.labeled-arithmetic-profiles.v1", "rows": historical}),
        "prospective_summary.json": core.stable_json(summary),
        "comparison.json": core.stable_json({"schema": "elliptic-curves.separate-cohort-comparison.v1", "groups": groups}),
        "comparison.csv": csv_text(groups),
        "SUMMARY.md": markdown(summary, groups, historical, plan),
    }
    outputs["REPORT.json"] = core.stable_json({
        "schema": "elliptic-curves.historical-external-report.v1", "status": "PASS_SEPARATE_COHORT_ANALYSIS",
        "prospective_count": len(prospective), "historical_count": len(historical),
        "historical_stage_counts": {mode: dict(sorted(Counter(r[f"{mode}_status"] for r in historical).items())) for mode in ("base", "local")},
        "prospective_files_unchanged": len(plan["prospective_files"]),
        "prospective_diagnostics_unchanged": True, "prospective_population_sha256": census_plan["population_sha256"],
        "class_workers": 0, "point_searches_launched": 0,
        "plan_sha256": sha(out / "plan.json"), "reporter_sha256": sha(Path(__file__)),
        "arithmetic_replay_sha256": sha(out / "ARITHMETIC_REPLAY.json"),
        "family_queue_sha256": sources[0]["sha256"],
        "output_sha256": {name: core.sha256_bytes(content.encode()) for name, content in outputs.items()},
    })
    return outputs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("report", "check"))
    ap.add_argument("--output", type=Path, default=lane.DEFAULT_OUTPUT)
    args = ap.parse_args()
    out = args.output.resolve()
    outputs = build(out)
    if args.command == "report":
        # Generated artifacts: refuse to overwrite any different prior result.
        for name, content in outputs.items():
            path = out / name
            require(not path.exists() or path.read_bytes() == content.encode(), f"refusing changed report: {name}")
        for name, content in outputs.items():
            path = out / name
            if not path.exists():
                path.write_bytes(content.encode())
    else:
        require(all((out / name).read_bytes() == content.encode() for name, content in outputs.items()), "deterministic report replay failed")
    print(f"HISTORICAL_ANALYSIS|{args.command}|PASS|prospective=2080|historical={len(read(out / 'plan.json')['rows'])}|no_cohort_mixing=PASS")


if __name__ == "__main__":
    main()
