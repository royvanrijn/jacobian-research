"""Positive-evidence-only tables for the Curve302 historical chart replay.

This analysis intentionally never interprets an absent chart hit as a negative.
It consumes the sealed replay ledger and historical trajectory/common-core inputs
and extracts four descriptive tables:

1. exposure lead times for L2/L1/L4 and the observed rank-5/rank-6 common cores;
2. pre-containment positive multiplicity lower bounds;
3. positive-positive co-exposure choice sets at historical gain stages;
4. saturation impact of every positively exposed candidate direction.

All hypothetical one-step effects use Sat_Z(S + Z v) = span_Q(S,v) intersect Z^14,
so containment of an integral target vector is equivalent to rational-span
containment. No chart-completeness flag is required for any promoted statistic.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from sympy import Matrix

from curve302_short_core_controls import integer, primitive, require
from curve302_chart_exposure import (
    contains_core_integrally,
    integer_contains,
    normalize_explicit_ledger,
    rationally_unknown,
    stage_prefixes,
    validate_ledger,
)


AXIS_INDEX = {"L1": 0, "L2": 1, "L4": 3}


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def hash_obj(obj) -> str:
    return sha256(canonical_json(obj).encode()).hexdigest()


def matrix_rank(rows: Sequence[Sequence[int]], n: int) -> int:
    if not rows:
        return 0
    return int(Matrix([[int(x) for x in r] for r in rows]).rank())


def span_intersection_rank(rows, target, n: int) -> int:
    """Rank of span_Q(rows) intersect span_Q(target)."""
    rr = matrix_rank(rows, n)
    rt = matrix_rank(target, n)
    both = matrix_rank(tuple(rows) + tuple(target), n)
    return rr + rt - both


def saturation_contains(rows, target, n: int) -> bool:
    """Whether Sat_Z(<rows>) contains the integral target lattice generators."""
    return span_intersection_rank(rows, target, n) == matrix_rank(target, n)


def exposure_index(ledger):
    out = {}
    for run in ledger["runs"]:
        for stage in run["stages"]:
            key = (run["seed"], integer(stage["epoch"]))
            charts = []
            for list_order, chart in enumerate(stage["charts"]):
                words = []
                for exp in chart["exposures"]:
                    word = tuple(integer(x) for x in exp["word"])
                    words.append(word)
                charts.append({
                    "chart_id": chart["chart_id"],
                    "list_order": list_order,
                    "order": integer(chart.get("order", list_order)),
                    "words": tuple(words),
                })
            out[key] = tuple(charts)
    return out


def stage_positive_words(charts, prefix, n: int):
    """Known-positive, rationally-new words with chart multiplicity metadata."""
    stats = {}
    for chart in charts:
        seen = set()
        for raw in chart["words"]:
            w = primitive(raw)
            if w in seen:
                continue
            seen.add(w)
            if not rationally_unknown(prefix, w, n):
                continue
            row = stats.setdefault(w, {
                "word": w,
                "chart_count": 0,
                "first_list_order": chart["list_order"],
                "first_order": chart["order"],
                "chart_ids": [],
            })
            row["chart_count"] += 1
            row["first_list_order"] = min(row["first_list_order"], chart["list_order"])
            row["first_order"] = min(row["first_order"], chart["order"])
            row["chart_ids"].append(chart["chart_id"])
    return stats


def target_roster(data, n: int):
    axes = {}
    for label, idx in AXIS_INDEX.items():
        e = tuple(int(i == idx) for i in range(n))
        axes[label] = {"label": label, "kind": "axis", "basis": (e,), "direct_word": e,
                       "rank": 1, "observed_dimension": None}
    rank_to_dim = {}
    for d in sorted(data["cores"]):
        r = len(data["cores"][d])
        rank_to_dim.setdefault(r, d)
    require(5 in rank_to_dim and 6 in rank_to_dim, "rank-5/rank-6 common cores unavailable")
    cores = {}
    for rank in (5, 6):
        d = rank_to_dim[rank]
        label = f"C{rank}"
        cores[label] = {"label": label, "kind": "common_core", "basis": tuple(data["cores"][d]),
                        "direct_word": None, "rank": rank, "observed_dimension": d}
    # Preserve the narrative order L2, L1, L4, C5, C6.
    return [axes["L2"], axes["L1"], axes["L4"], cores["C5"], cores["C6"]]


def _target_contained(prefix, target, n: int) -> bool:
    if target["kind"] == "axis":
        return integer_contains(prefix, target["basis"][0])
    return contains_core_integrally(prefix, target["basis"])


def _post_prefix(prefix, gains):
    return tuple(prefix) + tuple(tuple(g["word"]) for g in gains)


def _stage_map(data):
    result = {}
    for run in data["runs"]:
        for stage in run["stages"]:
            result[(run["seed"], integer(stage["epoch"]))] = stage
    return result


def _run_map(data):
    return {run["seed"]: run for run in data["runs"]}


def _median(values):
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    m = len(vals)//2
    return vals[m] if len(vals)%2 else (vals[m-1]+vals[m]) / 2


def extract_lead_times(data, ledger, n=14):
    prefixes = stage_prefixes(data["runs"], n)
    eidx = exposure_index(ledger)
    targets = target_roster(data, n)
    rows = []
    for run in data["runs"]:
        seed = run["seed"]
        for target in targets:
            first_direct = first_enable = first_progress = None
            containment = None
            if run["stages"]:
                first_epoch = integer(run["stages"][0]["epoch"])
                initial_prefix = prefixes[(seed, first_epoch)]
                if _target_contained(initial_prefix, target, n):
                    containment = {"epoch": None, "stage_index": -1, "post_dimension": len(initial_prefix), "initial_seed_contains": True}
            for stage_i, stage in enumerate(run["stages"]):
                epoch = integer(stage["epoch"])
                prefix = prefixes[(seed, epoch)]
                pre_dim = len(prefix)
                pre_contained = _target_contained(prefix, target, n)
                positives = stage_positive_words(eidx[(seed, epoch)], prefix, n)
                current_ir = span_intersection_rank(prefix, target["basis"], n)
                if not pre_contained:
                    if target["direct_word"] is not None and target["direct_word"] in positives and first_direct is None:
                        first_direct = {"epoch": epoch, "stage_index": stage_i, "pre_dimension": pre_dim,
                                        "chart_count": positives[target["direct_word"]]["chart_count"]}
                    for w, meta in positives.items():
                        ir = span_intersection_rank(tuple(prefix)+(w,), target["basis"], n)
                        if ir > current_ir and first_progress is None:
                            first_progress = {"epoch": epoch, "stage_index": stage_i, "pre_dimension": pre_dim,
                                              "word": w, "chart_count": meta["chart_count"],
                                              "intersection_rank_before": current_ir, "intersection_rank_after": ir}
                        if saturation_contains(tuple(prefix)+(w,), target["basis"], n) and first_enable is None:
                            first_enable = {"epoch": epoch, "stage_index": stage_i, "pre_dimension": pre_dim,
                                            "word": w, "chart_count": meta["chart_count"]}
                post = _post_prefix(prefix, stage.get("gains", ()))
                if containment is None and _target_contained(post, target, n):
                    containment = {"epoch": epoch, "stage_index": stage_i, "post_dimension": len(post)}
            def lead(first):
                if first is None or containment is None:
                    return None
                return containment["stage_index"] - first["stage_index"]
            rows.append({
                "seed": seed, "target": target["label"], "target_kind": target["kind"],
                "target_rank": target["rank"], "target_observed_common_dimension": target["observed_dimension"],
                "first_direct_positive": first_direct,
                "first_saturation_enabling_positive": first_enable,
                "first_progress_positive": first_progress,
                "first_integral_containment": containment,
                "direct_lead_stages": lead(first_direct),
                "enabling_lead_stages": lead(first_enable),
                "progress_lead_stages": lead(first_progress),
            })
    summary = []
    for target in targets:
        rr = [r for r in rows if r["target"] == target["label"]]
        summary.append({
            "target": target["label"], "kind": target["kind"], "rank": target["rank"],
            "observed_common_dimension": target["observed_dimension"], "runs": len(rr),
            "contained_runs": sum(r["first_integral_containment"] is not None for r in rr),
            "direct_preexposed_runs": sum(r["first_direct_positive"] is not None for r in rr),
            "saturation_enabling_preexposed_runs": sum(r["first_saturation_enabling_positive"] is not None for r in rr),
            "progress_preexposed_runs": sum(r["first_progress_positive"] is not None for r in rr),
            "median_direct_lead_stages": _median([r["direct_lead_stages"] for r in rr]),
            "median_enabling_lead_stages": _median([r["enabling_lead_stages"] for r in rr]),
            "median_progress_lead_stages": _median([r["progress_lead_stages"] for r in rr]),
        })
    return {"status": "PASS_POSITIVE_EXPOSURE_LEAD_TIMES", "summary": summary, "runs": rows,
            "boundary": "Only recorded positive exposures are used. Missing chart hits are never interpreted as non-exposure."}


def extract_precontainment_multiplicity(data, ledger, n=14):
    prefixes = stage_prefixes(data["runs"], n); eidx = exposure_index(ledger); targets = target_roster(data,n)
    rows=[]
    for run in data["runs"]:
        seed=run["seed"]
        for target in targets:
            direct_charts=set(); enable_charts=set(); progress_charts=set()
            direct_epochs=set(); enable_epochs=set(); progress_epochs=set()
            enable_words=set(); progress_words=set(); examined=0
            for stage in run["stages"]:
                epoch=integer(stage["epoch"]); prefix=prefixes[(seed,epoch)]
                if _target_contained(prefix,target,n):
                    break
                examined += 1
                current_ir=span_intersection_rank(prefix,target["basis"],n)
                positives=stage_positive_words(eidx[(seed,epoch)],prefix,n)
                for w,meta in positives.items():
                    if target["direct_word"] is not None and w == target["direct_word"]:
                        direct_epochs.add(epoch)
                        direct_charts.update((epoch,cid) for cid in meta["chart_ids"])
                    ir=span_intersection_rank(tuple(prefix)+(w,),target["basis"],n)
                    if ir > current_ir:
                        progress_epochs.add(epoch); progress_words.add(w)
                        progress_charts.update((epoch,cid) for cid in meta["chart_ids"])
                    if saturation_contains(tuple(prefix)+(w,),target["basis"],n):
                        enable_epochs.add(epoch); enable_words.add(w)
                        enable_charts.update((epoch,cid) for cid in meta["chart_ids"])
            rows.append({"seed":seed,"target":target["label"],"target_kind":target["kind"],"target_rank":target["rank"],
                         "stages_before_containment":examined,
                         "direct_chart_lower_bound":len(direct_charts),"direct_epoch_count":len(direct_epochs),
                         "progress_chart_lower_bound":len(progress_charts),"progress_epoch_count":len(progress_epochs),
                         "completion_enabling_chart_lower_bound":len(enable_charts),"completion_enabling_epoch_count":len(enable_epochs),
                         "distinct_progress_directions":len(progress_words),"distinct_completion_enabling_directions":len(enable_words)})
    summary=[]
    for target in targets:
        rr=[r for r in rows if r["target"]==target["label"]]
        summary.append({"target":target["label"],"kind":target["kind"],"rank":target["rank"],"runs":len(rr),
                        "direct_chart_lower_bound_total":sum(r["direct_chart_lower_bound"] for r in rr),
                        "progress_chart_lower_bound_total":sum(r["progress_chart_lower_bound"] for r in rr),
                        "completion_enabling_chart_lower_bound_total":sum(r["completion_enabling_chart_lower_bound"] for r in rr),
                        "runs_with_direct_positive":sum(r["direct_chart_lower_bound"]>0 for r in rr),
                        "runs_with_progress_positive":sum(r["progress_chart_lower_bound"]>0 for r in rr),
                        "runs_with_completion_enabling_positive":sum(r["completion_enabling_chart_lower_bound"]>0 for r in rr)})
    return {"status":"PASS_PRECONTAINMENT_POSITIVE_MULTIPLICITY_LOWER_BOUNDS","summary":summary,"runs":rows,
            "boundary":"Counts are lower bounds from recorded positive chart hits only; no missing hit contributes a zero."}


def extract_coexposure_choices(data, ledger, n=14):
    prefixes=stage_prefixes(data["runs"],n); eidx=exposure_index(ledger); rows=[]
    for run in data["runs"]:
        seed=run["seed"]
        for stage in run["stages"]:
            gains=stage.get("gains",())
            if not gains:
                continue
            epoch=integer(stage["epoch"]); prefix=prefixes[(seed,epoch)]
            pos=stage_positive_words(eidx[(seed,epoch)],prefix,n)
            actual={primitive(g["primitive"]) for g in gains}
            require(actual <= set(pos),f"actual batch not fully positively exposed for {seed} epoch {epoch}")
            alternatives=set(pos)-actual
            min_all=min(m["first_list_order"] for m in pos.values())
            actual_first=min(pos[w]["first_list_order"] for w in actual)
            alt_first=min((pos[w]["first_list_order"] for w in alternatives),default=None)
            max_actual=max(pos[w]["chart_count"] for w in actual)
            max_alt=max((pos[w]["chart_count"] for w in alternatives),default=0)
            earliest_words=sorted([w for w,m in pos.items() if m["first_list_order"]==min_all])
            positive_directions=[{"word":w,"actual_member":w in actual,"chart_count_lower_bound":pos[w]["chart_count"],
                                  "first_positive_chart_order":pos[w]["first_list_order"]} for w in sorted(pos)]
            rows.append({"seed":seed,"epoch":epoch,"pre_dimension":len(prefix),"post_dimension":len(prefix)+len(gains),
                         "batch_size":len(gains),"actual_words":sorted(actual),"positive_directions":positive_directions,
                         "positive_unknown_direction_count":len(pos),"positive_chart_direction_records":sum(m["chart_count"] for m in pos.values()),
                         "positive_alternative_direction_count":len(alternatives),
                         "earliest_positive_chart_order":min_all,"earliest_positive_words":earliest_words,
                         "actual_first_positive_chart_order":actual_first,"alternative_first_positive_chart_order":alt_first,
                         "actual_hits_earliest_positive":bool(actual & set(earliest_words)),
                         "alternative_precedes_all_actual":alt_first is not None and alt_first < actual_first,
                         "actual_chart_count_lower_bounds":{','.join(map(str,w)):pos[w]["chart_count"] for w in sorted(actual)},
                         "max_actual_chart_count_lower_bound":max_actual,"max_alternative_chart_count_lower_bound":max_alt,
                         "some_alternative_has_higher_recorded_multiplicity":max_alt>max_actual})
    summary={"gain_stages":len(rows),"single_gain_stages":sum(r["batch_size"]==1 for r in rows),"double_gain_stages":sum(r["batch_size"]==2 for r in rows),
             "stages_with_positive_alternatives":sum(r["positive_alternative_direction_count"]>0 for r in rows),
             "stages_actual_hits_earliest_positive":sum(r["actual_hits_earliest_positive"] for r in rows),
             "stages_alternative_precedes_all_actual":sum(r["alternative_precedes_all_actual"] for r in rows),
             "stages_some_alternative_has_higher_recorded_multiplicity":sum(r["some_alternative_has_higher_recorded_multiplicity"] for r in rows)}
    return {"status":"PASS_POSITIVE_COEXPOSURE_CHOICE_TABLE","summary":summary,"stages":rows,
            "boundary":"Comparisons are only among directions with explicit positive exposure evidence in the same historical stage. They are not complete choice-set claims."}


def extract_saturation_impact(data, ledger, n=14):
    prefixes=stage_prefixes(data["runs"],n); eidx=exposure_index(ledger); targets=target_roster(data,n)
    target_by={t["label"]:t for t in targets}; details=[]; stages=[]
    for run in data["runs"]:
        seed=run["seed"]
        for stage in run["stages"]:
            gains=stage.get("gains",())
            if not gains:
                continue
            epoch=integer(stage["epoch"]); prefix=prefixes[(seed,epoch)]; pos=stage_positive_words(eidx[(seed,epoch)],prefix,n)
            actual={primitive(g["primitive"]) for g in gains}; post_dim=len(prefix)+len(gains)
            next_core=data["cores"].get(post_dim,())
            current_next=span_intersection_rank(prefix,next_core,n) if next_core else 0
            candidate_rows=[]
            for w,meta in sorted(pos.items()):
                row={"seed":seed,"epoch":epoch,"pre_dimension":len(prefix),"post_dimension":post_dim,
                     "candidate_word":w,"actual_member":w in actual,"positive_chart_count_lower_bound":meta["chart_count"],
                     "first_positive_chart_order":meta["first_list_order"]}
                if next_core:
                    ir=span_intersection_rank(tuple(prefix)+(w,),next_core,n)
                    row.update(next_common_core_rank=len(next_core),current_next_common_core_intersection_rank=current_next,
                               candidate_next_common_core_intersection_rank=ir,next_common_core_delta=ir-current_next,
                               completes_next_common_core_after_saturation=(ir==len(next_core)))
                else:
                    row.update(next_common_core_rank=0,current_next_common_core_intersection_rank=0,
                               candidate_next_common_core_intersection_rank=0,next_common_core_delta=0,
                               completes_next_common_core_after_saturation=True)
                for label in ("L2","L1","L4","C5","C6"):
                    t=target_by[label]; before=span_intersection_rank(prefix,t["basis"],n); after=span_intersection_rank(tuple(prefix)+(w,),t["basis"],n)
                    row[f"{label}_delta"]=after-before
                    row[f"{label}_complete_after_saturation"]=after==t["rank"]
                candidate_rows.append(row); details.append(row)
            alt=[r for r in candidate_rows if not r["actual_member"]]; act=[r for r in candidate_rows if r["actual_member"]]
            post=_post_prefix(prefix,gains)
            actual_batch_ir=span_intersection_rank(post,next_core,n) if next_core else 0
            max_single=max((r["next_common_core_delta"] for r in candidate_rows),default=0)
            max_actual=max((r["next_common_core_delta"] for r in act),default=0)
            max_alt=max((r["next_common_core_delta"] for r in alt),default=0)
            stages.append({"seed":seed,"epoch":epoch,"pre_dimension":len(prefix),"post_dimension":post_dim,"batch_size":len(gains),
                           "positive_unknown_direction_count":len(candidate_rows),"positive_alternative_direction_count":len(alt),
                           "next_common_core_rank":len(next_core),"current_next_common_core_intersection_rank":current_next,
                           "actual_batch_next_common_core_intersection_rank":actual_batch_ir,
                           "actual_batch_next_common_core_delta":actual_batch_ir-current_next,
                           "best_positive_single_delta":max_single,"best_actual_single_delta":max_actual,"best_alternative_single_delta":max_alt,
                           "actual_single_tied_for_best_positive":bool(act) and max_actual==max_single,
                           "positive_single_improvers":sum(r["next_common_core_delta"]>0 for r in candidate_rows),
                           "alternative_single_improvers":sum(r["next_common_core_delta"]>0 for r in alt),
                           "positive_single_completers":sum(r["completes_next_common_core_after_saturation"] for r in candidate_rows),
                           "alternative_single_completers":sum(r["completes_next_common_core_after_saturation"] for r in alt)})
    summary={"gain_stages":len(stages),"positive_candidate_stage_pairs":len(details),
             "stages_with_positive_alternatives":sum(s["positive_alternative_direction_count"]>0 for s in stages),
             "stages_with_any_positive_single_core_improver":sum(s["positive_single_improvers"]>0 for s in stages),
             "stages_with_alternative_single_core_improver":sum(s["alternative_single_improvers"]>0 for s in stages),
             "stages_actual_single_tied_for_best_positive":sum(s["actual_single_tied_for_best_positive"] for s in stages),
             "stages_with_alternative_single_core_completer":sum(s["alternative_single_completers"]>0 for s in stages)}
    return {"status":"PASS_POSITIVE_SATURATION_IMPACT","summary":summary,"stages":stages,"candidates":details,
            "boundary":"Hypothetical impact is Sat_Z(S+Zv), evaluated exactly as rational-span intersection with integral target cores. Only explicitly positively exposed v are compared."}


def extract_all(data, ledger, n=14):
    return {
        "lead_times": extract_lead_times(data,ledger,n),
        "precontainment_multiplicity": extract_precontainment_multiplicity(data,ledger,n),
        "coexposure_choices": extract_coexposure_choices(data,ledger,n),
        "saturation_impact": extract_saturation_impact(data,ledger,n),
    }
