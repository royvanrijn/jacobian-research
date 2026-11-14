#!/usr/bin/env python3
"""Checkpointed determinant-1020 MW17 specialization scan.

This controller ranks a fixed rational-height population by a local Frobenius
score, interleaves ranked, hash-control, and score-diverse arms, then invokes
the exact full-basis pointed search for each selected fibre.  Scores schedule
work only.  New independent points are exact lower-bound witnesses; misses,
timeouts, and setup failures are retained without rank-exclusion meaning.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, log, lcm
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / "research/artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1"
MODEL = PACKET / "rr-weierstrass-invariants-QQ.json"
SECTIONS = PACKET / "rr-seventeen-sections-QQ.json"
RUNNER = Path(__file__).with_name("run_mw_search.py")
DEFAULT = ROOT / "research/artifacts/local/det1020-amplification-scan-v1"
PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def write_new(path: Path, value: object) -> None:
    data = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text() != data:
            raise RuntimeError(f"immutable file changed: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data)


def checkpoint(path: Path, value: object) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def valuation(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def character(value: int, prime: int) -> int:
    value %= prime
    if not value:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def local_trace(a: int, b: int, prime: int) -> int | None:
    scale = min(valuation(a, prime) // 4, valuation(b, prime) // 6)
    a //= prime ** (4 * scale)
    b //= prime ** (6 * scale)
    if (4 * pow(a, 3, prime) + 27 * pow(b, 2, prime)) % prime == 0:
        return None
    trace = -sum(character(x**3 + a*x + b, prime) for x in range(prime))
    if trace * trace > 4 * prime:
        raise ArithmeticError("Hasse bound failed")
    return trace


def evaluate(coefficients: list[str], parameter: Fraction) -> Fraction:
    return sum(Fraction(value) * parameter**i for i, value in enumerate(coefficients))


def homogeneous(coefficients: list[int], numerator: int, denominator: int, degree: int) -> int:
    return sum(value * numerator**i * denominator**(degree - i) for i, value in enumerate(coefficients))


def population(height: int, model: dict) -> list[dict]:
    a_coefficients = [int(Fraction(value)) for value in model["A"]]
    b_coefficients = [int(Fraction(value)) for value in model["B"]]
    rows = []
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if gcd(numerator, denominator) != 1:
                continue
            a = homogeneous(a_coefficients, numerator, denominator, 8)
            b = homogeneous(b_coefficients, numerator, denominator, 12)
            if 4*a**3 + 27*b**2 == 0:
                continue
            units = 0
            good = 0
            for prime in PRIMES:
                trace = local_trace(a, b, prime)
                if trace is None:
                    continue
                good += 1
                units += round(10**12 * (2 - trace) * log(prime) / (prime + 1 - trace))
            rows.append({"parameter": str(Fraction(numerator, denominator)), "numerator": numerator,
                         "denominator": denominator, "score_units": units, "good_primes": good,
                         "coefficient_bits": max(abs(a).bit_length(), abs(b).bit_length())})
    return rows


def select(rows: list[dict], ranked_count: int, control_count: int, diversity_count: int) -> list[dict]:
    if len(rows) < ranked_count + control_count + diversity_count:
        raise ValueError("population is too small for selected arms")
    ranked = sorted(rows, key=lambda row: (-row["score_units"], row["coefficient_bits"], row["numerator"], row["denominator"]))[:ranked_count]
    used = {row["parameter"] for row in ranked}
    remaining = [row for row in rows if row["parameter"] not in used]
    controls = sorted(remaining, key=lambda row: (sha256(("det1020-amplification-controls-v1/" + row["parameter"]).encode()).hexdigest(), row["parameter"]))[:control_count]
    used.update(row["parameter"] for row in controls)
    pool = sorted((row for row in rows if row["parameter"] not in used), key=lambda row: (row["score_units"], row["parameter"]))
    diversity = []
    for index in range(diversity_count):
        target = round(index * (len(pool) - 1) / max(1, diversity_count - 1))
        candidate = pool[target]
        if candidate["parameter"] not in {row["parameter"] for row in diversity}:
            diversity.append(candidate)
    if len(diversity) != diversity_count:
        raise RuntimeError("diversity selection collided")
    arms = {"ranked": ranked, "control": controls, "diversity": diversity}
    out = []
    maximum = max(map(len, arms.values()))
    for index in range(maximum):
        for arm, width in (("ranked", 4), ("control", 1), ("diversity", 1)):
            start = index * width
            for row in arms[arm][start:start + width]:
                out.append({**row, "arm": arm})
    if len({row["parameter"] for row in out}) != len(out):
        raise RuntimeError("selection arms overlap")
    return out


def section_scale(sections: dict) -> int:
    scale = 1
    for row in sections["records"]:
        for coordinate in ("x", "y"):
            for coefficient in row[coordinate]:
                scale = lcm(scale, Fraction(coefficient).denominator)
    return scale


def request_for(row: dict, model: dict, sections: dict, model_hash: str, sections_hash: str, config: dict) -> dict:
    parameter = Fraction(row["parameter"])
    denominator = parameter.denominator
    scale = section_scale(sections)
    a = evaluate(model["A"], parameter) * denominator**8 * scale**4
    b = evaluate(model["B"], parameter) * denominator**12 * scale**6
    if a.denominator != 1 or b.denominator != 1:
        raise ArithmeticError("specialized model did not clear denominators")
    points = []
    for source in sections["records"]:
        x = evaluate(source["x"], parameter) * denominator**4 * scale**2
        y = evaluate(source["y"], parameter) * denominator**6 * scale**3
        if x.denominator != 1 or y.denominator != 1 or y*y != x*x*x + a*x + b:
            raise ArithmeticError("specialized section failed exact identity")
        points.append([str(x.numerator), str(y.numerator)])
    return {
        "schema": "elliptic-curves.lazy-mw-search.v1",
        "source": {"model": str(MODEL.relative_to(ROOT)), "model_sha256": model_hash,
                   "sections": str(SECTIONS.relative_to(ROOT)), "sections_sha256": sections_hash,
                   "parameter": row["parameter"], "arm": row["arm"], "score_units": row["score_units"],
                   "denominator_clearing": "X=(s*d^2)^2*x, Y=(s*d^2)^3*y for u=n/d", "section_coordinate_scale": str(scale)},
        "curve": ["0", "0", "0", str(a.numerator), str(b.numerator)], "points": points,
        "metric_gram": sections["height_gram"], "next_holes": config["charts_per_fibre"],
        "height": config["height"], "seconds_per_chart": str(config["seconds_per_chart"]),
        "cvp_node_budget": config["cvp_node_budget"], "prime_bound": config["prime_bound"],
        "extra_primes": [1009, 1013, 1019, 1021, 1031, 1033],
        "policy": {"model_normalization": "raw", "chart_parameterization": "pointed-quartic",
                   "enumeration_backend": "gmp-pointed-sieve", "metric_name": "det1020-generic-height-gram",
                   "quotient_weight": "1", "chart_metric_kind": "raw", "chart_metric_weight": "1", "diversity_window": 1},
        "boundary": "Exact full-MW17 bounded point search. The score and metric order work only. A new independent point is a lower-bound gain; misses, timeouts and score strata do not bound rank. Any rank18 amplification continuation requires a separately declared metric.",
    }


def slug(index: int, parameter: str) -> str:
    return f"{index:04d}-" + parameter.replace("-", "minus-").replace("/", "over-")


def prepare(folder: Path, args: argparse.Namespace) -> dict:
    plan_path = folder / "plan.json"
    if plan_path.exists():
        return read(plan_path)
    model_bytes, section_bytes = MODEL.read_bytes(), SECTIONS.read_bytes()
    model, sections = json.loads(model_bytes), json.loads(section_bytes)
    if model["status"] != "PASS_EXACT_QQ_POINTED_CUBIC_JACOBIAN_INVARIANTS":
        raise RuntimeError("model packet is not exact-certified")
    if sections["status"] != "PASS_EXACT_QQ_SEVENTEEN_SECTION_IDENTITIES_HEIGHTS_AND_SATURATION":
        raise RuntimeError("section packet is not a saturated 17-point basis")
    excluded = []
    if args.exclude_completed_from is not None:
        prior_folder = args.exclude_completed_from.resolve()
        prior_plan_path = prior_folder / "plan.json"
        if not prior_plan_path.exists():
            raise FileNotFoundError(f"missing prior plan: {prior_plan_path}")
        prior_plan = read(prior_plan_path)
        for prior in prior_plan["selected"]:
            prior_case = ROOT / Path(prior["request"]).parent
            if (prior_case / "result.json").exists() or (prior_case / "terminal.json").exists():
                excluded.append(prior["parameter"])
        excluded = sorted(set(excluded))
    config = {"population_height": args.population_height, "ranked": args.ranked, "controls": args.controls,
              "diversity": args.diversity, "charts_per_fibre": args.charts_per_fibre, "height": args.height,
              "seconds_per_chart": args.seconds_per_chart, "cvp_node_budget": args.cvp_node_budget,
              "prime_bound": args.prime_bound, "hours": args.hours, "outer_wall_seconds": args.outer_wall_seconds,
              "rss_bytes": args.rss_bytes, "minimum_free_bytes": args.minimum_free_bytes,
              "exclude_completed_from": (str(args.exclude_completed_from.resolve().relative_to(ROOT)) if args.exclude_completed_from else None),
              "excluded_completed_parameters": excluded}
    candidate_rows = [row for row in population(args.population_height, model) if row["parameter"] not in set(excluded)]
    selected = select(candidate_rows, args.ranked, args.controls, args.diversity)
    items = []
    for index, row in enumerate(selected):
        request = request_for(row, model, sections, sha256(model_bytes).hexdigest(), sha256(section_bytes).hexdigest(), config)
        request_path = folder / "cases" / slug(index, row["parameter"]) / "request.json"
        write_new(request_path, request)
        items.append({**row, "index": index, "request": str(request_path.relative_to(ROOT)), "request_sha256": digest(request_path)})
    plan = {"schema": "elkies-k3.det1020-amplification-scan.v1", "status": "PREPARED_NOT_RUNNING",
            "created_unix": time.time(), "config": config, "selected": items,
            "inputs": {str(path.relative_to(ROOT)): digest(path) for path in (Path(__file__), MODEL, SECTIONS, RUNNER)},
            "score_definition": "sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) over good reductions after exact local scaling",
            "selection": "Ranked local-score arm, score-independent SHA256 control arm, and evenly spaced score-diverse arm, interleaved 4:1:1.",
            "checkpoint_semantics": "Requests and successful results are immutable. Each failed runner invocation receives a terminal receipt and is not retried. STOP prevents new dispatch.",
            "amplification_semantics": "A lower-bound gain is recorded as an amplification candidate but does not automatically alter the rank17 metric or dispatch rank18 work.",
            "boundary": "This is an explicitly bounded discovery search, not an exact-rank, conductor, minimal-model, or absence theorem."}
    write_new(plan_path, plan)
    return plan


def status(folder: Path, plan: dict, run_state: dict | None = None) -> dict:
    rows, gains = [], []
    for item in plan["selected"]:
        case = ROOT / Path(item["request"]).parent
        result_path, terminal_path = case / "result.json", case / "terminal.json"
        if result_path.exists():
            result = read(result_path)
            row = {"parameter": item["parameter"], "arm": item["arm"], "status": result["status"],
                   "charts": len(result["charts"]), "rank_gain": result["certified_rank_gain"]}
            if row["rank_gain"]:
                gains.append({"parameter": item["parameter"], "result": str(result_path.relative_to(ROOT)), "rank_gain": row["rank_gain"]})
        elif terminal_path.exists():
            terminal = read(terminal_path)
            row = {"parameter": item["parameter"], "arm": item["arm"], "status": terminal["status"], "charts": 0, "rank_gain": 0}
        else:
            row = {"parameter": item["parameter"], "arm": item["arm"], "status": "NOT_STARTED", "charts": 0, "rank_gain": 0}
        rows.append(row)
    current = {"schema": "elkies-k3.det1020-amplification-scan-status.v1", "updated_unix": time.time(),
               "rows": rows, "completed_cases": sum(row["status"] != "NOT_STARTED" for row in rows),
               "planned_cases": len(rows), "completed_charts": sum(row["charts"] for row in rows),
               "amplification_candidates": gains, "stopped": (folder / "STOP").exists(), "run_state": run_state}
    checkpoint(folder / "status.json", current)
    if gains:
        checkpoint(folder / "AMPLIFICATION_CANDIDATES.json", {"candidates": gains, "boundary": plan["amplification_semantics"]})
    return current


def run(folder: Path, args: argparse.Namespace) -> None:
    plan = prepare(folder, args)
    run_path = folder / "run-state.json"
    if run_path.exists():
        run_state = read(run_path)
    else:
        hours = plan["config"]["hours"]
        run_state = {"started_unix": time.time(),
                     "deadline_unix": (time.time() + 3600 * hours) if hours else None,
                     "status": "RUNNING"}
        checkpoint(run_path, run_state)
    for item in plan["selected"]:
        deadline = run_state.get("deadline_unix")
        if (deadline is not None and time.time() >= deadline) or (folder / "STOP").exists():
            break
        if shutil.disk_usage(folder).free < plan["config"]["minimum_free_bytes"]:
            run_state["status"] = "STOPPED_LOW_DISK"
            checkpoint(run_path, run_state)
            break
        case = ROOT / Path(item["request"]).parent
        result, terminal = case / "result.json", case / "terminal.json"
        if result.exists() or terminal.exists():
            continue
        command = [sys.executable, str(RUNNER), "--request", str(case / "request.json"), "--output", str(result),
                   "--wall-seconds", str(plan["config"]["outer_wall_seconds"]), "--rss-bytes", str(plan["config"]["rss_bytes"])]
        completed = subprocess.run(command, cwd=ROOT, env={**os.environ, "PYTHONINTMAXSTRDIGITS": "0"})
        if completed.returncode:
            checkpoint(terminal, {"status": "RUNNER_FAILED", "returncode": completed.returncode,
                                  "request_sha256": item["request_sha256"], "command": command,
                                  "boundary": "The runner failed before a complete result. This is neither a point-search miss nor a rank conclusion."})
        status(folder, plan, run_state)
    if run_state.get("status") != "STOPPED_LOW_DISK":
        deadline = run_state.get("deadline_unix")
        run_state["status"] = "STOPPED" if (folder / "STOP").exists() else (
            "DEADLINE_REACHED" if deadline is not None and time.time() >= deadline else "COMPLETE_QUEUE")
    checkpoint(run_path, run_state)
    print(json.dumps(status(folder, plan, run_state), indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("mode", choices=("prepare", "run", "status"))
    result.add_argument("--folder", type=Path, default=DEFAULT)
    result.add_argument("--hours", type=int, default=2,
                        help="Wall-clock deadline; 0 runs the declared queue without a time deadline.")
    result.add_argument("--population-height", type=int, default=64)
    result.add_argument("--ranked", type=int, default=128)
    result.add_argument("--controls", type=int, default=32)
    result.add_argument("--diversity", type=int, default=32)
    result.add_argument("--charts-per-fibre", type=int, default=8)
    result.add_argument("--height", type=int, default=8192)
    result.add_argument("--seconds-per-chart", type=int, default=45)
    result.add_argument("--cvp-node-budget", type=int, default=200000)
    result.add_argument("--prime-bound", type=int, default=997)
    result.add_argument("--outer-wall-seconds", type=int, default=420)
    result.add_argument("--rss-bytes", type=int, default=3221225472)
    result.add_argument("--minimum-free-bytes", type=int, default=6 * 1024**3,
                        help="Do not dispatch another fibre below this much free filesystem space.")
    result.add_argument("--exclude-completed-from", type=Path,
                        help="Prior campaign folder; completed parameters are excluded from this successor.")
    return result


if __name__ == "__main__":
    args = parser().parse_args()
    args.folder = args.folder.resolve()
    if args.hours < 0 or min(args.population_height, args.ranked, args.controls, args.diversity, args.charts_per_fibre, args.height, args.seconds_per_chart, args.cvp_node_budget, args.prime_bound, args.outer_wall_seconds, args.rss_bytes, args.minimum_free_bytes) < 1:
        raise SystemExit("all limits except --hours must be positive; --hours may be zero")
    plan = prepare(args.folder, args)
    if args.mode == "prepare":
        print(json.dumps({"status": plan["status"], "selected_cases": len(plan["selected"]), "folder": str(args.folder)}, indent=2))
    elif args.mode == "status":
        print(json.dumps(status(args.folder, plan, read(args.folder / "run-state.json") if (args.folder / "run-state.json").exists() else None), indent=2, sort_keys=True))
    else:
        run(args.folder, args)
