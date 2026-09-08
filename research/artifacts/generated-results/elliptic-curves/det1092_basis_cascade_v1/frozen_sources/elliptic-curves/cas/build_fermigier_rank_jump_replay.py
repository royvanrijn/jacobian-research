#!/usr/bin/env python3
"""Replay certified Fermigier quotient jumps in a complete frozen score box.

This is a retrospective retrieval benchmark.  It applies the four score
orderings declared by the archived leakage-controlled global scan to every
primitive ``T=a/b`` in its original box.  The two certified anchors are the
only labels; every other fibre is censored, never assigned quotient rank zero.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import log
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ARCHIVED_GLOBAL = (
    ROOT
    / "archive/elliptic-curves/artifacts/generated-results"
    / "elliptic_fermigier_global.json"
)
SCANNER = CAS / "scan_fermigier_rank_jump_replay.cpp"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "fermigier_rank_jump_replay_v1.json"
)
SCORE_SCALE = 10**9
COMPOSITE_POWER_DIVISOR = 16
EXPECTED_POPULATION = 60_815_684
ANCHORS = (
    {
        "id": "fermigier-E22",
        "canonical_parameter_u": "19754/39",
        "normalized_parameter_T": (39508, 39),
        "generic_rank": 12,
        "certified_rank_lower_bound": 22,
        "exceptional_quotient_rank_lower_bound": 10,
        "certificate": ROOT
        / "artifacts/generated-results/elliptic-curves"
        / "elliptic_fermigier_rank22_points.json",
    },
    {
        "id": "fermigier-rank20-near-miss",
        "canonical_parameter_u": "28917/20",
        "normalized_parameter_T": (28917, 10),
        "generic_rank": 12,
        "certified_rank_lower_bound": 20,
        "exceptional_quotient_rank_lower_bound": 8,
        "certificate": ROOT
        / "artifacts/generated-results/elliptic-curves"
        / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json",
    },
)
METRICS = (
    "discovery_rank",
    "discovery_composite",
    "held_rank",
    "held_composite",
)
BUDGETS = (10, 100, 1_000, 8_000, 16_133, 100_000)


import sys

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from ek_k3 import primes_up_to  # noqa: E402
from fermigier_mestre import (  # noqa: E402
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
)


@dataclass(frozen=True)
class ModularTable:
    prime: int
    rank_weights: tuple[int, ...]
    power_weights: tuple[int, ...]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def polynomial_mod(value: int, modulus: int) -> int:
    answer = 0
    for coefficient in reversed(DISCRIMINANT_FACTOR_COEFFICIENTS):
        answer = (answer * value + coefficient) % modulus
    return answer


def build_tables(primes: Sequence[int]) -> tuple[ModularTable, ...]:
    answer = []
    for prime in primes:
        rank_weights = []
        for residue in range(prime):
            local = FermigierMestreFamily.local_data(residue, prime)
            term = (
                0.0
                if not local.good_reduction
                else (2 - local.trace) / local.point_count * log(prime)
            )
            rank_weights.append(round(term * SCORE_SCALE))
        # The last entries represent projective infinity / a denominator prime.
        rank_weights.append(0)
        modulus = prime * prime
        power_weights = [0] * (modulus + 1)
        saving = round(log(prime) * SCORE_SCALE)
        for residue in range(prime):
            if polynomial_mod(residue, prime) != 0:
                continue
            for lift_digit in range(prime):
                lifted = residue + lift_digit * prime
                if polynomial_mod(lifted, modulus) == 0:
                    power_weights[lifted] = saving
        answer.append(
            ModularTable(prime, tuple(rank_weights), tuple(power_weights))
        )
    return tuple(answer)


def scanner_input(*bands: Sequence[ModularTable]) -> str:
    lines: list[str] = []
    for band in bands:
        lines.append(str(len(band)))
        for table in band:
            lines.append(
                " ".join([str(table.prime), *(str(v) for v in table.rank_weights)])
            )
            lines.append(" ".join(str(v) for v in table.power_weights))
    return "\n".join(lines) + "\n"


def validate_frozen_protocol(document: dict[str, object]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    population = document["population"]
    features = document["leakage_free_features"]
    if population["definition"] != (
        "every primitive 0<=a<=A, 1<=b<=B; T=a/b; T and -T identified"
    ):
        raise AssertionError("the archived population definition changed")
    if (
        population["a_max"],
        population["b_max"],
        population["primitive_pairs_enumerated"],
    ) != (100_000, 1_000, EXPECTED_POPULATION):
        raise AssertionError("the archived Fermigier box changed")
    if features["rank_formula"] != (
        "sum good ((2-a_p)/(p+1-a_p))*log(p); bad and denominator primes contribute zero"
    ):
        raise AssertionError("the frozen rank score changed")
    if features["power_formula"] != (
        "sum log(p) for exact affine congruences p^2 | H(a/b); denominator primes contribute zero"
    ):
        raise AssertionError("the frozen power score changed")
    if features["composite_power_divisor"] != COMPOSITE_POWER_DIVISOR:
        raise AssertionError("the frozen composite divisor changed")
    discovery = tuple(map(int, features["discovery_primes"]))
    held = tuple(map(int, features["held_primes"]))
    if set(discovery) & set(held) or features["bands_disjoint"] is not True:
        raise AssertionError("the score bands are not disjoint")
    return discovery, held


def run_scan(
    discovery: Sequence[ModularTable],
    held: Sequence[ModularTable],
    *,
    compile_timeout: float,
    scan_timeout: float,
) -> tuple[int, list[dict[str, int]], dict[str, str]]:
    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="fermigier-jump-replay-") as directory:
        executable = Path(directory) / "scan"
        subprocess.run(
            [compiler, "-std=c++17", "-O3", "-DNDEBUG", str(SCANNER), "-o", str(executable)],
            check=True,
            capture_output=True,
            text=True,
            timeout=compile_timeout,
        )
        command = [
            str(executable),
            "100000",
            "1000",
            "39508",
            "39",
            "28917",
            "10",
            str(EXPECTED_POPULATION),
        ]
        completed = subprocess.run(
            command,
            input=scanner_input(discovery, held),
            capture_output=True,
            text=True,
            check=True,
            timeout=scan_timeout,
        )
    lines = completed.stdout.splitlines()
    if not lines or not lines[0].startswith("SUMMARY\t"):
        raise AssertionError("the scanner omitted its summary")
    enumerated = int(lines[0].split("\t")[1])
    rows = []
    field_names = (
        "numerator",
        "denominator",
        "discovery_rank_scaled",
        "discovery_power_scaled",
        "held_rank_scaled",
        "held_power_scaled",
        "discovery_rank_position",
        "discovery_rank_equal_score_count",
        "discovery_composite_position",
        "discovery_composite_equal_score_count",
        "held_rank_position",
        "held_rank_equal_score_count",
        "held_composite_position",
        "held_composite_equal_score_count",
    )
    for line in lines[1:]:
        fields = line.split("\t")
        if fields[0] != "TARGET" or len(fields) != len(field_names) + 1:
            raise AssertionError("the scanner target row changed shape")
        rows.append(dict(zip(field_names, map(int, fields[1:]), strict=True)))
    if len(rows) != len(ANCHORS):
        raise AssertionError("the scanner returned the wrong number of anchors")
    return enumerated, rows, {
        "compiler_command": Path(compiler).name,
        "language_standard": "C++17",
        "optimization": "-O3 -DNDEBUG",
    }


def score_record(row: dict[str, int], metric: str, population_count: int) -> dict[str, object]:
    prefix = "discovery" if metric.startswith("discovery") else "held"
    rank_scaled = row[f"{prefix}_rank_scaled"]
    power_scaled = row[f"{prefix}_power_scaled"]
    score_scaled = (
        rank_scaled
        if metric.endswith("rank")
        else rank_scaled + power_scaled // COMPOSITE_POWER_DIVISOR
    )
    position = row[f"{metric}_position"]
    return {
        "score_scaled_by_1e9": score_scaled,
        "score": score_scaled / SCORE_SCALE,
        "rank_position_one_based": position,
        "equal_primary_score_count_before_tie_break": row[
            f"{metric}_equal_score_count"
        ],
        "top_fraction": position / population_count,
        "percentile_from_bottom": 1 - (position - 1) / population_count,
        "tie_break": "higher score, then lower projective height, denominator, numerator",
    }


def retrieval_metrics(anchors: Sequence[dict[str, object]]) -> dict[str, object]:
    thresholds = (8, 10)
    answer: dict[str, object] = {}
    for metric in METRICS:
        by_threshold = {}
        for threshold in thresholds:
            positives = [
                anchor
                for anchor in anchors
                if anchor["exceptional_quotient_rank_lower_bound"] >= threshold
            ]
            positions = [
                anchor["scores"][metric]["rank_position_one_based"]
                for anchor in positives
            ]
            by_threshold[f"gain_at_least_{threshold}"] = {
                "positive_count": len(positives),
                "positions": positions,
                "mean_reciprocal_rank": sum(1 / value for value in positions)
                / len(positions),
                "recall_at_budget": {
                    str(budget): sum(value <= budget for value in positions)
                    / len(positions)
                    for budget in BUDGETS
                },
            }
        answer[metric] = by_threshold
    return answer


def build_payload(*, compile_timeout: float, scan_timeout: float) -> dict[str, object]:
    archived = json.loads(ARCHIVED_GLOBAL.read_text())
    discovery_primes, held_primes = validate_frozen_protocol(archived)
    discovery = build_tables(discovery_primes)
    held = build_tables(held_primes)
    enumerated, rows, timing = run_scan(
        discovery,
        held,
        compile_timeout=compile_timeout,
        scan_timeout=scan_timeout,
    )
    anchors = []
    for specification, row in zip(ANCHORS, rows, strict=True):
        if (row["numerator"], row["denominator"]) != specification[
            "normalized_parameter_T"
        ]:
            raise AssertionError("scanner anchors returned out of order")
        anchors.append(
            {
                **{
                    key: value
                    for key, value in specification.items()
                    if key not in {"certificate", "normalized_parameter_T"}
                },
                "normalized_parameter_T": str(
                    Fraction(*specification["normalized_parameter_T"])
                ),
                "certificate": {
                    "path": display(specification["certificate"]),
                    "sha256": digest(specification["certificate"]),
                },
                "scores": {
                    metric: score_record(row, metric, enumerated)
                    for metric in METRICS
                },
            }
        )
    return {
        "schema": "elliptic-curves.fermigier-rank-jump-replay.v1",
        "status": "PASS_COMPLETE_RETROSPECTIVE_FROZEN_SCORE_REPLAY",
        "target": "retrieval of certified escape from the generic rank-12 lattice",
        "evaluation_role": "retrospective_frozen_rule_replay_not_prospective_holdout",
        "population": {
            "definition": archived["population"]["definition"],
            "a_max": 100_000,
            "b_max": 1_000,
            "primitive_parameter_count": enumerated,
            "certified_positive_count": len(anchors),
            "censored_unlabelled_count": enumerated - len(anchors),
            "negative_label_count": 0,
            "censoring_rule": (
                "unlabelled fibres are never assigned quotient rank zero and do not "
                "enter precision, specificity, ROC, or classification accuracy"
            ),
        },
        "frozen_feature_specification": {
            "score_scale": SCORE_SCALE,
            "discovery_primes": list(discovery_primes),
            "held_primes": list(held_primes),
            "bands_disjoint": True,
            "rank_formula": archived["leakage_free_features"]["rank_formula"],
            "power_formula": archived["leakage_free_features"]["power_formula"],
            "composite_formula": "rank_scaled + floor(power_scaled/16)",
            "fit_to_rank_jump_labels": False,
            "point_or_rank_results_used_by_scanner": False,
        },
        "anchors": anchors,
        "retrieval": retrieval_metrics(anchors),
        "interpretation": {
            "valid": (
                "exact within-box ranks and recall for four predeclared score orderings, "
                "using only the two certified positive anchors"
            ),
            "invalid": (
                "this is not a prospective holdout: E22 was an explicitly known "
                "calibration when the historical global score was designed"
            ),
            "next_gate": (
                "use only a score ordering with acceptable positive retrieval to rank "
                "new fibres, then require a residual Selmer/quotient certificate"
            ),
        },
        "inputs": {
            display(ARCHIVED_GLOBAL): digest(ARCHIVED_GLOBAL),
            display(SCANNER): digest(SCANNER),
            "elliptic-curves/cas/fermigier_mestre.py": digest(
                CAS / "fermigier_mestre.py"
            ),
        },
        "software": timing,
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 "
            "elliptic-curves/cas/build_fermigier_rank_jump_replay.py"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--compile-timeout", type=float, default=60.0)
    parser.add_argument("--scan-timeout", type=float, default=180.0)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload(
        compile_timeout=args.compile_timeout, scan_timeout=args.scan_timeout
    )
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing replay artifact: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    positions = ",".join(
        f"{anchor['id']}:{anchor['scores']['held_rank']['rank_position_one_based']}"
        for anchor in payload["anchors"]
    )
    print(
        "FERMIGIERRANKJUMPREPLAY|metric=held_rank|positions="
        f"{positions}|population={payload['population']['primitive_parameter_count']}"
        f"|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
