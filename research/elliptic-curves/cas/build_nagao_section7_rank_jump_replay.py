#!/usr/bin/env python3
"""Replay Nagao's rank-20 escape anchor in the complete frozen global box.

This is a retrospective external-positive retrieval benchmark.  It applies
the two label-blind score bands from the archived global section-7 scan to
every primitive positive ``T=a/b`` in that scan's rectangle.  The known
rank-20 fibre is the sole positive label; every other fibre remains censored.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
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
    / "elliptic_nagao_section7_global.json"
)
SCANNER = CAS / "scan_nagao_section7_rank_jump_replay.cpp"
CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_nagao_rank20_t5081_rank20_certificate.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "nagao_section7_rank_jump_replay_v1.json"
)
SCORE_SCALE = 10**12
EXPECTED_POPULATION = 18_244_819
TARGET = (5081, 47)
BUDGETS = (10, 100, 1_000, 8_000, 16_133, 100_000)

import sys

sys.path[:0] = [str(CAS)]

from ek_k3 import primes_up_to  # noqa: E402
from search_nagao_rank20_t5081_neighborhood import build_residue_tables  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def frozen_prime_bands() -> tuple[tuple[int, ...], tuple[int, ...]]:
    primes = tuple(primes_up_to(397))
    return (
        tuple(prime for prime in primes if 5 <= prime <= 199),
        tuple(prime for prime in primes if 211 <= prime <= 397),
    )


def validate_protocol(document: dict[str, object]) -> None:
    population = document["population"]
    scores = document["leakage_free_scoring"]
    if document.get("schema_version") != 1:
        raise AssertionError("the archived Nagao global schema changed")
    if population["definition"] != (
        "every positive primitive (a,b) with 1<=a<=A_MAX and 1<=b<=B_MAX, "
        "before bounded training-only frontier retention"
    ):
        raise AssertionError("the archived population definition changed")
    if (
        population["a_max"],
        population["b_max"],
        population["primitive_pairs_enumerated"],
        population["score_scale"],
    ) != (30_000, 1_000, EXPECTED_POPULATION, SCORE_SCALE):
        raise AssertionError("the archived Nagao population changed")
    if scores["training_prime_band"] != [5, 199]:
        raise AssertionError("the frozen training band changed")
    if scores["heldout_validation_prime_band"] != [211, 397]:
        raise AssertionError("the frozen validation band changed")
    if scores["bands_pairwise_disjoint"] is not True:
        raise AssertionError("the frozen bands stopped being disjoint")


def scanner_input(bands: Sequence[Sequence[object]]) -> str:
    lines: list[str] = []
    for tables in bands:
        lines.append(str(len(tables)))
        for table in tables:
            weights = [
                round(float(symbol.contribution) * SCORE_SCALE) for symbol in table
            ]
            prime = len(table) - 1
            lines.append(" ".join([str(prime), *(str(value) for value in weights)]))
    return "\n".join(lines) + "\n"


def run_scan(*, compile_timeout: float, scan_timeout: float) -> tuple[dict[str, int], dict[str, str]]:
    training_primes, validation_primes = frozen_prime_bands()
    tables = build_residue_tables(max(validation_primes))
    training = tuple(tables[prime] for prime in training_primes)
    validation = tuple(tables[prime] for prime in validation_primes)
    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="nagao-section7-jump-replay-") as directory:
        executable = Path(directory) / "scan"
        subprocess.run(
            [compiler, "-std=c++17", "-O3", "-DNDEBUG", str(SCANNER), "-o", str(executable)],
            check=True,
            capture_output=True,
            text=True,
            timeout=compile_timeout,
        )
        completed = subprocess.run(
            [str(executable), "30000", "1000", "5081", "47", str(EXPECTED_POPULATION)],
            input=scanner_input((training, validation)),
            check=True,
            capture_output=True,
            text=True,
            timeout=scan_timeout,
        )
    lines = completed.stdout.splitlines()
    if len(lines) != 2 or lines[0] != f"SUMMARY\t{EXPECTED_POPULATION}":
        raise AssertionError("the replay scanner summary changed")
    fields = lines[1].split("\t")
    names = (
        "numerator",
        "denominator",
        "training_scaled",
        "validation_scaled",
        "training_position",
        "training_equal_score_count",
        "validation_position",
        "validation_equal_score_count",
    )
    if fields[0] != "TARGET" or len(fields) != len(names) + 1:
        raise AssertionError("the replay scanner target row changed")
    row = dict(zip(names, map(int, fields[1:]), strict=True))
    if (row["numerator"], row["denominator"]) != TARGET:
        raise AssertionError("the replay scanner returned the wrong target")
    return row, {
        "compiler_command": Path(compiler).name,
        "language_standard": "C++17",
        "optimization": "-O3 -DNDEBUG",
    }


def score_record(row: dict[str, int], band: str) -> dict[str, object]:
    position = row[f"{band}_position"]
    score = row[f"{band}_scaled"]
    return {
        "score_scaled_by_1e12": score,
        "score": score / SCORE_SCALE,
        "rank_position_one_based": position,
        "equal_primary_score_count_before_tie_break": row[
            f"{band}_equal_score_count"
        ],
        "top_fraction": position / EXPECTED_POPULATION,
        "percentile_from_bottom": 1 - (position - 1) / EXPECTED_POPULATION,
        "tie_break": "higher score, then lower projective height, numerator, denominator",
    }


def build_payload(*, compile_timeout: float, scan_timeout: float) -> dict[str, object]:
    archived = json.loads(ARCHIVED_GLOBAL.read_text())
    validate_protocol(archived)
    row, software = run_scan(
        compile_timeout=compile_timeout, scan_timeout=scan_timeout
    )
    training_primes, validation_primes = frozen_prime_bands()
    scores = {band: score_record(row, band) for band in ("training", "validation")}
    retrieval = {}
    for band, record in scores.items():
        position = record["rank_position_one_based"]
        retrieval[band] = {
            "gain_at_least_8": {
                "positive_count": 1,
                "positions": [position],
                "mean_reciprocal_rank": 1 / position,
                "recall_at_budget": {
                    str(budget): float(position <= budget) for budget in BUDGETS
                },
            }
        }
    return {
        "schema": "elliptic-curves.nagao-section7-rank-jump-replay.v1",
        "status": "PASS_COMPLETE_RETROSPECTIVE_EXTERNAL_POSITIVE_REPLAY",
        "target": "retrieval of certified escape from the generic rank-12 lattice",
        "evaluation_role": "retrospective_external_positive_replay_not_prospective_holdout",
        "population": {
            "definition": archived["population"]["definition"],
            "a_max": 30_000,
            "b_max": 1_000,
            "primitive_parameter_count": EXPECTED_POPULATION,
            "certified_positive_count": 1,
            "censored_unlabelled_count": EXPECTED_POPULATION - 1,
            "negative_label_count": 0,
            "censoring_rule": (
                "unlabelled fibres are never assigned quotient rank zero and do not "
                "enter precision, specificity, ROC, or classification accuracy"
            ),
        },
        "frozen_feature_specification": {
            "score_scale": SCORE_SCALE,
            "training_primes": list(training_primes),
            "validation_primes": list(validation_primes),
            "bands_disjoint": True,
            "rank_formula": (
                "sum good ((2-a_p)/(p+1-a_p))*log(p); bad and denominator "
                "primes contribute zero"
            ),
            "fit_to_rank_jump_labels": False,
            "point_or_rank_results_used_by_scanner": False,
        },
        "anchors": [
            {
                "id": "nagao-section7-rank20-t5081",
                "constructor_parameter_T": "5081/47",
                "paper_parameter_t": "5081/94",
                "generic_rank": 12,
                "certified_rank_lower_bound": 20,
                "exceptional_quotient_rank_lower_bound": 8,
                "certificate": {"path": display(CERTIFICATE), "sha256": digest(CERTIFICATE)},
                "scores": scores,
            }
        ],
        "retrieval": retrieval,
        "interpretation": {
            "valid": (
                "exact within-box ranks and recall for the two frozen label-blind "
                "score bands, using one externally known certified positive"
            ),
            "invalid": (
                "this is not a prospective holdout and supplies no negative labels; "
                "the fibre was historically known before this replay"
            ),
            "selection_boundary": (
                "the archived global search explicitly excluded 5081/47 before "
                "frontier selection, but this full-population diagnostic ranks it"
            ),
        },
        "inputs": {
            display(ARCHIVED_GLOBAL): digest(ARCHIVED_GLOBAL),
            display(SCANNER): digest(SCANNER),
            display(CERTIFICATE): digest(CERTIFICATE),
            "elliptic-curves/cas/search_nagao_rank20_t5081_neighborhood.py": digest(
                CAS / "search_nagao_rank20_t5081_neighborhood.py"
            ),
        },
        "software": software,
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 "
            "elliptic-curves/cas/build_nagao_section7_rank_jump_replay.py"
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
    anchor = payload["anchors"][0]
    print(
        "NAGAOSECTION7RANKJUMPREPLAY|training="
        f"{anchor['scores']['training']['rank_position_one_based']}|validation="
        f"{anchor['scores']['validation']['rank_position_one_based']}|population="
        f"{payload['population']['primitive_parameter_count']}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
