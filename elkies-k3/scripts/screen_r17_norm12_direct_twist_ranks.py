#!/usr/bin/env python3
"""Vectorized Nagao-style singleton census on a direct norm-12 fibration.

This is a target-selection heuristic only.  It scores every exact quadratic
character already constructed for the direct alternate-Q80 or hidden-103b2
equation.  The rational constant squareclass is retained.  A high score does
not prove a second twist section or a rank jump.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, log
from pathlib import Path
import shlex
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from screen_elkies_2026_quadratic_twist_ranks import (  # noqa: E402
    Candidate,
    DEFAULT_PRIME_BLOCKS,
    build_prime_data,
    parse_prime_block,
    reduce_twist_polynomial,
    square_equivalent_integer_polynomial,
    verify_trace_identity,
)
from search_h92_q12o5867_rootless_nagao import FamilyModel  # noqa: E402


CONFIG = {
    "norm12-orbit-11952": {
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
        "covers": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-singleton-twist-census-v1.json",
        "direct_status": "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS",
        "cover_status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "count": 39147,
        "vector_key": "direct_alternate_w",
    },
    "norm12-orbit-103b2": {
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json",
        "covers": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisections-full-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-singleton-twist-census-v1.json",
        "direct_status": "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS",
        "cover_status": "PASS_EXACT_COMPLETE_103B2_HIDDEN_BISECTION_EQUATIONS",
        "count": 39120,
        "vector_key": "direct_hidden_w",
    },
}


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def direct_model(path: Path, expected_status: str) -> FamilyModel:
    payload = path.read_bytes()
    document = json.loads(payload)
    if document.get("status") != expected_status:
        raise ValueError("direct fibration has the wrong exact status")
    model = document["weierstrass_model"]
    if model.get("degrees_A_B_Delta") != [8, 12, 24]:
        raise ValueError("direct fibration does not have degree pattern (8,12,24)")
    a = tuple(Fraction(value) for value in model["A_coefficients_low_to_high"])
    b = tuple(Fraction(value) for value in model["B_coefficients_low_to_high"])
    if len(a) != 9 or len(b) != 13:
        raise ValueError("direct coefficient arrays have the wrong lengths")
    return FamilyModel(
        source=path.resolve(),
        source_sha256=sha256(payload).hexdigest(),
        a_coefficients=a,
        b_coefficients=b,
        a_degree=8,
        b_degree=12,
        coordinate="direct_norm12_u",
        coefficient_source_keys=(
            "weierstrass_model.A_coefficients_low_to_high",
            "weierstrass_model.B_coefficients_low_to_high",
        ),
    )


def load_candidates(path: Path, configuration: dict) -> list[Candidate]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("cover artifact has the wrong generic schema")
    if payload.get("status") != configuration["cover_status"]:
        raise ValueError("cover artifact has the wrong exact status")
    records = payload.get("bisections")
    if not isinstance(records, list) or len(records) != configuration["count"]:
        raise ValueError("cover artifact has the wrong complete record count")
    candidates = []
    labels = set()
    for record in records:
        label = str(record["label"])
        if label in labels:
            raise ValueError(f"duplicate label {label}")
        labels.add(label)
        branch = record.get("branch")
        if not isinstance(branch, dict) or branch.get("denominator_coefficients") != ["1"]:
            raise ValueError(f"{label}: expected a polynomial branch")
        coefficients = square_equivalent_integer_polynomial(
            branch["numerator_coefficients"]
        )
        if len(coefficients) != 3:
            raise ValueError(f"{label}: branch is not quadratic")
        candidates.append(Candidate(
            kind="singleton",
            key=label,
            masks=(int(record["lattice_orbit_mask"]),),
            coefficients=coefficients,
            forced_twist_rank=1,
            metadata={
                "priority_rank": int(record["priority_rank"]),
                "lattice_orbit_mask": int(record["lattice_orbit_mask"]),
                "orbit_hex": f"0x{int(record['lattice_orbit_mask']):05x}",
                "group_addition_upper_bound": int(
                    record["equation_complexity"]["group_addition_upper_bound"]
                ),
            },
        ))
    return candidates


def score_one_prime(
    candidates: list[Candidate], model: FamilyModel, prime: int, batch_size: int
) -> tuple[np.ndarray, dict]:
    data = build_prime_data(model, prime)
    character = np.asarray(data.characters, dtype=np.int64)
    traces = np.asarray(data.traces, dtype=np.int64)
    parameter = np.arange(prime, dtype=np.int64)
    parameter_squared = parameter * parameter % prime
    result = np.full(len(candidates), np.nan, dtype=np.float64)
    usable_count = 0
    for start in range(0, len(candidates), batch_size):
        stop = min(start + batch_size, len(candidates))
        rows = []
        positions = []
        for position in range(start, stop):
            reduced = reduce_twist_polynomial(candidates[position].coefficients, prime)
            if reduced is None:
                continue
            padded = (*reduced, *(0 for _ in range(3 - len(reduced))))
            rows.append(padded)
            positions.append(position)
        if not rows:
            continue
        coefficients = np.asarray(rows, dtype=np.int64)
        values = (
            coefficients[:, 0, None]
            + coefficients[:, 1, None] * parameter[None, :]
            + coefficients[:, 2, None] * parameter_squared[None, :]
        ) % prime
        trace_sums = character[values] @ traces
        result[np.asarray(positions, dtype=np.int64)] = -trace_sums / prime
        usable_count += len(positions)
    return result, {
        "prime": prime,
        "usable_candidate_count": usable_count,
        "bad_candidate_count": len(candidates) - usable_count,
        "singular_fibre_count": data.singular_fibre_count,
        "trace_sha256": sha256(",".join(map(str, data.traces)).encode("ascii")).hexdigest(),
    }


def quantile(values: np.ndarray, fraction: float) -> float | None:
    if not len(values):
        return None
    return float(np.quantile(values, fraction))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", choices=tuple(CONFIG), default="norm12-orbit-11952")
    parser.add_argument("--direct", type=Path)
    parser.add_argument("--covers", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--prime-block", type=parse_prime_block, action="append")
    parser.add_argument("--top", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--label", action="append", default=[])
    arguments = parser.parse_args()
    if arguments.top < 1 or arguments.batch_size < 1:
        raise ValueError("top and batch-size must be positive")
    configuration = CONFIG[arguments.source_label]
    direct_path = arguments.direct or configuration["direct"]
    covers_path = arguments.covers or configuration["covers"]
    output = arguments.output or configuration["output"]
    blocks = tuple(arguments.prime_block or DEFAULT_PRIME_BLOCKS)
    primes = [prime for block in blocks for prime in block]
    if len(primes) != len(set(primes)):
        raise ValueError("prime blocks must be disjoint")

    started = perf_counter()
    model = direct_model(direct_path, configuration["direct_status"])
    all_candidates = load_candidates(covers_path, configuration)
    requested = set(arguments.label)
    candidates = [item for item in all_candidates if not requested or item.key in requested]
    if requested != {item.key for item in candidates} and requested:
        missing = sorted(requested - {item.key for item in candidates})
        raise ValueError(f"unknown labels: {missing}")

    local_columns = []
    prime_records = []
    usable_primes = []
    skipped_primes = []
    for prime in primes:
        try:
            values, record = score_one_prime(
                candidates, model, prime, arguments.batch_size
            )
        except ZeroDivisionError as error:
            skipped_primes.append({"prime": prime, "reason": str(error)})
            continue
        local_columns.append(values)
        prime_records.append(record)
        usable_primes.append(prime)
        print(
            f"DIRECTTWISTCENSUS|source={arguments.source_label}|p={prime}"
            f"|usable={record['usable_candidate_count']}/{len(candidates)}",
            flush=True,
        )
    if not local_columns:
        raise ArithmeticError("no usable prime remains")
    local = np.column_stack(local_columns)
    prime_to_column = {prime: index for index, prime in enumerate(usable_primes)}
    block_scores = np.full((len(candidates), len(blocks)), np.nan, dtype=np.float64)
    block_usable_counts = np.zeros((len(candidates), len(blocks)), dtype=np.int64)
    for block_index, block in enumerate(blocks):
        columns = [prime_to_column[prime] for prime in block if prime in prime_to_column]
        if not columns:
            continue
        weights = np.asarray([log(usable_primes[column]) for column in columns])
        values = local[:, columns]
        finite = np.isfinite(values)
        weighted = np.where(finite, values, 0.0) @ weights
        weight_sums = finite @ weights
        block_scores[:, block_index] = np.divide(
            weighted,
            weight_sums,
            out=np.full(len(candidates), np.nan),
            where=weight_sums > 0,
        )
        block_usable_counts[:, block_index] = finite.sum(axis=1)
    complete = np.isfinite(block_scores).all(axis=1)
    weakest = np.where(complete, np.nanmin(block_scores, axis=1), np.nan)
    means = np.where(complete, np.nanmean(block_scores, axis=1), np.nan)
    order = sorted(
        range(len(candidates)),
        key=lambda index: (
            0 if isfinite(float(weakest[index])) else 1,
            -float(weakest[index]) if isfinite(float(weakest[index])) else 0.0,
            -float(means[index]) if isfinite(float(means[index])) else 0.0,
            candidates[index].key,
        ),
    )
    top = []
    for rank, index in enumerate(order[: arguments.top], start=1):
        candidate = candidates[index]
        top.append({
            "census_rank": rank,
            "label": candidate.key,
            **candidate.metadata,
            "block_scores": [float(value) for value in block_scores[index]],
            "weakest_block_score": float(weakest[index]),
            "mean_block_score": float(means[index]),
            "weakest_block_excess_over_forced_rank": float(weakest[index] - 1.0),
            "usable_primes_by_block": [int(value) for value in block_usable_counts[index]],
        })
    finite_weakest = weakest[np.isfinite(weakest)]
    finite_means = means[np.isfinite(means)]
    ledger = sha256()
    for index in sorted(range(len(candidates)), key=lambda item: candidates[item].key):
        ledger.update(candidates[index].key.encode("utf-8"))
        ledger.update(b"|")
        ledger.update("|".join(format(float(value), ".17g") for value in block_scores[index]).encode("ascii"))
        ledger.update(b"\n")

    self_test = None
    for prime in usable_primes:
        data = build_prime_data(model, prime)
        for candidate in candidates:
            if reduce_twist_polynomial(candidate.coefficients, prime) is not None:
                self_test = verify_trace_identity(model, candidate, data)
                break
        if self_test is not None:
            break
    if self_test is None:
        raise ArithmeticError("trace identity self-test found no good candidate")

    result = {
        "schema": "elkies-k3.r17-norm12-direct-singleton-twist-census.v1",
        "status": "PASS_BOUNDED_HEURISTIC_DIRECT_SINGLETON_TWIST_CENSUS",
        "source_label": arguments.source_label,
        "candidate_count": len(candidates),
        "complete_frame_candidate_count": len(all_candidates),
        "forced_twist_rank": 1,
        "prime_blocks_requested": [list(block) for block in blocks],
        "usable_primes": usable_primes,
        "skipped_model_primes": skipped_primes,
        "prime_records": prime_records,
        "trace_identity_self_test": self_test,
        "complete_score_ledger_sha256": ledger.hexdigest(),
        "score_distribution": {
            "finite_candidate_count": int(len(finite_weakest)),
            "weakest_block_quantiles": {
                label: quantile(finite_weakest, fraction)
                for label, fraction in (
                    ("q00", 0.0), ("q25", 0.25), ("q50", 0.5),
                    ("q75", 0.75), ("q90", 0.9), ("q99", 0.99), ("q100", 1.0),
                )
            },
            "mean_block_quantiles": {
                label: quantile(finite_means, fraction)
                for label, fraction in (
                    ("q00", 0.0), ("q25", 0.25), ("q50", 0.5),
                    ("q75", 0.75), ("q90", 0.9), ("q99", 0.99), ("q100", 1.0),
                )
            },
            "weakest_threshold_counts": {
                format(threshold, ".2f"): int(np.sum(finite_weakest >= threshold))
                for threshold in (0.5, 0.75, 1.0, 1.25, 1.5)
            },
        },
        "top": top,
        "method": {
            "trace_identity": "a_p(E^q_u)=chi_p(q(u))*a_p(E_u)",
            "local_score": "-(1/p)*sum_{u in F_p} a_p(E^q_u)",
            "ranking": "descending weakest logarithmically weighted prime-block score",
            "constant_squareclass_preserved": True,
            "numpy_version": np.__version__,
        },
        "inputs": {
            relative(direct_path): digest(direct_path),
            relative(covers_path): digest(covers_path),
        },
        "proof_boundary": (
            "This complete finite-prime census is a target-selection heuristic. It proves "
            "neither a second twist section nor a Mordell--Weil rank lower bound beyond the "
            "one exact section already attached to each quadratic character."
        ),
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    leader = top[0]
    print(
        f"DIRECTTWISTCENSUS|status=PASS_HEURISTIC|source={arguments.source_label}"
        f"|candidates={len(candidates)}|leader={leader['label']}"
        f"|weakest={leader['weakest_block_score']:.6f}|mean={leader['mean_block_score']:.6f}"
        f"|seconds={result['runtime_seconds']:.3f}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
