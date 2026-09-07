#!/usr/bin/env python3
"""Rank product twists among a fixed prefix of direct singleton candidates.

The input singleton census supplies a deterministic candidate prefix.  Every
unordered product in that prefix is scored from exact finite-field trace and
quadratic-character tables on disjoint prime blocks.  This is only a Nagao-
style target selector: it neither constructs a product-character section nor
proves a Mordell--Weil rank gain.
"""

from __future__ import annotations

import argparse
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
    multiply_integer_polynomials,
    parse_prime_block,
    reduce_twist_polynomial,
)
from screen_r17_norm12_direct_twist_ranks import (  # noqa: E402
    CONFIG,
    digest,
    direct_model,
    load_candidates,
    quantile,
    relative,
)


DEFAULT_SINGLETON_CENSUS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-singleton-twist-census-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-product-twist-census-top200-v1.json"
)


def score_one_prime(candidates, model, prime, batch_size):
    data = build_prime_data(model, prime)
    character = np.asarray(data.characters, dtype=np.int64)
    traces = np.asarray(data.traces, dtype=np.int64)
    parameter = np.arange(prime, dtype=np.int64)
    powers = np.vstack(
        [np.ones(prime, dtype=np.int64), parameter]
        + [np.power(parameter, degree, dtype=np.int64) % prime for degree in range(2, 5)]
    )
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
            rows.append((*reduced, *(0 for unused in range(5 - len(reduced)))))
            positions.append(position)
        if not rows:
            continue
        coefficient_array = np.asarray(rows, dtype=np.int64)
        values = coefficient_array @ powers % prime
        trace_sums = character[values] @ traces
        result[np.asarray(positions, dtype=np.int64)] = -trace_sums / prime
        usable_count += len(positions)
    return result, {
        "prime": prime,
        "usable_candidate_count": usable_count,
        "bad_candidate_count": len(candidates) - usable_count,
        "singular_fibre_count": data.singular_fibre_count,
        "trace_sha256": sha256(
            ",".join(map(str, data.traces)).encode("ascii")
        ).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", choices=tuple(CONFIG), default="norm12-orbit-11952")
    parser.add_argument("--direct", type=Path)
    parser.add_argument("--covers", type=Path)
    parser.add_argument("--singleton-census", type=Path, default=DEFAULT_SINGLETON_CENSUS)
    parser.add_argument("--prefix", type=int, default=200)
    parser.add_argument("--prime-block", type=parse_prime_block, action="append")
    parser.add_argument("--top", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    if arguments.prefix < 2 or arguments.top < 1 or arguments.batch_size < 1:
        parser.error("prefix must be at least two; top and batch size must be positive")

    started = perf_counter()
    configuration = CONFIG[arguments.source_label]
    direct_path = arguments.direct or configuration["direct"]
    covers_path = arguments.covers or configuration["covers"]
    model = direct_model(direct_path, configuration["direct_status"])
    singleton_candidates = load_candidates(covers_path, configuration)
    singleton_by_label = {candidate.key: candidate for candidate in singleton_candidates}
    census = json.loads(arguments.singleton_census.read_text())
    if census.get("schema") != "elkies-k3.r17-norm12-direct-singleton-twist-census.v1":
        raise ValueError("unexpected singleton census schema")
    if census.get("source_label") != arguments.source_label:
        raise ValueError("singleton census source disagrees with this product run")
    prefix_rows = census["top"][: arguments.prefix]
    if len(prefix_rows) != arguments.prefix:
        raise ValueError("singleton census did not retain the requested prefix")
    prefix = [singleton_by_label[row["label"]] for row in prefix_rows]

    products = []
    for left_index, left in enumerate(prefix):
        for right in prefix[left_index + 1 :]:
            products.append(
                Candidate(
                    kind="direct_product",
                    key=f"{left.key}:{right.key}",
                    masks=(left.masks[0], right.masks[0]),
                    coefficients=multiply_integer_polynomials(
                        left.coefficients, right.coefficients
                    ),
                    forced_twist_rank=0,
                    metadata={
                        "left_label": left.key,
                        "right_label": right.key,
                        "left_mask": left.masks[0],
                        "right_mask": right.masks[0],
                    },
                )
            )

    blocks = tuple(arguments.prime_block or DEFAULT_PRIME_BLOCKS)
    primes = [prime for block in blocks for prime in block]
    if len(primes) != len(set(primes)):
        raise ValueError("prime blocks must be disjoint")
    local_columns = []
    prime_records = []
    usable_primes = []
    for prime in primes:
        values, record = score_one_prime(
            products, model, prime, arguments.batch_size
        )
        local_columns.append(values)
        prime_records.append(record)
        usable_primes.append(prime)
        print(
            f"DIRECTPRODUCTCENSUS|p={prime}|usable={record['usable_candidate_count']}/{len(products)}",
            flush=True,
        )

    local = np.column_stack(local_columns)
    prime_to_column = {prime: index for index, prime in enumerate(usable_primes)}
    block_scores = np.full((len(products), len(blocks)), np.nan, dtype=np.float64)
    block_usable_counts = np.zeros((len(products), len(blocks)), dtype=np.int64)
    for block_index, block in enumerate(blocks):
        columns = [prime_to_column[prime] for prime in block]
        weights = np.asarray([log(prime) for prime in block])
        values = local[:, columns]
        finite = np.isfinite(values)
        weighted = np.where(finite, values, 0.0) @ weights
        weight_sums = finite @ weights
        block_scores[:, block_index] = np.divide(
            weighted,
            weight_sums,
            out=np.full(len(products), np.nan),
            where=weight_sums > 0,
        )
        block_usable_counts[:, block_index] = finite.sum(axis=1)
    complete = np.isfinite(block_scores).all(axis=1)
    weakest = np.where(complete, np.nanmin(block_scores, axis=1), np.nan)
    means = np.where(complete, np.nanmean(block_scores, axis=1), np.nan)
    order = sorted(
        range(len(products)),
        key=lambda index: (
            0 if isfinite(float(weakest[index])) else 1,
            -float(weakest[index]) if isfinite(float(weakest[index])) else 0,
            -float(means[index]) if isfinite(float(means[index])) else 0,
            products[index].key,
        ),
    )
    top = []
    for rank, index in enumerate(order[: arguments.top], start=1):
        product_candidate = products[index]
        top.append(
            {
                "census_rank": rank,
                "key": product_candidate.key,
                **product_candidate.metadata,
                "block_scores": [float(value) for value in block_scores[index]],
                "weakest_block_score": float(weakest[index]),
                "mean_block_score": float(means[index]),
                "usable_primes_by_block": [
                    int(value) for value in block_usable_counts[index]
                ],
            }
        )
    finite_weakest = weakest[np.isfinite(weakest)]
    finite_means = means[np.isfinite(means)]
    result = {
        "schema": "elkies-k3.r17-norm12-direct-product-twist-census.v1",
        "status": "PASS_BOUNDED_HEURISTIC_DIRECT_PRODUCT_TWIST_CENSUS",
        "source_label": arguments.source_label,
        "singleton_prefix_count": len(prefix),
        "product_candidate_count": len(products),
        "prime_blocks": [list(block) for block in blocks],
        "prime_records": prime_records,
        "top": top,
        "score_distribution": {
            "finite_candidate_count": len(finite_weakest),
            "weakest_block_quantiles": {
                label: quantile(finite_weakest, fraction)
                for label, fraction in (
                    ("q00", 0.0),
                    ("q50", 0.5),
                    ("q90", 0.9),
                    ("q99", 0.99),
                    ("q100", 1.0),
                )
            },
            "mean_block_quantiles": {
                label: quantile(finite_means, fraction)
                for label, fraction in (
                    ("q00", 0.0),
                    ("q50", 0.5),
                    ("q90", 0.9),
                    ("q99", 0.99),
                    ("q100", 1.0),
                )
            },
            "weakest_threshold_counts": {
                format(threshold, ".2f"): int(np.sum(finite_weakest >= threshold))
                for threshold in (0.25, 0.5, 0.75, 1.0)
            },
        },
        "method": {
            "trace_identity": "a_p(E^(q_i*q_j)_u)=chi_p(q_i(u)q_j(u))*a_p(E_u)",
            "local_score": "-(1/p)*sum_{u in F_p} a_p(E^(q_i*q_j)_u)",
            "ranking": "descending weakest logarithmically weighted prime-block score",
            "constant_squareclasses_preserved": True,
            "numpy_version": np.__version__,
        },
        "inputs": {
            relative(direct_path): digest(direct_path),
            relative(covers_path): digest(covers_path),
            relative(arguments.singleton_census): digest(arguments.singleton_census),
        },
        "runtime_seconds": perf_counter() - started,
        "proof_boundary": (
            "This finite-prime census is a heuristic over the displayed singleton prefix. "
            "It does not construct a product-twist section or prove a third nonzero V4 character."
        ),
        "reproducing_command": shlex.join(sys.argv),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    leader = top[0]
    print(
        f"DIRECTPRODUCTCENSUS|status=PASS_HEURISTIC|pairs={len(products)}"
        f"|leader={leader['key']}|weakest={leader['weakest_block_score']:.6f}"
        f"|mean={leader['mean_block_score']:.6f}|seconds={result['runtime_seconds']:.3f}"
        f"|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
