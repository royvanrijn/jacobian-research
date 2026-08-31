#!/usr/bin/env python3
"""Nagao sieve on the rational parameter of Elkies's first rank-18 cover.

The hot path reuses the certified published-family local tables.  For each
prime it pulls those tables back along

    t = (289444-r^2)/(130*r-38636),

so the staged scan enumerates primitive rational ``r`` directly.  This is a
heuristic candidate sieve on a generic-rank-at-least-18 family, not a rank
certificate for any specialization.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import importlib.util
import json
from math import gcd
from pathlib import Path
import shlex
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
COMMON_PATH = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
COVER = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank18_first_cover.json"


def load_common():
    spec = importlib.util.spec_from_file_location("elkies_r17_nagao_common", COMMON_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {COMMON_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def conic_t_pair(r_numerator: int, r_denominator: int) -> tuple[int, int]:
    if r_denominator == 0:
        return 1, 0
    t_numerator = 289444 * r_denominator**2 - r_numerator**2
    t_denominator = (
        130 * r_numerator * r_denominator - 38636 * r_denominator**2
    )
    common = gcd(abs(t_numerator), abs(t_denominator))
    t_numerator //= common
    t_denominator //= common
    if t_denominator < 0 or (t_denominator == 0 and t_numerator < 0):
        t_numerator, t_denominator = -t_numerator, -t_denominator
    return t_numerator, t_denominator


def pull_back_table(common, prime, table):
    pulled = []
    for r_index in range(prime + 1):
        if r_index == prime:
            t_index = prime
        else:
            t_numerator = (289444 - r_index**2) % prime
            t_denominator = (130 * r_index - 38636) % prime
            t_index = common.projective_index(t_numerator, t_denominator, prime)
        pulled.append(replace(table[t_index], projective_index=r_index))
    return tuple(pulled)


def main() -> None:
    common = load_common()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--cover", type=Path, default=COVER)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--numerator-bound", type=int, default=1000)
    parser.add_argument("--denominator-bound", type=int, default=1000)
    parser.add_argument("--prime-blocks", default=common.default_prime_blocks_text())
    parser.add_argument("--keep-per-bucket", default="32,16,8")
    parser.add_argument("--height-bucket-width", type=int, default=100)
    parser.add_argument("--finalists", type=int, default=1000)
    args = parser.parse_args()
    if args.finalists < 1:
        parser.error("--finalists must be positive")

    cover = json.loads(args.cover.read_text())
    if cover.get("status") != "PASS_EXACT_ELKIES_2026_FIRST_RANK18_COVER_SECTION":
        raise ValueError("the first-cover exact certificate status is missing")
    started = perf_counter()
    model = common.load_family_model(args.model)
    prime_blocks = common.parse_prime_blocks(args.prime_blocks)
    keep_per_bucket = tuple(int(value) for value in args.keep_per_bucket.split(","))
    if len(keep_per_bucket) != len(prime_blocks):
        parser.error("--keep-per-bucket must have one entry per prime block")
    base_blocks, rejected = common.build_residue_tables(model, prime_blocks)
    pulled_blocks_list = []
    rejected = list(rejected)
    for block_number, block in enumerate(base_blocks, start=1):
        pulled_block = {}
        for prime, table in block.items():
            try:
                pulled_block[prime] = pull_back_table(common, prime, table)
            except ValueError:
                rejected.append(
                    {
                        "prime": prime,
                        "block": block_number,
                        "reason": "conic_parameter_map_indeterminate_mod_prime",
                    }
                )
        if not pulled_block:
            raise ValueError(f"prime block {block_number} has no usable pullback primes")
        pulled_blocks_list.append(pulled_block)
    pulled_blocks = tuple(pulled_blocks_list)
    survivors, stages = common.run_staged_sieve(
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
        table_blocks=pulled_blocks,
        keep_per_bucket=keep_per_bucket,
        bucket_width=args.height_bucket_width,
    )
    finalists = sorted(survivors, key=common.candidate_sort_key)[: args.finalists]
    finalist_records = []
    for candidate in finalists:
        record = common.candidate_record(candidate)
        t_numerator, t_denominator = conic_t_pair(
            candidate.numerator, candidate.denominator
        )
        record["base_t_projective_pair"] = [t_numerator, t_denominator]
        record["base_t"] = (
            "infinity"
            if t_denominator == 0
            else f"{t_numerator}/{t_denominator}"
        )
        finalist_records.append(record)

    payload = {
        "schema": "elkies-k3.elkies-2026-rank18-conic-nagao-sieve.v1",
        "status": "PASS_BOUNDED_HEURISTIC_RANK18_CONIC_NAGAO_SIEVE",
        "model": {
            "source": str(model.source),
            "source_sha256": model.source_sha256,
            "coordinate": model.coordinate,
        },
        "cover": {
            "source": str(args.cover.resolve()),
            "parameter": "r",
            "map_to_t": "t=(289444-r^2)/(130*r-38636)",
            "generic_rank_lower_bound": 18,
        },
        "search": {
            "numerator_bound": args.numerator_bound,
            "denominator_bound": args.denominator_bound,
            "prime_blocks": [list(block) for block in prime_blocks],
            "keep_per_bucket": list(keep_per_bucket),
            "height_bucket_width": args.height_bucket_width,
            "rejected_primes": list(rejected),
        },
        "stages": stages,
        "final_survivor_count": len(survivors),
        "finalists": finalist_records,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The conic map and generic eighteenth section are exact. Nagao scores and "
            "staged survival are heuristics; finalists require exact specialization and "
            "independence certificates before promotion."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS rank18_conic survivors={len(survivors)} finalists={len(finalists)} "
        f"seconds={payload['runtime_seconds']:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
