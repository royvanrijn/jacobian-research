#!/usr/bin/env python3
"""Nagao sieve on an exactly certified extreme-anchored R17 bisection cover.

The cover certificate supplies a rational map ``t=t(r)`` with ``r=0`` equal
to a known high-jump fibre.  This script pulls the native R17 local tables
back along that map and scans primitive rational ``r``.  It is a bounded
heuristic candidate sieve, not a specialization-rank certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import gcd, lcm
from pathlib import Path
import shlex
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
COMMON_PATH = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json"
CERTIFICATE_STATUS = "PASS_EXACT_EXTREME_ANCHORED_MW18_COVERS"


def load_common():
    spec = importlib.util.spec_from_file_location("r17_anchored_nagao_common", COMMON_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {COMMON_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def coefficients(record: dict, key: str) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in record[key])


def map_data(cover: dict) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...], int]:
    record = cover["anchor_line_parameterization"]["t_of_r"]
    numerator = coefficients(record, "numerator_coefficients_low_to_high")
    denominator = coefficients(record, "denominator_coefficients_low_to_high")
    degree = max(len(numerator), len(denominator)) - 1
    if degree < 1:
        raise ValueError("the anchored cover map is constant")
    return numerator, denominator, degree


def homogeneous_fraction_value(
    values: tuple[Fraction, ...], degree: int, numerator: int, denominator: int
) -> Fraction:
    return sum(
        (value * numerator**power * denominator ** (degree - power) for power, value in enumerate(values)),
        Fraction(0),
    )


def normalize_pair(numerator: int, denominator: int) -> tuple[int, int]:
    if numerator == 0 and denominator == 0:
        raise ValueError("the rational cover map is indeterminate")
    common = gcd(abs(numerator), abs(denominator))
    numerator //= common
    denominator //= common
    if denominator < 0 or (denominator == 0 and numerator < 0):
        numerator, denominator = -numerator, -denominator
    return numerator, denominator


def exact_map_pair(
    r_numerator: int,
    r_denominator: int,
    numerator_coefficients: tuple[Fraction, ...],
    denominator_coefficients: tuple[Fraction, ...],
    degree: int,
) -> tuple[int, int]:
    numerator = homogeneous_fraction_value(
        numerator_coefficients, degree, r_numerator, r_denominator
    )
    denominator = homogeneous_fraction_value(
        denominator_coefficients, degree, r_numerator, r_denominator
    )
    multiplier = lcm(numerator.denominator, denominator.denominator)
    return normalize_pair(int(numerator * multiplier), int(denominator * multiplier))


def coefficient_mod(value: Fraction, prime: int) -> int:
    denominator = value.denominator % prime
    if denominator == 0:
        raise ValueError("cover-map coefficient denominator vanishes")
    return value.numerator % prime * pow(denominator, -1, prime) % prime


def homogeneous_mod_value(
    values: tuple[Fraction, ...], degree: int, numerator: int, denominator: int, prime: int
) -> int:
    return sum(
        coefficient_mod(value, prime)
        * pow(numerator, power, prime)
        * pow(denominator, degree - power, prime)
        for power, value in enumerate(values)
    ) % prime


def pull_back_table(common, prime, table, numerator_coefficients, denominator_coefficients, degree):
    pulled = []
    for r_index in range(prime + 1):
        r_numerator, r_denominator = ((r_index, 1) if r_index < prime else (1, 0))
        t_numerator = homogeneous_mod_value(
            numerator_coefficients, degree, r_numerator, r_denominator, prime
        )
        t_denominator = homogeneous_mod_value(
            denominator_coefficients, degree, r_numerator, r_denominator, prime
        )
        t_index = common.projective_index(t_numerator, t_denominator, prime)
        pulled.append(replace(table[t_index], projective_index=r_index))
    return tuple(pulled)


def select_cover(document: dict, curve_id: int, cover_label: str | None):
    matches = []
    for chart in document["charts"]:
        for fibre in chart["fibres"]:
            if int(fibre["curve_id"]) != curve_id:
                continue
            for cover in fibre["covers"]:
                if cover.get("extreme_anchored") is not True:
                    continue
                if cover_label is None or cover["label"] == cover_label:
                    matches.append((chart, fibre, cover))
    if not matches:
        suffix = "" if cover_label is None else f" with label {cover_label}"
        raise ValueError(f"curve {curve_id} has no certified extreme-anchored cover{suffix}")
    if cover_label is not None and len(matches) != 1:
        raise ValueError("the requested cover label is not unique")
    return matches[0]


def main() -> None:
    common = load_common()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    anchor_group = parser.add_mutually_exclusive_group(required=True)
    anchor_group.add_argument("--curve-id", type=int)
    anchor_group.add_argument("--historical-rank28", action="store_true")
    parser.add_argument("--cover-label")
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

    certificate_bytes = args.certificate.read_bytes()
    certificate = json.loads(certificate_bytes)
    if certificate.get("status") != CERTIFICATE_STATUS:
        raise ValueError("the exact extreme-anchored cover status is missing")
    if args.historical_rank28:
        cover = certificate["historical_rank28_anchor"]
        if args.cover_label is not None and args.cover_label != cover["label"]:
            raise ValueError("the requested label is not the historical rank-28 cover")
        chart = {"direct_model": cover["direct_model"]}
        fibre = cover
        anchor_id = "historical-rank28"
    else:
        chart, fibre, cover = select_cover(certificate, args.curve_id, args.cover_label)
        anchor_id = f"curve-{args.curve_id}"
    if cover["anchor_line_parameterization"].get("anchor_parameter") != "0":
        raise ValueError("the selected cover is not centered at the extreme anchor")
    numerator_coefficients, denominator_coefficients, map_degree = map_data(cover)
    anchor_pair = exact_map_pair(
        0, 1, numerator_coefficients, denominator_coefficients, map_degree
    )
    expected_anchor = Fraction(fibre["native_parameter"])
    if anchor_pair[1] == 0 or Fraction(*anchor_pair) != expected_anchor:
        raise ValueError("r=0 does not map to the certified extreme fibre")

    model_path = ROOT / chart["direct_model"]
    started = perf_counter()
    model = common.load_family_model(model_path)
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
                pulled_block[prime] = pull_back_table(
                    common,
                    prime,
                    table,
                    numerator_coefficients,
                    denominator_coefficients,
                    map_degree,
                )
            except ValueError:
                rejected.append(
                    {
                        "prime": prime,
                        "block": block_number,
                        "reason": "anchored_cover_map_bad_or_indeterminate_mod_prime",
                    }
                )
        if not pulled_block:
            raise ValueError(f"prime block {block_number} has no usable pullback primes")
        pulled_blocks_list.append(pulled_block)
    survivors, stages = common.run_staged_sieve(
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
        table_blocks=tuple(pulled_blocks_list),
        keep_per_bucket=keep_per_bucket,
        bucket_width=args.height_bucket_width,
    )
    finalists = sorted(survivors, key=common.candidate_sort_key)[: args.finalists]
    finalist_records = []
    for candidate in finalists:
        record = common.candidate_record(candidate)
        t_numerator, t_denominator = exact_map_pair(
            candidate.numerator,
            candidate.denominator,
            numerator_coefficients,
            denominator_coefficients,
            map_degree,
        )
        record["base_t_projective_pair"] = [t_numerator, t_denominator]
        record["base_t"] = (
            "infinity" if t_denominator == 0 else f"{t_numerator}/{t_denominator}"
        )
        record["is_certified_anchor"] = candidate.numerator == 0 and candidate.denominator == 1
        finalist_records.append(record)

    payload = {
        "schema": "elkies-k3.r17-extreme-anchored-mw18-nagao-sieve.v1",
        "status": "PASS_BOUNDED_HEURISTIC_EXTREME_ANCHORED_MW18_NAGAO_SIEVE",
        "model": {
            "source": str(model.source),
            "source_sha256": model.source_sha256,
            "coordinate": model.coordinate,
        },
        "cover": {
            "certificate_source": str(args.certificate.resolve()),
            "certificate_sha256": sha256(certificate_bytes).hexdigest(),
            "anchor_id": anchor_id,
            "curve_id": args.curve_id,
            "label": cover["label"],
            "generic_rank_lower_bound": 18,
            "displayed_anchor_jump_over_MW17": fibre["displayed_jump_over_MW17"],
            "anchor_parameter_r": "0",
            "anchor_base_t": fibre["native_parameter"],
            "map_degree": map_degree,
            "t_of_r": cover["anchor_line_parameterization"]["t_of_r"],
        },
        "search": {
            "numerator_bound": args.numerator_bound,
            "denominator_bound": args.denominator_bound,
            "prime_blocks": [list(block) for block in prime_blocks],
            "keep_per_bucket": list(keep_per_bucket),
            "height_bucket_width": args.height_bucket_width,
            "rejected_primes": rejected,
        },
        "stages": stages,
        "final_survivor_count": len(survivors),
        "finalists": finalist_records,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The conic map, r=0 anchor, generic eighteenth section, and nonzero anchor "
            "class are exact in the source certificate. Nagao scores and staged survival "
            "are heuristics; finalists require exact specialization and independence "
            "certificates before promotion."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "PASS extreme_anchored_mw18_nagao "
        f"anchor={anchor_id} cover={cover['label']} survivors={len(survivors)} "
        f"finalists={len(finalists)} seconds={payload['runtime_seconds']:.3f} "
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
