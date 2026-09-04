#!/usr/bin/env sage
"""Collect exact relations in Kummer-shaped mod-two S-class quotients.

Status: ACTIVE_SEARCH
Claim: exact relations in bounded materialized ideal-class presentations only.
Inputs: pinned public-fibre and Kummer/class-group pressure certificates.
Outputs: a user-selected local checkpoint; no generated theorem artifact.

The collector maintains two presentations from every exact relation:

``generic``
    quotient by the first seventeen generic MW half-ideals, retaining each
    known exceptional half-ideal as a formal class coordinate;

``full-known``
    quotient by every certified known half-ideal.

Candidate reductions cycle through ``alpha*I_i``, ``alpha*I_i*I_j``, and
short products ``alpha*prod(I_k^e_k)`` (optionally with signed exponents).
Exceptional half-ideals are
sampled with a declared extra weight.  They therefore shape the norm lattice
even in the full-known lane, where their classes are killed in the answer.
The active target is a cyclically chosen nonzero product of free columns in
the requested quotient, and each closed relation records its actual rank gain
in both projections.

Every accepted row stores enough ideal data to replay its exact equality.
Until a separate factor-base generation proof is supplied, either displayed
dimension is a relation-collection fingerprint only: it is not a class-group,
Selmer, or Mordell--Weil rank upper bound.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random
import time
from typing import Any

from sage.all import QQ, ZZ, PolynomialRing, pari, prime_range
from sage.version import version as sage_version

from r17_kummer_quotient_search import (
    BinaryRows,
    GENERIC_POINT_COUNT,
    parse_strategies,
    projection_masks,
    select_companion_terms,
)


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
PRESSURE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
PUBLIC_STATUS = "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES"
PRESSURE_STATUS = "PROVED_KUMMER_FORCED_CUBIC_CLASS_GROUP_2RANK_LOWER_BOUNDS"
SCHEMA = "elliptic-curves.r17-kummer-quotient-sclass-collector.v2"
PROTOCOL = "R17KUMMERQSCLASS"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ideal_key(nf, ideal) -> str:
    return str(pari.idealhnf(nf, ideal))


def mask_hex(mask: int) -> str:
    return format(mask, "x")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--curve-id", type=int, choices=(351, 356, 376, 377, 385), required=True
    )
    parser.add_argument(
        "--objective",
        choices=("generic", "full-known"),
        default="full-known",
        help="presentation whose free columns determine the next target",
    )
    parser.add_argument(
        "--kill-point-count",
        type=int,
        default=None,
        help=(
            "compatibility override: 17 selects --objective generic and the "
            "full supplied point count selects --objective full-known"
        ),
    )
    parser.add_argument("--generic-point-count", type=int, default=GENERIC_POINT_COUNT)
    parser.add_argument("--factor-base-bound", type=int, default=5000)
    parser.add_argument("--attempts", type=int, default=10000)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument(
        "--max-target-columns",
        type=int,
        default=1,
        help="cycle through unresolved target products of widths 1 through this value",
    )
    parser.add_argument(
        "--companion-strategies",
        default="single,pair,sparse",
        help="comma-separated cycle drawn from single,pair,sparse",
    )
    parser.add_argument(
        "--sparse-min",
        "--min-kummer-companions",
        dest="sparse_min",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--sparse-max",
        "--max-kummer-companions",
        dest="sparse_max",
        type=int,
        default=6,
    )
    parser.add_argument("--companion-exponent-radius", type=int, default=1)
    parser.add_argument("--signed-companion-exponents", action="store_true")
    parser.add_argument("--exceptional-companion-weight", type=int, default=4)
    parser.add_argument("--max-s-companions", type=int, default=2)
    parser.add_argument("--direction-radius", type=int, default=16)
    parser.add_argument(
        "--reduction-engine",
        choices=("idealred", "idealredmodpower2"),
        default="idealred",
        help="ordinary directed reduction or PARI reduction modulo ideal squares",
    )
    parser.add_argument(
        "--large-prime-bound",
        type=int,
        default=None,
        help="large-prime ideal norm bound (default: factor-base bound squared)",
    )
    parser.add_argument("--max-large-primes", type=int, choices=(1, 2), default=1)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--checkpoint-every-attempts", type=int, default=1000)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.factor_base_bound < 2:
        parser.error("the factor-base bound must be at least 2")
    if args.attempts <= 0 or args.timeout_seconds <= 0:
        parser.error("attempt and time bounds must be positive")
    if args.max_target_columns <= 0:
        parser.error("the maximum target width must be positive")
    if not (0 <= args.generic_point_count <= 29):
        parser.error("the generic point count must lie between 0 and 29")
    if args.sparse_min < 0 or args.sparse_min > args.sparse_max:
        parser.error("invalid sparse Kummer companion interval")
    if args.companion_exponent_radius <= 0:
        parser.error("the companion exponent radius must be positive")
    if args.exceptional_companion_weight <= 0:
        parser.error("the exceptional companion weight must be positive")
    if args.max_s_companions < 0 or args.direction_radius <= 0:
        parser.error("invalid S-companion or direction bound")
    if args.checkpoint_every_attempts <= 0:
        parser.error("the checkpoint interval must be positive")
    try:
        args.companion_strategies = parse_strategies(args.companion_strategies)
    except ValueError as error:
        parser.error(str(error))
    if args.large_prime_bound is None:
        args.large_prime_bound = args.factor_base_bound**2
    if args.large_prime_bound < args.factor_base_bound:
        parser.error("the large-prime bound cannot be below the factor-base bound")
    return args


def load_curve(curve_id: int):
    public = json.loads(PUBLIC.read_text())
    pressure = json.loads(PRESSURE.read_text())
    if public.get("status") != PUBLIC_STATUS:
        raise ArithmeticError("the public-fibre certificate is not passing")
    if pressure.get("status") != PRESSURE_STATUS:
        raise ArithmeticError("the Kummer half-ideal certificate is not passing")
    public_record = next(
        record for record in public["records"] if int(record["id"]) == curve_id
    )
    pressure_record = next(
        record for record in pressure["curves"] if int(record["curve_id"]) == curve_id
    )
    if len(public_record["points"]) != len(pressure_record["point_half_ideals"]):
        raise ArithmeticError("the point and half-ideal counts differ")
    return public_record, pressure_record


def factor_base(nf, bound: int):
    ideals = []
    decompositions = {}
    for rational_prime in prime_range(bound + 1):
        prime = int(rational_prime)
        decomposition = list(pari.idealprimedec(nf, prime))
        decompositions[prime] = decomposition
        for ideal in decomposition:
            if prime ** int(ideal[3]) <= bound:
                ideals.append(ideal)
    index = {ideal_key(nf, ideal): offset for offset, ideal in enumerate(ideals)}
    if len(index) != len(ideals):
        raise ArithmeticError("the factor base contains duplicate prime ideals")
    return ideals, decompositions, index


def checkpoint(path: Path, document: dict[str, Any]):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def ideal_power(nf, ideal, exponent: int):
    if exponent == 0:
        return pari.idealhnf(nf, 1)
    if exponent > 0:
        return pari.idealpow(nf, ideal, exponent)
    return pari.idealinv(nf, pari.idealpow(nf, ideal, -exponent))


def main():
    args = parse_args()
    output = args.checkpoint.resolve()
    if output.exists() and not args.overwrite:
        raise FileExistsError(output)
    rng = random.Random(args.seed)
    public_record, pressure_record = load_curve(args.curve_id)

    available_half_ideal_records = pressure_record["point_half_ideals"]
    point_count = len(available_half_ideal_records)
    if args.generic_point_count > point_count:
        raise ArithmeticError("the generic point prefix exceeds the certified point list")
    if args.kill_point_count is not None:
        if args.kill_point_count == args.generic_point_count:
            args.objective = "generic"
        elif args.kill_point_count == point_count:
            args.objective = "full-known"
        else:
            raise ArithmeticError(
                "--kill-point-count must equal the generic prefix or all known points"
            )
    objective_killed_point_count = (
        args.generic_point_count if args.objective == "generic" else point_count
    )

    ainvs = tuple(QQ(value) for value in public_record["ainvs"])
    a1, a2, a3, a4, a6 = ainvs
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    polynomial = (
        z**3
        + ZZ(a1**2 + 4 * a2) * z**2
        + ZZ(8 * (a1 * a3 + 2 * a4)) * z
        + ZZ(16 * (a3**2 + 4 * a6))
    )
    bad_primes = [ZZ(value) for value in public_record["bad_primes"]]
    pari.addprimes(bad_primes)
    nf = pari.nfinit([pari(polynomial), bad_primes])
    if list(pari.nfcertify(nf)):
        raise ArithmeticError("the cubic maximal order failed certification")
    signature_length = sum(int(value) for value in nf.nf_get_sign())

    ideals, decompositions, column_index = factor_base(nf, args.factor_base_bound)
    factor_base_width = len(ideals)
    columns = [
        {
            "column": index,
            "prime_ideal": str(ideal),
            "prime_ideal_hnf": ideal_key(nf, ideal),
            "norm": str(pari.idealnorm(nf, ideal)),
        }
        for index, ideal in enumerate(ideals)
    ]
    exceptional_half_ideal_records = available_half_ideal_records[
        args.generic_point_count :
    ]
    exceptional_count = len(exceptional_half_ideal_records)
    generic_width = factor_base_width + exceptional_count
    full_known_width = factor_base_width
    generic_rows = BinaryRows()
    full_known_rows = BinaryRows()
    initial_rows = []

    def add_initial_base_row(mask: int) -> tuple[bool, bool]:
        return generic_rows.add(mask), full_known_rows.add(mask)

    for prime, decomposition in decompositions.items():
        keys = [ideal_key(nf, ideal) for ideal in decomposition]
        if not all(key in column_index for key in keys):
            continue
        mask = 0
        for ideal, key in zip(decomposition, keys):
            if int(ideal[2]) & 1:
                mask ^= 1 << column_index[key]
        if mask:
            generic_gain, full_known_gain = add_initial_base_row(mask)
            initial_rows.append(
                {
                    "kind": "principal_rational_prime",
                    "rational_prime": prime,
                    "base_row_mask_hex": mask_hex(mask),
                    "generic_increased_rank": generic_gain,
                    "full_known_increased_rank": full_known_gain,
                }
            )

    s_ideals = [
        ideal for prime in bad_primes for ideal in pari.idealprimedec(nf, prime)
    ]
    for ideal in s_ideals:
        key = ideal_key(nf, ideal)
        if key not in column_index:
            continue
        mask = 1 << column_index[key]
        generic_gain, full_known_gain = add_initial_base_row(mask)
        initial_rows.append(
            {
                "kind": "killed_s_prime_ideal",
                "prime_ideal": str(ideal),
                "base_row_mask_hex": mask_hex(mask),
                "generic_increased_rank": generic_gain,
                "full_known_increased_rank": full_known_gain,
            }
        )

    half_ideals = [
        pari(record["half_ideal_hnf"]) for record in available_half_ideal_records
    ]
    half_ideal_labels = [record["label"] for record in available_half_ideal_records]
    initial_generic_rank = generic_rows.rank
    initial_full_known_rank = full_known_rows.rank
    accepted_relations = []
    large_prime_partials = {}
    large_prime_collisions = []
    strategy_stats = {
        strategy: {
            "attempts": 0,
            "fully_factored_reductions": 0,
            "factor_base_smooth_relations": 0,
            "large_prime_partials": 0,
            "large_prime_cycles": 0,
            "generic_rank_gains": 0,
            "full_known_rank_gains": 0,
            "generic_companion_terms": 0,
            "exceptional_companion_terms": 0,
            "negative_exponent_terms": 0,
            "even_exponent_terms": 0,
        }
        for strategy in args.companion_strategies
    }
    smooth_reduction_count = 0
    fully_factored_reduction_count = 0
    started = time.monotonic()
    status = "attempt_limit"
    attempts_completed = 0

    def projection_record(base_mask: int, exceptional_mask: int, add: bool):
        generic_mask, full_known_mask = projection_masks(
            base_mask=base_mask,
            exceptional_parity_mask=exceptional_mask,
            factor_base_width=factor_base_width,
        )
        generic_residual = generic_rows.reduce(generic_mask)
        full_known_residual = full_known_rows.reduce(full_known_mask)
        record = {
            "generic_row_mask_hex": mask_hex(generic_mask),
            "full_known_row_mask_hex": mask_hex(full_known_mask),
            "generic_residual_mask_before_hex": mask_hex(generic_residual),
            "full_known_residual_mask_before_hex": mask_hex(full_known_residual),
            "generic_residual_weight_before": generic_residual.bit_count(),
            "full_known_residual_weight_before": full_known_residual.bit_count(),
            "generic_increased_rank": bool(generic_residual),
            "full_known_increased_rank": bool(full_known_residual),
        }
        if add:
            if generic_rows.add(generic_mask) != bool(generic_residual):
                raise ArithmeticError("generic projection rank prediction failed")
            if full_known_rows.add(full_known_mask) != bool(full_known_residual):
                raise ArithmeticError("full-known projection rank prediction failed")
        return record

    def projection_summary(rows: BinaryRows, width: int, initial_rank: int):
        return {
            "presentation_width": width,
            "initial_mod2_rank": initial_rank,
            "current_mod2_rank": rows.rank,
            "relation_rank_gain": rows.rank - initial_rank,
            "materialized_quotient_dimension": width - rows.rank,
        }

    def settings_record():
        return {
            "curve_id": args.curve_id,
            "objective": args.objective,
            "objective_killed_point_count": objective_killed_point_count,
            "generic_point_count": args.generic_point_count,
            "known_point_count": point_count,
            "factor_base_bound": args.factor_base_bound,
            "attempts": args.attempts,
            "timeout_seconds": args.timeout_seconds,
            "max_target_columns": args.max_target_columns,
            "companion_strategies": list(args.companion_strategies),
            "sparse_min": args.sparse_min,
            "sparse_max": args.sparse_max,
            "companion_exponent_radius": args.companion_exponent_radius,
            "signed_companion_exponents": args.signed_companion_exponents,
            "exceptional_companion_weight": args.exceptional_companion_weight,
            "max_s_companions": args.max_s_companions,
            "direction_radius": args.direction_radius,
            "reduction_engine": args.reduction_engine,
            "large_prime_bound": args.large_prime_bound,
            "max_large_primes": args.max_large_primes,
            "seed": args.seed,
            "checkpoint_every_attempts": args.checkpoint_every_attempts,
            "checkpoint": str(args.checkpoint),
        }

    def build_document(document_status: str):
        projections = {
            "generic": projection_summary(
                generic_rows, generic_width, initial_generic_rank
            ),
            "full-known": projection_summary(
                full_known_rows, full_known_width, initial_full_known_rank
            ),
        }
        objective = projections[args.objective]
        return {
            "schema": SCHEMA,
            "status": document_status,
            "curve_id": args.curve_id,
            "certified_displayed_mw_gain_over_generic_mw17": int(
                pressure_record["residual_gain_over_mw17"]
            ),
            "settings": settings_record(),
            "field_polynomial": str(polynomial),
            "field_signature": list(map(int, nf.nf_get_sign())),
            "factor_base_columns": columns,
            "exceptional_half_ideal_columns_in_generic_projection": [
                {
                    "column": factor_base_width + index,
                    "point_index": args.generic_point_count + index,
                    "label": record["label"],
                    "half_ideal_hnf": record["half_ideal_hnf"],
                }
                for index, record in enumerate(exceptional_half_ideal_records)
            ],
            "initial_exact_rows": initial_rows,
            "objective": args.objective,
            "projections": projections,
            # Compatibility aliases always refer to the declared objective.
            "initial_mod2_rank": objective["initial_mod2_rank"],
            "current_mod2_rank": objective["current_mod2_rank"],
            "materialized_quotient_dimension": objective[
                "materialized_quotient_dimension"
            ],
            "accepted_quotient_relations": accepted_relations,
            "large_prime_collisions": large_prime_collisions,
            "unmatched_large_prime_partial_count": len(large_prime_partials),
            "fully_factored_reduction_count": fully_factored_reduction_count,
            "smooth_reduction_count": smooth_reduction_count,
            "strategy_statistics": strategy_stats,
            "attempts_completed": attempts_completed,
            "elapsed_seconds": time.monotonic() - started,
            "claim_boundary": [
                "Every closed stored relation is an exact ideal-class relation modulo the declared S and generic or full-known Kummer half-ideal classes.",
                "Exceptional half-ideal columns in the generic presentation are formal exact class generators, not an assertion that they lie in the materialized factor base.",
                "Neither materialized factor-base presentation is a global class-group upper bound until factor-base generation is separately proved.",
                "The cross-curve dimensions and correlations are bounded relation fingerprints, not S-class, Selmer, or rank theorems.",
            ],
            "inputs": {
                str(PUBLIC.relative_to(ROOT)): digest(PUBLIC),
                str(PRESSURE.relative_to(ROOT)): digest(PRESSURE),
            },
        }

    for attempt in range(1, args.attempts + 1):
        attempts_completed = attempt - 1
        elapsed = time.monotonic() - started
        if elapsed >= args.timeout_seconds:
            status = "strict_wall_timeout"
            break
        attempts_completed = attempt
        objective_rows = generic_rows if args.objective == "generic" else full_known_rows
        objective_width = generic_width if args.objective == "generic" else full_known_width
        # Rotate through unresolved column products.  Repeatedly asking for
        # the first free column trapped the v1 pilot on a single hard ideal;
        # allowing short products is the quotient analogue of a multi-special-q
        # relation search.
        requested_target_width = 1 + (
            (attempt - 1) % args.max_target_columns
        )
        target_start = (
            (attempt - 1) // args.max_target_columns
        ) % objective_width
        target_columns = objective_rows.free_column_combination(
            objective_width,
            start=target_start,
            count=requested_target_width,
        )
        if not target_columns:
            status = "materialized_objective_quotient_closed"
            break

        target = pari.idealhnf(nf, 1)
        base_mask = 0
        exceptional_mask = 0
        target_terms = []
        for target_column in target_columns:
            if target_column < factor_base_width:
                target_factor = ideals[target_column]
                base_mask ^= 1 << target_column
                target_terms.append(
                    {
                        "kind": "factor_base_prime_ideal",
                        "column": target_column,
                        "prime_ideal": str(target_factor),
                        "prime_ideal_hnf": ideal_key(nf, target_factor),
                    }
                )
            else:
                exceptional_offset = target_column - factor_base_width
                point_index = args.generic_point_count + exceptional_offset
                target_factor = half_ideals[point_index]
                exceptional_mask ^= 1 << exceptional_offset
                target_terms.append(
                    {
                        "kind": "known_exceptional_half_ideal",
                        "column": target_column,
                        "point_index": point_index,
                        "label": available_half_ideal_records[point_index]["label"],
                        "half_ideal_hnf": available_half_ideal_records[point_index][
                            "half_ideal_hnf"
                        ],
                    }
                )
            target = pari.idealmul(nf, target, target_factor)

        if len(target_terms) == 1:
            target_record = target_terms[0]
        else:
            target_record = {
                "kind": "multi_unresolved_column_product",
                "requested_column_count": requested_target_width,
                "actual_column_count": len(target_terms),
                "terms": target_terms,
            }

        strategy, companion_terms = select_companion_terms(
            rng=rng,
            attempt=attempt,
            labels=half_ideal_labels,
            generic_point_count=args.generic_point_count,
            strategies=args.companion_strategies,
            sparse_min=args.sparse_min,
            sparse_max=args.sparse_max,
            exponent_radius=args.companion_exponent_radius,
            exceptional_weight=args.exceptional_companion_weight,
            signed_exponents=args.signed_companion_exponents,
        )
        strategy_stats[strategy]["attempts"] += 1
        strategy_stats[strategy]["generic_companion_terms"] += sum(
            term.role == "generic_MW17" for term in companion_terms
        )
        strategy_stats[strategy]["exceptional_companion_terms"] += sum(
            term.role == "known_exceptional" for term in companion_terms
        )
        strategy_stats[strategy]["negative_exponent_terms"] += sum(
            term.exponent < 0 for term in companion_terms
        )
        strategy_stats[strategy]["even_exponent_terms"] += sum(
            not term.parity for term in companion_terms
        )
        search_ideal = target
        companion_records = []
        for term in companion_terms:
            search_ideal = pari.idealmul(
                nf,
                search_ideal,
                ideal_power(nf, half_ideals[term.index], term.exponent),
            )
            if term.index >= args.generic_point_count and term.parity:
                exceptional_mask ^= 1 << (term.index - args.generic_point_count)
            companion_records.append(
                {
                    "point_index": term.index,
                    "label": term.label,
                    "role": term.role,
                    "exponent": term.exponent,
                }
            )

        s_labels = []
        if s_ideals:
            for _ in range(rng.randint(0, args.max_s_companions)):
                choice = rng.randrange(len(s_ideals))
                search_ideal = pari.idealmul(nf, search_ideal, s_ideals[choice])
                s_labels.append(str(s_ideals[choice]))
        if args.reduction_engine == "idealred":
            direction = [
                rng.randrange(args.direction_radius)
                for _ in range(signature_length)
            ]
            reduced = pari.idealred(nf, [search_ideal, 1], direction)
            reduced_ideal = reduced[0]
            multiplier = reduced[1]
            expected = pari.idealmul(nf, reduced_ideal, multiplier)
            if ideal_key(nf, search_ideal) != ideal_key(nf, expected):
                raise ArithmeticError("idealred multiplier identity failed")
            reduction_witness = {
                "engine": "idealred",
                "archimedean_direction": direction,
                "multiplier": str(multiplier),
                "verified_identity": "search_ideal = reduced_ideal * (multiplier)",
            }
        else:
            direction = None
            multiplier = pari.idealredmodpower(
                nf, search_ideal, 2, args.factor_base_bound
            )
            multiplier_square = pari.nfeltpow(nf, multiplier, 2)
            reduced_ideal = pari.idealmul(
                nf, search_ideal, multiplier_square
            )
            expected = pari.idealmul(
                nf, search_ideal, multiplier_square
            )
            if ideal_key(nf, reduced_ideal) != ideal_key(nf, expected):
                raise ArithmeticError(
                    "idealredmodpower square-multiplier identity failed"
                )
            reduction_witness = {
                "engine": "idealredmodpower2",
                "factor_limit": args.factor_base_bound,
                "multiplier": str(multiplier),
                "multiplier_square": str(multiplier_square),
                "verified_identity": "reduced_ideal = search_ideal * (multiplier)^2",
            }

        reduced_norm = ZZ(pari.idealnorm(nf, reduced_ideal))
        factors = pari.idealfactor(nf, reduced_ideal, args.factor_base_bound + 1)
        factored_norm = ZZ(1)
        for index in range(int(factors.nrows())):
            factor_ideal = factors[index, 0]
            exponent = int(factors[index, 1])
            factored_norm *= ZZ(pari.idealnorm(nf, factor_ideal)) ** exponent
        if reduced_norm % factored_norm:
            raise ArithmeticError("partial ideal factorization has the wrong norm")
        cofactor = reduced_norm // factored_norm
        if cofactor != 1:
            if cofactor > args.large_prime_bound**args.max_large_primes:
                if attempt % args.checkpoint_every_attempts == 0:
                    attempts_completed = attempt
                    checkpoint(output, build_document("collecting"))
                continue
            rational_large_factors = list(cofactor.factor())
            if (
                len(rational_large_factors) > args.max_large_primes
                or any(
                    prime > args.large_prime_bound
                    for prime, _exponent in rational_large_factors
                )
            ):
                if attempt % args.checkpoint_every_attempts == 0:
                    attempts_completed = attempt
                    checkpoint(output, build_document("collecting"))
                continue
            factors = pari.idealfactor(nf, reduced_ideal)

        fully_factored_reduction_count += 1
        strategy_stats[strategy]["fully_factored_reductions"] += 1
        factor_records = []
        odd_large_ideals = []
        factored_norm = ZZ(1)
        for index in range(int(factors.nrows())):
            factor_ideal = factors[index, 0]
            exponent = int(factors[index, 1])
            key = ideal_key(nf, factor_ideal)
            factor_norm = ZZ(pari.idealnorm(nf, factor_ideal))
            factored_norm *= factor_norm**exponent
            if key in column_index and exponent & 1:
                base_mask ^= 1 << column_index[key]
            elif key not in column_index and exponent & 1:
                if factor_norm > args.large_prime_bound:
                    odd_large_ideals.append("OUT_OF_RANGE:" + key)
                else:
                    odd_large_ideals.append(key)
            factor_records.append(
                {
                    "prime_ideal": str(factor_ideal),
                    "prime_ideal_hnf": key,
                    "exponent": exponent,
                }
            )
        if (
            factored_norm != reduced_norm
            or len(odd_large_ideals) > args.max_large_primes
            or any(key.startswith("OUT_OF_RANGE:") for key in odd_large_ideals)
        ):
            continue

        smooth_reduction_count += 1
        base_record = {
            "attempt": attempt,
            "objective": args.objective,
            "target": target_record,
            "companion_strategy": strategy,
            "kummer_companions": companion_records,
            "killed_s_companion_ideals": s_labels,
            "archimedean_direction": direction,
            "search_ideal_hnf": ideal_key(nf, search_ideal),
            "reduction_witness": reduction_witness,
            "reduced_ideal_hnf": ideal_key(nf, reduced_ideal),
            "reduced_ideal_norm": str(reduced_norm),
            "reduced_ideal_factorization": factor_records,
            "base_row_mask_hex": mask_hex(base_mask),
            "exceptional_parity_mask_hex": mask_hex(exceptional_mask),
        }
        objective_gain = False
        if not odd_large_ideals:
            metrics = projection_record(base_mask, exceptional_mask, add=True)
            base_record.update(metrics)
            base_record["relation_kind"] = "factor_base_smooth"
            accepted_relations.append(base_record)
            strategy_stats[strategy]["factor_base_smooth_relations"] += 1
            strategy_stats[strategy]["generic_rank_gains"] += int(
                metrics["generic_increased_rank"]
            )
            strategy_stats[strategy]["full_known_rank_gains"] += int(
                metrics["full_known_increased_rank"]
            )
            objective_gain = metrics[
                f"{args.objective.replace('-', '_')}_increased_rank"
            ]
        else:
            strategy_stats[strategy]["large_prime_partials"] += 1
            large_keys = set(odd_large_ideals)
            partial_records = [base_record]
            combined_base_mask = base_mask
            combined_exceptional_mask = exceptional_mask
            while large_keys:
                pivot = max(large_keys)
                previous = large_prime_partials.get(pivot)
                if previous is None:
                    base_record["relation_kind"] = (
                        f"{len(odd_large_ideals)}_large_prime_partial"
                    )
                    base_record["odd_large_prime_ideal_hnfs"] = sorted(large_keys)
                    large_prime_partials[pivot] = {
                        "large_keys": sorted(large_keys),
                        "base_row_mask_hex": mask_hex(combined_base_mask),
                        "exceptional_parity_mask_hex": mask_hex(
                            combined_exceptional_mask
                        ),
                        "partial_records": partial_records,
                    }
                    break
                large_keys.symmetric_difference_update(previous["large_keys"])
                combined_base_mask ^= int(previous["base_row_mask_hex"], 16)
                combined_exceptional_mask ^= int(
                    previous["exceptional_parity_mask_hex"], 16
                )
                partial_records.extend(previous["partial_records"])
            else:
                metrics = projection_record(
                    combined_base_mask, combined_exceptional_mask, add=True
                )
                collision = {
                    "relation_kind": "large_prime_cycle",
                    "closing_attempt": attempt,
                    "partial_records": partial_records,
                    "combined_base_row_mask_hex": mask_hex(combined_base_mask),
                    "combined_exceptional_parity_mask_hex": mask_hex(
                        combined_exceptional_mask
                    ),
                    **metrics,
                }
                large_prime_collisions.append(collision)
                strategy_stats[strategy]["large_prime_cycles"] += 1
                strategy_stats[strategy]["generic_rank_gains"] += int(
                    metrics["generic_increased_rank"]
                )
                strategy_stats[strategy]["full_known_rank_gains"] += int(
                    metrics["full_known_increased_rank"]
                )
                objective_gain = metrics[
                    f"{args.objective.replace('-', '_')}_increased_rank"
                ]

        attempts_completed = attempt
        if objective_gain:
            objective_rows = generic_rows if args.objective == "generic" else full_known_rows
            objective_width = generic_width if args.objective == "generic" else full_known_width
            print(
                f"{PROTOCOL}|curve={args.curve_id}|objective={args.objective}|"
                f"attempt={attempt}|rank={objective_rows.rank}|"
                f"dimension={objective_width-objective_rows.rank}|status=drop",
                flush=True,
            )
            checkpoint(output, build_document("collecting"))
        elif attempt % args.checkpoint_every_attempts == 0:
            checkpoint(output, build_document("collecting"))

    elapsed = time.monotonic() - started
    objective_rows = generic_rows if args.objective == "generic" else full_known_rows
    objective_width = generic_width if args.objective == "generic" else full_known_width
    if objective_rows.free_column(objective_width) is None:
        status = "materialized_objective_quotient_closed"
    final = build_document(status)
    final["elapsed_seconds"] = elapsed
    final["software_assumptions"] = {
        "sage": str(sage_version),
        "pari": ".".join(str(part) for part in pari.version()),
    }
    checkpoint(output, final)
    try:
        output_label = output.relative_to(ROOT)
    except ValueError:
        output_label = output
    print(
        f"{PROTOCOL}|curve={args.curve_id}|objective={args.objective}|"
        f"attempts={attempts_completed}|rank={objective_rows.rank}|"
        f"dimension={objective_width-objective_rows.rank}|status={status}|"
        f"seconds={elapsed:.6f}|output={output_label}",
        flush=True,
    )


if __name__ == "__main__":
    main()
