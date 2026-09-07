#!/usr/bin/env python3
"""Pilot Hensel--CRT--lattice search in Fermigier's rank-12 family.

This is a bounded experiment.  It produces candidates and exact local and
conductor data; it does not infer rank from a Mestre--Nagao score.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable
from decimal import Decimal

from crt_lattice import ConstraintChoice, beam_combine, short_rational_representatives
from ek_k3 import primes_up_to, rational_to_string
from fermigier_mestre import FermigierMestreFamily, ROOTS
from pari_bridge import minimal_curve_data, pari_version


TARGET_LOG_CONDUCTOR = Decimal("182.72")


def parse_power_primes(value: str) -> tuple[tuple[int, int], ...]:
    if not value:
        return ()
    result: list[tuple[int, int]] = []
    for item in value.split(","):
        prime_text, separator, exponent_text = item.partition(":")
        if not separator:
            raise argparse.ArgumentTypeError("power primes use p:k syntax")
        result.append((int(prime_text), int(exponent_text)))
    return tuple(result)


def parse_primes(value: str) -> tuple[int, ...]:
    if not value:
        return ()
    return tuple(int(item) for item in value.split(","))


def residue_of_rational(numerator: int, denominator: int, prime: int) -> int | None:
    if gcd(denominator, prime) != 1:
        return None
    return numerator * pow(denominator, -1, prime) % prime


def score_term(trace: int, point_count: int, prime: int, score: str) -> float:
    if score == "fermigier-good":
        return (2 - trace) / point_count * log(prime)
    if score == "nagao-log":
        return log(point_count / prime)
    raise ValueError(f"unknown score {score}")


def rank_choices(
    prime: int, *, count: int, score: str
) -> tuple[ConstraintChoice, ...]:
    choices: list[ConstraintChoice] = []
    for residue in range(prime):
        local = FermigierMestreFamily.local_data(residue, prime)
        if not local.good_reduction:
            continue
        value = score_term(local.trace, local.point_count, prime, score)
        choices.append(
            ConstraintChoice(
                prime=prime,
                modulus=prime,
                residue=residue,
                kind="rank",
                label=f"good residue {residue} mod {prime}, a_p={local.trace}",
                local_score=value,
            )
        )
    return tuple(
        sorted(choices, key=lambda item: (-item.local_score, item.residue))[:count]
    )


def power_choices(
    prime: int, exponent: int, *, split_only: bool
) -> tuple[ConstraintChoice, ...]:
    choices = []
    for root in FermigierMestreFamily.power_roots(
        prime, exponent, split_only=split_only
    ):
        # Use the actual bad Euler factor as a small tie-breaker.  The main
        # objective remains height; discriminant roots are never mislabeled
        # as good-reduction a_p observations.
        local_trace = 1 if root.split_multiplicative else -1
        choices.append(
            ConstraintChoice(
                prime=prime,
                modulus=root.modulus,
                residue=root.residue,
                kind="power",
                label=(
                    f"v_{prime}(H(T))>={exponent}; "
                    + ("split" if root.split_multiplicative else "nonsplit")
                    + " multiplicative"
                ),
                local_score=-log(1 - local_trace / prime),
            )
        )
    return tuple(choices)


def full_score(
    numerator: int, denominator: int, *, bound: int, score: str
) -> dict[str, Any]:
    total = 0.0
    used = 0
    skipped_denominator = 0
    skipped_bad = 0
    traces: list[dict[str, int]] = []
    for prime in primes_up_to(bound):
        if prime < 5:
            continue
        residue = residue_of_rational(numerator, denominator, prime)
        if residue is None:
            skipped_denominator += 1
            continue
        local = FermigierMestreFamily.local_data(residue, prime)
        if not local.good_reduction:
            skipped_bad += 1
            continue
        term = score_term(local.trace, local.point_count, prime, score)
        total += term
        used += 1
        traces.append({"prime": prime, "trace": local.trace})
    return {
        "value": total,
        "primes_used": used,
        "skipped_denominator_primes": skipped_denominator,
        "skipped_bad_primes": skipped_bad,
        "traces": traces,
    }


def candidate_records(
    states: Iterable[Any],
    *,
    representatives_per_state: int,
    coefficient_radius: int,
    score_bound: int,
    score: str,
) -> list[dict[str, Any]]:
    records: dict[tuple[int, int], dict[str, Any]] = {}
    for state_index, state in enumerate(states):
        representatives = short_rational_representatives(
            state.residue,
            state.modulus,
            coefficient_radius=coefficient_radius,
            limit=representatives_per_state,
        )
        for representative in representatives:
            numerator = representative.numerator
            denominator = representative.denominator
            t = Fraction(numerator, denominator)
            if not FermigierMestreFamily.discriminant_factor(t):
                continue
            valuations: list[dict[str, int]] = []
            for choice in state.choices:
                residue = residue_of_rational(numerator, denominator, choice.prime)
                if residue is None or residue != choice.residue % choice.prime:
                    raise AssertionError("a rational representative lost a CRT constraint")
                if choice.kind == "power":
                    exponent = 0
                    modulus = choice.modulus
                    while modulus % choice.prime == 0:
                        exponent += 1
                        modulus //= choice.prime
                    actual = FermigierMestreFamily.verify_power_constraint(
                        numerator, denominator, choice.prime, exponent
                    )
                    valuations.append(
                        {
                            "prime": choice.prime,
                            "forced": exponent,
                            "actual": actual,
                        }
                    )
            score_data = full_score(
                numerator,
                denominator,
                bound=score_bound,
                score=score,
            )
            record = {
                "t": rational_to_string(t),
                "numerator": numerator,
                "denominator": denominator,
                "height": representative.height,
                "crt_residue": state.residue,
                "crt_modulus": state.modulus,
                "state_index": state_index,
                "local_constraint_score": state.local_score,
                "score": score_data,
                "forced_valuations": valuations,
                "constraints": [asdict(choice) for choice in state.choices],
            }
            records[(numerator, denominator)] = record
    return sorted(
        records.values(),
        key=lambda item: (-item["score"]["value"], item["height"], item["t"]),
    )


def select_for_pari(records: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    """Select a deterministic mix of high-score and low-height candidates."""

    if count <= 0:
        return []
    by_score = records[:count]
    by_height = sorted(records, key=lambda item: (item["height"], -item["score"]["value"]))[
        :count
    ]
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in by_score + by_height:
        if record["t"] not in seen:
            selected.append(record)
            seen.add(record["t"])
        if len(selected) == count:
            break
    return selected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power-primes", type=parse_power_primes, default=((89, 2), (131, 2)))
    parser.add_argument("--rank-primes", type=parse_primes, default=(7, 11))
    parser.add_argument("--rank-options", type=int, default=2)
    parser.add_argument("--split-only", action="store_true")
    parser.add_argument("--beam-width", type=int, default=128)
    parser.add_argument("--height-weight", type=float, default=0.08)
    parser.add_argument("--coefficient-radius", type=int, default=12)
    parser.add_argument("--representatives-per-state", type=int, default=8)
    parser.add_argument("--score-bound", type=int, default=200)
    parser.add_argument(
        "--score",
        choices=("fermigier-good", "nagao-log"),
        default="fermigier-good",
        help=(
            "good-reduction-only score; 'fermigier-good' uses Fermigier's "
            "summand but a numerical prime bound"
        ),
    )
    parser.add_argument("--keep", type=int, default=30)
    parser.add_argument("--pari-count", type=int, default=8)
    parser.add_argument("--pari-timeout", type=float, default=30.0)
    parser.add_argument("--pari-rank-count", type=int, default=0)
    parser.add_argument("--pari-rank-effort", type=int, default=0)
    default_output = (
        Path(__file__).resolve().parents[2]
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_crt_lattice_pilot.json"
    )
    parser.add_argument("--output", type=Path, default=default_output)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    power_primes = tuple(args.power_primes)
    rank_primes = tuple(args.rank_primes)
    all_primes = [prime for prime, _ in power_primes] + list(rank_primes)
    if len(all_primes) != len(set(all_primes)):
        raise SystemExit("power and rank constraints must use distinct primes")

    groups: list[tuple[ConstraintChoice, ...]] = []
    group_summary: list[dict[str, Any]] = []
    for prime, exponent in power_primes:
        group = power_choices(prime, exponent, split_only=args.split_only)
        if not group:
            raise SystemExit(f"no clean power roots found at p={prime}")
        groups.append(group)
        group_summary.append(
            {"kind": "power", "prime": prime, "exponent": exponent, "choices": len(group)}
        )
    for prime in rank_primes:
        group = rank_choices(prime, count=args.rank_options, score=args.score)
        if not group:
            raise SystemExit(f"no good rank residues found at p={prime}")
        groups.append(group)
        group_summary.append({"kind": "rank", "prime": prime, "choices": len(group)})
    groups.sort(key=lambda group: group[0].modulus)

    states = beam_combine(
        groups,
        beam_width=args.beam_width,
        height_weight=args.height_weight,
        coefficient_radius=args.coefficient_radius,
    )
    records = candidate_records(
        states,
        representatives_per_state=args.representatives_per_state,
        coefficient_radius=args.coefficient_radius,
        score_bound=args.score_bound,
        score=args.score,
    )
    records = records[: args.keep]

    pari_errors: list[dict[str, str]] = []
    pari_selected = select_for_pari(records, args.pari_count)
    rank_budget = args.pari_rank_count
    for index, record in enumerate(pari_selected):
        t = Fraction(record["numerator"], record["denominator"])
        rank_effort = args.pari_rank_effort if index < rank_budget else None
        try:
            known_points = None
            if rank_effort is not None:
                # The first twelve visible images have a known relation.
                # Omitting the first and retaining Mestre's extra section is
                # generically a rank-12 set; PARI rechecks every specialization.
                known_points = FermigierMestreFamily.known_jacobian_points(t)[1:]
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(t),
                timeout=args.pari_timeout,
                rank_effort=rank_effort,
                known_points=known_points,
            )
            log_n = Decimal(record["pari"]["log_conductor"])
            record["below_log_conductor_target"] = log_n < TARGET_LOG_CONDUCTOR
        except Exception as error:  # Preserve a bounded search even if factoring times out.
            pari_errors.append({"t": record["t"], "error": str(error)})

    target_hits = []
    for record in records:
        lower_bound = record.get("pari", {}).get("pari_ellrank", {}).get(
            "lower_bound", 0
        )
        small_conductor_hit = (
            lower_bound >= 21 and record.get("below_log_conductor_target", False)
        )
        record_rank_hit = lower_bound >= 30
        if small_conductor_hit or record_rank_hit:
            target_hits.append(record)
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": "bounded experiment; scores are heuristics and no rank follows from them",
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "family": {
            "name": "normalized Fermigier--Mestre rank-at-least-12 family",
            "root_tuple": list(ROOTS),
            "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        },
        "method": {
            "power_constraints": "simple roots of H(T)=disc_X(R_T)/16, Hensel lifted",
            "globalization": "CRT followed by exact two-dimensional Gauss reduction",
            "rank_point_seed": (
                "when ellrank is requested, supply Jacobian images 2--13 of the "
                "thirteen visible quartic sections"
            ),
            "score": args.score,
            "score_scope": (
                "good reduction at primes p <= score_bound; the historical "
                "Fermigier table instead uses prime ordinal M and includes bad-prime ellap"
            ),
            "score_bound": args.score_bound,
            "group_summary": group_summary,
        },
        "parameters": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(args).items()
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
        "states_retained": len(states),
        "candidates_retained": len(records),
        "pari_errors": pari_errors,
        "target_hits": target_hits,
        "candidates": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(f"states={len(states)} candidates={len(records)} target_hits={len(target_hits)}")
    for record in pari_selected:
        if "pari" in record:
            print(
                f"t={record['t']} height={record['height']} "
                f"score={record['score']['value']:.6f} "
                f"logN={record['pari']['log_conductor']}"
            )


if __name__ == "__main__":
    main()
