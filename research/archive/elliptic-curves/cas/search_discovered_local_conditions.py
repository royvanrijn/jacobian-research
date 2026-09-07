#!/usr/bin/env python3
"""Discover, optimize, and search Fermigier local power conditions.

The declared prime range is scanned rather than supplied as a hand-picked
list.  Exact p-adic root lifting and ball compression expose congruence
classes with more forced discriminant valuation than their modulus exponent.
Clean split-multiplicative unions are ranked by radical-saving per logarithmic
congruence cost.  A bounded number is selected under a CRT-class cap, after
which the exact CRT/Gauss-lattice experiment follows the same bounded status
discipline as the earlier hand-picked pilot.

No score, conductor computation, or numerical point calculation is promoted
to a Mordell--Weil rank statement.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import prod
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any

from ek_k3 import primes_up_to
from fermigier_mestre import FermigierMestreFamily, ROOTS
from local_condition_discovery import (
    ConditionGroup,
    discover_prime,
    select_condition_groups,
    validate_invariant_polynomials,
)
from pari_bridge import minimal_curve_data, pari_version
from search_multiple_root_crt import (
    LocalConstraintGroup,
    crt_classes,
    enumerate_candidates,
    exact_group_certificate,
    retained_records,
    score_height_pool,
    select_for_pari,
    verify_pari_local_data,
)
from search_record_residue_class import build_score_tables, score_rational


TARGET_LOG_CONDUCTOR = Decimal("182.72")


def convert_group(group: ConditionGroup) -> LocalConstraintGroup:
    scaling_note = (
        f" after {group.presented_model_scaling} minimalizing p-scaling"
        if group.presented_model_scaling
        else ""
    )
    return LocalConstraintGroup(
        prime=group.prime,
        modulus=group.modulus,
        residues=group.residues,
        forced_h_valuation=group.forced_h_valuation,
        reduction=f"{group.reduction}{scaling_note}",
        presented_model_scaling=group.presented_model_scaling,
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-min", type=int, default=5)
    parser.add_argument("--prime-max", type=int, default=199)
    parser.add_argument("--lift-exponent", type=int, default=4)
    parser.add_argument("--classification-exponent", type=int, default=2)
    parser.add_argument("--max-roots", type=int, default=200_000)
    parser.add_argument("--group-count", type=int, default=5)
    parser.add_argument("--maximum-crt-classes", type=int, default=500)
    parser.add_argument("--coefficient-radius", type=int, default=12)
    parser.add_argument("--representatives-per-class", type=int, default=12)
    parser.add_argument("--height-pool", type=int, default=512)
    parser.add_argument("--score-bound", type=int, default=200)
    parser.add_argument("--keep", type=int, default=24)
    parser.add_argument("--pari-count", type=int, default=4)
    parser.add_argument("--pari-timeout", type=float, default=30.0)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_discovered_local_conditions.json"
        ),
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> tuple[int, ...]:
    if args.prime_min < 5 or args.prime_max < args.prime_min:
        raise SystemExit("require 5 <= --prime-min <= --prime-max")
    if args.lift_exponent < 2:
        raise SystemExit("--lift-exponent must be at least 2")
    if args.classification_exponent < 1:
        raise SystemExit("--classification-exponent must be positive")
    positive = {
        "--max-roots": args.max_roots,
        "--group-count": args.group_count,
        "--maximum-crt-classes": args.maximum_crt_classes,
        "--coefficient-radius": args.coefficient_radius,
        "--representatives-per-class": args.representatives_per_class,
        "--height-pool": args.height_pool,
        "--score-bound": args.score_bound,
        "--keep": args.keep,
        "--pari-timeout": args.pari_timeout,
        "--pari-stack-bytes": args.pari_stack_bytes,
    }
    for option, value in positive.items():
        if value <= 0:
            raise SystemExit(f"{option} must be positive")
    if args.score_bound < 5:
        raise SystemExit("--score-bound must be at least 5")
    if args.pari_count < 0:
        raise SystemExit("--pari-count must be nonnegative")
    primes = tuple(
        prime
        for prime in primes_up_to(args.prime_max)
        if prime >= args.prime_min
    )
    if not primes:
        raise SystemExit("the declared interval contains no primes")
    return primes


def main() -> None:
    args = build_parser().parse_args()
    primes = validate_arguments(args)
    validate_invariant_polynomials()

    discovery_start = time.monotonic()
    prime_records: list[dict[str, Any]] = []
    all_groups: list[ConditionGroup] = []
    for prime in primes:
        record = discover_prime(
            prime,
            lift_exponent=args.lift_exponent,
            classification_exponent=args.classification_exponent,
            max_roots=args.max_roots,
        )
        all_groups.extend(record.pop("_groups"))
        prime_records.append(record)
    discovery_seconds = time.monotonic() - discovery_start

    selected_discovered = select_condition_groups(
        all_groups,
        count=args.group_count,
        maximum_crt_classes=args.maximum_crt_classes,
    )
    selected_groups = tuple(convert_group(group) for group in selected_discovered)
    group_certificates = [
        exact_group_certificate(group) for group in selected_groups
    ]

    search_start = time.monotonic()
    classes = crt_classes(selected_groups)
    height_ranked, enumeration = enumerate_candidates(
        classes,
        selected_groups,
        coefficient_radius=args.coefficient_radius,
        representatives_per_class=args.representatives_per_class,
    )
    score_ranked = score_height_pool(
        height_ranked,
        height_pool=args.height_pool,
        score_bound=args.score_bound,
        score="fermigier-good",
    )
    search_seconds = time.monotonic() - search_start

    pari_errors: list[dict[str, str]] = []
    pari_start = time.monotonic()
    pari_selected = select_for_pari(score_ranked, args.pari_count)
    for record, reason in pari_selected:
        parameter = Fraction(record["numerator"], record["denominator"])
        record["pari_selection_reason"] = reason
        try:
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.pari_timeout,
                local_primes=tuple(group.prime for group in selected_groups),
                stack_bytes=args.pari_stack_bytes,
            )
            record["below_log_conductor_target"] = (
                Decimal(record["pari"]["log_conductor"])
                < TARGET_LOG_CONDUCTOR
            )
            record["target_status"] = (
                "conductor inequality only; rank was not computed"
                if record["below_log_conductor_target"]
                else "not a target hit; rank was not computed"
            )
            verify_pari_local_data(record, selected_groups)
        except Exception as error:
            pari_errors.append({"t": record["t"], "error": str(error)})
    pari_seconds = time.monotonic() - pari_start

    retained = retained_records(score_ranked, args.keep)
    retained_keys = {
        (record["numerator"], record["denominator"]) for record in retained
    }
    for record, _ in pari_selected:
        key = (record["numerator"], record["denominator"])
        if key not in retained_keys:
            retained.append(record)
            retained_keys.add(key)
    retained.sort(
        key=lambda record: (
            record["score_rank_within_height_pool"],
            record["height_rank"],
        )
    )
    score_tables = build_score_tables(args.score_bound, "fermigier-good")
    for record in retained:
        record["score"] = score_rational(
            record["numerator"],
            record["denominator"],
            score_tables,
            include_traces=True,
        )

    completed = [record for record, _ in pari_selected if "pari" in record]
    completed.sort(key=lambda record: Decimal(record["pari"]["log_conductor"]))
    classifications = Counter(
        ball["reduction"]
        for prime_record in prime_records
        for ball in prime_record["classified_balls"]
    )
    eligible_groups = [
        group
        for group in all_groups
        if group.reduction == "split multiplicative"
        and group.forced_minimal_discriminant_valuation >= 2
        and group.reciprocal_density > 1
    ]

    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact local-condition discovery and bounded CRT/lattice "
            "experiment; scores are heuristic, conductor-only successes are not "
            "target hits, and no Mordell--Weil rank is computed or inferred"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "family": {
            "name": "normalized Fermigier--Mestre family",
            "root_tuple": list(ROOTS),
            "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        },
        "parameters": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "discovery_method": {
            "prime_interval": [args.prime_min, args.prime_max],
            "primes_scanned": list(primes),
            "root_enumeration": (
                "complete digit lifting of H(T)=0 mod p^k for every 2<=k<=lift_exponent; "
                "the max_roots limit raises rather than truncates"
            ),
            "compression": (
                "complete p-sibling sets are compressed to disjoint maximal balls"
            ),
            "power_efficiency_filter": (
                "retain coarsest balls with forced v_p(H) strictly larger than "
                "their congruence exponent"
            ),
            "classification": (
                "exact invariant-polynomial content, uniform unit tests modulo p, "
                "and split tangent symbols; mixed balls are refined to the declared "
                "classification exponent"
            ),
            "group_cost": (
                "reciprocal density p^e/number_of_residues and its natural logarithm"
            ),
            "benefit": (
                "(forced minimal-discriminant valuation - "
                "max(1, unconditional fixed H-valuation))*log(p), reflecting "
                "one radical/conductor copy at clean multiplicative reduction "
                "without crediting a fixed divisor already present everywhere"
            ),
            "optimization": (
                "greedy descending benefit/log-cost, clean split groups only, "
                "distinct primes, under the declared total CRT-class cap"
            ),
        },
        "discovery_summary": {
            "primes_scanned": len(primes),
            "classified_ball_count": sum(
                len(record["classified_balls"]) for record in prime_records
            ),
            "condition_group_count": len(all_groups),
            "eligible_clean_split_group_count": len(eligible_groups),
            "classification_counts": dict(sorted(classifications.items())),
            "elapsed_seconds": round(discovery_seconds, 6),
        },
        "prime_discovery": prime_records,
        "selected_groups_by_automatic_objective": [
            group.serializable() for group in selected_discovered
        ],
        "selected_group_certificates": group_certificates,
        "search_method": {
            "crt": "all combinations of automatically selected residue unions",
            "lattice": (
                "bounded coefficient box in each exact Gauss-reduced CRT lattice; "
                "T and -T are deduplicated"
            ),
            "score": (
                "good-reduction Fermigier score at numerical primes "
                "5<=p<=score_bound, excluding denominator and bad primes"
            ),
            "pari": "bounded conductor and selected-prime local checks; never ellrank",
        },
        "search_summary": {
            "selected_primes": [group.prime for group in selected_groups],
            "expected_crt_classes": prod(
                len(group.residues) for group in selected_groups
            ),
            "crt_modulus": classes[0]["crt_modulus"],
            **enumeration,
            "height_pool_scored": len(score_ranked),
            "records_output": len(retained),
            "search_elapsed_seconds": round(search_seconds, 6),
            "pari_elapsed_seconds": round(pari_seconds, 6),
            "pari_calls_requested": args.pari_count,
            "pari_calls_completed": len(completed),
            "conductor_only_below_target": sum(
                record["below_log_conductor_target"] for record in completed
            ),
            "best_checked_conductor": (
                {
                    "t": completed[0]["t"],
                    "log_conductor": completed[0]["pari"]["log_conductor"],
                    "rank_status": "not computed",
                }
                if completed
                else None
            ),
        },
        "pari_errors": pari_errors,
        "candidates": retained,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    print(f"wrote {args.output}")
    print(
        f"scanned_primes={len(primes)} groups={len(all_groups)} "
        f"selected={[group.prime for group in selected_groups]} "
        f"crt_classes={len(classes)} candidates={len(height_ranked)}"
    )
    for record, reason in pari_selected:
        print(
            f"t={record['t']} reason={reason} "
            f"logN={record.get('pari', {}).get('log_conductor')} "
            f"below_target={record.get('below_log_conductor_target')}"
        )


if __name__ == "__main__":
    main()
