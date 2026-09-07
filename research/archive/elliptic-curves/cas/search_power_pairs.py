#!/usr/bin/env python3
"""Bounded two-prime discriminant-power search in the Fermigier family.

For every pair of distinct primes in a configured interval, this experiment
chooses one clean simple root of the Fermigier discriminant factor modulo
``p^k`` at each prime, combines the roots by CRT, and searches the resulting
two-dimensional lattice for short rational parameters ``T=a/b``.

Selection deliberately has two stages: retain the lowest-height rational
parameters first, then rank only that bounded pool with a good-reduction
Mestre--Nagao score.  Optional PARI/GP calls compute conductors for a fixed
mixture of the lowest-height and highest-score candidates.  No rank is
computed or inferred by this script.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable

from crt_lattice import crt_pair, short_rational_representatives
from ek_k3 import primes_up_to, rational_to_string
from fermigier_mestre import FermigierMestreFamily, FermigierPowerRoot, ROOTS
from pari_bridge import minimal_curve_data, pari_version
from search_crt_lattice import full_score


TARGET_LOG_CONDUCTOR = Decimal("182.72")


def root_record(root: FermigierPowerRoot) -> dict[str, Any]:
    """Serialize a lifted root without assigning it a rank interpretation."""

    record = asdict(root)
    record["reduction"] = (
        "split multiplicative"
        if root.split_multiplicative
        else "nonsplit multiplicative"
    )
    return record


def enumerate_root_groups(
    *, prime_min: int, prime_max: int, exponent: int, split_only: bool
) -> tuple[tuple[FermigierPowerRoot, ...], ...]:
    """Return the nonempty groups of clean lifted roots, one group per prime."""

    groups: list[tuple[FermigierPowerRoot, ...]] = []
    for prime in primes_up_to(prime_max):
        if prime < prime_min or prime < 5:
            continue
        roots = FermigierMestreFamily.power_roots(
            prime, exponent, split_only=split_only
        )
        if roots:
            groups.append(roots)
    return tuple(groups)


def pair_paths(
    groups: tuple[tuple[FermigierPowerRoot, ...], ...]
) -> Iterable[tuple[FermigierPowerRoot, FermigierPowerRoot]]:
    """Yield every root choice at every unordered pair of distinct primes."""

    for left_index, left_group in enumerate(groups):
        for right_group in groups[left_index + 1 :]:
            for left_root in left_group:
                for right_root in right_group:
                    yield left_root, right_root


def enumerate_candidates(
    groups: tuple[tuple[FermigierPowerRoot, ...], ...],
    *,
    coefficient_radius: int,
    representatives_per_pair: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Enumerate and deduplicate nonsingular rational representatives."""

    candidates: dict[tuple[int, int], dict[str, Any]] = {}
    root_choice_pairs = 0
    representatives_generated = 0
    singular_representatives = 0
    symmetry_or_duplicate_merges = 0
    for left_root, right_root in pair_paths(groups):
        root_choice_pairs += 1
        residue, modulus = crt_pair(
            left_root.residue,
            left_root.modulus,
            right_root.residue,
            right_root.modulus,
        )
        representatives = short_rational_representatives(
            residue,
            modulus,
            coefficient_radius=coefficient_radius,
            limit=representatives_per_pair,
        )
        for representative in representatives:
            representatives_generated += 1
            numerator = representative.numerator
            path_left = left_root
            path_right = right_root
            path_residue = residue
            # Every coefficient of this Jacobian is even in T.  Store one
            # representative of the isomorphic T <-> -T pair, while also
            # negating the lifted roots so the retained congruences replay.
            if numerator < 0:
                numerator = -numerator
                path_left = FermigierPowerRoot(
                    left_root.prime,
                    left_root.exponent,
                    left_root.modulus,
                    (-left_root.residue) % left_root.modulus,
                    left_root.split_multiplicative,
                )
                path_right = FermigierPowerRoot(
                    right_root.prime,
                    right_root.exponent,
                    right_root.modulus,
                    (-right_root.residue) % right_root.modulus,
                    right_root.split_multiplicative,
                )
                path_residue = (-residue) % modulus
            parameter = Fraction(numerator, representative.denominator)
            if FermigierMestreFamily.discriminant_factor(parameter) == 0:
                singular_representatives += 1
                continue
            key = (numerator, representative.denominator)
            path = {
                "roots": [root_record(path_left), root_record(path_right)],
                "crt_residue": path_residue,
                "crt_modulus": modulus,
            }
            previous = candidates.get(key)
            if previous is None:
                candidates[key] = {
                    "t": rational_to_string(parameter),
                    "numerator": numerator,
                    "denominator": representative.denominator,
                    "height": representative.height,
                    "forcing_paths": [path],
                }
            else:
                symmetry_or_duplicate_merges += 1
                if path not in previous["forcing_paths"]:
                    previous["forcing_paths"].append(path)

    ordered = sorted(
        candidates.values(),
        key=lambda item: (
            item["height"],
            abs(item["numerator"]),
            item["denominator"],
            item["t"],
        ),
    )
    counts = {
        "root_choice_pairs": root_choice_pairs,
        "representatives_generated": representatives_generated,
        "singular_representatives_rejected": singular_representatives,
        "sign_symmetry_or_duplicate_merges": symmetry_or_duplicate_merges,
        "unique_nonsingular_representatives": len(ordered),
    }
    return ordered, counts


def add_exact_constraints(record: dict[str, Any], exponent: int) -> None:
    """Replay the first forcing path and record exact achieved valuations."""

    roots = record["forcing_paths"][0]["roots"]
    valuations = []
    for root in roots:
        actual = FermigierMestreFamily.verify_power_constraint(
            record["numerator"],
            record["denominator"],
            root["prime"],
            exponent,
        )
        valuations.append(
            {
                "prime": root["prime"],
                "forced_at_least": exponent,
                "actual": actual,
            }
        )
    record["forced_valuations"] = valuations


def score_height_pool(
    height_ranked: list[dict[str, Any]],
    *,
    height_pool: int,
    exponent: int,
    score_bound: int,
    score: str,
) -> list[dict[str, Any]]:
    """Score only the bounded lowest-height prefix and assign both ranks."""

    pool = height_ranked[:height_pool]
    for height_rank, record in enumerate(pool, start=1):
        record["height_rank"] = height_rank
        add_exact_constraints(record, exponent)
        record["score"] = full_score(
            record["numerator"],
            record["denominator"],
            bound=score_bound,
            score=score,
        )
    score_ranked = sorted(
        pool,
        key=lambda item: (
            -item["score"]["value"],
            item["height"],
            abs(item["numerator"]),
            item["denominator"],
        ),
    )
    for score_rank, record in enumerate(score_ranked, start=1):
        record["score_rank_within_height_pool"] = score_rank
    return score_ranked


def select_for_pari(
    score_ranked: list[dict[str, Any]], count: int
) -> list[tuple[dict[str, Any], str]]:
    """Select a deterministic half-height, half-score bounded PARI subset."""

    if count <= 0:
        return []
    height_ranked = sorted(
        score_ranked,
        key=lambda item: (
            item["height"],
            abs(item["numerator"]),
            item["denominator"],
        ),
    )
    selected: list[tuple[dict[str, Any], str]] = []
    seen: set[tuple[int, int]] = set()
    sources = ((height_ranked, "lowest height"), (score_ranked, "highest score"))
    index = 0
    while len(selected) < min(count, len(score_ranked)):
        made_progress = False
        for records, reason in sources:
            if index >= len(records):
                continue
            record = records[index]
            key = (record["numerator"], record["denominator"])
            if key not in seen:
                selected.append((record, reason))
                seen.add(key)
                made_progress = True
                if len(selected) == count:
                    break
        index += 1
        if not made_progress and index >= len(score_ranked):
            break
    return selected


def retained_records(
    score_ranked: list[dict[str, Any]], keep: int
) -> list[dict[str, Any]]:
    """Retain the union of the best height and score prefixes."""

    if keep <= 0:
        return []
    height_ranked = sorted(score_ranked, key=lambda item: item["height_rank"])
    keys = {
        (record["numerator"], record["denominator"])
        for record in score_ranked[:keep] + height_ranked[:keep]
    }
    return [
        record
        for record in score_ranked
        if (record["numerator"], record["denominator"]) in keys
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-min", type=int, default=5)
    parser.add_argument("--prime-max", type=int, default=500)
    parser.add_argument("--exponent", type=int, default=2)
    parser.add_argument("--split-only", action="store_true")
    parser.add_argument("--coefficient-radius", type=int, default=8)
    parser.add_argument("--representatives-per-pair", type=int, default=4)
    parser.add_argument(
        "--height-pool",
        type=int,
        default=256,
        help="number of lowest-height candidates admitted to score ranking",
    )
    parser.add_argument("--score-bound", type=int, default=200)
    parser.add_argument(
        "--score",
        choices=("fermigier-good", "nagao-log"),
        default="fermigier-good",
    )
    parser.add_argument(
        "--keep",
        type=int,
        default=40,
        help="retain the union of this many lowest-height and highest-score records",
    )
    parser.add_argument(
        "--pari-count",
        type=int,
        default=0,
        help=(
            "optional bounded conductor calls, alternating lowest height and "
            "highest score; disabled by default"
        ),
    )
    parser.add_argument("--pari-timeout", type=float, default=20.0)
    default_output = (
        Path(__file__).resolve().parents[2]
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_power_pairs.json"
    )
    parser.add_argument("--output", type=Path, default=default_output)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.prime_min < 5:
        raise SystemExit("--prime-min must be at least 5 for the normalized model")
    if args.prime_max < args.prime_min:
        raise SystemExit("--prime-max must be at least --prime-min")
    if args.exponent < 2:
        raise SystemExit("--exponent must be at least 2 in a power-pair search")
    positive = {
        "--coefficient-radius": args.coefficient_radius,
        "--representatives-per-pair": args.representatives_per_pair,
        "--height-pool": args.height_pool,
        "--score-bound": args.score_bound,
        "--keep": args.keep,
    }
    for name, value in positive.items():
        if value < 1:
            raise SystemExit(f"{name} must be positive")
    if args.pari_count < 0:
        raise SystemExit("--pari-count must be nonnegative")
    if args.pari_timeout <= 0:
        raise SystemExit("--pari-timeout must be positive")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)

    groups = enumerate_root_groups(
        prime_min=args.prime_min,
        prime_max=args.prime_max,
        exponent=args.exponent,
        split_only=args.split_only,
    )
    if len(groups) < 2:
        raise SystemExit("fewer than two primes have clean roots in this interval")

    height_ranked, enumeration = enumerate_candidates(
        groups,
        coefficient_radius=args.coefficient_radius,
        representatives_per_pair=args.representatives_per_pair,
    )
    if not height_ranked:
        raise SystemExit("the bounded lattice enumeration found no candidate")
    score_ranked = score_height_pool(
        height_ranked,
        height_pool=args.height_pool,
        exponent=args.exponent,
        score_bound=args.score_bound,
        score=args.score,
    )

    pari_errors: list[dict[str, str]] = []
    pari_selected = select_for_pari(score_ranked, args.pari_count)
    for record, reason in pari_selected:
        roots = record["forcing_paths"][0]["roots"]
        forced_primes = tuple(root["prime"] for root in roots)
        record["pari_selection_reason"] = reason
        try:
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(
                    Fraction(record["numerator"], record["denominator"])
                ),
                timeout=args.pari_timeout,
                local_primes=forced_primes,
            )
            log_conductor = Decimal(record["pari"]["log_conductor"])
            record["below_log_conductor_bound"] = (
                log_conductor < TARGET_LOG_CONDUCTOR
            )
            record["target_status"] = (
                "conductor inequality only; no rank computation or claim"
                if record["below_log_conductor_bound"]
                else "not a target hit; no rank computation or claim"
            )
        except Exception as error:  # Keep the bounded run if factoring times out.
            pari_errors.append({"t": record["t"], "error": str(error)})

    records = retained_records(score_ranked, args.keep)
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded experiment; the score is heuristic and no Mordell-Weil "
            "rank is computed or claimed"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
            "explanation": (
                "a conductor below the bound is not a hit without 21 "
                "independent rational points"
            ),
        },
        "family": {
            "name": "normalized Fermigier--Mestre rank-at-least-12 family",
            "root_tuple": list(ROOTS),
            "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        },
        "method": {
            "power_constraints": (
                "all clean simple roots of H(T)=disc_X(R_T)/16, Hensel "
                "lifted to the configured exponent"
            ),
            "globalization": (
                "every distinct-prime root pair, CRT, then exact "
                "two-dimensional Gauss reduction"
            ),
            "parameter_symmetry": (
                "the Jacobian coefficients are even in T, so T and -T are "
                "deduplicated by retaining a nonnegative numerator"
            ),
            "selection_order": (
                "retain the lowest-height pool first; rank only that pool "
                "by the configured good-reduction score"
            ),
            "score": args.score,
            "score_bound": args.score_bound,
            "pari_subset": (
                "deterministic alternation between lowest height and highest "
                "score; conductor only, never ellrank"
            ),
        },
        "parameters": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in vars(args).items()
        },
        "root_groups": [
            {
                "prime": group[0].prime,
                "root_count": len(group),
                "roots": [root_record(root) for root in group],
            }
            for group in groups
        ],
        "enumeration": enumeration,
        "height_pool_scored": len(score_ranked),
        "records_retained": len(records),
        "pari_calls_requested": args.pari_count,
        "pari_calls_completed": sum("pari" in record for record, _ in pari_selected),
        "pari_errors": pari_errors,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "candidates": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    print(f"wrote {args.output}")
    print(
        f"root_primes={len(groups)} root_choice_pairs={enumeration['root_choice_pairs']} "
        f"unique={enumeration['unique_nonsingular_representatives']} "
        f"height_pool={len(score_ranked)} retained={len(records)}"
    )
    for record, reason in pari_selected:
        summary = (
            f"t={record['t']} height={record['height']} "
            f"height_rank={record['height_rank']} "
            f"score_rank={record['score_rank_within_height_pool']} "
            f"selected={reason!r}"
        )
        if "pari" in record:
            summary += f" logN={record['pari']['log_conductor']}"
        else:
            summary += " logN=unavailable"
        print(summary)


if __name__ == "__main__":
    main()
