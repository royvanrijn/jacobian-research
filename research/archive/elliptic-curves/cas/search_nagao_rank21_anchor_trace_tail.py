#!/usr/bin/env python3
"""Preserve the trace/root tail omitted by the broad anchor-union cutoff.

The parent anchor-union experiment generated several thousand candidates, so
its global B=200 retention was dominated by smooth/even-denominator controls.
This complementary lane regenerates only CRT strata containing exact trace
conditions, retains every proxy-feasible member through exact B=2000, and
then gives each trace/root stratum a conductor quota.  Every fiber point-
searched by the parent artifact is excluded exactly.

Subtarget root-number -1 curves have point-search priority.  Numerical ranks
remain triage; stable rank at least 18 triggers exact finite reductions.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Sequence

from ek_k3 import rational_to_string
from pari_bridge import pari_version
from search_nagao_rank21_anchor_union import (
    ANCHORS,
    DESIGN_PRIMES,
    PROXY_LIMIT,
    ROOT_TARGETS,
    TARGET_LOG_CONDUCTOR,
    TRACE_BEAM_WIDTH,
    GeneratedCandidate,
    build_residue_tables,
    build_strata,
    build_trace_beam,
    conductor_radical_proxy,
    exact_checkpoints,
    exact_decontaminated_scores,
    gauss_shell,
    generated_record,
    learn_trace_fingerprint,
    outside_old_box,
    parallel_conductors,
    prior_parameter_exclusions,
    projective_index,
    root_ball_union,
    select_point_population,
    sha256_file,
    staged_points,
)


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
BROAD_ARTIFACT = (
    REPOSITORY / "artifacts/generated-results/elliptic_nagao_rank21_anchor_union.json"
)
BROAD_ARTIFACT_SHA256 = "c518ee794118c74c94c9621f702d18047c11d6c37bc773a0b5e7134d1616801a"
CONDUCTOR_KEEP = 96
POINT_KEEP = 32
TARGET_RANK = 21
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_anchor_trace_tail.py"
)


def parent_point_exclusions() -> tuple[frozenset[Fraction], dict[str, Any]]:
    if sha256_file(BROAD_ARTIFACT) != BROAD_ARTIFACT_SHA256:
        raise AssertionError("the pinned broad anchor-union artifact changed")
    data = json.loads(BROAD_ARTIFACT.read_text(encoding="utf-8"))
    excluded = {
        abs(Q(record["constructor_parameter"]))
        for stage in data["point_stages"]
        for record in stage["ranked_population"]
    }
    if len(excluded) != data["point_population_selection"]["count"]:
        raise AssertionError("the parent point population no longer deduplicates exactly")
    return frozenset(excluded), {
        "path": str(BROAD_ARTIFACT),
        "sha256": BROAD_ARTIFACT_SHA256,
        "exact_point_searched_fibers_excluded": len(excluded),
        "constructor_parameters": sorted(rational_to_string(value) for value in excluded),
    }


def trace_only_candidates(
    strata: dict[str, tuple[Any, ...]],
    exclusions: frozenset[Fraction],
    *,
    proxy_limit: Decimal,
) -> tuple[tuple[GeneratedCandidate, ...], dict[str, Any]]:
    retained: dict[Fraction, GeneratedCandidate] = {}
    raw = old_box = prior = proxy_rejected = singular = 0
    for stratum, states in strata.items():
        if "trace" not in stratum:
            continue
        for state in states:
            shell = gauss_shell(state.residue, state.modulus, radius=4, limit=6)
            indices = (0, 2, 5) if "trace-full" in stratum and "root" not in stratum else (0, 3)
            for index in indices:
                if index >= len(shell):
                    continue
                raw += 1
                parameter = abs(shell[index][0])
                if not outside_old_box(parameter):
                    old_box += 1
                    continue
                if parameter in exclusions:
                    prior += 1
                    continue
                try:
                    proxy = conductor_radical_proxy(parameter)
                except ValueError:
                    singular += 1
                    continue
                if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
                    proxy_rejected += 1
                    continue
                traces = tuple(
                    (
                        choice.prime,
                        projective_index(parameter.numerator, parameter.denominator, choice.prime),
                        choice.threshold,
                    )
                    for choice in state.trace_choices
                )
                roots = tuple(
                    (
                        ball.prime,
                        parameter.numerator * pow(parameter.denominator, -1, ball.modulus) % ball.modulus,
                        ball.modulus,
                        ball.forced_valuation,
                    )
                    for ball in state.root_balls
                )
                candidate = GeneratedCandidate(parameter, stratum, traces, roots, proxy)
                existing = retained.get(parameter)
                if existing is None or (
                    len(candidate.root_conditions),
                    len(candidate.trace_conditions),
                    candidate.stratum,
                ) > (
                    len(existing.root_conditions),
                    len(existing.trace_conditions),
                    existing.stratum,
                ):
                    retained[parameter] = candidate
    answer = tuple(
        sorted(
            retained.values(),
            key=lambda candidate: (
                candidate.proxy["log_radical_upper_proxy"],
                candidate.height,
                candidate.identifier,
            ),
        )
    )
    digest = hashlib.sha256()
    for candidate in answer:
        digest.update(
            f"{candidate.parameter}|{candidate.stratum}|{candidate.proxy['log_radical_upper_proxy']!r}\n".encode()
        )
    return answer, {
        "raw_trace_stratum_records": raw,
        "excluded_old_box": old_box,
        "excluded_prior_or_parent_point_fibers": prior,
        "proxy_rejected": proxy_rejected,
        "singular": singular,
        "exactly_deduplicated_proxy_survivors": len(answer),
        "stratum_counts": {
            stratum: sum(candidate.stratum == stratum for candidate in answer)
            for stratum in sorted({candidate.stratum for candidate in answer})
        },
        "survivor_stream_sha256": digest.hexdigest(),
    }


def select_conductor_population(candidates: Sequence[Any], *, keep: int) -> tuple[Any, ...]:
    selected = {candidate.identifier: candidate for candidate in candidates[:48]}
    by_proxy = sorted(
        candidates,
        key=lambda candidate: (
            candidate.generated.proxy["log_radical_upper_proxy"],
            -Decimal(candidate.exact_score_b2000),
            candidate.identifier,
        ),
    )
    for candidate in by_proxy[:24]:
        selected[candidate.identifier] = candidate
    for stratum in sorted({candidate.generated.stratum for candidate in candidates}):
        quota = 0
        for candidate in candidates:
            if candidate.generated.stratum != stratum:
                continue
            selected[candidate.identifier] = candidate
            quota += 1
            if quota == 2:
                break
    # Stratum leaders are inserted first in the final ordering, then the
    # score/proxy union fills the declared cap.
    stratum_ids = {
        candidate.identifier
        for stratum in {candidate.generated.stratum for candidate in candidates}
        for candidate in [
            item
            for item in candidates
            if item.generated.stratum == stratum
        ][:2]
    }
    ordered = sorted(
        selected.values(),
        key=lambda candidate: (
            candidate.identifier not in stratum_ids,
            -Decimal(candidate.exact_score_b2000),
            candidate.generated.proxy["log_radical_upper_proxy"],
            candidate.identifier,
        ),
    )
    return tuple(ordered[: min(keep, len(ordered))])


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    prior, prior_audit = prior_parameter_exclusions()
    parent, parent_audit = parent_point_exclusions()
    exclusions = frozenset(set(prior) | set(parent))
    tables = build_residue_tables(200)
    trace_data = {}
    trace_audit = {}
    for anchor in ANCHORS:
        fingerprint = learn_trace_fingerprint(anchor, tables)
        trace_data[anchor.label] = build_trace_beam(
            fingerprint,
            tables,
            width=args.trace_beam_width,
        )
        trace_audit[anchor.label] = {
            "constructor_parameter": rational_to_string(anchor.parameter),
            "fingerprint": [list(value) for value in fingerprint],
            "beam_stages": list(trace_data[anchor.label][2]),
        }
    roots = {prime: root_ball_union(prime, target) for prime, target in ROOT_TARGETS}
    strata = build_strata(trace_data, roots)
    generated, generation = trace_only_candidates(
        strata,
        exclusions,
        proxy_limit=args.proxy_limit,
    )
    scoring_tables = {
        prime: table for prime, table in tables.items() if prime not in DESIGN_PRIMES
    }
    # The trace tail is small enough that no second heuristic cutoff is made.
    from search_nagao_rank21_anchor_union import prefilter_candidates

    prefiltered = prefilter_candidates(
        generated,
        scoring_tables,
        keep=len(generated),
    )
    exact = exact_decontaminated_scores(
        prefiltered,
        cutoff=2_000,
        batch_size=args.score_batch_size,
        timeout=args.score_timeout,
        stack_bytes=args.stack_bytes,
    )
    conductor_candidates = select_conductor_population(
        exact,
        keep=args.conductor_keep,
    )
    conductors = parallel_conductors(
        conductor_candidates,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.workers,
    )
    conductor_by_id = {replay.candidate.identifier: replay for replay in conductors}
    point_candidates = select_point_population(conductors, keep=args.point_keep)
    stages, best = staged_points(point_candidates, args)
    checkpoints = exact_checkpoints(best, conductor_by_id, args)
    hits = [record for record in checkpoints if record["target_rank21_under_log_conductor_hit"]]
    return {
        "schema_version": 1,
        "status": "bounded complementary trace/root-tail search; numerical ranks are triage only",
        "target": {
            "rank_at_least": TARGET_RANK,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "certified_hits": hits,
        },
        "parent_broad_search": parent_audit,
        "prior_population_exclusion": prior_audit,
        "trace_design": trace_audit,
        "design_primes_omitted_from_scores": list(DESIGN_PRIMES),
        "generation": generation,
        "exact_b2000_population": [generated_record(candidate) for candidate in exact],
        "conductor_population": [
            {
                **generated_record(replay.candidate),
                "status": replay.status,
                "conductor": replay.data if replay.status == "completed" else None,
                "error": replay.error,
            }
            for replay in conductors
        ],
        "point_population_selection": {
            "count": len(point_candidates),
            "constructor_parameters": [
                rational_to_string(candidate.parameter) for candidate in point_candidates
            ],
            "subtarget_root_minus_one_count": sum(
                conductor_by_id[candidate.identifier].status == "completed"
                and conductor_by_id[candidate.identifier].data.get("below_strict_log_conductor_target") is True
                and conductor_by_id[candidate.identifier].data.get("root_number") == -1
                for candidate in point_candidates
            ),
        },
        "point_stages": stages,
        "exact_checkpoints_stable_numerical_rank_at_least_18": checkpoints,
        "bounds_and_caveats": {
            "proxy_limit": str(args.proxy_limit),
            "no_b200_candidate_cut_after_trace_generation": True,
            "conductor_keep": args.conductor_keep,
            "point_keep": args.point_keep,
            "bounded_search_is_not_a_rank_upper_bound": True,
            "root_number_priority_is_a_parity_heuristic": True,
            "all_subprocesses_synchronous_with_finite_timeouts": True,
        },
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-beam-width", type=int, default=TRACE_BEAM_WIDTH)
    parser.add_argument("--proxy-limit", type=Decimal, default=PROXY_LIMIT)
    parser.add_argument("--score-batch-size", type=int, default=40)
    parser.add_argument("--score-timeout", type=float, default=40.0)
    parser.add_argument("--conductor-keep", type=int, default=CONDUCTOR_KEEP)
    parser.add_argument("--conductor-timeout", type=float, default=60.0)
    parser.add_argument("--point-keep", type=int, default=POINT_KEEP)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic_nagao_rank21_anchor_trace_tail.json"
        ),
    )
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        raise SystemExit("--workers must lie in [1,4]")
    if any(
        not 0 < value <= 120
        for value in (
            args.score_timeout,
            args.conductor_timeout,
            args.height_timeout,
            args.saturation_timeout,
        )
    ):
        raise SystemExit("all subprocess timeouts must lie in (0,120]")
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
