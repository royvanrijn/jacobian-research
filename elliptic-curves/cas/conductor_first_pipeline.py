"""Strict stage and Pareto bookkeeping for conductor-first family searches.

This module does not compute conductors, Selmer groups, covers, or points.  It
prevents later-stage records from being promoted without the required exact
earlier gates and computes the rank/conductor Pareto frontier without a blended
heuristic score.
"""

from __future__ import annotations

from math import prod
from typing import Any, Iterable, Mapping, Sequence


def discriminant_sieve_record(
    factorization: Sequence[Sequence[int]], *, complete: bool
) -> dict[str, Any]:
    rows = tuple((int(prime), int(exponent)) for prime, exponent in factorization)
    if not rows or any(prime < 2 or exponent < 1 for prime, exponent in rows):
        raise ValueError("factorization rows require a prime and positive exponent")
    if len({prime for prime, _ in rows}) != len(rows):
        raise ValueError("factorization primes must be distinct")
    radical = prod(prime for prime, _ in rows)
    absolute_value = prod(prime**exponent for prime, exponent in rows)
    squareful_quotient = absolute_value // radical
    repeated = tuple((prime, exponent) for prime, exponent in rows if exponent > 1)
    return {
        "status": "pass",
        "factorization": [[prime, exponent] for prime, exponent in rows],
        "factorization_complete": bool(complete),
        "radical": str(radical),
        "squareful_quotient": str(squareful_quotient),
        "repeated_prime_count": len(repeated),
        "priority_key": [
            radical.bit_length(),
            radical,
            -squareful_quotient.bit_length(),
            -len(repeated),
        ],
        "interpretation": (
            "cheap lexicographic work ordering only; neither a conductor nor a rank score"
        ),
    }


def validate_candidate(record: Mapping[str, Any]) -> str:
    if not record.get("id"):
        raise ValueError("candidate id is required")
    sieve = record.get("cheap_sieve")
    tate = record.get("tate")
    selmer = record.get("residual_selmer")
    covers = record.get("residual_covers")
    recovery = record.get("point_recovery")
    if sieve is None:
        if any(value is not None for value in (tate, selmer, covers, recovery)):
            raise ValueError("later stages require the cheap sieve")
        return "cheap_sieve"
    if sieve.get("status") == "reject":
        if any(value is not None for value in (tate, selmer, covers, recovery)):
            raise ValueError("a sieve rejection cannot have later stages")
        return "rejected"
    if sieve.get("status") != "pass":
        raise ValueError("cheap_sieve status must be pass or reject")
    if tate is None:
        if any(value is not None for value in (selmer, covers, recovery)):
            raise ValueError("Selmer/cover/point stages require exact Tate data")
        return "tate"
    if tate.get("status") != "complete" or not tate.get("global_minimal"):
        if any(value is not None for value in (selmer, covers, recovery)):
            raise ValueError("later stages require completed global minimal Tate data")
        return "tate"
    if not tate.get("conductor") or not tate.get("local_reductions"):
        raise ValueError("completed Tate data require conductor and local reductions")
    if selmer is None:
        if any(value is not None for value in (covers, recovery)):
            raise ValueError("cover/point stages require a residual Selmer result")
        return "residual_selmer"
    if selmer.get("status") != "complete":
        if any(value is not None for value in (covers, recovery)):
            raise ValueError("later stages require a complete residual Selmer result")
        return "residual_selmer"
    residual_dimension = int(selmer["residual_dimension"])
    if residual_dimension < 0:
        raise ValueError("residual dimension cannot be negative")
    if residual_dimension == 0:
        if covers is not None or recovery is not None:
            raise ValueError("a zero residual quotient closes the fibre")
        return "closed"
    if covers is None:
        if recovery is not None:
            raise ValueError("point recovery requires explicit residual covers")
        return "residual_covers"
    if covers.get("status") != "complete":
        if recovery is not None:
            raise ValueError("point recovery requires completed cover construction")
        return "residual_covers"
    if not covers.get("locally_surviving_cover_ids"):
        if recovery is not None:
            raise ValueError("no locally surviving cover is available for point recovery")
        return "locally_closed"
    if recovery is None:
        return "point_recovery"
    if recovery.get("status") not in {"bounded", "certified"}:
        raise ValueError("point_recovery status must be bounded or certified")
    return "certified" if recovery.get("status") == "certified" else "point_recovery"


def work_queues(records: Iterable[Mapping[str, Any]]) -> dict[str, list[str]]:
    queues: dict[str, list[str]] = {}
    for record in records:
        stage = validate_candidate(record)
        queues.setdefault(stage, []).append(str(record["id"]))
    return {stage: sorted(ids) for stage, ids in sorted(queues.items())}


def pareto_frontier(records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    certified = []
    for record in records:
        if validate_candidate(record) != "certified":
            continue
        rank = int(record["point_recovery"]["certified_rank_lower_bound"])
        conductor = int(record["tate"]["conductor"])
        certified.append((record, rank, conductor))
    frontier = []
    for record, rank, conductor in certified:
        dominated = any(
            other_rank >= rank
            and other_conductor <= conductor
            and (other_rank > rank or other_conductor < conductor)
            for other, other_rank, other_conductor in certified
            if other is not record
        )
        if not dominated:
            frontier.append(record)
    return sorted(
        frontier,
        key=lambda record: (
            -int(record["point_recovery"]["certified_rank_lower_bound"]),
            int(record["tate"]["conductor"]),
            str(record["id"]),
        ),
    )


def rank_first_order(records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Order certified candidates by rank first and exact conductor second."""
    certified = [record for record in records if validate_candidate(record) == "certified"]
    return sorted(
        certified,
        key=lambda record: (
            -int(record["point_recovery"]["certified_rank_lower_bound"]),
            int(record["tate"]["conductor"]),
            str(record["id"]),
        ),
    )
