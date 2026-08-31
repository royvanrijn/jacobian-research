#!/usr/bin/env sage-python
"""Resumably search and certify ranks of pointed pair-base Jacobians.

The input may be the 5,566-row zero/infinity immediate-point catalogue or the
300-row t=3/8 control-selected catalogue. Sage/mwrank supplies candidate
generators. This script then certifies their independence by exact finite
quotients; an uncertified search result is retained as such and never promoted
to a rank lower bound.

Repeated runs append new arithmetic-complexity ranks to the same output.  Use
``--start`` and ``--limit`` to split a long scan into reproducible intervals.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from tempfile import NamedTemporaryFile

from sage.all import EllipticCurve, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/elkies-2026-immediate-point-pair-catalogue-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-immediate-point-pair-rank-ledger.json"
INPUT_SCHEMA = "elkies-k3.elkies-2026-immediate-point-pair-catalogue.v1"
OUTPUT_SCHEMA = "elkies-k3.elkies-2026-immediate-point-pair-rank-ledger.v1"
CONTROL_INPUT_SCHEMA = "elkies-k3.elkies-2026-control-pair-base-catalogue.v1"
CONTROL_OUTPUT_SCHEMA = "elkies-k3.elkies-2026-control-pair-base-rank-ledger.v1"
FINITE_QUOTIENT_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
SHORT_MODEL_HELPER = ROOT / "elliptic-curves/ecsearch/q12o5867_specialization.py"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def fraction(value) -> Fraction:
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def atomic_write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", dir=path.parent, prefix=path.name + ".", delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def timeout_call(seconds, function):
    """Run a Sage arithmetic call with a process-local real-time bound."""

    if seconds <= 0:
        return function()
    previous_handler = signal.getsignal(signal.SIGALRM)

    def timeout_handler(_signum, _frame):
        raise TimeoutError(f"arithmetic call exceeded {seconds} seconds")

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return function()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]
from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    source_point_to_target,
    verify_finite_quotient_certificate,
)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--start", type=int, default=1, help="first one-based arithmetic complexity rank")
parser.add_argument("--limit", type=int, default=100)
parser.add_argument("--prime-bound", type=int, default=1000)
parser.add_argument("--pari-effort", type=int, default=3)
parser.add_argument("--search-timeout", type=int, default=30)
parser.add_argument(
    "--backend",
    choices=("mwrank-first", "pari-only"),
    default="mwrank-first",
    help="use PARI directly for a fast bounded baseline, or try mwrank first",
)
parser.add_argument("--checkpoint-every", type=int, default=10)
parser.add_argument("--retry-failures", action="store_true")
parser.add_argument(
    "--rebind-compatible-input",
    action="store_true",
    help="accept a new catalogue hash only after every stored row invariant is unchanged",
)
args = parser.parse_args()

if (
    args.start < 1
    or args.limit < 0
    or args.prime_bound < 3
    or args.pari_effort < 0
    or args.search_timeout < 0
    or args.checkpoint_every < 1
):
    raise ValueError("invalid nonnegative interval or certificate bound")

catalogue_sha = digest(args.input)
catalogue = json.loads(args.input.read_text())
if catalogue.get("schema") not in (INPUT_SCHEMA, CONTROL_INPUT_SCHEMA):
    raise ValueError(f"unexpected catalogue schema: {catalogue.get('schema')!r}")
output_schema = CONTROL_OUTPUT_SCHEMA if catalogue["schema"] == CONTROL_INPUT_SCHEMA else OUTPUT_SCHEMA

input_fingerprints = {
    display_path(args.input): catalogue_sha,
    display_path(FINITE_QUOTIENT_HELPER): digest(FINITE_QUOTIENT_HELPER),
    display_path(SHORT_MODEL_HELPER): digest(SHORT_MODEL_HELPER),
}
if args.output.exists():
    ledger = json.loads(args.output.read_text())
    if ledger.get("schema") != output_schema:
        raise ValueError("existing ledger has an incompatible schema")
    if ledger.get("inputs") != input_fingerprints:
        if not args.rebind_compatible_input:
            raise ValueError("existing ledger has incompatible input hashes")
        old_inputs = ledger.get("inputs", {})
        for helper in (FINITE_QUOTIENT_HELPER, SHORT_MODEL_HELPER):
            name = display_path(helper)
            if old_inputs.get(name) != input_fingerprints[name]:
                raise ValueError(f"cannot rebind after helper change: {name}")
        new_rows = {row["pair_key"]: row for row in catalogue["pairs"]}
        for key, result in ledger.get("results", {}).items():
            row = new_rows.get(key)
            if row is None:
                raise ValueError(f"stored pair disappeared from new catalogue: {key}")
            invariants = (
                ("arithmetic_complexity_rank", int(row["arithmetic_complexity_rank"])),
                ("global_root_number", int(row["global_root_number"])),
                ("minimal_jacobian_a1_a2_a3_a4_a6", row["minimal_jacobian_a1_a2_a3_a4_a6"]),
            )
            for name, expected in invariants:
                if result.get(name) != expected:
                    raise ValueError(f"stored result changed under catalogue rebind: {key}:{name}")
        ledger["inputs"] = input_fingerprints
        ledger.setdefault("compatible_input_rebindings", []).append(
            {
                "old_catalogue_sha256": old_inputs.get(display_path(args.input)),
                "new_catalogue_sha256": catalogue_sha,
                "verified_stored_result_count": len(ledger.get("results", {})),
                "invariants": [
                    "pair_key",
                    "arithmetic_complexity_rank",
                    "global_root_number",
                    "minimal_jacobian_a1_a2_a3_a4_a6",
                    "finite_quotient_certificates unchanged because helper hashes are unchanged",
                ],
            }
        )
else:
    ledger = {
        "schema": output_schema,
        "status": "PARTIAL_EXACT_RANK_LOWER_BOUND_LEDGER",
        "inputs": input_fingerprints,
        "method": {
            "generator_search": (
                "Sage EllipticCurve.gens(proof=False), mwrank_lib first and lower-level PARI "
                "ellrank partial-point fallback after a search exception"
            ),
            "certification": "exact finite-quotient infinite descent",
            "relation_primes_tried": [3, 5, 7],
            "proof_boundary": (
                "Only finite-quotient-certified independent points contribute to a rank lower "
                "bound. Search failures and uncertified generators are not rank bounds; no upper "
                "rank is asserted."
            ),
        },
        "results": {},
    }

# Keep the human-readable method description current across resumptions; the
# mathematical compatibility gate is supplied by the pinned input/helper
# hashes above.
ledger["method"]["generator_search"] = (
    "Sage EllipticCurve.gens(proof=False), mwrank_lib first and lower-level PARI "
    "ellrank partial-point fallback after a search exception"
)

rows = catalogue["pairs"]
stop = min(len(rows), args.start - 1 + args.limit)
selected = rows[args.start - 1 : stop]
completed_this_run = 0

for row in selected:
    key = row["pair_key"]
    previous = ledger["results"].get(key)
    if previous is not None and not (args.retry_failures and previous["status"] != "CERTIFIED"):
        continue

    model = tuple(Fraction(value) for value in row["minimal_jacobian_a1_a2_a3_a4_a6"])
    curve = EllipticCurve(QQ, list(model))
    try:
        primary_error = None
        pari_rank_bounds = None
        if args.backend == "pari-only":
            pari_result = timeout_call(
                args.search_timeout,
                lambda: curve.pari_curve().ellrank(args.pari_effort),
            )
            pari_rank_bounds = [int(pari_result[0]), int(pari_result[1])]
            generators = tuple(
                curve(QQ(str(point[0])), QQ(str(point[1]))) for point in pari_result[3]
            )
            search_algorithm = "pari_ellrank_bounded_baseline"
        else:
            try:
                generators = tuple(
                    timeout_call(
                        args.search_timeout,
                        lambda: curve.gens(proof=False, algorithm="mwrank_lib"),
                    )
                )
                search_algorithm = "mwrank_lib"
            except Exception as error:
                primary_error = {"type": type(error).__name__, "message": str(error)}
                pari_result = timeout_call(
                    args.search_timeout,
                    lambda: curve.pari_curve().ellrank(args.pari_effort),
                )
                pari_rank_bounds = [int(pari_result[0]), int(pari_result[1])]
                generators = tuple(
                    curve(QQ(str(point[0])), QQ(str(point[1]))) for point in pari_result[3]
                )
                search_algorithm = "pari_ellrank_partial_fallback"
        points = tuple((fraction(point[0]), fraction(point[1])) for point in generators)
        certificate = None
        if points:
            short_model, change = short_certificate_model(model)
            short_points = tuple(source_point_to_target(point, change) for point in points)
            for relation_prime in (3, 5, 7):
                candidate = build_finite_quotient_certificate(
                    short_model,
                    short_points,
                    relation_prime=relation_prime,
                    prime_bound=args.prime_bound,
                )
                if candidate["certified_independent"]:
                    verify_finite_quotient_certificate(short_model, short_points, candidate)
                    certificate = candidate
                    break
        else:
            certificate = {
                "certificate_type": "empty-point-lower-bound",
                "point_count": 0,
                "certified_independent": True,
                "certified_rank_lower_bound": 0,
            }

        certified = certificate is not None and certificate["certified_independent"]
        result = {
            "status": "CERTIFIED" if certified else "GENERATORS_UNCERTIFIED",
            "pair_key": key,
            "arithmetic_complexity_rank": int(row["arithmetic_complexity_rank"]),
            "global_root_number": int(row["global_root_number"]),
            "minimal_jacobian_a1_a2_a3_a4_a6": row["minimal_jacobian_a1_a2_a3_a4_a6"],
            "generator_search_algorithm": search_algorithm,
            "primary_search_error": primary_error,
            "pari_rank_bounds_if_used": pari_rank_bounds,
            "generator_search_count": len(points),
            "generators": [[rational_text(x), rational_text(y)] for x, y in points],
            "independence_certificate": certificate,
            "certified_rank_lower_bound": len(points) if certified else None,
        }
    except Exception as error:
        result = {
            "status": "SEARCH_ERROR",
            "pair_key": key,
            "arithmetic_complexity_rank": int(row["arithmetic_complexity_rank"]),
            "global_root_number": int(row["global_root_number"]),
            "minimal_jacobian_a1_a2_a3_a4_a6": row["minimal_jacobian_a1_a2_a3_a4_a6"],
            "error_type": type(error).__name__,
            "error": str(error),
            "certified_rank_lower_bound": None,
        }
    ledger["results"][key] = result
    completed_this_run += 1
    if completed_this_run % args.checkpoint_every == 0:
        atomic_write(args.output, ledger)
        print(
            "ELKIES2026PAIRRANKS|checkpoint=true|"
            f"completed_this_run={completed_this_run}|last_complexity_rank={row['arithmetic_complexity_rank']}|"
            f"last_status={result['status']}|last_lower_bound={result['certified_rank_lower_bound']}",
            flush=True,
        )

result_values = list(ledger["results"].values())
status_counts = Counter(result["status"] for result in result_values)
rank_counts = Counter(
    int(result["certified_rank_lower_bound"])
    for result in result_values
    if result["certified_rank_lower_bound"] is not None
)
max_rank = max(rank_counts, default=None)
leaders = sorted(
    (
        {
            "pair_key": result["pair_key"],
            "arithmetic_complexity_rank": result["arithmetic_complexity_rank"],
            "certified_rank_lower_bound": result["certified_rank_lower_bound"],
        }
        for result in result_values
        if result["certified_rank_lower_bound"] == max_rank
    ),
    key=lambda item: (item["arithmetic_complexity_rank"], item["pair_key"]),
)
ledger["summary"] = {
    "catalogue_pair_count": len(rows),
    "ledger_pair_count": len(result_values),
    "complete": len(result_values) == len(rows),
    "status_counts": dict(sorted(status_counts.items())),
    "certified_rank_lower_bound_counts": {
        str(rank): count for rank, count in sorted(rank_counts.items())
    },
    "maximum_certified_rank_lower_bound": max_rank,
    "leaders": leaders,
}
ledger["status"] = (
    "PASS_COMPLETE_EXACT_RANK_LOWER_BOUND_LEDGER"
    if ledger["summary"]["complete"] and status_counts.get("SEARCH_ERROR", 0) == 0
    else "PARTIAL_EXACT_RANK_LOWER_BOUND_LEDGER"
)
atomic_write(args.output, ledger)
print(
    "ELKIES2026PAIRRANKS|"
    f"interval={args.start}:{stop}|completed_this_run={completed_this_run}|"
    f"ledger_pairs={len(result_values)}|max_certified_lower_bound={max_rank}|"
    f"status={ledger['status']}|output={display_path(args.output)}"
)
