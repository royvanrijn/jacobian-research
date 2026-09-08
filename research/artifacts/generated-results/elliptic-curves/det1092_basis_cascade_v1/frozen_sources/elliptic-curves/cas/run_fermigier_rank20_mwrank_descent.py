#!/usr/bin/env python3
"""Run an eclib/mwrank 2-descent for the Fermigier rank-20 curve.

Run with SageMath, not plain CPython::

    sage -python elliptic-curves/cas/run_fermigier_rank20_mwrank_descent.py

This is a second, independent open-source backend for the relative-rank
experiment.  mwrank does not accept the repository's twenty known generators
as input to its Selmer-only descent, so the program combines two separate exact
facts:

* the repository certificate proves ``rank(E) >= 20``;
* mwrank supplies a rigorous 2-descent upper bound.

The default is Selmer-only mode: do not spend time searching for points that the
repository already knows.  The raw mwrank lower bound is retained for provenance
but is not expected to reproduce rank 20 in this mode.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys
import time


CAS_ROOT = Path(__file__).resolve().parent
PROGRAM_ROOT = CAS_ROOT.parent
sys.path.insert(0, str(CAS_ROOT))
sys.path.insert(0, str(PROGRAM_ROOT))

from run_fermigier_rank20_pari_descent import (  # noqa: E402
    KNOWN_RANK,
    load_descent_basis,
    sage_q,
)


PROTOCOL = "R20MWRANK"
PROTOCOL_VERSION = 1


@dataclass(frozen=True)
class MwrankResult:
    selmer_rank: int
    mwrank_lower: int
    mwrank_upper: int
    effective_lower: int
    certain_internal: bool
    two_torsion_rank: int
    elapsed_seconds: float
    classification: str


def classify_bounds(
    mwrank_lower: int,
    mwrank_upper: int,
    known_rank: int = KNOWN_RANK,
) -> tuple[int, str]:
    """Combine mwrank's interval with the independent certified lower bound."""

    if mwrank_lower < 0 or mwrank_upper < mwrank_lower:
        raise ValueError("invalid mwrank rank interval")
    if mwrank_upper < known_rank:
        raise ArithmeticError(
            f"mwrank upper bound {mwrank_upper} contradicts certified rank >= {known_rank}"
        )
    effective_lower = max(known_rank, mwrank_lower)
    if effective_lower > mwrank_upper:
        raise ArithmeticError("combined mwrank/repository interval is empty")
    if effective_lower == mwrank_upper == known_rank:
        return effective_lower, "M0_exact_rank20"
    if effective_lower >= known_rank + 1:
        return effective_lower, "M3_rank_at_least21"
    return effective_lower, "M2_residual_rank_interval"


def collect_result(mwrank_curve, two_torsion_rank: int, elapsed_seconds: float) -> MwrankResult:
    """Read completed eclib descent state and classify it."""

    selmer_rank = int(mwrank_curve.selmer_rank())
    mwrank_lower = int(mwrank_curve.rank())
    mwrank_upper = int(mwrank_curve.rank_bound())
    certain_internal = bool(mwrank_curve.certain())
    effective_lower, classification = classify_bounds(mwrank_lower, mwrank_upper)
    return MwrankResult(
        selmer_rank=selmer_rank,
        mwrank_lower=mwrank_lower,
        mwrank_upper=mwrank_upper,
        effective_lower=effective_lower,
        certain_internal=certain_internal,
        two_torsion_rank=int(two_torsion_rank),
        elapsed_seconds=float(elapsed_seconds),
        classification=classification,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/"
            "fermigier_rank20_near_miss_v1.json"
        ),
    )
    parser.add_argument(
        "--candidate-record",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/"
            "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier_rank20_mwrank_descent.json"
        ),
    )
    parser.add_argument(
        "--first-limit",
        type=int,
        default=20,
        help="mwrank first quartic-search bound (normally irrelevant in Selmer-only mode)",
    )
    parser.add_argument(
        "--second-limit",
        type=int,
        default=8,
        help="mwrank logarithmic second quartic-search bound",
    )
    parser.add_argument(
        "--n-aux",
        type=int,
        default=-1,
        help="number of auxiliary sieve primes; -1 uses eclib's default",
    )
    parser.add_argument(
        "--search-points",
        action="store_true",
        help="disable Selmer-only mode and let mwrank search homogeneous spaces for points",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress eclib's detailed 2-descent progress output",
    )
    args = parser.parse_args()

    if args.first_limit < 0 or args.second_limit < 0:
        raise ValueError("search limits must be nonnegative")
    if args.n_aux == 0 or args.n_aux < -1:
        raise ValueError("--n-aux must be -1 or a positive integer")

    try:
        from sage.all import EllipticCurve, QQ, version as sage_version
    except ImportError as exc:
        raise SystemExit(
            "This program must run inside SageMath: sage -python " + str(Path(__file__))
        ) from exc

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != KNOWN_RANK or not basis.mod2_certified:
        raise ArithmeticError("refusing descent without the exact full mod-2 basis certificate")

    E = EllipticCurve(QQ, [sage_q(value) for value in basis.model])
    known_points = [E(sage_q(x), sage_q(y)) for x, y in basis.points]
    if len(known_points) != KNOWN_RANK or len(set(known_points)) != KNOWN_RANK:
        raise ArithmeticError("the pinned descent basis did not replay as 20 distinct Sage points")

    two_torsion_rank = int(E.two_torsion_rank())
    selmer_only = not args.search_points
    verbose = not args.quiet

    print(
        f"{PROTOCOL}|version={PROTOCOL_VERSION}|stage=input|sage={sage_version()}"
        f"|known_external={KNOWN_RANK}|basis={basis.basis_id}"
        f"|basis_sha256={basis.basis_sha256}|two_torsion_rank={two_torsion_rank}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=two_descent|status=start|selmer_only={str(selmer_only).lower()}"
        f"|first_limit={args.first_limit}|second_limit={args.second_limit}"
        f"|n_aux={args.n_aux}|verbose={str(verbose).lower()}",
        flush=True,
    )

    mwrank_curve = E.mwrank_curve(verbose=verbose)
    started = time.monotonic()
    try:
        mwrank_curve.two_descent(
            verbose=verbose,
            selmer_only=selmer_only,
            first_limit=args.first_limit,
            second_limit=args.second_limit,
            n_aux=args.n_aux,
            second_descent=True,
        )
    except Exception as exc:
        elapsed = time.monotonic() - started
        print(
            f"{PROTOCOL}|stage=two_descent|status=error|seconds={elapsed:.6f}"
            f"|error_type={type(exc).__name__}|error={exc}",
            flush=True,
        )
        raise

    elapsed = time.monotonic() - started
    result = collect_result(mwrank_curve, two_torsion_rank, elapsed)
    print(
        f"{PROTOCOL}|stage=two_descent|status=complete|seconds={elapsed:.6f}"
        f"|selmer_rank={result.selmer_rank}|mwrank_lower={result.mwrank_lower}"
        f"|mwrank_upper={result.mwrank_upper}"
        f"|effective_lower={result.effective_lower}"
        f"|certain_internal={str(result.certain_internal).lower()}"
        f"|classification={result.classification}",
        flush=True,
    )

    payload = {
        "schema": "elliptic-curves.fermigier-rank20-mwrank-descent.v1",
        "engine": "SageMath/eclib mwrank two_descent",
        "mathematical_status": "external_cas_rank_upper_bound",
        "basis": {
            "id": basis.basis_id,
            "count": KNOWN_RANK,
            "sha256": basis.basis_sha256,
            "mod2_rank": basis.mod2_rank,
            "mod2_certified": basis.mod2_certified,
            "candidate_record_sha256": basis.candidate_record_sha256,
            "manifest_sha256": basis.manifest_sha256,
        },
        "options": {
            "selmer_only": selmer_only,
            "first_limit": args.first_limit,
            "second_limit": args.second_limit,
            "n_aux": args.n_aux,
            "verbose": verbose,
        },
        "result": asdict(result),
        "rank_interval": [result.effective_lower, result.mwrank_upper],
        "classification": result.classification,
        "notes": [
            "The lower endpoint comes from the repository's exact rank>=20 certificate.",
            "Selmer-only mwrank does not consume the twenty known points as generators.",
            "certain_internal describes only mwrank's own lower/upper bounds, not the combined external certificate.",
            "For curves without rational 2-torsion, mwrank's rank upper bound equals its 2-Selmer rank.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=summary|rank_lower={result.effective_lower}"
        f"|rank_upper={result.mwrank_upper}|classification={result.classification}"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
