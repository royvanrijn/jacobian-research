#!/usr/bin/env python3
"""Run Denis Simon's older GP 2-descent on the Fermigier rank-20 anchor.

This intentionally bypasses PARI ellrank/ellrankinit.  For curves over QQ,
Sage's wrapper calls Simon's ellQ_ellrank directly and does not construct a
BNF for a cubic number field.

Run from jacobian-research root:

  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_simon_descent.py \
    --verbose 2 \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_simon_descent.log
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

CAS_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS_ROOT))

from run_fermigier_rank20_pari_descent import (
    KNOWN_RANK,
    load_descent_basis,
    sage_q,
)

PROTOCOL = "R20SIMON"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/"
            "fermigier_rank20_near_miss_v1.json"
        ),
    )
    ap.add_argument(
        "--candidate-record",
        type=Path,
        default=Path(
            "artifacts/generated-results/"
            "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
        ),
    )
    ap.add_argument("--verbose", type=int, default=2, choices=(0, 1, 2, 3))
    ap.add_argument("--lim1", type=int, default=5)
    ap.add_argument("--lim3", type=int, default=50)
    ap.add_argument("--limtriv", type=int, default=3)
    ap.add_argument("--maxprob", type=int, default=20)
    ap.add_argument("--limbigprime", type=int, default=30)
    args = ap.parse_args()

    from sage.all import EllipticCurve, QQ, version as sage_version
    from sage.schemes.elliptic_curves.gp_simon import simon_two_descent

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != KNOWN_RANK or not basis.mod2_certified:
        raise ArithmeticError("missing exact rank-20 mod-2 certificate")

    E = EllipticCurve(QQ, [sage_q(v) for v in basis.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis.points]

    print(
        f"{PROTOCOL}|stage=input|sage={sage_version()}|known={len(known)}"
        f"|basis_sha256={basis.basis_sha256}"
        f"|verbose={args.verbose}|lim1={args.lim1}|lim3={args.lim3}"
        f"|limtriv={args.limtriv}|maxprob={args.maxprob}"
        f"|limbigprime={args.limbigprime}",
        flush=True,
    )

    t0 = time.monotonic()
    print(f"{PROTOCOL}|stage=two_descent|status=start", flush=True)
    try:
        lower, upper, points = simon_two_descent(
            E,
            verbose=args.verbose,
            lim1=args.lim1,
            lim3=args.lim3,
            limtriv=args.limtriv,
            maxprob=args.maxprob,
            limbigprime=args.limbigprime,
            known_points=known,
        )
    except Exception as exc:
        print(
            f"{PROTOCOL}|stage=two_descent|status=error"
            f"|seconds={time.monotonic()-t0:.6f}"
            f"|error_type={type(exc).__name__}|error={exc}",
            flush=True,
        )
        raise

    elapsed = time.monotonic() - t0
    effective_lower = max(KNOWN_RANK, int(lower))
    upper = int(upper)

    if upper < KNOWN_RANK:
        classification = "CONTRADICTION"
    elif upper == KNOWN_RANK:
        classification = "S0_EXACT_RANK20"
    elif effective_lower >= KNOWN_RANK + 1:
        classification = "S3_RANK_AT_LEAST21"
    else:
        classification = "S2_RESIDUAL_SELMER_ROOM"

    print(
        f"{PROTOCOL}|stage=two_descent|status=complete"
        f"|seconds={elapsed:.6f}|simon_lower={int(lower)}"
        f"|simon_upper={upper}|effective_lower={effective_lower}"
        f"|returned_points={len(points)}|classification={classification}",
        flush=True,
    )

    for i, P in enumerate(points):
        print(f"{PROTOCOL}|stage=returned_point|i={i}|point={P}", flush=True)


if __name__ == "__main__":
    main()
