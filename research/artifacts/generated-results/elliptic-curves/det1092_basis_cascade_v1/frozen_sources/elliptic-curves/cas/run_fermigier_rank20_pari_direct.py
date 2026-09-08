#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

CAS_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS_ROOT))

from run_fermigier_rank20_pari_descent import (
    KNOWN_RANK,
    _pari_int,
    _pari_points,
    classify_bounds,
    load_descent_basis,
    sage_q,
)

PROTOCOL = "R20PARIDIRECT"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=Path(
        "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"))
    ap.add_argument("--candidate-record", type=Path, default=Path(
        "artifacts/generated-results/elliptic-curves/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"))
    ap.add_argument("--effort", type=int, default=0)
    ap.add_argument("--pari-debug", type=int, default=0)
    args = ap.parse_args()

    from sage.all import EllipticCurve, QQ, pari, version as sage_version

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != KNOWN_RANK or not basis.mod2_certified:
        raise ArithmeticError("missing exact rank-20 mod-2 certificate")

    E = EllipticCurve(QQ, [sage_q(v) for v in basis.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis.points]

    if args.pari_debug:
        try:
            pari.default("debug", args.pari_debug)
        except Exception:
            pass

    print(f"{PROTOCOL}|stage=input|sage={sage_version()}|known={len(known)}|effort={args.effort}|pari_debug={args.pari_debug}", flush=True)
    PE = E.pari_curve()
    Pknown = _pari_points(known)

    print(f"{PROTOCOL}|stage=ellrank_direct|status=start", flush=True)
    t0 = time.monotonic()
    try:
        try:
            result = PE.ellrank(args.effort, Pknown)
        except (AttributeError, TypeError):
            result = pari.ellrank(PE, args.effort, Pknown)
    except Exception as exc:
        print(f"{PROTOCOL}|stage=ellrank_direct|status=error|seconds={time.monotonic()-t0:.6f}|error={exc}", flush=True)
        raise

    elapsed = time.monotonic() - t0
    lo = _pari_int(result[0])
    hi = _pari_int(result[1])
    sha = _pari_int(result[2])
    pts = result[3]
    eff_lo, cls = classify_bounds(lo, hi)

    print(f"{PROTOCOL}|stage=ellrank_direct|status=complete|seconds={elapsed:.6f}|pari_lower={lo}|pari_upper={hi}|effective_lower={eff_lo}|sha_pairing_rank={sha}|returned_points={len(pts)}|classification={cls}", flush=True)
    for i, p in enumerate(pts):
        print(f"{PROTOCOL}|stage=returned_point|i={i}|point={p}", flush=True)

if __name__ == "__main__":
    main()
