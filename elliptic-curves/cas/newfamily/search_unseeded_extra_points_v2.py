#!/usr/bin/env sage -python
"""Unseeded extra-point discovery with robust height-span triage.

This is the successor to search_unseeded_extra_points.py.  It fixes two
problems exposed by the T=11 calibration:

* eclib discovery uses pp=0, because we want raw points rather than subgroup
  processing/saturation;
* augmented-matrix numerical rank is not used to classify a new direction.
  Instead we compute the Schur-complement residual of Q against the known
  rank-11 height lattice at two precisions.

For G = (<Ui,Uj>) and v = (<Ui,Q>), the residual

    r(Q) = <Q,Q> - v^T G^{-1} v

is the squared canonical-height component orthogonal to span(U0,...,U10).
It is zero for a dependent point in exact arithmetic.  A TRIAGE hit requires
a positive residual that is large relative to height scale and stable when
precision is increased.  Any hit still requires exact eclib verification.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import QQ, RR, RealField, ZZ, matrix, vector
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from search_unseeded_extra_points import (
    build_integral_specialization,
    mw_point_to_sage,
)
from screen_seeded_rational_candidates_fast import qbits


def schur_residual(E, known, Q, precision):
    RF = RealField(precision)
    M = matrix(RF, E.height_pairing_matrix(known + [Q], precision=precision))
    G = M[:11, :11]
    v = vector(RF, [M[i, 11] for i in range(11)])
    hq = RF(M[11, 11])
    coeff = G.solve_right(v)
    projection = v * coeff
    residual = hq - projection
    scale = max(RF(1), abs(hq), abs(projection))
    relative = residual / scale
    return {
        "precision": int(precision),
        "height": hq,
        "projection": projection,
        "residual": residual,
        "relative": relative,
    }


def classify_residual(E, known, Q, p1, p2, relative_threshold, stability_factor):
    lo = schur_residual(E, known, Q, p1)
    hi = schur_residual(E, known, Q, p2)

    RF = RealField(p2)
    threshold = RF(relative_threshold)
    rlo = RF(lo["relative"])
    rhi = RF(hi["relative"])

    # A true orthogonal component must be positive and visible at both
    # precisions.  We also require the two estimates to agree multiplicatively
    # to avoid promoting numerical cancellation as a new rank direction.
    positive_visible = rlo > threshold and rhi > threshold
    if positive_visible:
        ratio = max(rlo / rhi, rhi / rlo)
        stable = ratio <= RF(stability_factor)
    else:
        ratio = RF(0)
        stable = False

    return {
        "low_precision": p1,
        "high_precision": p2,
        "low_residual": str(lo["residual"]),
        "high_residual": str(hi["residual"]),
        "low_relative_residual": str(lo["relative"]),
        "high_relative_residual": str(hi["relative"]),
        "stability_ratio": str(ratio),
        "relative_threshold": str(relative_threshold),
        "stable_positive_residual": bool(positive_visible and stable),
    }


def run_single(args):
    started = time.monotonic()
    print("PHASE build_start", flush=True)
    E, known, scale = build_integral_specialization(
        args.numerator, args.denominator, args.sections_sobj
    )
    print("PHASE build_done", flush=True)

    # Search on a global minimal model only.  Parent timeout intentionally
    # treats pathological minimalization as a cheap-search skip.
    print("PHASE minimal_start", flush=True)
    ms = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - ms
    iso_to_min = E.isomorphism_to(Emin)
    iso_from_min = ~iso_to_min
    print(f"PHASE minimal_done seconds={minimal_seconds:.6f}", flush=True)

    print("PHASE eclib_init_start", flush=True)
    es = time.monotonic()
    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    # pp=0: retain raw points found; do not spend time processing a subgroup.
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=0, maxr=args.maxr)
    init_seconds = time.monotonic() - es
    print(f"PHASE eclib_init_done seconds={init_seconds:.6f}", flush=True)

    print(f"PHASE search_start height={args.height}", flush=True)
    ss = time.monotonic()
    mw.search(args.height, verbose=args.verbose_search)
    search_seconds = time.monotonic() - ss
    raw = list(mw.points())
    print(f"PHASE search_done seconds={search_seconds:.6f} raw={len(raw)}", flush=True)

    known_signless = set(known) | set(-P for P in known)
    candidates = []
    seen = set()
    for triple in raw:
        Qmin = mw_point_to_sage(Emin, triple)
        if Qmin.is_zero():
            continue
        Q = iso_from_min(Qmin)
        if Q in known_signless or -Q in known_signless:
            continue
        key = min(str(Q), str(-Q))
        if key in seen:
            continue
        seen.add(key)

        test = classify_residual(
            E, known, Q,
            args.precision,
            args.verify_precision,
            args.relative_threshold,
            args.stability_factor,
        )
        candidates.append({
            "point": [str(Q[0]), str(Q[1])],
            "point_minimal": [str(Qmin[0]), str(Qmin[1])],
            "point_bits": max(qbits(Q[0]), qbits(Q[1])),
            "point_minimal_bits": max(qbits(Qmin[0]), qbits(Qmin[1])),
            **test,
            "numerical_new_direction": test["stable_positive_residual"],
        })

    hits = [c for c in candidates if c["numerical_new_direction"]]
    absrels = []
    for c in candidates:
        try:
            absrels.append(abs(RR(c["high_relative_residual"])))
        except Exception:
            pass

    result = {
        "status": "completed",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "height_limit": args.height,
        "extra_scale": str(scale),
        "minimal_seconds": minimal_seconds,
        "eclib_init_seconds": init_seconds,
        "search_seconds": search_seconds,
        "mwrank_points_found": len(raw),
        "nonknown_candidates": len(candidates),
        "numerical_new_direction_hits": len(hits),
        "max_abs_high_relative_residual": str(max(absrels)) if absrels else None,
        "candidate_points": candidates,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0


def priority(r):
    return (
        r.get("minimum", min(r.get("discovery", 0), r.get("held", 0))),
        r.get("held", 0), r.get("discovery", 0),
        -r.get("generator_bits_max", 10**9),
    )


def run_parent(args):
    source = json.loads(Path(args.screen_json).read_text())
    eligible = [r for r in source if r.get("status") == "completed"
                and r.get("height_rank") == 11
                and r.get("positive_definite") is True]
    eligible.sort(key=priority, reverse=True)
    if args.limit is not None:
        eligible = eligible[:args.limit]

    print(
        f"eligible={len(eligible)} height={args.height} timeout={args.timeout} "
        f"precision={args.precision}/{args.verify_precision} "
        f"rel_threshold={args.relative_threshold}", flush=True
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []

    for pos, r in enumerate(eligible, 1):
        a, b = int(r["numerator"]), int(r["denominator"])
        print(f"[{pos}/{len(eligible)}] T={a}/{b} D={r.get('discovery',0):.4f} H={r.get('held',0):.4f}", flush=True)
        cmd = [
            sys.executable, str(Path(__file__).resolve()), "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(a), "--denominator", str(b),
            "--discovery", repr(r.get("discovery", 0.0)),
            "--held", repr(r.get("held", 0.0)),
            "--height", str(args.height),
            "--precision", str(args.precision),
            "--verify-precision", str(args.verify_precision),
            "--relative-threshold", str(args.relative_threshold),
            "--stability-factor", str(args.stability_factor),
            "--maxr", str(args.maxr),
        ]
        if args.verbose_search:
            cmd.append("--verbose-search")
        try:
            completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT,
                                       timeout=args.timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            partial = exc.stdout or ""
            if isinstance(partial, bytes):
                partial = partial.decode(errors="replace")
            rec = {"status":"timeout", "parameter":f"{a}/{b}",
                   "numerator":a, "denominator":b,
                   "timeout_seconds":args.timeout,
                   "output_tail":"\n".join(partial.splitlines()[-20:])}
            records.append(rec)
            output.write_text(json.dumps(records, indent=2, sort_keys=True)+"\n")
            print("  TIMEOUT", flush=True)
            continue

        lines = [x for x in completed.stdout.splitlines() if x.strip()]
        try:
            rec = json.loads(lines[-1])
        except Exception:
            rec = {"status":"error", "parameter":f"{a}/{b}",
                   "returncode":completed.returncode,
                   "output_tail":"\n".join(lines[-30:])}
        records.append(rec)
        output.write_text(json.dumps(records, indent=2, sort_keys=True)+"\n")
        print(
            "  status=%s found=%s candidates=%s newdir=%s maxrel=%s min_s=%s init_s=%s search_s=%s wall=%s" % (
                rec.get("status"), rec.get("mwrank_points_found"),
                rec.get("nonknown_candidates"), rec.get("numerical_new_direction_hits"),
                rec.get("max_abs_high_relative_residual"), rec.get("minimal_seconds"),
                rec.get("eclib_init_seconds"), rec.get("search_seconds"),
                rec.get("wall_seconds")), flush=True)
        if rec.get("numerical_new_direction_hits", 0):
            print("  *** STABLE SCHUR NEW-DIRECTION HIT -- EXACT VERIFY NEXT ***", flush=True)

    done = [r for r in records if r.get("status") == "completed"]
    hits = [r for r in done if r.get("numerical_new_direction_hits",0)]
    print(json.dumps({
        "output":str(output), "attempted":len(records), "completed":len(done),
        "timeouts":sum(r.get("status")=="timeout" for r in records),
        "stable_schur_hit_specializations":len(hits),
    }, sort_keys=True), flush=True)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--screen-json", default="artifacts/local/elliptic-curves/newfamily/height_rank_screen16.json")
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--height", type=float, default=14.0)
    p.add_argument("--precision", type=int, default=180)
    p.add_argument("--verify-precision", type=int, default=260)
    p.add_argument("--relative-threshold", default="1e-8")
    p.add_argument("--stability-factor", type=float, default=2.0)
    p.add_argument("--maxr", type=int, default=64)
    p.add_argument("--limit", type=int)
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--verbose-search", action="store_true")
    p.add_argument("--output", default="artifacts/local/elliptic-curves/newfamily/unseeded_extra_points_v2.json")
    p.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--numerator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--denominator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--discovery", type=float, default=0.0, help=argparse.SUPPRESS)
    p.add_argument("--held", type=float, default=0.0, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.single:
        if args.numerator is None or args.denominator is None:
            p.error("single mode requires numerator and denominator")
        return run_single(args)
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
