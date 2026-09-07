#!/usr/bin/env sage -python
"""Cheap extra-point discovery for rational newfamily specializations.

This stage deliberately does NOT seed eclib with the known rank-11 subgroup.
For each healthy specialization from the height-rank screen it:

1. builds the homogeneous integral short model and the 11 hidden sections;
2. checks the known section lattice numerically by canonical heights;
3. converts only the SEARCH curve to a global minimal model;
4. runs an empty mwrank_MordellWeil.search(height) on that minimal model;
5. maps points found by eclib back to the fixed homogeneous model;
6. rejects exact duplicates of the known sections;
7. tests each remaining point by augmenting the 11x11 canonical-height Gram
   matrix to 12x12.

The hybrid model choice matters: height arithmetic on the fixed short model is
fast, while eclib initialization can be extremely slow on the large nonminimal
homogeneous integral model. Search therefore happens on a global minimal model
without ever processing the known rank-11 subgroup through eclib.

Discovery uses pp=0 intentionally: eclib retains raw points found by the sieve
instead of processing/saturating them into a Mordell-Weil subgroup.  Exact
processing is deferred until a numerical new-direction hit exists.

A numerical height-rank jump to 12 is only a TRIAGE HIT. It is not promoted as
an exact rank statement. Exact eclib verification is intentionally deferred
until such a hit exists.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, QQ, RR, RealField, ZZ, matrix
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from screen_seeded_rational_candidates_fast import (
    ROOTS,
    extra_integral_scale,
    load_builder,
    load_sections,
    qbits,
)


def numerical_height_rank(E, points, precision: int, rank_bits: int):
    G = E.height_pairing_matrix(points, precision=precision)
    RF = RealField(precision)
    GR = matrix(RF, G)
    eigenvalues = sorted(RF(x) for x in GR.eigenvalues())
    largest = max(abs(x) for x in eigenvalues) if eigenvalues else RF(0)
    tolerance = max(RF(2) ** (-rank_bits), largest * RF(2) ** (-rank_bits))
    rank = sum(abs(x) > tolerance for x in eigenvalues)
    positive = all(x > tolerance for x in eigenvalues)
    return int(rank), bool(positive), str(tolerance), str(eigenvalues[0]), str(eigenvalues[-1])


def mw_point_to_sage(E, triple):
    X, Y, Z = [ZZ(v) for v in triple]
    if Z == 0:
        return E(0)
    return E([QQ(X) / QQ(Z), QQ(Y) / QQ(Z)])


def build_integral_specialization(numerator: int, denominator: int, sections_sobj: str):
    a = ZZ(numerator)
    b = ZZ(denominator)
    t = QQ(a) / QQ(b)
    sections = load_sections(sections_sobj)
    family = load_builder()(ROOTS)

    A = QQ(family["Amin"](t))
    B = QQ(family["Bmin"](t))
    Ah = A * b**8
    Bh = B * b**12
    c = extra_integral_scale(Ah, Bh)
    Aint_q = Ah * c**4
    Bint_q = Bh * c**6
    if Aint_q.denominator() != 1 or Bint_q.denominator() != 1:
        raise RuntimeError("failed to integralize homogeneous short model")
    E = EllipticCurve(QQ, [0, 0, 0, ZZ(Aint_q), ZZ(Bint_q)])
    if E.discriminant() == 0:
        raise RuntimeError("singular specialization")

    xscale = b**4 * c**2
    yscale = b**6 * c**3
    known = []
    for xf, yf in sections:
        known.append(E([QQ(xf(t)) * xscale, QQ(yf(t)) * yscale]))
    if len(set(known)) != 11:
        raise RuntimeError("hidden sections collide at specialization")
    return E, known, c


def run_single(args):
    started = time.monotonic()
    print("PHASE build_start", flush=True)
    E, known, scale = build_integral_specialization(
        args.numerator, args.denominator, args.sections_sobj
    )
    print("PHASE build_done", flush=True)

    height_started = time.monotonic()
    baseline_rank, baseline_pd, tolerance, _, _ = numerical_height_rank(
        E, known, args.precision, args.rank_bits
    )
    baseline_height_seconds = time.monotonic() - height_started
    print(
        f"PHASE baseline_height_done seconds={baseline_height_seconds:.6f} "
        f"rank={baseline_rank} pd={baseline_pd}",
        flush=True,
    )
    if baseline_rank != 11 or not baseline_pd:
        raise RuntimeError(
            f"unhealthy baseline height lattice rank={baseline_rank} pd={baseline_pd}"
        )

    print("PHASE minimal_start", flush=True)
    minimal_started = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - minimal_started
    iso_to_min = E.isomorphism_to(Emin)
    iso_from_min = ~iso_to_min
    print(
        f"PHASE minimal_done seconds={minimal_seconds:.6f} "
        f"disc_bits={ZZ(abs(Emin.discriminant())).nbits()}",
        flush=True,
    )

    print("PHASE eclib_init_start", flush=True)
    init_started = time.monotonic()
    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=0, maxr=args.maxr)
    init_seconds = time.monotonic() - init_started
    print(f"PHASE eclib_init_done seconds={init_seconds:.6f} pp=0", flush=True)

    print(
        f"PHASE search_start height={args.height} baseline_hrank=11 "
        f"known_max_bits={max(max(qbits(P[0]), qbits(P[1])) for P in known)}",
        flush=True,
    )
    search_started = time.monotonic()
    mw.search(args.height, verbose=args.verbose_search)
    search_seconds = time.monotonic() - search_started
    raw = list(mw.points())
    print(
        f"PHASE search_done seconds={search_seconds:.6f} mw_points={len(raw)}",
        flush=True,
    )

    known_signless = set()
    for P in known:
        known_signless.add(P)
        known_signless.add(-P)

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
        hrank, pd, tol12, emin, emax = numerical_height_rank(
            E, known + [Q], args.precision, args.rank_bits
        )
        candidates.append({
            "point": [str(Q[0]), str(Q[1])],
            "point_minimal": [str(Qmin[0]), str(Qmin[1])],
            "point_bits": max(qbits(Q[0]), qbits(Q[1])),
            "point_minimal_bits": max(qbits(Qmin[0]), qbits(Qmin[1])),
            "augmented_height_rank": hrank,
            "augmented_positive_definite": pd,
            "tolerance": tol12,
            "smallest_eigenvalue": emin,
            "largest_eigenvalue": emax,
            "numerical_new_direction": bool(hrank == 12 and pd),
        })

    hits = [c for c in candidates if c["numerical_new_direction"]]
    result = {
        "status": "completed",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "height_limit": args.height,
        "baseline_height_rank": baseline_rank,
        "baseline_positive_definite": baseline_pd,
        "baseline_tolerance": tolerance,
        "extra_scale": str(scale),
        "minimal_seconds": minimal_seconds,
        "eclib_init_seconds": init_seconds,
        "mwrank_points_found": len(raw),
        "nonknown_candidates": len(candidates),
        "numerical_new_direction_hits": len(hits),
        "candidate_points": candidates,
        "baseline_height_seconds": baseline_height_seconds,
        "search_seconds": search_seconds,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0


def priority(record):
    return (
        record.get("minimum", min(record.get("discovery", 0), record.get("held", 0))),
        record.get("held", 0),
        record.get("discovery", 0),
        -record.get("generator_bits_max", 10**9),
    )


def run_parent(args):
    source = json.loads(Path(args.screen_json).read_text())
    eligible = [
        r for r in source
        if r.get("status") == "completed"
        and r.get("height_rank") == 11
        and r.get("positive_definite") is True
    ]
    eligible.sort(key=priority, reverse=True)
    if args.limit is not None:
        eligible = eligible[:args.limit]

    print(
        f"eligible={len(eligible)} height={args.height} timeout={args.timeout} "
        f"precision={args.precision}", flush=True
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []

    for pos, r in enumerate(eligible, 1):
        a = int(r["numerator"])
        b = int(r["denominator"])
        print(
            f"[{pos}/{len(eligible)}] T={a}/{b} "
            f"D={r.get('discovery',0):.4f} H={r.get('held',0):.4f} "
            f"bits={r.get('generator_bits_max')}", flush=True
        )
        cmd = [
            sys.executable, str(Path(__file__).resolve()), "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(a), "--denominator", str(b),
            "--discovery", repr(r.get("discovery", 0.0)),
            "--held", repr(r.get("held", 0.0)),
            "--height", str(args.height),
            "--precision", str(args.precision),
            "--rank-bits", str(args.rank_bits),
            "--maxr", str(args.maxr),
        ]
        if args.verbose_search:
            cmd.append("--verbose-search")
        try:
            completed = subprocess.run(
                cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=args.timeout, check=False,
            )
        except subprocess.TimeoutExpired as exc:
            partial = exc.stdout or ""
            if isinstance(partial, bytes):
                partial = partial.decode(errors="replace")
            record = {
                "status": "timeout", "parameter": f"{a}/{b}",
                "numerator": a, "denominator": b,
                "timeout_seconds": args.timeout,
                "output_tail": "\n".join(partial.splitlines()[-30:]),
            }
            records.append(record)
            output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
            print("  TIMEOUT", flush=True)
            for line in partial.splitlines()[-8:]:
                print("   ", line, flush=True)
            continue

        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        try:
            record = json.loads(lines[-1])
        except Exception:
            record = {
                "status": "error", "parameter": f"{a}/{b}",
                "returncode": completed.returncode,
                "output_tail": "\n".join(lines[-30:]),
            }
        records.append(record)
        output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
        print(
            "  status=%s found=%s candidates=%s newdir=%s min_s=%s init_s=%s "
            "search_s=%s wall=%s" % (
                record.get("status"), record.get("mwrank_points_found"),
                record.get("nonknown_candidates"),
                record.get("numerical_new_direction_hits"),
                record.get("minimal_seconds"), record.get("eclib_init_seconds"),
                record.get("search_seconds"), record.get("wall_seconds"),
            ), flush=True
        )
        if record.get("numerical_new_direction_hits", 0):
            print("  *** NUMERICAL NEW-DIRECTION HIT -- EXACT VERIFY NEXT ***", flush=True)

    done = [r for r in records if r.get("status") == "completed"]
    hits = [r for r in done if r.get("numerical_new_direction_hits", 0) > 0]
    print(json.dumps({
        "output": str(output),
        "attempted": len(records),
        "completed": len(done),
        "timeouts": sum(r.get("status") == "timeout" for r in records),
        "numerical_new_direction_specializations": len(hits),
    }, sort_keys=True), flush=True)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--screen-json",
        default="artifacts/local/elliptic-curves/newfamily/height_rank_screen16.json",
    )
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--height", type=float, default=12.0)
    p.add_argument("--precision", type=int, default=180)
    p.add_argument("--rank-bits", type=int, default=80)
    p.add_argument("--maxr", type=int, default=32)
    p.add_argument("--limit", type=int)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--verbose-search", action="store_true")
    p.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/unseeded_extra_points.json",
    )
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
