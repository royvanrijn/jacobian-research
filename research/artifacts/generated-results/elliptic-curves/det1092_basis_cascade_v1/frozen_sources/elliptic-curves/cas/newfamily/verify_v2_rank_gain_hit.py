#!/usr/bin/env sage -python
"""Exact verification of a v2 Schur-complement rank-gain hit.

Consumes search_unseeded_extra_points_v2.py JSON, selects a specialization and
one flagged candidate point, rebuilds the specialization, maps the known 11
hidden sections to the same global minimal model, and processes all 12 exact
rational points through eclib/mwrank with saturation_bound=0.

If the processed subgroup rank reaches 12, this proves the specialization has
Mordell--Weil rank at least 12.  It does not prove the full rank is exactly 12.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from sage.all import QQ, ZZ
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from search_unseeded_extra_points import build_integral_specialization
from screen_seeded_rational_candidates_fast import qbits, to_mwrank_triple


def parse_parameter(text: str):
    if "/" in text:
        a, b = text.split("/", 1)
        return int(a), int(b)
    return int(text), 1


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-json", required=True)
    p.add_argument("--parameter", required=True, help="e.g. 58/1")
    p.add_argument("--hit-index", type=int, default=0,
                   help="index among numerical_new_direction candidates")
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--output", default="artifacts/local/elliptic-curves/newfamily/exact_rank_gain_verify.json")
    args = p.parse_args()

    rows = json.loads(Path(args.input_json).read_text())
    matches = [r for r in rows if r.get("parameter") == args.parameter]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one record for {args.parameter}, got {len(matches)}")
    record = matches[0]
    hits = [c for c in record.get("candidate_points", []) if c.get("numerical_new_direction")]
    if not hits:
        raise SystemExit(f"no numerical_new_direction hit for {args.parameter}")
    if not (0 <= args.hit_index < len(hits)):
        raise SystemExit(f"hit-index {args.hit_index} out of range 0..{len(hits)-1}")
    hit = hits[args.hit_index]

    a, b = parse_parameter(args.parameter)
    started = time.monotonic()
    print(f"STAGE build T={a}/{b}", flush=True)
    E, known, _ = build_integral_specialization(a, b, args.sections_sobj)

    ms = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - ms
    iso = E.isomorphism_to(Emin)
    known_min = [iso(P) for P in known]

    xq = QQ(hit["point_minimal"][0])
    yq = QQ(hit["point_minimal"][1])
    Qmin = Emin([xq, yq])
    if Qmin.is_zero():
        raise RuntimeError("flagged candidate is zero")

    labels = [f"U{i}" for i in range(11)] + ["Q"]
    points = known_min + [Qmin]
    bits = [max(qbits(P[0]), qbits(P[1])) for P in points]
    order = sorted(range(12), key=lambda i: bits[i])

    print(
        f"STAGE minimal seconds={minimal_seconds:.6f} "
        f"disc_bits={ZZ(abs(Emin.discriminant())).nbits()} "
        f"Q_bits={bits[11]} schur_rel={hit.get('high_relative_residual')}",
        flush=True,
    )

    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=32)

    steps = []
    for position, index in enumerate(order, 1):
        before = len(mw.points())
        ss = time.monotonic()
        print(
            f"STAGE process position={position}/12 label={labels[index]} "
            f"bits={bits[index]} rank_before={before}",
            flush=True,
        )
        mw.process([to_mwrank_triple(points[index])], saturation_bound=0)
        after = len(mw.points())
        seconds = time.monotonic() - ss
        print(
            f"STAGE processed label={labels[index]} rank={before}->{after} "
            f"seconds={seconds:.6f}",
            flush=True,
        )
        steps.append({
            "label": labels[index],
            "bits": bits[index],
            "rank_before": before,
            "rank_after": after,
            "seconds": seconds,
        })

    final_rank = len(mw.points())
    result = {
        "status": "completed",
        "parameter": args.parameter,
        "input_json": args.input_json,
        "hit_index": args.hit_index,
        "schur_high_relative_residual": hit.get("high_relative_residual"),
        "schur_low_relative_residual": hit.get("low_relative_residual"),
        "candidate_point_minimal": hit.get("point_minimal"),
        "candidate_point_fixed": hit.get("point"),
        "processed_subgroup_rank": final_rank,
        "proved_rank_at_least_12": bool(final_rank >= 12),
        "minimal_seconds": minimal_seconds,
        "steps": steps,
        "wall_seconds": time.monotonic() - started,
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True), flush=True)
    if final_rank >= 12:
        print("*** EXACT SUBGROUP RANK >= 12 VERIFIED ***", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
