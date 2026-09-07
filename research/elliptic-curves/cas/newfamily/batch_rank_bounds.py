#!/usr/bin/env sage -python
"""Batch PARI rank-bound classification for pinned newfamily specializations.

Loads the versioned rank-gain batch, reconstructs each selected specialization
from committed family/section data, replays the pinned exact subgroup with eclib,
and then calls PARI ellrank with that subgroup as known points.

This is intended as a breadth-first classifier before expensive H18/H20 point
searches.  Exact-rank claims are made only when the PARI interval closes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import time

from sage.all import QQ, ZZ, pari, version as sage_version
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from search_unseeded_extra_points_v3 import build_integral_specialization_git
from screen_seeded_rational_candidates_fast import qbits, to_mwrank_triple


@dataclass(frozen=True)
class PariPass:
    effort: int
    pari_lower: int
    pari_upper: int
    effective_lower: int
    sha_pairing_rank: int
    returned_points: int
    seconds: float


def parse_parameter(text: str):
    if "/" in text:
        a, b = text.split("/", 1)
        return int(a), int(b)
    return int(text), 1


def parse_efforts(text: str):
    vals = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if not vals or any(x < 0 for x in vals):
        raise ValueError("--efforts must be nonempty nonnegative integers")
    return tuple(dict.fromkeys(vals))


def _pari_int(value):
    try:
        return int(value)
    except TypeError:
        return int(value.sage())


def _ellrankinit(pari_curve):
    try:
        return pari_curve.ellrankinit()
    except (AttributeError, TypeError):
        return pari.ellrankinit(pari_curve)


def _ellrank(ctx, effort, points):
    pp = pari([[P[0], P[1]] for P in points])
    try:
        return ctx.ellrank(effort, pp)
    except (AttributeError, TypeError):
        return pari.ellrank(ctx, effort, pp)


def choose_rows(payload, args):
    rows = [r for r in payload.get("records", []) if r.get("processed_subgroup_rank") is not None]
    if args.parameter:
        wanted = set(args.parameter)
        rows = [r for r in rows if r.get("parameter") in wanted]
    else:
        rows = [r for r in rows if int(r.get("processed_subgroup_rank", -1)) >= args.min_rank]
    if args.exclude_parameter:
        excluded = set(args.exclude_parameter)
        rows = [r for r in rows if r.get("parameter") not in excluded]
    rows.sort(key=lambda r: (-int(r.get("processed_subgroup_rank", -1)), -int(r.get("unique_hit_count", 0)), r.get("parameter", "")))
    return rows


def replay_row(row):
    parameter = row["parameter"]
    a, b = parse_parameter(parameter)
    E, known_integral, _ = build_integral_specialization_git(a, b)
    Emin = E.global_minimal_model()
    iso = E.isomorphism_to(Emin)
    known = [iso(P) for P in known_integral]

    extras = []
    for hit in row.get("hit_rank_increases", []):
        x, y = hit["point_minimal"]
        extras.append(Emin([QQ(x), QQ(y)]))

    expected_baseline = int(row["known_subgroup_rank"])
    expected_lower = int(row["processed_subgroup_rank"])

    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=max(32, expected_lower + 8))
    known_order = sorted(range(11), key=lambda i: max(qbits(known[i][0]), qbits(known[i][1])))
    growth = []
    for i in known_order:
        before = len(mw.points())
        mw.process([to_mwrank_triple(known[i])], saturation_bound=0)
        after = len(mw.points())
        growth.append({"label": f"U{i}", "before": before, "after": after})
    baseline = len(mw.points())
    if baseline != expected_baseline:
        raise RuntimeError(f"{parameter}: baseline {baseline} != pinned {expected_baseline}")

    for j, P in enumerate(extras):
        before = len(mw.points())
        mw.process([to_mwrank_triple(P)], saturation_bound=0)
        after = len(mw.points())
        growth.append({"label": f"Q{j}", "before": before, "after": after})
        if after != before + 1:
            raise RuntimeError(f"{parameter}: pinned extra Q{j} did not increase rank")
    lower = len(mw.points())
    if lower != expected_lower:
        raise RuntimeError(f"{parameter}: replay lower {lower} != pinned {expected_lower}")
    return Emin, known + extras, growth


def classify(row, efforts):
    parameter = row["parameter"]
    pinned_lower = int(row["processed_subgroup_rank"])
    print(f"NFRANK|T={parameter}|stage=eclib|status=start|pinned_lower={pinned_lower}", flush=True)
    Emin, points, growth = replay_row(row)
    print(f"NFRANK|T={parameter}|stage=eclib|status=complete|lower={len(points)}", flush=True)

    started = time.monotonic()
    ctx = _ellrankinit(Emin.pari_curve())
    init_seconds = time.monotonic() - started
    print(f"NFRANK|T={parameter}|stage=ellrankinit|status=complete|seconds={init_seconds:.6f}", flush=True)

    passes = []
    for effort in efforts:
        print(f"NFRANK|T={parameter}|stage=ellrank|status=start|effort={effort}", flush=True)
        started = time.monotonic()
        result = _ellrank(ctx, effort, points)
        seconds = time.monotonic() - started
        if len(result) != 4:
            raise RuntimeError(f"{parameter}: unexpected ellrank result {result}")
        pl, pu, sha = _pari_int(result[0]), _pari_int(result[1]), _pari_int(result[2])
        returned = len(result[3])
        if pu < pinned_lower:
            raise RuntimeError(f"{parameter}: PARI upper {pu} contradicts pinned lower {pinned_lower}")
        effective = max(pinned_lower, pl)
        p = PariPass(effort, pl, pu, effective, sha, returned, seconds)
        passes.append(p)
        print(
            f"NFRANK|T={parameter}|stage=ellrank|status=complete|effort={effort}"
            f"|pari_lower={pl}|pari_upper={pu}|effective_lower={effective}"
            f"|sha_pairing_rank={sha}|returned_points={returned}|seconds={seconds:.6f}",
            flush=True,
        )
        if effective == pu:
            break

    lower = max(p.effective_lower for p in passes)
    upper = min(p.pari_upper for p in passes)
    if lower > upper:
        raise RuntimeError(f"{parameter}: inconsistent PARI intervals")
    classification = f"exact_rank_{lower}" if lower == upper else "residual_rank_interval"
    return {
        "parameter": parameter,
        "pinned_lower_bound": pinned_lower,
        "rank_interval": [lower, upper],
        "classification": classification,
        "root_number": int(Emin.root_number()),
        "minimal_model": [str(v) for v in Emin.ainvs()],
        "minimal_discriminant_bits": ZZ(abs(Emin.discriminant())).nbits(),
        "ellrankinit_seconds": init_seconds,
        "passes": [asdict(p) for p in passes],
        "eclib_growth": growth,
        "known_point_count": len(points),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=Path, default=Path("artifacts/generated-results/elliptic-curves/newfamily_exact_subgroup_rank_gain_batch_v1.json"))
    ap.add_argument("--min-rank", type=int, default=13)
    ap.add_argument("--parameter", action="append")
    ap.add_argument("--exclude-parameter", action="append")
    ap.add_argument("--efforts", default="0,1,2")
    ap.add_argument("--output", type=Path, default=Path("artifacts/local/elliptic-curves/newfamily/rank_bounds_rank13plus.json"))
    args = ap.parse_args()

    payload = json.loads(args.batch.read_text())
    rows = choose_rows(payload, args)
    efforts = parse_efforts(args.efforts)
    print(f"NFRANK|stage=input|sage={sage_version()}|targets={len(rows)}|efforts={','.join(map(str, efforts))}", flush=True)

    records = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for row in rows:
        try:
            rec = classify(row, efforts)
            rec["status"] = "completed"
        except Exception as exc:
            rec = {"status": "error", "parameter": row.get("parameter"), "error": repr(exc)}
            print(f"NFRANK|T={row.get('parameter')}|status=error|error={exc!r}", flush=True)
        records.append(rec)
        args.output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")

    exact = [r for r in records if r.get("status") == "completed" and str(r.get("classification", "")).startswith("exact_rank_")]
    unresolved = [r for r in records if r.get("status") == "completed" and r.get("classification") == "residual_rank_interval"]
    print(f"NFRANK|stage=done|completed={sum(r.get('status') == 'completed' for r in records)}|exact={len(exact)}|unresolved={len(unresolved)}|output={args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
