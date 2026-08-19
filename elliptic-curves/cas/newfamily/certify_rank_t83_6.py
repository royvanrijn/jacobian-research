#!/usr/bin/env sage -python
"""Replay the pinned T=83/6 rank-14 subgroup and ask PARI for rank bounds.

This driver is git-only: the family, hidden generic sections, and three extra
rank-increasing points are loaded from committed source/artifacts.  It first
replays the exact baseline-first eclib rank growth 11 -> 14, then calls PARI
``ellrank`` with the fourteen known independent points.

A result [14,14] proves exact rank 14.  Any wider interval is recorded without
being promoted to an exact-rank claim.
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

ROOTS = [-47, -43, -31, 30, 45, 46]
PARAMETER = "83/6"
KNOWN_GENERIC_RANK = 11
EXPECTED_LOWER_BOUND = 14


@dataclass(frozen=True)
class PariPass:
    effort: int
    pari_lower: int
    pari_upper: int
    effective_lower: int
    sha_pairing_rank: int
    returned_points: int
    seconds: float


def parse_efforts(text: str):
    values = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if not values or any(x < 0 for x in values):
        raise ValueError("--efforts must be nonempty nonnegative integers")
    return tuple(dict.fromkeys(values))


def load_record(path: Path):
    payload = json.loads(path.read_text())
    if payload.get("family") != "newfamily six-root quartic":
        raise RuntimeError("unexpected family in pinned batch")
    rows = [r for r in payload.get("records", []) if r.get("parameter") == PARAMETER]
    if len(rows) != 1:
        raise RuntimeError(f"expected one {PARAMETER} record, found {len(rows)}")
    row = rows[0]
    if row.get("known_subgroup_rank") != KNOWN_GENERIC_RANK:
        raise RuntimeError("pinned baseline is not rank 11")
    if row.get("processed_subgroup_rank") != EXPECTED_LOWER_BOUND:
        raise RuntimeError("pinned processed subgroup is not rank 14")
    if row.get("exact_rank_gain_over_known") != 3:
        raise RuntimeError("pinned T=83/6 record does not contain three extra directions")
    if len(row.get("hit_rank_increases", [])) != 3:
        raise RuntimeError("expected exactly three rank-increasing hit records")
    return row


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--batch",
        type=Path,
        default=Path("artifacts/generated-results/elliptic-curves/newfamily_rank_gain_batch_v1.json"),
    )
    ap.add_argument("--efforts", default="0,1,2")
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/local/elliptic-curves/newfamily/rank_bounds_t83_6.json"),
    )
    args = ap.parse_args()

    row = load_record(args.batch)
    efforts = parse_efforts(args.efforts)

    E, known_integral, _ = build_integral_specialization_git(83, 6)
    Emin = E.global_minimal_model()
    iso = E.isomorphism_to(Emin)
    known = [iso(P) for P in known_integral]

    extras = []
    for hit in row["hit_rank_increases"]:
        x, y = hit["point_minimal"]
        P = Emin([QQ(x), QQ(y)])
        if P.is_zero():
            raise RuntimeError("pinned extra point is zero")
        extras.append(P)
    if len(set(extras)) != 3:
        raise RuntimeError("pinned extra points are not distinct")

    # Exact baseline-first replay with eclib.
    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=32)
    known_order = sorted(range(11), key=lambda i: max(qbits(known[i][0]), qbits(known[i][1])))
    growth = []
    print(f"T83RANK|stage=eclib|status=start|sage={sage_version()}", flush=True)
    for i in known_order:
        before = len(mw.points())
        mw.process([to_mwrank_triple(known[i])], saturation_bound=0)
        after = len(mw.points())
        growth.append({"label": f"U{i}", "before": before, "after": after})
    baseline = len(mw.points())
    if baseline != 11:
        raise RuntimeError(f"git-only known baseline replayed as {baseline}, expected 11")
    for j, P in enumerate(extras):
        before = len(mw.points())
        mw.process([to_mwrank_triple(P)], saturation_bound=0)
        after = len(mw.points())
        growth.append({"label": f"Q{j}", "before": before, "after": after})
        if after != before + 1:
            raise RuntimeError(f"pinned extra Q{j} did not increase rank")
    lower = len(mw.points())
    if lower != 14:
        raise RuntimeError(f"exact subgroup replay ended at {lower}, expected 14")
    print("T83RANK|stage=eclib|status=complete|baseline=11|lower=14", flush=True)

    points14 = known + extras
    print("T83RANK|stage=ellrankinit|status=start", flush=True)
    started = time.monotonic()
    ctx = _ellrankinit(Emin.pari_curve())
    init_seconds = time.monotonic() - started
    print(f"T83RANK|stage=ellrankinit|status=complete|seconds={init_seconds:.6f}", flush=True)

    passes = []
    for effort in efforts:
        print(f"T83RANK|stage=ellrank|status=start|effort={effort}", flush=True)
        started = time.monotonic()
        result = _ellrank(ctx, effort, points14)
        seconds = time.monotonic() - started
        if len(result) != 4:
            raise RuntimeError(f"unexpected ellrank result: {result}")
        pl = _pari_int(result[0])
        pu = _pari_int(result[1])
        sha = _pari_int(result[2])
        returned = len(result[3])
        if pu < 14:
            raise RuntimeError(f"PARI upper bound {pu} contradicts exact rank>=14")
        effective = max(14, pl)
        p = PariPass(effort, pl, pu, effective, sha, returned, seconds)
        passes.append(p)
        print(
            f"T83RANK|stage=ellrank|status=complete|effort={effort}"
            f"|pari_lower={pl}|pari_upper={pu}|effective_lower={effective}"
            f"|sha_pairing_rank={sha}|returned_points={returned}|seconds={seconds:.6f}",
            flush=True,
        )
        if effective == pu:
            break

    best_lower = max(p.effective_lower for p in passes)
    best_upper = min(p.pari_upper for p in passes)
    if best_lower > best_upper:
        raise RuntimeError("PARI passes produced inconsistent intervals")
    classification = "exact_rank_14" if best_lower == best_upper == 14 else "residual_rank_interval"

    payload = {
        "schema": "newfamily_t83_6_rank_bounds_v1",
        "family": "newfamily six-root quartic",
        "roots": ROOTS,
        "parameter": PARAMETER,
        "engine": "Sage/eclib + PARI ellrank",
        "exact_subgroup_lower_bound": 14,
        "rank_interval": [best_lower, best_upper],
        "classification": classification,
        "root_number": int(Emin.root_number()),
        "minimal_model": [str(v) for v in Emin.ainvs()],
        "minimal_discriminant_bits": ZZ(abs(Emin.discriminant())).nbits(),
        "eclib_growth": growth,
        "extra_points_minimal": [[str(P[0]), str(P[1])] for P in extras],
        "ellrankinit_seconds": init_seconds,
        "passes": [asdict(p) for p in passes],
        "claim_boundary": (
            "exact rank only if rank_interval endpoints agree; otherwise the exact statement is rank>=14"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"T83RANK|stage=done|interval={best_lower},{best_upper}|classification={classification}"
        f"|output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
