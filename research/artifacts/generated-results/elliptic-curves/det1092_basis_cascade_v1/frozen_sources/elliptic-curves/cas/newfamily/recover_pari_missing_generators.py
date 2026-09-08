#!/usr/bin/env sage -python
"""Recover explicit PARI generators missing from the pinned newfamily subgroup.

For selected specializations, replay the pinned subgroup from git, ask PARI ellrank
for a rank interval and returned points, then exact-test those returned points with
eclib against the pinned subgroup.  Any point that increases the exact processed
subgroup rank is recorded as a recovered missing generator.

The driver first tries ellrank on an ellrankinit context and falls back to direct
ellrank on the PARI curve.  This is useful for fibers where ellrankinit raised a
PARI exception in the breadth-first classifier.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from sage.all import QQ, ZZ, pari, version as sage_version
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from batch_rank_bounds import parse_parameter, replay_row
from screen_seeded_rational_candidates_fast import to_mwrank_triple


def _pari_int(value):
    try:
        return int(value)
    except TypeError:
        return int(value.sage())


def _sage_point(E, item):
    if len(item) < 2:
        raise ValueError(f"unexpected PARI point {item}")
    return E([QQ(item[0]), QQ(item[1])])


def call_ellrank(E, known_points, effort):
    pp = pari([[P[0], P[1]] for P in known_points])
    errors = []
    started = time.monotonic()
    try:
        ctx = E.pari_curve().ellrankinit()
        init_seconds = time.monotonic() - started
        rs = time.monotonic()
        result = ctx.ellrank(effort, pp)
        return "ellrankinit", init_seconds, time.monotonic() - rs, result, errors
    except Exception as exc:
        errors.append(f"ellrankinit:{exc!r}")

    started = time.monotonic()
    try:
        result = pari.ellrank(E.pari_curve(), effort, pp)
        return "direct", 0.0, time.monotonic() - started, result, errors
    except Exception as exc:
        errors.append(f"direct:{exc!r}")
        raise RuntimeError("; ".join(errors)) from exc


def recover(row, effort):
    parameter = row["parameter"]
    pinned_lower = int(row["processed_subgroup_rank"])
    E, known_points, growth = replay_row(row)

    method, init_seconds, rank_seconds, result, errors = call_ellrank(E, known_points, effort)
    if len(result) != 4:
        raise RuntimeError(f"unexpected ellrank result {result}")
    pl, pu, sha = _pari_int(result[0]), _pari_int(result[1]), _pari_int(result[2])
    effective_lower = max(pinned_lower, pl)
    returned = [_sage_point(E, item) for item in result[3]]

    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in E.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1,
                            maxr=max(32, pu + len(returned) + 8))
    for P in known_points:
        mw.process([to_mwrank_triple(P)], saturation_bound=0)
    if len(mw.points()) != pinned_lower:
        raise RuntimeError(f"{parameter}: replay subgroup rank changed")

    recovered = []
    seen = set()
    for idx, P in enumerate(returned):
        key = (str(P[0]), min(str(P[1]), str(-P[1])))
        if key in seen:
            continue
        seen.add(key)
        before = len(mw.points())
        mw.process([to_mwrank_triple(P)], saturation_bound=0)
        after = len(mw.points())
        if after > before:
            recovered.append({
                "pari_index": idx,
                "point_minimal": [str(P[0]), str(P[1])],
                "rank_before": before,
                "rank_after": after,
            })
            print(
                f"NFRECOVER|T={parameter}|point={idx}|rank={before}->{after}"
                f"|x={P[0]}|y={P[1]}", flush=True,
            )

    final_rank = len(mw.points())
    classification = f"exact_rank_{effective_lower}" if effective_lower == pu else "residual_rank_interval"
    return {
        "status": "completed",
        "parameter": parameter,
        "pinned_lower_bound": pinned_lower,
        "rank_interval": [effective_lower, pu],
        "classification": classification,
        "pari_lower": pl,
        "pari_upper": pu,
        "sha_pairing_rank": sha,
        "pari_method": method,
        "ellrankinit_seconds": init_seconds,
        "ellrank_seconds": rank_seconds,
        "fallback_errors": errors,
        "pari_returned_point_count": len(returned),
        "recovered_missing_generators": recovered,
        "recovered_generator_count": len(recovered),
        "processed_rank_after_returned_points": final_rank,
        "root_number": int(E.root_number()),
        "minimal_model": [str(v) for v in E.ainvs()],
        "minimal_discriminant_bits": ZZ(abs(E.discriminant())).nbits(),
        "pinned_growth": growth,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--batch", type=Path,
        default=Path("artifacts/generated-results/elliptic-curves/newfamily_exact_subgroup_rank_gain_batch_v1.json"),
    )
    ap.add_argument("--parameter", action="append", required=True)
    ap.add_argument("--effort", type=int, default=0)
    ap.add_argument(
        "--output", type=Path,
        default=Path("artifacts/local/elliptic-curves/newfamily/recovered_pari_generators.json"),
    )
    args = ap.parse_args()

    payload = json.loads(args.batch.read_text())
    by_parameter = {r["parameter"]: r for r in payload.get("records", [])}
    print(
        f"NFRECOVER|stage=input|sage={sage_version()}|targets={len(args.parameter)}"
        f"|effort={args.effort}", flush=True,
    )

    records = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for parameter in args.parameter:
        row = by_parameter.get(parameter)
        if row is None:
            rec = {"status": "error", "parameter": parameter, "error": "not in pinned batch"}
        else:
            print(f"NFRECOVER|T={parameter}|stage=start|pinned={row.get('processed_subgroup_rank')}", flush=True)
            try:
                rec = recover(row, args.effort)
                print(
                    f"NFRECOVER|T={parameter}|stage=done|interval={rec['rank_interval'][0]},{rec['rank_interval'][1]}"
                    f"|recovered={rec['recovered_generator_count']}|method={rec['pari_method']}",
                    flush=True,
                )
            except Exception as exc:
                rec = {"status": "error", "parameter": parameter, "error": repr(exc)}
                print(f"NFRECOVER|T={parameter}|stage=error|error={exc!r}", flush=True)
        records.append(rec)
        args.output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")

    print(
        f"NFRECOVER|stage=done|completed={sum(r.get('status') == 'completed' for r in records)}"
        f"|output={args.output}", flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
