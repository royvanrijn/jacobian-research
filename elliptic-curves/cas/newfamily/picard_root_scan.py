#!/usr/bin/env sage -python
"""
Staged Picard-jump search over the exact max-root-100 six-root Mestre census.

Pipeline:
  exact max-root-100 census (235 nonsingular nonreflection families)
    -> specialize the 12 automatic sections at a small admissible T
    -> exact eclib subgroup rank
    -> retain only automatic-rank-11 families
    -> clean-prime Picard probes using F_p, F_(p^2), F_(p^3)
    -> reject immediately when rho(reduction)=16
    -> retain families surviving several independent clean primes

A rank-11 specialization proves the 11 generic section directions independent.
For those families the trivial lattice rank is 5 (20 I1 + split I4), so the
known divisor rank is 16 and the residual H^2 Frobenius factor has degree 6.
Three traces reconstruct it exactly.

A finite-field rho > 16 is NOT a characteristic-zero Picard jump.  A family is
only a search survivor until an extra characteristic-zero divisor/section is
constructed or a separate argument proves the jump.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from fractions import Fraction

from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CAS = REPO / "elliptic-curves" / "cas"
sys.path.insert(0, str(CAS))
sys.path.insert(0, str(HERE))

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    classify_nonreflection,
    primitive_visible_points,
    quartic_point_to_jacobian,
    tuple_digest,
    verify_enumerator_records,
)
from search_mestre_root_tuple_scale_max100 import (
    EXPECTED_MAX100_COUNTS,
    EXPECTED_MAX100_NONSINGULAR_SHA256,
    compiled_enumeration_max100,
)
from screen_seeded_rational_candidates_fast import to_mwrank_triple, qbits
from picard_probe import (
    clean_model,
    count_extension,
    residual_polynomial,
    cyclotomic_part,
    artin_tate_square_class,
)

Q = Fraction

DEFAULT_SPECIALIZATIONS = (1,2,3,5,7,11,13,17,19)
DEFAULT_PRIMES = tuple(int(p) for p in prime_range(43, 211))
TRIVIAL_RANK = 5
TARGET_AUTOMATIC_RANK = 11
KNOWN_DIVISOR_RANK = TRIVIAL_RANK + TARGET_AUTOMATIC_RANK

def py(v):
    if isinstance(v, (Integer, Rational)):
        return int(v) if v in ZZ else str(v)
    if isinstance(v, dict):
        return {str(k): py(x) for k,x in v.items()}
    if isinstance(v, (list, tuple)):
        return [py(x) for x in v]
    return v

def root_id(roots):
    return "r" + "_".join(map(str, roots))

def eclib_rank_of_automatic_sections(roots):
    C = SixRootMestreConstruction(tuple(Q(int(r)) for r in roots))

    last_error = None
    for t in DEFAULT_SPECIALIZATIONS:
        tq = Q(t)
        try:
            if C.quartic_discriminant(tq) == 0:
                continue
            deg = C.visible_point_degeneracy(tq)
            if deg.collision_loss or deg.zero_ordinates:
                continue

            qpts = primitive_visible_points(C, tq)
            if len(qpts) != 12:
                continue
            jpts = tuple(quartic_point_to_jacobian(C, tq, P) for P in qpts)

            coeff = C.primitive_jacobian_coefficients(tq)
            E = EllipticCurve(QQ, list(map(QQ, coeff)))
            Emin = E.global_minimal_model()
            iso = E.isomorphism_to(Emin)

            points = []
            for x,y in jpts:
                P = iso(E([QQ(x), QQ(y)]))
                if not P.is_zero():
                    points.append(P)

            # deterministic small-coordinate order improves mwrank processing
            points.sort(key=lambda P: max(qbits(P[0]), qbits(P[1])))

            mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
            mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=16)
            growth = []
            for i,P in enumerate(points):
                before = len(mw.points())
                mw.process([to_mwrank_triple(P)], saturation_bound=0)
                after = len(mw.points())
                if after > before:
                    growth.append(i)

            rank = len(mw.points())
            return {
                "status": "ok",
                "specialization": t,
                "rank": int(rank),
                "rank_increase_indices": growth,
                "point_count": len(points),
            }
        except Exception as e:
            last_error = repr(e)
            continue

    return {
        "status": "failed",
        "rank": None,
        "error": last_error,
    }

def probe_family(roots, primes, max_probes, progress_every):
    probes = []
    for p in primes:
        if len(probes) >= max_probes:
            break
        try:
            D = clean_model(tuple(roots), int(p))
        except Exception:
            continue

        print(f"PICARD|roots={','.join(map(str,roots))}|p={p}|stage=start", flush=True)
        counts = {}
        try:
            for n in (1,2,3):
                counts[n] = count_extension(D, int(p), n, progress_every)
            Pres, residual_traces = residual_polynomial(int(p), counts)
            norm, rem, cyc, cycdegree, orders = cyclotomic_part(int(p), Pres)
            rho = int(KNOWN_DIVISOR_RANK + cycdegree)

            R = Pres.parent()
            X = R.gen()
            fullP = (X-p)**KNOWN_DIVISOR_RANK * Pres
            at = artin_tate_square_class(int(p), fullP, rho, orders)

            rec = {
                "prime": int(p),
                "rho_reduction": rho,
                "mw_reduction": rho - TRIVIAL_RANK,
                "residual_polynomial": str(Pres),
                "normalized_residual_polynomial": str(norm),
                "cyclotomic_factors": cyc,
                "noncyclotomic_remainder": str(rem),
                "residual_traces": residual_traces,
                "artin_tate": at,
                "counts": {str(k): counts[k] for k in counts},
            }
            probes.append(rec)

            print(
                f"PICARD|roots={','.join(map(str,roots))}|p={p}|"
                f"rho={rho}|disc={at['signed_NS_discriminant_square_class']}|"
                f"Pres={Pres}",
                flush=True,
            )

            # One clean rho=16 reduction proves char-0 rho=16 exactly.
            if rho == 16:
                return probes, "closed_rho16"

        except Exception as e:
            print(
                f"PICARD|roots={','.join(map(str,roots))}|p={p}|"
                f"stage=error|error={type(e).__name__}:{e}",
                flush=True,
            )
            continue

    if len(probes) >= max_probes:
        return probes, "survived_all_probes"
    return probes, "insufficient_clean_probes"

def invariants(roots):
    roots = [ZZ(r) for r in roots]
    # translate to sum zero over Q; elementary symmetric invariants there
    mean = QQ(sum(roots)) / 6
    centered = [QQ(r)-mean for r in roots]
    R = PolynomialRing(QQ, "z")
    z = R.gen()
    f = prod(z-r for r in centered)
    # f=z^6-e1 z^5+e2 z^4-e3 z^3+...
    e = {}
    for k in range(2,7):
        coeff = f[6-k]
        e[k] = ((-1)**k) * coeff
    return {
        "centered_e2": str(e[2]),
        "centered_e3": str(e[3]),
        "centered_e4": str(e[4]),
        "centered_e5": str(e[5]),
        "centered_e6": str(e[6]),
        "quartic_relation_2e5_minus_e2e3": str(2*e[5]-e[2]*e[3]),
    }

def build_census(args):
    source = CAS / "enumerate_mestre_root_tuples_scale.cpp"
    print("CENSUS|stage=enumerate", flush=True)
    enumeration, timings = compiled_enumeration_max100(
        source,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    got = (
        enumeration.normalized_count,
        enumeration.obstruction_count,
        enumeration.reflection_count,
        enumeration.nonreflection_count,
    )
    if got != EXPECTED_MAX100_COUNTS:
        raise RuntimeError(f"max100 census changed: {got} != {EXPECTED_MAX100_COUNTS}")
    verify_enumerator_records(enumeration)
    nonsingular, singular, witnesses = classify_nonreflection(enumeration)
    if tuple_digest(nonsingular) != EXPECTED_MAX100_NONSINGULAR_SHA256:
        raise RuntimeError("max100 nonsingular-family digest changed")
    print(
        f"CENSUS|status=complete|nonsingular={len(nonsingular)}|"
        f"enumeration_seconds={timings['enumeration_wall_seconds']:.3f}",
        flush=True,
    )
    return tuple(nonsingular)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compile-timeout", type=float, default=30.0)
    ap.add_argument("--enumeration-timeout", type=float, default=30.0)
    ap.add_argument("--max-probes", type=int, default=3)
    ap.add_argument("--prime-min", type=int, default=43)
    ap.add_argument("--prime-max", type=int, default=200)
    ap.add_argument("--progress-every", type=int, default=0)
    ap.add_argument("--shard-count", type=int, default=1)
    ap.add_argument("--shard-index", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0,
                    help="debug: limit number of shard families, 0=all")
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
    )
    args = ap.parse_args()

    if args.shard_count < 1 or not (0 <= args.shard_index < args.shard_count):
        raise SystemExit("invalid shard")
    if args.max_probes < 1:
        raise SystemExit("--max-probes must be >=1")

    primes = tuple(int(p) for p in prime_range(args.prime_min, args.prime_max+1))
    families = build_census(args)
    shard = [r for i,r in enumerate(families) if i % args.shard_count == args.shard_index]
    if args.limit:
        shard = shard[:args.limit]

    print(
        f"SCAN|families_total={len(families)}|shard={args.shard_index}/"
        f"{args.shard_count}|families_here={len(shard)}",
        flush=True,
    )

    records = []
    rank11 = []
    for idx, roots in enumerate(shard, 1):
        rid = root_id(roots)
        rr = eclib_rank_of_automatic_sections(roots)
        rec = {
            "id": rid,
            "roots": list(map(int, roots)),
            "automatic_rank": rr,
            "root_invariants": invariants(roots),
        }
        records.append(rec)
        print(
            f"RANK|{idx}/{len(shard)}|roots={','.join(map(str,roots))}|"
            f"status={rr['status']}|T={rr.get('specialization')}|rank={rr.get('rank')}",
            flush=True,
        )
        if rr.get("rank") == TARGET_AUTOMATIC_RANK:
            rank11.append(rec)

    print(f"RANK11|count={len(rank11)}", flush=True)

    survivors = []
    closed = []
    for i,rec in enumerate(rank11, 1):
        roots = rec["roots"]
        print(
            f"FROBENIUS|{i}/{len(rank11)}|roots={','.join(map(str,roots))}",
            flush=True,
        )
        probes, status = probe_family(
            roots, primes, args.max_probes, args.progress_every
        )
        rec["picard_probes"] = probes
        rec["picard_status"] = status
        rec["minimum_reduction_rho"] = (
            min(p["rho_reduction"] for p in probes) if probes else None
        )

        if status == "closed_rho16":
            closed.append(rec)
        else:
            survivors.append(rec)

    result = {
        "schema": "newfamily.picard-root-scan.v1",
        "claim_boundary": (
            "rank-11 specialization proves the displayed generic subgroup has "
            "rank 11; rho=16 at one good reduction proves geometric generic MW "
            "rank exactly 11. Surviving several rho>16 reductions is only a "
            "search signal, not proof of a characteristic-zero Picard jump."
        ),
        "max_root": 100,
        "nonsingular_family_count": len(families),
        "shard_count": args.shard_count,
        "shard_index": args.shard_index,
        "families_scanned_this_shard": len(shard),
        "automatic_rank11_count": len(rank11),
        "closed_at_rho16_count": len(closed),
        "survivor_count": len(survivors),
        "survivors": survivors,
        "rank11_families": rank11,
        "all_rank_screen_records": records,
    }

    if args.output is None:
        suffix = (
            f".shard_{args.shard_index:03d}_of_{args.shard_count:03d}"
            if args.shard_count > 1 else ""
        )
        out = (
            REPO / "artifacts/local/elliptic-curves/newfamily"
            / f"picard_root_scan_max100{suffix}.json"
        )
    else:
        out = args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(py(result), indent=2, sort_keys=True, default=str) + "\n")

    print(
        f"DONE|rank11={len(rank11)}|closed={len(closed)}|"
        f"survivors={len(survivors)}|output={out}",
        flush=True,
    )
    for rec in survivors:
        print(
            "SURVIVOR|roots=" + ",".join(map(str, rec["roots"]))
            + "|rhos=" + ",".join(str(p["rho_reduction"]) for p in rec["picard_probes"]),
            flush=True,
        )

if __name__ == "__main__":
    main()
