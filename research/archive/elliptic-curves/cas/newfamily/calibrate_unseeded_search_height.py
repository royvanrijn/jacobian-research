#!/usr/bin/env sage -python
"""Calibrate empty eclib point-search height on the known T=11 specialization.

The purpose is diagnostic: measure the actual naive-log heights of the eleven
known hidden sections on a global minimal model, compare them with canonical
heights, then run fresh unseeded eclib searches at requested height limits.

This tells us whether H=12 (or another bound) is even large enough to
rediscover known points in this family before launching broad rational scans.
"""

from __future__ import annotations

import argparse
import math
import time

from sage.all import EllipticCurve, QQ, ZZ
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from screen_seeded_rational_candidates_fast import ROOTS, load_builder, load_sections


def x_naive_log_height(P):
    if P.is_zero():
        return float("-inf")
    x = QQ(P[0])
    return math.log(max(abs(int(x.numerator())), int(x.denominator())))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--T", type=int, default=11)
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--heights", default="8,10,12,14,16")
    p.add_argument("--maxr", type=int, default=64)
    p.add_argument("--verbose-search", action="store_true")
    args = p.parse_args()

    t = QQ(args.T)
    family = load_builder()(ROOTS)
    sections = load_sections(args.sections_sobj)
    A = QQ(family["Amin"](t))
    B = QQ(family["Bmin"](t))
    E = EllipticCurve(QQ, [0, 0, 0, A, B])
    known = [E([QQ(xf(t)), QQ(yf(t))]) for xf, yf in sections]

    print(f"T={args.T} discriminant_bits={ZZ(abs(E.discriminant())).nbits()}", flush=True)
    started = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_s = time.monotonic() - started
    iso = E.isomorphism_to(Emin)
    known_min = [iso(P) for P in known]
    print(
        f"minimal_s={minimal_s:.6f} minimal_disc_bits={ZZ(abs(Emin.discriminant())).nbits()} "
        f"ainvs={list(Emin.ainvs())}",
        flush=True,
    )

    rows = []
    print("KNOWN SECTION HEIGHTS ON MINIMAL MODEL", flush=True)
    for i, P in enumerate(known_min):
        naive = x_naive_log_height(P)
        canonical = float(P.height())
        x = QQ(P[0])
        row = (i, naive, canonical, ZZ(abs(x.numerator())).nbits(), ZZ(x.denominator()).nbits())
        rows.append(row)
        print(
            f"U{i:02d} naive_log_x={naive:.6f} canonical={canonical:.6f} "
            f"x_num_bits={row[3]} x_den_bits={row[4]}",
            flush=True,
        )

    print(
        "SUMMARY "
        f"min_naive={min(r[1] for r in rows):.6f} "
        f"median_naive={sorted(r[1] for r in rows)[5]:.6f} "
        f"max_naive={max(r[1] for r in rows):.6f} "
        f"min_canonical={min(r[2] for r in rows):.6f} "
        f"max_canonical={max(r[2] for r in rows):.6f}",
        flush=True,
    )

    known_signless = set(known_min) | set(-P for P in known_min)
    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])

    for H in [float(x.strip()) for x in args.heights.split(",") if x.strip()]:
        mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=0, maxr=args.maxr)
        s = time.monotonic()
        mw.search(H, verbose=args.verbose_search)
        elapsed = time.monotonic() - s
        raw = list(mw.points())

        exact_known = 0
        other = 0
        for triple in raw:
            X, Y, Z = [ZZ(v) for v in triple]
            if Z == 0:
                continue
            Q = Emin([QQ(X) / QQ(Z), QQ(Y) / QQ(Z)])
            if Q in known_signless or -Q in known_signless:
                exact_known += 1
            else:
                other += 1

        print(
            f"SEARCH H={H:g} seconds={elapsed:.6f} raw={len(raw)} "
            f"exact_hidden={exact_known} other={other}",
            flush=True,
        )


if __name__ == "__main__":
    main()
