#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import random
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients


PROTOCOL = "R30REL2"

# Exact bad rational-prime support already recovered from the
# 2-division polynomial discriminant.
BAD_RATIONAL_PRIMES = (
    2,
    3,
    5,
    7,
    13,
    31,
    41,
    47,
    53,
    67,
    379,
    4349,
    25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
)


def sage_q(QQ, ZZ, value):
    n = value.numerator
    d = value.denominator

    if callable(n):
        n = n()
    if callable(d):
        d = d()

    return QQ(ZZ(n)) / QQ(ZZ(d))


def insert_row(pivots, row):
    """
    Insert a packed F2 row.

    Returns True iff rank increased.
    """
    v = int(row)

    while v:
        k = v.bit_length() - 1

        if k in pivots:
            v ^= pivots[k]
        else:
            pivots[k] = v
            return True

    return False


def project(row, cols):
    out = 0

    for j, c in enumerate(cols):
        if (row >> c) & 1:
            out |= 1 << j

    return out


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--factor-base-bound",
        type=int,
        default=500,
    )

    ap.add_argument(
        "--box",
        type=int,
        default=12,
    )

    ap.add_argument(
        "--random-trials",
        type=int,
        default=20000,
    )

    ap.add_argument(
        "--stall",
        type=int,
        default=5000,
    )

    ap.add_argument(
        "--seed",
        type=int,
        default=20260820,
    )

    args = ap.parse_args()

    from sage.all import (
        QQ,
        ZZ,
        NumberField,
        PolynomialRing,
        prime_range,
    )

    coeffs = short_coefficients()

    A = ZZ(sage_q(QQ, ZZ, coeffs[3]))
    B = ZZ(sage_q(QQ, ZZ, coeffs[4]))

    x = PolynomialRing(QQ, "x").gen()

    f = x**3 + A*x + B

    print(
        f"{PROTOCOL}|stage=input"
        f"|A_bits={abs(A).nbits()}"
        f"|B_bits={abs(B).nbits()}"
        f"|disc_bits={abs(ZZ(f.discriminant())).nbits()}"
        f"|fb_bound={args.factor_base_bound}"
        f"|box={args.box}",
        flush=True,
    )

    # --------------------------------------------------------
    # Cubic field + maximal order
    # --------------------------------------------------------

    t0 = time.monotonic()

    K = NumberField(f, "a")

    print(
        f"{PROTOCOL}|stage=number_field"
        f"|seconds={time.monotonic()-t0:.6f}",
        flush=True,
    )

    t0 = time.monotonic()

    OK = K.maximal_order()
    basis = OK.basis()

    print(
        f"{PROTOCOL}|stage=maximal_order"
        f"|seconds={time.monotonic()-t0:.6f}"
        f"|basis={basis}",
        flush=True,
    )

    bad = tuple(
        ZZ(p)
        for p in BAD_RATIONAL_PRIMES
    )

    # --------------------------------------------------------
    # Factor base
    #
    # Includes every S-prime, even when the rational prime is
    # much larger than factor-base-bound.
    # --------------------------------------------------------

    fb = []
    fb_rational = []
    s_cols = []

    def idx_of(P):
        for i, Q in enumerate(fb):
            if P == Q:
                return i
        return None

    def add(P, rational_prime, is_s):
        i = idx_of(P)

        if i is None:
            i = len(fb)
            fb.append(P)
            fb_rational.append(int(rational_prime))

        if is_s and i not in s_cols:
            s_cols.append(i)

        return i

    t0 = time.monotonic()

    for p in bad:
        places = list(
            K.primes_above(p)
        )

        print(
            f"{PROTOCOL}|stage=S_prime"
            f"|p={p}"
            f"|places={len(places)}",
            flush=True,
        )

        for P in places:
            add(P, p, True)

    for q in prime_range(
        3,
        args.factor_base_bound + 1,
    ):
        q = ZZ(q)

        if q in bad:
            continue

        for P in K.primes_above(q):
            add(P, q, False)

    print(
        f"{PROTOCOL}|stage=factor_base"
        f"|columns={len(fb)}"
        f"|S_columns={len(s_cols)}"
        f"|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )

    # Map rational prime -> prime ideals in the factor base.
    fb_by_p = {}

    for idx, (P, q) in enumerate(
        zip(fb, fb_rational)
    ):
        fb_by_p.setdefault(
            int(q),
            [],
        ).append(
            (idx, P)
        )

    fb_rational_primes = sorted(
        fb_by_p
    )

    # --------------------------------------------------------
    # Relation spaces
    # --------------------------------------------------------

    pivots = {}
    s_pivots = {}

    sampled = 0
    smooth = 0
    accepted = 0

    last_gain = 0

    def process(alpha, label):
        nonlocal sampled
        nonlocal smooth
        nonlocal accepted
        nonlocal last_gain

        sampled += 1

        if alpha == 0:
            return

        N = abs(
            ZZ(alpha.norm())
        )

        if N == 0:
            return

        # ----------------------------------------------------
        # Factor-base smoothness test.
        #
        # Deliberately DO NOT factor N.
        #
        # Divide by our known rational-prime factor base only.
        # ----------------------------------------------------

        cofactor = N
        used_q = []

        for q in fb_rational_primes:

            if cofactor % q != 0:
                continue

            used_q.append(q)

            while cofactor % q == 0:
                cofactor //= q

            if cofactor == 1:
                break

        if cofactor != 1:
            return

        smooth += 1

        # ----------------------------------------------------
        # Exact prime-ideal valuation parity.
        # ----------------------------------------------------

        row = 0

        for q in used_q:

            for index, P in fb_by_p[q]:

                exponent = int(
                    alpha.valuation(P)
                )

                if exponent & 1:
                    row ^= 1 << index

        if row == 0:
            return

        accepted += 1

        gained = insert_row(
            pivots,
            row,
        )

        srow = project(
            row,
            s_cols,
        )

        if srow:
            insert_row(
                s_pivots,
                srow,
            )

        if gained:
            last_gain = sampled

            print(
                f"{PROTOCOL}|stage=relation"
                f"|status=RANK_GAIN"
                f"|label={label}"
                f"|sampled={sampled}"
                f"|smooth={smooth}"
                f"|accepted={accepted}"
                f"|rank={len(pivots)}"
                f"|fb_qdim={len(fb)-len(pivots)}"
                f"|S_proj_rank={len(s_pivots)}",
                flush=True,
            )

    # --------------------------------------------------------
    # Deterministic coefficient box
    # --------------------------------------------------------

    box = args.box

    print(
        f"{PROTOCOL}|stage=deterministic"
        f"|status=START",
        flush=True,
    )

    t0 = time.monotonic()

    for a in range(-box, box + 1):

        for b in range(-box, box + 1):

            for c in range(-box, box + 1):

                if a == b == c == 0:
                    continue

                # Global sign canonicalization.
                first = (
                    a
                    if a
                    else (
                        b
                        if b
                        else c
                    )
                )

                if first < 0:
                    continue

                # Skip obvious integer multiples.
                if math.gcd(
                    abs(a),
                    math.gcd(
                        abs(b),
                        abs(c),
                    ),
                ) > 1:
                    continue

                alpha = (
                    a*basis[0]
                    + b*basis[1]
                    + c*basis[2]
                )

                process(
                    alpha,
                    f"box:{a},{b},{c}",
                )

    print(
        f"{PROTOCOL}|stage=deterministic"
        f"|status=COMPLETE"
        f"|seconds={time.monotonic()-t0:.3f}"
        f"|sampled={sampled}"
        f"|smooth={smooth}"
        f"|accepted={accepted}"
        f"|rank={len(pivots)}",
        flush=True,
    )

    # --------------------------------------------------------
    # Random extension
    # --------------------------------------------------------

    rng = random.Random(
        args.seed
    )

    print(
        f"{PROTOCOL}|stage=random"
        f"|status=START"
        f"|trials={args.random_trials}"
        f"|seed={args.seed}",
        flush=True,
    )

    for trial in range(
        args.random_trials
    ):

        # Every fifth sample explores a larger box.
        scale = (
            box
            if trial % 5
            else 5*box
        )

        a, b, c = [
            rng.randint(
                -scale,
                scale,
            )
            for _ in range(3)
        ]

        if a == b == c == 0:
            continue

        alpha = (
            a*basis[0]
            + b*basis[1]
            + c*basis[2]
        )

        process(
            alpha,
            f"rnd:{a},{b},{c}",
        )

        if (
            sampled - last_gain
            >= args.stall
        ):
            print(
                f"{PROTOCOL}|stage=random"
                f"|status=STALLED"
                f"|since_gain={sampled-last_gain}",
                flush=True,
            )
            break

    # --------------------------------------------------------
    # Quotient after killing S primes.
    #
    # This is a FACTOR-BASE MODEL only until completeness of
    # the principal relation set has been certified.
    # --------------------------------------------------------

    kill = dict(
        pivots
    )

    for index in s_cols:
        insert_row(
            kill,
            1 << index,
        )

    relation_rank = len(pivots)

    fb_qdim = (
        len(fb)
        - relation_rank
    )

    after_killing_S = (
        len(fb)
        - len(kill)
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|fb_columns={len(fb)}"
        f"|S_columns={len(s_cols)}"
        f"|relation_rank={relation_rank}"
        f"|fb_qdim={fb_qdim}"
        f"|S_projected_rank={len(s_pivots)}"
        f"|after_killing_S={after_killing_S}"
        f"|sampled={sampled}"
        f"|smooth={smooth}"
        f"|accepted={accepted}",
        flush=True,
    )

    if after_killing_S == 0:

        print(
            f"{PROTOCOL}|result="
            "HEURISTIC_S_CLASS_MOD2_ZERO"
            "|exact_relations=true"
            "|completeness_proved=false",
            flush=True,
        )

    else:

        print(
            f"{PROTOCOL}|result="
            "RESIDUAL_S_CLASS_DIRECTIONS_REMAIN"
            f"|dimension_model={after_killing_S}"
            "|exact_relations=true"
            "|completeness_proved=false",
            flush=True,
        )

    print(
        f"{PROTOCOL}|warning="
        "factor_base_model_is_not_a_selmer_upper_bound;"
        "auxiliary_good_primes_are_not_selmer_conditions;"
        "relation_completeness_not_yet_certified",
        flush=True,
    )


if __name__ == "__main__":
    main()
