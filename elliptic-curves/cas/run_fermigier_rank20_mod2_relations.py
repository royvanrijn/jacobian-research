#!/usr/bin/env python3
from __future__ import annotations

import argparse, math, random, time
from pathlib import Path

PROTOCOL = "R20REL2"

def insert_row(pivots, row):
    v = row
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
    ap.add_argument("--factor-base-bound", type=int, default=500)
    ap.add_argument("--box", type=int, default=12)
    ap.add_argument("--random-trials", type=int, default=10000)
    ap.add_argument("--stall", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=20260820)
    args = ap.parse_args()

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    f = x**3 - ZZ(5750886029903523759416717668139307)*x + ZZ(167347710468055045100164888198438918505621536951206)
    K = NumberField(f, "a")
    OK = K.maximal_order()
    B = OK.basis()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))

    print(f"{PROTOCOL}|stage=input|bad={bad}|fb_bound={args.factor_base_bound}|box={args.box}", flush=True)

    fb = []
    s_cols = []

    def idx_of(P):
        for i, Q in enumerate(fb):
            if P == Q:
                return i
        return None

    def add(P, is_s):
        i = idx_of(P)
        if i is None:
            i = len(fb)
            fb.append(P)
        if is_s and i not in s_cols:
            s_cols.append(i)

    for p in bad:
        for P in K.primes_above(p):
            add(P, True)

    for q in prime_range(3, args.factor_base_bound + 1):
        q = ZZ(q)
        if q in bad:
            continue
        for P in K.primes_above(q):
            add(P, False)

    print(f"{PROTOCOL}|stage=factor_base|columns={len(fb)}|S_columns={len(s_cols)}", flush=True)

    # Rational-prime factor-base map.  Relation collection should NEVER fully
    # factor a huge norm: trial-divide only by this known factor base.
    fb_by_p = {}
    for idx, P in enumerate(fb):
        q = int(P.smallest_integer())
        fb_by_p.setdefault(q, []).append((idx, P))
    fb_rational_primes = sorted(fb_by_p)

    pivots = {}
    s_pivots = {}
    sampled = smooth = accepted = 0
    last_gain = 0

    def process(alpha, label):
        nonlocal sampled, smooth, accepted, last_gain
        sampled += 1
        if alpha == 0:
            return
        # Fast factor-base smoothness test.  Do NOT call factor(N): these norms
        # can be enormous.  Divide only by rational primes represented in the
        # factor base.  If anything remains, reject this sample immediately.
        N = abs(ZZ(alpha.norm()))
        if N == 0:
            return

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

        # The norm is rational-factor-base smooth.  Now compute the exact
        # prime-ideal valuation parity only above rational primes that occurred.
        row = 0
        for q in used_q:
            for i, P in fb_by_p[q]:
                e = int(alpha.valuation(P))
                if e & 1:
                    row ^= 1 << i
        smooth += 1
        if row == 0:
            return
        accepted += 1
        gained = insert_row(pivots, row)
        sr = project(row, s_cols)
        if sr:
            insert_row(s_pivots, sr)
        if gained:
            last_gain = sampled
            print(
                f"{PROTOCOL}|stage=relation|status=rank_gain|label={label}"
                f"|sampled={sampled}|smooth={smooth}|rank={len(pivots)}"
                f"|fb_qdim={len(fb)-len(pivots)}|S_proj_rank={len(s_pivots)}",
                flush=True,
            )

    box = args.box
    print(f"{PROTOCOL}|stage=deterministic|status=start", flush=True)
    t0 = time.monotonic()
    for a in range(-box, box + 1):
        for b in range(-box, box + 1):
            for c in range(-box, box + 1):
                if a == b == c == 0:
                    continue
                first = a if a else (b if b else c)
                if first < 0:
                    continue
                if math.gcd(abs(a), math.gcd(abs(b), abs(c))) > 1:
                    continue
                process(a*B[0] + b*B[1] + c*B[2], f"box:{a},{b},{c}")
    print(f"{PROTOCOL}|stage=deterministic|status=complete|seconds={time.monotonic()-t0:.3f}|sampled={sampled}|smooth={smooth}|rank={len(pivots)}", flush=True)

    rng = random.Random(args.seed)
    print(f"{PROTOCOL}|stage=random|status=start|trials={args.random_trials}", flush=True)
    for t in range(args.random_trials):
        scale = box if t % 5 else 5*box
        a,b,c = [rng.randint(-scale, scale) for _ in range(3)]
        if a == b == c == 0:
            continue
        process(a*B[0] + b*B[1] + c*B[2], f"rnd:{a},{b},{c}")
        if sampled - last_gain >= args.stall:
            print(f"{PROTOCOL}|stage=random|status=stalled|since_gain={sampled-last_gain}", flush=True)
            break

    kill = dict(pivots)
    for i in s_cols:
        insert_row(kill, 1 << i)

    print(
        f"{PROTOCOL}|stage=summary|fb_columns={len(fb)}|S_columns={len(s_cols)}"
        f"|relation_rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}"
        f"|S_projected_rank={len(s_pivots)}|after_killing_S={len(fb)-len(kill)}"
        f"|sampled={sampled}|smooth={smooth}|accepted={accepted}",
        flush=True,
    )
    print(f"{PROTOCOL}|warning=heuristic_factor_base_model|exact_relations=true|completeness_proved=false", flush=True)

if __name__ == "__main__":
    main()
