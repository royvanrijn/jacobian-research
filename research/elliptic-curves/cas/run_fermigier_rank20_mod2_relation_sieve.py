#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import time

PROTOCOL = "R20REL2SIEVE"

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
    ap.add_argument("--factor-base-bound", type=int, default=2000)
    ap.add_argument("--m-bound", type=int, default=1_000_000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 40)
    ap.add_argument("--progress-every", type=int, default=10000)
    args = ap.parse_args()

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    A = ZZ(5750886029903523759416717668139307)
    C = ZZ(167347710468055045100164888198438918505621536951206)
    f = x**3 - A*x + C

    K = NumberField(f, "a")
    th = K.gen()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))
    print(f"{PROTOCOL}|stage=input|fb_bound={args.factor_base_bound}|m_bound={args.m_bound}|large_prime_bound={args.large_prime_bound}|bad={bad}", flush=True)

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

    fb_rational = set()
    for q in prime_range(2, args.factor_base_bound + 1):
        q = ZZ(q)
        fb_rational.add(q)
        if q in bad:
            continue
        for P in K.primes_above(q):
            add(P, False)

    fb_by_q = {}
    for q in sorted(fb_rational | set(bad)):
        rows = []
        for P in K.primes_above(q):
            i = idx_of(P)
            if i is not None:
                rows.append((i, P))
        if rows:
            fb_by_q[q] = rows

    trial_primes = sorted(q for q in fb_by_q if q <= args.factor_base_bound)
    print(f"{PROTOCOL}|stage=factor_base|columns={len(fb)}|S_columns={len(s_cols)}|rational_primes={len(trial_primes)}", flush=True)

    pivots = {}
    s_pivots = {}
    partials = {}
    full = partial_seen = matched = 0
    t0 = time.monotonic()

    def add_relation(row, source, m):
        nonlocal full, matched
        if row == 0:
            return
        gained = insert_row(pivots, row)
        sr = project(row, s_cols)
        if sr:
            insert_row(s_pivots, sr)
        if source == "full":
            full += 1
        else:
            matched += 1
        if gained:
            kill = dict(pivots)
            for i in s_cols:
                insert_row(kill, 1 << i)
            print(f"{PROTOCOL}|stage=relation|source={source}|m={m}|rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}|after_killing_S={len(fb)-len(kill)}|partials_open={len(partials)}", flush=True)

    def exact_fb_row(alpha, used_q):
        row = 0
        for q in used_q:
            for i, P in fb_by_q.get(q, ()):
                v = int(alpha.valuation(P))
                if v & 1:
                    row ^= 1 << i
        return row

    def large_sig(alpha, ell):
        support = []
        for P in K.primes_above(ell):
            v = int(alpha.valuation(P))
            if v & 1:
                support.append((int(P.residue_class_degree()), str(P.pari_hnf())))
        if len(support) != 1:
            return None
        return (int(ell), support[0])

    total = 2 * args.m_bound + 1
    for n, m0 in enumerate(range(-args.m_bound, args.m_bound + 1), start=1):
        m = ZZ(m0)
        N = abs(m**3 - A*m + C)
        if N == 0:
            continue

        co = N
        used_q = []
        for q in trial_primes:
            if co % q:
                continue
            used_q.append(q)
            while co % q == 0:
                co //= q
            if co == 1:
                break

        if co == 1:
            alpha = K(m) - th
            add_relation(exact_fb_row(alpha, used_q), "full", m0)
        elif co <= args.large_prime_bound and co.is_prime(proof=False):
            alpha = K(m) - th
            row = exact_fb_row(alpha, used_q)
            sig = large_sig(alpha, co)
            if sig is not None:
                partial_seen += 1
                old = partials.pop(sig, None)
                if old is None:
                    partials[sig] = row
                else:
                    add_relation(old ^ row, "matched_large_prime", m0)

        if args.progress_every and n % args.progress_every == 0:
            kill = dict(pivots)
            for i in s_cols:
                insert_row(kill, 1 << i)
            print(f"{PROTOCOL}|stage=progress|processed={n}/{total}|m={m0}|rank={len(pivots)}|after_killing_S={len(fb)-len(kill)}|full={full}|partial_seen={partial_seen}|matched={matched}|partials_open={len(partials)}|seconds={time.monotonic()-t0:.3f}", flush=True)

    kill = dict(pivots)
    for i in s_cols:
        insert_row(kill, 1 << i)

    print(f"{PROTOCOL}|stage=summary|factor_base_columns={len(fb)}|S_columns={len(s_cols)}|relation_rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}|S_projected_rank={len(s_pivots)}|after_killing_S={len(fb)-len(kill)}|full_relations={full}|partial_seen={partial_seen}|matched_relations={matched}|partials_open={len(partials)}|seconds={time.monotonic()-t0:.3f}", flush=True)
    print(f"{PROTOCOL}|warning=heuristic_relation_lattice|exact_relations=true|completeness_proved=false|linear_family=m_minus_theta", flush=True)

if __name__ == "__main__":
    main()
