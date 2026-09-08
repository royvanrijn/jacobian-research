#!/usr/bin/env python3
"""Special-q lattice sieve for mod-2 ideal-class relations in the Fermigier cubic field.

We use alpha = a + b*theta and force divisibility by a chosen degree-1 prime
ideal q = (q, theta-r) via
    a + b*r == 0 (mod q).

For each special q:
  * enumerate small b and choose the nearest a satisfying the congruence;
  * compute N = |Norm(alpha)| and divide out q-adic rational contribution;
  * trial-divide by the rational factor base;
  * retain full relations and one-large-prime partials;
  * verify exact prime-ideal valuation parity before adding any relation.

No bnfinit/class_group/regulator computation is used.

Run:
  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_specialq_sieve.py \
    --factor-base-bound 5000 \
    --special-q-min 5000 \
    --special-q-max 100000 \
    --b-bound 1000 \
    --large-prime-bound 1125899906842624 \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_specialq_sieve.log
"""

from __future__ import annotations

import argparse
from pathlib import Path
import time

PROTOCOL = "R20SPECIALQ"


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
    ap.add_argument("--factor-base-bound", type=int, default=5000)
    ap.add_argument("--special-q-min", type=int, default=5000)
    ap.add_argument("--special-q-max", type=int, default=100000)
    ap.add_argument("--b-bound", type=int, default=1000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 50)
    ap.add_argument("--max-special-q", type=int, default=200)
    ap.add_argument("--progress-every-q", type=int, default=10)
    args = ap.parse_args()

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, RealField, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    A = ZZ(5750886029903523759416717668139307)
    C = ZZ(167347710468055045100164888198438918505621536951206)
    f = x**3 - A*x + C

    K = NumberField(f, "a")
    th = K.gen()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))

    print(
        f"{PROTOCOL}|stage=input|fb_bound={args.factor_base_bound}"
        f"|special_q_range={args.special_q_min}:{args.special_q_max}"
        f"|b_bound={args.b_bound}|large_prime_bound={args.large_prime_bound}"
        f"|bad={bad}",
        flush=True,
    )

    # Factor base.
    fb = []
    s_cols = []

    def idx_of(P):
        for i, Q in enumerate(fb):
            if P == Q:
                return i
        return None

    def add(P, is_s=False):
        i = idx_of(P)
        if i is None:
            i = len(fb)
            fb.append(P)
        if is_s and i not in s_cols:
            s_cols.append(i)
        return i

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
        lst = []
        for P in K.primes_above(q):
            i = idx_of(P)
            if i is not None:
                lst.append((i, P))
        if lst:
            fb_by_q[q] = lst

    trial_primes = sorted(q for q in fb_by_q if q <= args.factor_base_bound)

    print(
        f"{PROTOCOL}|stage=factor_base|columns={len(fb)}|S_columns={len(s_cols)}"
        f"|rational_primes={len(trial_primes)}",
        flush=True,
    )

    # Degree-1 special-q candidates.
    specials = []
    for q in prime_range(args.special_q_min, args.special_q_max + 1):
        q = ZZ(q)
        if q in bad:
            continue
        # roots of f mod q correspond to degree-1 prime ideals.
        roots = f.change_ring(ZZ.quotient(q)).roots(multiplicities=False)
        for r in roots:
            r = ZZ(int(r))
            specials.append((q, r))
            if len(specials) >= args.max_special_q:
                break
        if len(specials) >= args.max_special_q:
            break

    print(
        f"{PROTOCOL}|stage=special_q_setup|count={len(specials)}"
        f"|first={specials[:10]}",
        flush=True,
    )

    # Real roots used to reduce the special-q lattice.
    # For alpha=a+b*theta choose a in its congruence class mod q so that
    # a+b*rho is tiny at one real embedding rho.
    RF = RealField(200)
    real_roots = sorted(r for r, mult in f.change_ring(RF).roots())
    print(
        f"{PROTOCOL}|stage=real_roots|roots={real_roots}",
        flush=True,
    )

    pivots = {}
    s_pivots = {}
    partials = {}
    full = partial_seen = matched = 0
    candidates = 0
    t0 = time.monotonic()

    def add_relation(row, source, q, r, a, b):
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
            print(
                f"{PROTOCOL}|stage=relation|source={source}|q={q}|r={r}|a={a}|b={b}"
                f"|rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}"
                f"|after_killing_S={len(fb)-len(kill)}|partials_open={len(partials)}",
                flush=True,
            )

    def exact_fb_row(alpha, used_q):
        row = 0
        for q in used_q:
            for i, P in fb_by_q.get(q, ()):
                v = int(alpha.valuation(P))
                if v & 1:
                    row ^= 1 << i
        return row

    def special_prime_ideal(q, r):
        # Find the degree-1 prime above q corresponding to theta = r mod q.
        for P in K.primes_above(q):
            if int(P.residue_class_degree()) != 1:
                continue
            if int((th - K(r)).valuation(P)) > 0:
                return P
        return None

    def large_sig(alpha, ell):
        support = []
        for P in K.primes_above(ell):
            v = int(alpha.valuation(P))
            if v & 1:
                support.append((int(P.residue_class_degree()), str(P.pari_hnf())))
        if len(support) != 1:
            return None
        return (int(ell), support[0])

    for sq_index, (q, r) in enumerate(specials, start=1):
        Pq = special_prime_ideal(q, r)
        if Pq is None:
            print(f"{PROTOCOL}|stage=special_q|q={q}|r={r}|status=no_matching_prime", flush=True)
            continue

        q_candidates = 0
        q_full0, q_partial0, q_matched0 = full, partial_seen, matched

        # For each b, choose a congruent to -b*r mod q and centered near 0.
        seen_ab = set()

        for b0 in range(-args.b_bound, args.b_bound + 1):
            if b0 == 0:
                continue
            b = ZZ(b0)

            # All solutions are a = residue + k*q.
            residue = ZZ((-b * r) % q)

            # Instead of taking the small centered representative, choose k
            # separately for each real embedding so a ~= -b*rho.
            # This makes one conjugate of alpha=a+b*theta very small.
            for rho in real_roots:
                target = -RF(b) * rho
                k = ZZ(((target - RF(residue)) / RF(q)).round())
                a = residue + k*q

                key = (int(a), int(b))
                if key in seen_ab:
                    continue
                seen_ab.add(key)

                alpha = K(a) + K(b) * th
                vq = int(alpha.valuation(Pq))
                if vq <= 0:
                    continue

                N = abs(ZZ(alpha.norm()))
                if N == 0:
                    continue
                candidates += 1
                q_candidates += 1

                # Remove the forced rational q contribution at least once.
                co = N
                used_q = []

                # Remove the forced special-q rational factor FIRST.
                # q may lie above the ordinary factor-base bound.
                if co % q != 0:
                    raise ArithmeticError(
                        f"special-q divisibility failed: q={q}, a={a}, b={b}"
                    )

                while co % q == 0:
                    co //= q
                used_q.append(q)

                # Now remove the ordinary rational factor-base primes.
                for p in trial_primes:
                    if p == q or co % p:
                        continue
                    used_q.append(p)
                    while co % p == 0:
                        co //= p
                    if co == 1:
                        break

                # Include special q in exact valuation reconstruction even if q > FB bound.
                if q not in fb_by_q:
                    # Add all prime ideals above this q lazily to the factor base.
                    lst = []
                    for P in K.primes_above(q):
                        i = add(P, False)
                        lst.append((i, P))
                    fb_by_q[q] = lst
                if q not in used_q:
                    used_q.append(q)

                if co == 1:
                    row = exact_fb_row(alpha, used_q)
                    add_relation(row, "full", q, r, a, b)
                elif co <= args.large_prime_bound and co.is_prime(proof=False):
                    row = exact_fb_row(alpha, used_q)
                    sig = large_sig(alpha, co)
                    if sig is not None:
                        partial_seen += 1
                        old = partials.pop(sig, None)
                        if old is None:
                            partials[sig] = row
                        else:
                            add_relation(old ^ row, "matched_large_prime", q, r, a, b)

        if args.progress_every_q and sq_index % args.progress_every_q == 0:
            kill = dict(pivots)
            for i in s_cols:
                insert_row(kill, 1 << i)
            print(
                f"{PROTOCOL}|stage=progress|special_q_done={sq_index}/{len(specials)}"
                f"|q={q}|r={r}|q_candidates={q_candidates}"
                f"|rank={len(pivots)}|after_killing_S={len(fb)-len(kill)}"
                f"|full_delta={full-q_full0}|partial_delta={partial_seen-q_partial0}"
                f"|matched_delta={matched-q_matched0}"
                f"|partials_open={len(partials)}|seconds={time.monotonic()-t0:.3f}",
                flush=True,
            )

    kill = dict(pivots)
    for i in s_cols:
        insert_row(kill, 1 << i)

    print(
        f"{PROTOCOL}|stage=summary|factor_base_columns={len(fb)}"
        f"|S_columns={len(s_cols)}|relation_rank={len(pivots)}"
        f"|fb_qdim={len(fb)-len(pivots)}"
        f"|S_projected_rank={len(s_pivots)}"
        f"|after_killing_S={len(fb)-len(kill)}"
        f"|full_relations={full}|partial_seen={partial_seen}"
        f"|matched_relations={matched}|partials_open={len(partials)}"
        f"|candidates={candidates}|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=heuristic_relation_lattice"
        f"|exact_relations=true|completeness_proved=false|family=a_plus_b_theta",
        flush=True,
    )


if __name__ == "__main__":
    main()
