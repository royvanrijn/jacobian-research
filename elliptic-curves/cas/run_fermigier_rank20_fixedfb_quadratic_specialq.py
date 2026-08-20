#!/usr/bin/env python3
"""Quadratic-in-theta special-q collector for the Fermigier cubic field.

Search family:
    alpha = a + b*theta + c*theta^2

Special-q condition at degree-1 prime (q, theta-r):
    a + b*r + c*r^2 == 0 mod q.

For each sampled (b,c), and for each real embedding rho, choose the integer
a in that congruence class which makes
    a + b*rho + c*rho^2
as small as possible.

This is deliberately a different norm geometry from the saturated a+b*theta
family.  It reuses the SAME fixed factor base and SAME checkpointed GF(2)
relation basis from R20FIXEDQ2.

No bnfinit/class_group/regulator computation is used.

Example:
  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_fixedfb_quadratic_specialq.py \
    --factor-base-bound 5000 \
    --special-q-min 50000 \
    --special-q-max 500000 \
    --pairs-per-q 20000 \
    --coeff-bound 5000 \
    --max-special-q 100 \
    --checkpoint artifacts/local/elliptic-curves/r20_fixedfb5000_checkpoint.json \
    2>&1 | tee artifacts/local/elliptic-curves/r20_fixedfb_quadratic_specialq.log
"""

from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict, deque
from pathlib import Path
import time

PROTOCOL = "R20QUADQ"
CHECKPOINT_SCHEMA = "r20-fixedfb-specialq-checkpoint-v1"


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


def rows_from_pivots(pivots):
    return [pivots[k] for k in sorted(pivots, reverse=True)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--factor-base-bound", type=int, default=5000)
    ap.add_argument("--special-q-min", type=int, default=50000)
    ap.add_argument("--special-q-max", type=int, default=500000)
    ap.add_argument("--coeff-bound", type=int, default=5000)
    ap.add_argument("--pairs-per-q", type=int, default=20000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 70)
    ap.add_argument("--max-special-q", type=int, default=100)
    ap.add_argument("--progress-every-q", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260820)
    ap.add_argument("--seed-specials", default="")
    ap.add_argument("--only-seeds", action="store_true")
    ap.add_argument("--checkpoint", type=Path, default=Path(
        "artifacts/local/elliptic-curves/r20_fixedfb5000_checkpoint.json"))
    args = ap.parse_args()

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, RealField, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    A = ZZ(5750886029903523759416717668139307)
    C0 = ZZ(167347710468055045100164888198438918505621536951206)
    f = x**3 - A*x + C0
    K = NumberField(f, "a")
    th = K.gen()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))

    # ---- identical fixed factor base to R20FIXEDQ2 -------------------------
    fb, s_cols = [], []

    def idx_of(P):
        for i, Q in enumerate(fb):
            if P == Q:
                return i
        return None

    def add_fb(P, is_s=False):
        i = idx_of(P)
        if i is None:
            i = len(fb)
            fb.append(P)
        if is_s and i not in s_cols:
            s_cols.append(i)
        return i

    for p in bad:
        for P in K.primes_above(p):
            add_fb(P, True)

    for q in prime_range(2, args.factor_base_bound + 1):
        q = ZZ(q)
        if q in bad:
            continue
        for P in K.primes_above(q):
            add_fb(P, False)

    fb_by_q = {}
    rationals = sorted(
        set(bad) |
        set(ZZ(q) for q in prime_range(2, args.factor_base_bound + 1))
    )
    for q in rationals:
        arr = []
        for P in K.primes_above(q):
            i = idx_of(P)
            if i is not None:
                arr.append((i, P))
        if arr:
            fb_by_q[q] = arr

    trial_primes = sorted(q for q in fb_by_q if q <= args.factor_base_bound)

    def fb_sig():
        return [
            {
                "hnf": str(P.pari_hnf()),
                "norm": int(P.norm()),
                "f": int(P.residue_class_degree()),
            }
            for P in fb
        ]

    # ---- load existing a+b*theta checkpoint -------------------------------
    pivots, s_pivots = {}, {}
    history = []

    if args.checkpoint.exists():
        cp = json.loads(args.checkpoint.read_text())
        if cp.get("schema") != CHECKPOINT_SCHEMA:
            raise ValueError("checkpoint schema mismatch")
        if int(cp["factor_base_bound"]) != args.factor_base_bound:
            raise ValueError("factor-base bound mismatch")
        if cp["factor_base"] != fb_sig():
            raise ValueError("factor-base structure mismatch")
        if cp["S_columns"] != s_cols:
            raise ValueError("S columns mismatch")
        for h in cp.get("basis_rows_hex", []):
            insert_row(pivots, int(h, 16))
        history = list(cp.get("history", []))
        print(
            f"{PROTOCOL}|stage=checkpoint_load|rank={len(pivots)}"
            f"|history={len(history)}",
            flush=True,
        )

    for row in rows_from_pivots(pivots):
        sr = project(row, s_cols)
        if sr:
            insert_row(s_pivots, sr)

    def afterS():
        qbasis = dict(pivots)
        for i in s_cols:
            insert_row(qbasis, 1 << i)
        return len(fb) - len(qbasis)

    def save_checkpoint():
        args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": CHECKPOINT_SCHEMA,
            "factor_base_bound": args.factor_base_bound,
            "factor_base": fb_sig(),
            "S_columns": s_cols,
            "basis_rows_hex": [hex(r) for r in rows_from_pivots(pivots)],
            "relation_rank": len(pivots),
            "after_killing_S": afterS(),
            "history": history,
        }
        tmp = args.checkpoint.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        tmp.replace(args.checkpoint)

    # ---- real roots / specials ---------------------------------------------
    RF = RealField(220)
    real_roots = sorted(r for r, mult in f.change_ring(RF).roots())

    specials, seen = [], set()
    for item in args.seed_specials.split(","):
        item = item.strip()
        if not item:
            continue
        q0, r0 = item.split(":")
        qr = (ZZ(q0), ZZ(r0))
        specials.append(qr)
        seen.add(qr)

    if not args.only_seeds:
        for q in prime_range(args.special_q_min, args.special_q_max + 1):
            q = ZZ(q)
            if q in bad:
                continue
            for r in f.change_ring(ZZ.quotient(q)).roots(multiplicities=False):
                qr = (q, ZZ(int(r)))
                if qr in seen:
                    continue
                specials.append(qr)
                seen.add(qr)
                if len(specials) >= args.max_special_q:
                    break
            if len(specials) >= args.max_special_q:
                break

    print(
        f"{PROTOCOL}|stage=input|fb_columns={len(fb)}|S_columns={len(s_cols)}"
        f"|start_rank={len(pivots)}|start_afterS={afterS()}"
        f"|specials={len(specials)}|pairs_per_q={args.pairs_per_q}"
        f"|coeff_bound={args.coeff_bound}",
        flush=True,
    )

    # ---- shared partial relation graph -------------------------------------
    ROOT = ("ROOT",)
    graph = defaultdict(list)
    rng = random.Random(args.seed)
    t0 = time.monotonic()

    def exact_fb_row(alpha, used_q):
        row = 0
        for q in used_q:
            for i, P in fb_by_q.get(q, ()):
                if int(alpha.valuation(P)) & 1:
                    row ^= 1 << i
        return row

    def prime_sig(alpha, ell):
        out = []
        for P in K.primes_above(ell):
            if int(alpha.valuation(P)) & 1:
                out.append(
                    (int(P.residue_class_degree()), str(P.pari_hnf()))
                )
        if len(out) != 1:
            return None
        return (int(ell), out[0])

    def special_prime(q, r):
        for P in K.primes_above(q):
            if (
                int(P.residue_class_degree()) == 1
                and int((th - K(r)).valuation(P)) > 0
            ):
                return P
        return None

    def find_path(src, dst):
        if src == dst:
            return 0
        dq = deque([src])
        parent = {src: None}
        edge = {}
        while dq:
            v = dq.popleft()
            for w, row in graph[v]:
                if w in parent:
                    continue
                parent[w] = v
                edge[w] = row
                if w == dst:
                    acc = 0
                    cur = dst
                    while parent[cur] is not None:
                        acc ^= edge[cur]
                        cur = parent[cur]
                    return acc
                dq.append(w)
        return None

    def add_partial(v1, v2, row):
        path = find_path(v1, v2)
        graph[v1].append((v2, row))
        graph[v2].append((v1, row))
        if path is None:
            return None
        return path ^ row

    for sq_index, (q, r) in enumerate(specials, 1):
        Pq = special_prime(q, r)
        if Pq is None:
            continue

        before_rank = len(pivots)
        before_after = afterS()

        candidates = p1 = p2 = cycles = indep = after_gain = 0
        accepted_residual_prime = 0
        norm_bits_min = None
        norm_bits_sum = 0
        seen_abc = set()

        # Random (b,c), plus some structured sparse pairs.
        pairs = []
        B = args.coeff_bound

        # deterministic sparse directions first
        sparse = min(500, args.pairs_per_q // 4)
        for t in range(1, sparse + 1):
            v = 1 + (t % B)
            pairs.extend([(v, 1), (-v, 1), (1, v), (1, -v)])

        # random remainder
        while len(pairs) < args.pairs_per_q:
            b = rng.randint(-B, B)
            c = rng.randint(-B, B)
            if b == 0 and c == 0:
                continue
            pairs.append((b, c))

        pairs = pairs[:args.pairs_per_q]

        for b0, c0 in pairs:
            b = ZZ(b0)
            c = ZZ(c0)

            # congruence: a == -(b*r + c*r^2) mod q
            residue = ZZ((-(b*r + c*r*r)) % q)

            for rho in real_roots:
                target = -(RF(b)*rho + RF(c)*rho*rho)
                k = ZZ(((target - RF(residue))/RF(q)).round())
                a = residue + k*q

                key = (int(a), int(b), int(c))
                if key in seen_abc:
                    continue
                seen_abc.add(key)

                alpha = K(a) + K(b)*th + K(c)*(th**2)
                vq = int(alpha.valuation(Pq))
                if vq <= 0:
                    continue

                N = abs(ZZ(alpha.norm()))
                if N == 0:
                    continue

                candidates += 1
                nb = int(N.nbits())
                norm_bits_sum += nb
                norm_bits_min = nb if norm_bits_min is None else min(norm_bits_min, nb)

                co = N
                while co % q == 0:
                    co //= q

                used = []
                for p in trial_primes:
                    if co % p:
                        continue
                    used.append(p)
                    while co % p == 0:
                        co //= p
                    if co == 1:
                        break

                row = exact_fb_row(alpha, used)
                vertices = []

                if vq & 1:
                    vertices.append(
                        (int(q), (1, str(Pq.pari_hnf())))
                    )

                # Same fast mode as successful linear-family checkpoint:
                # one prime residual cofactor at most.
                if co != 1:
                    if (
                        co > args.large_prime_bound
                        or not co.is_prime(proof=False)
                    ):
                        continue
                    accepted_residual_prime += 1
                    sig = prime_sig(alpha, co)
                    if sig is None:
                        continue
                    vertices.append(sig)

                cycle = None
                if len(vertices) == 0:
                    cycle = row
                elif len(vertices) == 1:
                    p1 += 1
                    cycle = add_partial(ROOT, vertices[0], row)
                elif len(vertices) == 2:
                    p2 += 1
                    cycle = add_partial(vertices[0], vertices[1], row)
                else:
                    continue

                if cycle:
                    cycles += 1
                    old_after = afterS()
                    if insert_row(pivots, cycle):
                        indep += 1
                        sr = project(cycle, s_cols)
                        if sr:
                            insert_row(s_pivots, sr)
                        if afterS() < old_after:
                            after_gain += 1

        stat = {
            "family": "a+b*theta+c*theta^2",
            "q": int(q),
            "r": int(r),
            "coeff_bound": args.coeff_bound,
            "pairs_per_q": args.pairs_per_q,
            "candidates": candidates,
            "partial1": p1,
            "partial2": p2,
            "cycles": cycles,
            "rank_gain": len(pivots) - before_rank,
            "afterS_gain": before_after - afterS(),
            "independent": indep,
            "independent_afterS": after_gain,
            "min_norm_bits": norm_bits_min,
            "avg_norm_bits": (
                norm_bits_sum / candidates if candidates else None
            ),
        }
        history.append(stat)
        save_checkpoint()

        print(
            f"{PROTOCOL}|stage=special_q_done|index={sq_index}/{len(specials)}"
            f"|q={q}|r={r}|rank_gain={stat['rank_gain']}"
            f"|afterS_gain={stat['afterS_gain']}"
            f"|rank={len(pivots)}|afterS={afterS()}"
            f"|candidates={candidates}|p1={p1}|p2={p2}|cycles={cycles}"
            f"|min_norm_bits={norm_bits_min}"
            f"|avg_norm_bits={stat['avg_norm_bits']}",
            flush=True,
        )

        if args.progress_every_q and sq_index % args.progress_every_q == 0:
            print(
                f"{PROTOCOL}|stage=progress|done={sq_index}/{len(specials)}"
                f"|rank={len(pivots)}|afterS={afterS()}"
                f"|seconds={time.monotonic()-t0:.3f}",
                flush=True,
            )

    # Rank just the quadratic-family runs for recommendations.
    qhist = [
        st for st in history
        if st.get("family") == "a+b*theta+c*theta^2"
    ]
    ranked = sorted(
        qhist,
        key=lambda st: (
            st["afterS_gain"],
            st["rank_gain"],
            st["cycles"],
        ),
        reverse=True,
    )

    print(
        f"{PROTOCOL}|stage=summary|final_rank={len(pivots)}"
        f"|afterS={afterS()}|checkpoint={args.checkpoint}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=recommended_seeds|top={ranked[:15]}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=heuristic_relation_lattice"
        f"|exact_relations=true|completeness_proved=false"
        f"|family=a_plus_btheta_plus_ctheta2|checkpointed=true",
        flush=True,
    )


if __name__ == "__main__":
    main()
