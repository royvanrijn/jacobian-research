#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, time
from collections import defaultdict, deque
from pathlib import Path

PROTOCOL = "R20FIXEDQ2"
SCHEMA = "r20-fixedfb-specialq-checkpoint-v1"

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
    ap.add_argument("--special-q-min", type=int, default=5000)
    ap.add_argument("--special-q-max", type=int, default=50000)
    ap.add_argument("--b-bound", type=int, default=5000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 70)
    ap.add_argument("--max-special-q", type=int, default=100)
    ap.add_argument("--progress-every-q", type=int, default=5)
    ap.add_argument("--seed-specials", default="5689:5096,5237:5128,5023:3099")
    ap.add_argument("--only-seeds", action="store_true")
    ap.add_argument("--checkpoint", type=Path, default=Path(
        "artifacts/local/elliptic-curves/r20_fixedfb5000_checkpoint.json"))
    args = ap.parse_args()

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, RealField, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    A = ZZ(5750886029903523759416717668139307)
    C = ZZ(167347710468055045100164888198438918505621536951206)
    f = x**3 - A*x + C
    K = NumberField(f, "a")
    th = K.gen()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))

    fb, s_cols = [], []
    def idx_of(P):
        for i, Q in enumerate(fb):
            if P == Q: return i
        return None
    def add_fb(P, is_s=False):
        i = idx_of(P)
        if i is None:
            i = len(fb); fb.append(P)
        if is_s and i not in s_cols: s_cols.append(i)
        return i

    for p in bad:
        for P in K.primes_above(p): add_fb(P, True)
    for q in prime_range(2, args.factor_base_bound + 1):
        q = ZZ(q)
        if q in bad: continue
        for P in K.primes_above(q): add_fb(P, False)

    fb_by_q = {}
    rationals = sorted(set(bad) | set(ZZ(q) for q in prime_range(2, args.factor_base_bound + 1)))
    for q in rationals:
        arr = []
        for P in K.primes_above(q):
            i = idx_of(P)
            if i is not None: arr.append((i, P))
        if arr: fb_by_q[q] = arr
    trial_primes = sorted(q for q in fb_by_q if q <= args.factor_base_bound)

    def fb_sig():
        return [{"hnf": str(P.pari_hnf()), "norm": int(P.norm()), "f": int(P.residue_class_degree())} for P in fb]

    pivots, s_pivots = {}, {}
    history = []
    if args.checkpoint.exists():
        cp = json.loads(args.checkpoint.read_text())
        if cp["schema"] != SCHEMA: raise ValueError("checkpoint schema mismatch")
        if cp["factor_base_bound"] != args.factor_base_bound: raise ValueError("factor-base mismatch")
        if cp["factor_base"] != fb_sig(): raise ValueError("factor-base structure mismatch")
        if cp["S_columns"] != s_cols: raise ValueError("S-column mismatch")
        for h in cp.get("basis_rows_hex", []): insert_row(pivots, int(h, 16))
        history = cp.get("history", [])
        print(f"{PROTOCOL}|stage=checkpoint_load|rank={len(pivots)}|history={len(history)}", flush=True)

    for row in rows_from_pivots(pivots):
        sr = project(row, s_cols)
        if sr: insert_row(s_pivots, sr)

    def afterS():
        q = dict(pivots)
        for i in s_cols: insert_row(q, 1 << i)
        return len(fb) - len(q)

    def save():
        args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema": SCHEMA,
            "factor_base_bound": args.factor_base_bound,
            "factor_base": fb_sig(),
            "S_columns": s_cols,
            "basis_rows_hex": [hex(r) for r in rows_from_pivots(pivots)],
            "relation_rank": len(pivots),
            "after_killing_S": afterS(),
            "history": history,
        }
        tmp = args.checkpoint.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
        tmp.replace(args.checkpoint)

    RF = RealField(200)
    real_roots = sorted(r for r, mult in f.change_ring(RF).roots())

    seeds = []
    for item in args.seed_specials.split(","):
        if item.strip():
            q, r = item.strip().split(":")
            seeds.append((ZZ(q), ZZ(r)))

    specials, seen = [], set()
    for qr in seeds:
        if qr not in seen: specials.append(qr); seen.add(qr)

    if not args.only_seeds:
        for q in prime_range(args.special_q_min, args.special_q_max + 1):
            q = ZZ(q)
            if q in bad: continue
            for r in f.change_ring(ZZ.quotient(q)).roots(multiplicities=False):
                qr = (q, ZZ(int(r)))
                if qr in seen: continue
                specials.append(qr); seen.add(qr)
                if len(specials) >= args.max_special_q: break
            if len(specials) >= args.max_special_q: break

    print(f"{PROTOCOL}|stage=input|fb_columns={len(fb)}|S_columns={len(s_cols)}|start_rank={len(pivots)}|start_afterS={afterS()}|specials={len(specials)}", flush=True)

    ROOT = ("ROOT",)
    graph = defaultdict(list)
    t0 = time.monotonic()

    def exact_fb_row(alpha, used_q):
        row = 0
        for q in used_q:
            for i, P in fb_by_q.get(q, ()):
                if int(alpha.valuation(P)) & 1: row ^= 1 << i
        return row

    def prime_sig(alpha, ell):
        out = []
        for P in K.primes_above(ell):
            if int(alpha.valuation(P)) & 1:
                out.append((int(P.residue_class_degree()), str(P.pari_hnf())))
        return None if len(out) != 1 else (int(ell), out[0])

    def special_prime(q, r):
        for P in K.primes_above(q):
            if int(P.residue_class_degree()) == 1 and int((th - K(r)).valuation(P)) > 0:
                return P
        return None

    def find_path(src, dst):
        if src == dst: return 0
        dq = deque([src]); parent = {src: None}; edge = {}
        while dq:
            v = dq.popleft()
            for w, row in graph[v]:
                if w in parent: continue
                parent[w] = v; edge[w] = row
                if w == dst:
                    acc = 0; cur = dst
                    while parent[cur] is not None:
                        acc ^= edge[cur]; cur = parent[cur]
                    return acc
                dq.append(w)
        return None

    def add_partial(v1, v2, row):
        path = find_path(v1, v2)
        graph[v1].append((v2, row)); graph[v2].append((v1, row))
        return None if path is None else (path ^ row)

    for idx, (q, r) in enumerate(specials, 1):
        Pq = special_prime(q, r)
        if Pq is None: continue

        before_rank, before_afterS = len(pivots), afterS()
        candidates = p1 = p2 = cycles = indep = after_gain = 0
        seen_ab = set()

        for b0 in range(-args.b_bound, args.b_bound + 1):
            if b0 == 0: continue
            b = ZZ(b0); residue = ZZ((-b*r) % q)
            for rho in real_roots:
                k = ZZ(((-RF(b)*rho - RF(residue))/RF(q)).round())
                a = residue + k*q
                key = (int(a), int(b))
                if key in seen_ab: continue
                seen_ab.add(key)

                alpha = K(a) + K(b)*th
                vq = int(alpha.valuation(Pq))
                if vq <= 0: continue
                N = abs(ZZ(alpha.norm()))
                if N == 0: continue
                candidates += 1

                co = N
                while co % q == 0: co //= q

                used = []
                for p in trial_primes:
                    if co % p: continue
                    used.append(p)
                    while co % p == 0: co //= p
                    if co == 1: break

                row = exact_fb_row(alpha, used)
                vertices = []
                if vq & 1: vertices.append((int(q), (1, str(Pq.pari_hnf()))))

                if co != 1:
                    if co > args.large_prime_bound or not co.is_prime(proof=False): continue
                    sig = prime_sig(alpha, co)
                    if sig is None: continue
                    vertices.append(sig)

                cyc = None
                if len(vertices) == 0:
                    cyc = row
                elif len(vertices) == 1:
                    p1 += 1; cyc = add_partial(ROOT, vertices[0], row)
                elif len(vertices) == 2:
                    p2 += 1; cyc = add_partial(vertices[0], vertices[1], row)
                else:
                    continue

                if cyc:
                    cycles += 1
                    old_after = afterS()
                    if insert_row(pivots, cyc):
                        indep += 1
                        sr = project(cyc, s_cols)
                        if sr: insert_row(s_pivots, sr)
                        if afterS() < old_after: after_gain += 1

        st = {
            "q": int(q), "r": int(r), "b_bound": args.b_bound,
            "candidates": candidates, "partial1": p1, "partial2": p2,
            "cycles": cycles, "rank_gain": len(pivots)-before_rank,
            "afterS_gain": before_afterS-afterS(), "independent": indep,
            "independent_afterS": after_gain,
        }
        history.append(st); save()

        print(f"{PROTOCOL}|stage=special_q_done|index={idx}/{len(specials)}|q={q}|r={r}|rank_gain={st['rank_gain']}|afterS_gain={st['afterS_gain']}|rank={len(pivots)}|afterS={afterS()}|candidates={candidates}|p1={p1}|p2={p2}|cycles={cycles}", flush=True)

        if args.progress_every_q and idx % args.progress_every_q == 0:
            print(f"{PROTOCOL}|stage=progress|done={idx}/{len(specials)}|rank={len(pivots)}|afterS={afterS()}|seconds={time.monotonic()-t0:.3f}", flush=True)

    agg = {}
    for st in history:
        key = (st["q"], st["r"])
        a = agg.setdefault(key, {"q":st["q"],"r":st["r"],"runs":0,"candidates":0,"rank_gain":0,"afterS_gain":0,"cycles":0})
        a["runs"] += 1; a["candidates"] += st["candidates"]; a["rank_gain"] += st["rank_gain"]; a["afterS_gain"] += st["afterS_gain"]; a["cycles"] += st["cycles"]

    ranked = sorted(agg.values(), key=lambda a:(a["afterS_gain"],a["rank_gain"],a["cycles"]), reverse=True)
    print(f"{PROTOCOL}|stage=summary|final_rank={len(pivots)}|afterS={afterS()}|checkpoint={args.checkpoint}", flush=True)
    print(f"{PROTOCOL}|stage=recommended_seeds|top={ranked[:15]}", flush=True)
    print(f"{PROTOCOL}|warning=heuristic_relation_lattice|exact_relations=true|completeness_proved=false|checkpointed=true", flush=True)

if __name__ == "__main__":
    main()
