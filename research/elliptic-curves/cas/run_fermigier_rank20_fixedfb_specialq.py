#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict, deque
import json
from pathlib import Path
import time

PROTOCOL = "R20FIXEDQ"

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
    ap.add_argument("--factor-base-bound", type=int, default=1000)
    ap.add_argument("--special-q-min", type=int, default=5000)
    ap.add_argument("--special-q-max", type=int, default=50000)
    ap.add_argument("--b-bound", type=int, default=5000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 70)
    ap.add_argument("--max-special-q", type=int, default=250)
    ap.add_argument("--progress-every-q", type=int, default=5)
    ap.add_argument("--seed-specials", default="5023:3099,5237:5128")
    ap.add_argument(
        "--norm-factor-mode",
        choices=("trial", "hybrid", "exact"),
        default="trial",
        help=(
            "trial divides by every factor-base rational prime; hybrid first "
            "trial-divides then exactly factors only a bounded residual; exact "
            "factors the norm once"
        ),
    )
    ap.add_argument(
        "--trial-prime-bound",
        type=int,
        help=(
            "trial division limit in hybrid mode (default: min(factor-base "
            "bound, 10000)); trial mode always uses the full factor base"
        ),
    )
    ap.add_argument(
        "--residual-factor-limit",
        type=int,
        default=1 << 40,
        help="hybrid mode factors a residual only through this exact integer limit",
    )
    ap.add_argument("--relation-ledger", type=Path)
    args = ap.parse_args()
    if args.trial_prime_bound is not None and args.trial_prime_bound < 2:
        raise ValueError("--trial-prime-bound must be at least 2")
    if args.residual_factor_limit < 2:
        raise ValueError("--residual-factor-limit must be at least 2")

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, RealField, factor, prime_range

    x = PolynomialRing(QQ, "x").gen()
    A = ZZ(5750886029903523759416717668139307)
    C = ZZ(167347710468055045100164888198438918505621536951206)
    f = x**3 - A*x + C
    K = NumberField(f, "a")
    th = K.gen()

    bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))]))
    print(f"{PROTOCOL}|stage=input|fb_bound={args.factor_base_bound}|special_q_range={args.special_q_min}:{args.special_q_max}|b_bound={args.b_bound}|large_prime_bound={args.large_prime_bound}|norm_factor_mode={args.norm_factor_mode}|bad={bad}", flush=True)

    fb = []
    fb_index_by_hnf = {}
    s_cols = []

    def idx_of(P):
        return fb_index_by_hnf.get(str(P.pari_hnf()))

    def add_fb(P, is_s=False):
        i = idx_of(P)
        if i is None:
            i = len(fb)
            fb.append(P)
            fb_index_by_hnf[str(P.pari_hnf())] = i
        if is_s and i not in s_cols:
            s_cols.append(i)
        return i

    def fb_sig():
        return [
            {
                "hnf": str(P.pari_hnf()),
                "norm": int(P.norm()),
                "residue_degree": int(P.residue_class_degree()),
                "rational_prime": int(P.smallest_integer()),
            }
            for P in fb
        ]

    for p in bad:
        for P in K.primes_above(p):
            add_fb(P, True)

    fb_rational = set()
    for q in prime_range(2, args.factor_base_bound + 1):
        q = ZZ(q)
        fb_rational.add(q)
        if q in bad:
            continue
        for P in K.primes_above(q):
            add_fb(P, False)

    fb_by_q = {}
    for q in sorted(fb_rational | set(bad)):
        lst = []
        for P in K.primes_above(q):
            i = idx_of(P)
            if i is not None:
                lst.append((i, P))
        if lst:
            fb_by_q[q] = lst

    if args.norm_factor_mode == "trial":
        trial_division_bound = args.factor_base_bound
    elif args.trial_prime_bound is None:
        trial_division_bound = min(args.factor_base_bound, 10_000)
    else:
        trial_division_bound = min(args.factor_base_bound, args.trial_prime_bound)
    trial_primes = sorted(q for q in fb_by_q if q <= trial_division_bound)
    print(f"{PROTOCOL}|stage=factor_base|columns={len(fb)}|S_columns={len(s_cols)}|factor_base_rational_primes={len(fb_rational)}|trial_division_bound={trial_division_bound}|trial_primes={len(trial_primes)}", flush=True)

    RF = RealField(200)
    real_roots = sorted(r for r, mult in f.change_ring(RF).roots())
    print(f"{PROTOCOL}|stage=real_roots|roots={real_roots}", flush=True)

    seeded = []
    for item in args.seed_specials.split(","):
        item = item.strip()
        if item:
            q0, r0 = item.split(":")
            seeded.append((ZZ(q0), ZZ(r0)))

    specials = []
    seen = set()
    for qr in seeded:
        if qr not in seen:
            specials.append(qr); seen.add(qr)

    for q in prime_range(args.special_q_min, args.special_q_max + 1):
        q = ZZ(q)
        if q in bad:
            continue
        for r in f.change_ring(ZZ.quotient(q)).roots(multiplicities=False):
            qr = (q, ZZ(int(r)))
            if qr in seen:
                continue
            specials.append(qr); seen.add(qr)
            if len(specials) >= args.max_special_q:
                break
        if len(specials) >= args.max_special_q:
            break

    print(f"{PROTOCOL}|stage=special_q_setup|count={len(specials)}|first={specials[:10]}", flush=True)

    pivots = {}
    s_pivots = {}
    full_relations = 0
    cycle_relations = 0
    partial1 = 0
    partial2 = 0
    candidates = 0
    ROOT = ("ROOT",)
    graph = defaultdict(list)
    generators = []
    closed_relations = []
    t0 = time.monotonic()

    def after_killing_s():
        kill = dict(pivots)
        for i in s_cols:
            insert_row(kill, 1 << i)
        return len(fb) - len(kill)

    def add_full(row, source, meta, provenance):
        nonlocal full_relations, cycle_relations
        closed_relations.append(
            {
                "fb_parity_mask_hex": hex(row),
                "generator_indices": sorted(provenance),
                "kind": "unit_dependency" if row == 0 else "fb_relation",
                "source": source,
            }
        )
        if row == 0:
            return
        gained = insert_row(pivots, row)
        sr = project(row, s_cols)
        if sr:
            insert_row(s_pivots, sr)
        if source == "full":
            full_relations += 1
        else:
            cycle_relations += 1
        if gained:
            print(f"{PROTOCOL}|stage=relation|source={source}|rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}|after_killing_S={after_killing_s()}|meta={meta}", flush=True)

    def exact_fb_row(alpha, used_q):
        row = 0
        for q in used_q:
            for i, P in fb_by_q.get(q, ()):
                if int(alpha.valuation(P)) & 1:
                    row ^= 1 << i
        return row

    def prime_sig(alpha, ell):
        support = []
        for P in K.primes_above(ell):
            if int(alpha.valuation(P)) & 1:
                support.append((int(P.residue_class_degree()), str(P.pari_hnf())))
        if len(support) != 1:
            return None
        return (int(ell), support[0])

    def find_path(src, dst):
        if src == dst:
            return 0, set()
        dq = deque([src])
        parent = {src: None}
        edge_row = {}
        while dq:
            v = dq.popleft()
            for w, row, generator_index in graph[v]:
                if w in parent:
                    continue
                parent[w] = v
                edge_row[w] = (row, generator_index)
                if w == dst:
                    acc = 0
                    provenance = set()
                    cur = dst
                    while parent[cur] is not None:
                        edge, edge_generator = edge_row[cur]
                        acc ^= edge
                        if edge_generator in provenance:
                            provenance.remove(edge_generator)
                        else:
                            provenance.add(edge_generator)
                        cur = parent[cur]
                    return acc, provenance
                dq.append(w)
        return None, None

    def add_partial(v1, v2, row, meta, generator_index):
        path, provenance = find_path(v1, v2)
        if path is not None:
            if generator_index in provenance:
                provenance.remove(generator_index)
            else:
                provenance.add(generator_index)
            add_full(path ^ row, "partial_cycle", meta, provenance)
        graph[v1].append((v2, row, generator_index))
        graph[v2].append((v1, row, generator_index))

    def special_prime(q, r):
        for P in K.primes_above(q):
            if int(P.residue_class_degree()) == 1 and int((th - K(r)).valuation(P)) > 0:
                return P
        return None

    for sq_index, (q, r) in enumerate(specials, start=1):
        Pq = special_prime(q, r)
        if Pq is None:
            continue

        seen_ab = set()
        q_candidates = 0

        for b0 in range(-args.b_bound, args.b_bound + 1):
            if b0 == 0:
                continue
            b = ZZ(b0)
            residue = ZZ((-b * r) % q)

            for rho in real_roots:
                target = -RF(b) * rho
                k = ZZ(((target - RF(residue)) / RF(q)).round())
                a = residue + k*q
                key = (int(a), int(b))
                if key in seen_ab:
                    continue
                seen_ab.add(key)

                alpha = K(a) + K(b)*th
                vq = int(alpha.valuation(Pq))
                if vq <= 0:
                    continue

                N = abs(ZZ(alpha.norm()))
                if N == 0:
                    continue
                candidates += 1
                q_candidates += 1
                if args.norm_factor_mode in {"trial", "hybrid"}:
                    co = N
                    while co % q == 0:
                        co //= q
                    used_q = []
                    for p in trial_primes:
                        if co % p:
                            continue
                        used_q.append(p)
                        while co % p == 0:
                            co //= p
                        if co == 1:
                            break
                else:
                    co = ZZ(1)
                    used_q = []
                    for p, exponent in factor(N):
                        p = ZZ(p)
                        exponent = ZZ(exponent)
                        if p == q:
                            continue
                        if p <= args.factor_base_bound:
                            used_q.append(p)
                        else:
                            co *= p**exponent

                if args.norm_factor_mode == "hybrid" and co != 1 and co <= args.residual_factor_limit:
                    residual = ZZ(1)
                    for p, exponent in factor(co):
                        p = ZZ(p)
                        exponent = ZZ(exponent)
                        if p <= args.factor_base_bound:
                            used_q.append(p)
                        else:
                            residual *= p**exponent
                    co = residual

                row = exact_fb_row(alpha, used_q)
                generator_index = len(generators)
                generators.append(
                    {
                        "power_basis": [str(a), str(b), "0"],
                        "source_special_q": [int(q), int(r)],
                        "norm": str(N),
                    }
                )
                vertices = []
                if vq & 1:
                    vertices.append((int(q), (1, str(Pq.pari_hnf()))))

                if co != 1:
                    # Fast 1-large-prime mode. The forced special-q is already
                    # one outside graph vertex, so a single prime cofactor gives
                    # the useful two-vertex partial relation. Avoid factor(co):
                    # arbitrary cofactor factorization dominated the runtime.
                    if co > args.large_prime_bound or not co.is_prime(proof=False):
                        continue

                    sig = prime_sig(alpha, co)
                    if sig is None:
                        continue
                    vertices.append(sig)

                if len(vertices) == 0:
                    add_full(
                        row,
                        "full",
                        ("smooth", int(q), int(r), int(a), int(b)),
                        {generator_index},
                    )
                elif len(vertices) == 1:
                    partial1 += 1
                    add_partial(
                        ROOT, vertices[0], row,
                        ("1LP", int(q), int(r), int(a), int(b)), generator_index,
                    )
                elif len(vertices) == 2:
                    partial2 += 1
                    add_partial(
                        vertices[0], vertices[1], row,
                        ("2LP", int(q), int(r), int(a), int(b)), generator_index,
                    )

        if args.progress_every_q and sq_index % args.progress_every_q == 0:
            print(f"{PROTOCOL}|stage=progress|special_q_done={sq_index}/{len(specials)}|q={q}|r={r}|q_candidates={q_candidates}|rank={len(pivots)}|after_killing_S={after_killing_s()}|partial1={partial1}|partial2={partial2}|cycle_relations={cycle_relations}|graph_vertices={len(graph)}|seconds={time.monotonic()-t0:.3f}", flush=True)

    print(f"{PROTOCOL}|stage=summary|factor_base_columns={len(fb)}|S_columns={len(s_cols)}|relation_rank={len(pivots)}|fb_qdim={len(fb)-len(pivots)}|S_projected_rank={len(s_pivots)}|after_killing_S={after_killing_s()}|full_relations={full_relations}|cycle_relations={cycle_relations}|partial1={partial1}|partial2={partial2}|graph_vertices={len(graph)}|candidates={candidates}|seconds={time.monotonic()-t0:.3f}", flush=True)
    if args.relation_ledger:
        ledger = {
            "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
            "status": "exact_new_relations_not_class_group_completion",
            "field_polynomial": str(f),
            "defining_polynomial_ascending": [str(C), str(-A), "0", "1"],
            "field_discriminant": str(K.discriminant()),
            "generator_coordinate_order": ["1", "theta", "theta^2"],
            "factor_base_bound": args.factor_base_bound,
            "norm_factor_mode": args.norm_factor_mode,
            "trial_division_bound": trial_division_bound,
            "residual_factor_limit": args.residual_factor_limit,
            "factor_base_completion": {
                "all_prime_ideals_above_rational_primes_through": (
                    args.factor_base_bound
                ),
                "materialized_complete_factor_base": True,
                "extra_declared_S_rational_primes": [int(p) for p in bad],
            },
            "selmer_rational_primes": [int(p) for p in bad],
            "factor_base": fb_sig(),
            "S_columns": s_cols,
            "generators": generators,
            "closed_relations": closed_relations,
        }
        args.relation_ledger.parent.mkdir(parents=True, exist_ok=True)
        args.relation_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
        print(f"{PROTOCOL}|stage=write_relation_ledger|path={args.relation_ledger}|generators={len(generators)}|closed_relations={len(closed_relations)}", flush=True)
    print(f"{PROTOCOL}|warning=heuristic_relation_lattice|exact_relations=true|completeness_proved=false|fixed_factor_base=true|double_large_prime=true", flush=True)

if __name__ == "__main__":
    main()
