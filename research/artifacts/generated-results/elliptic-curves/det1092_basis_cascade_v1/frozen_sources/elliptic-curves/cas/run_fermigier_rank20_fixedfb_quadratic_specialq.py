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
family.  By default it reuses the SAME Fermigier factor base and checkpointed
GF(2) relation basis from R20FIXEDQ2.  Explicit cubic/S-prime arguments make
the exact collector reusable for another fixed no-rational-2-torsion curve;
use a target-specific checkpoint in that mode.

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


class SparseLargePrimeEliminator:
    """Eliminate arbitrary large-prime hyperedges with exact witnesses.

    Each input carries a factor-base row and a set of large-prime ideal
    vertices.  A dependency in the large-prime columns returns the resulting
    factor-base relation and the symmetric-difference set of principal
    generators which proves it.
    """

    def __init__(self):
        self.vertex_columns = {}
        self.pivots = {}
        self.edge_count = 0
        self.dependency_count = 0

    def _mask(self, vertices):
        mask = 0
        parity = set()
        for vertex in vertices:
            if vertex in parity:
                parity.remove(vertex)
            else:
                parity.add(vertex)
        for vertex in sorted(parity, key=repr):
            column = self.vertex_columns.setdefault(vertex, len(self.vertex_columns))
            mask ^= 1 << column
        return mask

    def add(self, vertices, row, generator_index):
        self.edge_count += 1
        mask = self._mask(vertices)
        provenance = {generator_index}
        while mask:
            pivot = mask.bit_length() - 1
            previous = self.pivots.get(pivot)
            if previous is None:
                self.pivots[pivot] = (mask, row, provenance)
                return None, None
            previous_mask, previous_row, previous_provenance = previous
            mask ^= previous_mask
            row ^= previous_row
            provenance.symmetric_difference_update(previous_provenance)
        self.dependency_count += 1
        return row, provenance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--factor-base-bound", type=int, default=5000)
    ap.add_argument(
        "--field-polynomial-ascending",
        help="comma-separated integral coefficients c0,c1,c2,1 for a monic cubic",
    )
    ap.add_argument(
        "--selmer-rational-primes",
        help="comma-separated rational primes whose ideals form S",
    )
    ap.add_argument("--special-q-min", type=int, default=50000)
    ap.add_argument("--special-q-max", type=int, default=500000)
    ap.add_argument("--coeff-bound", type=int, default=5000)
    ap.add_argument("--pairs-per-q", type=int, default=20000)
    ap.add_argument("--large-prime-bound", type=int, default=1 << 70)
    ap.add_argument(
        "--large-prime-merge-engine",
        choices=("graph", "sparse-hypergraph"),
        default="graph",
    )
    ap.add_argument(
        "--max-residual-primes",
        type=int,
        default=1,
        help="maximum odd outside-factor primes retained in exact-factor mode",
    )
    ap.add_argument(
        "--adaptive-residual-specials",
        type=int,
        default=0,
        help="maximum degree-one residual ideals appended as new special ideals",
    )
    ap.add_argument(
        "--adaptive-residual-prime-bound",
        type=int,
        default=1 << 40,
    )
    ap.add_argument(
        "--adaptive-residual-depth",
        type=int,
        default=1,
        help="maximum number of residual-special generations after the initial list",
    )
    ap.add_argument("--max-special-q", type=int, default=100)
    ap.add_argument("--progress-every-q", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260820)
    ap.add_argument("--seed-specials", default="")
    ap.add_argument("--only-seeds", action="store_true")
    ap.add_argument(
        "--norm-factor-mode",
        choices=("trial", "hybrid", "exact"),
        default="trial",
    )
    ap.add_argument("--trial-prime-bound", type=int)
    ap.add_argument("--residual-factor-limit", type=int, default=1 << 40)
    ap.add_argument("--checkpoint", type=Path, default=Path(
        "artifacts/local/elliptic-curves/r20_fixedfb5000_checkpoint.json"))
    ap.add_argument(
        "--relation-ledger",
        type=Path,
        help=(
            "write exact principal generators for every newly closed relation "
            "and unit dependency in this run"
        ),
    )
    args = ap.parse_args()
    if args.trial_prime_bound is not None and args.trial_prime_bound < 2:
        raise ValueError("--trial-prime-bound must be at least 2")
    if args.residual_factor_limit < 2:
        raise ValueError("--residual-factor-limit must be at least 2")
    if args.max_residual_primes < 1:
        raise ValueError("--max-residual-primes must be positive")
    if args.adaptive_residual_specials < 0 or args.adaptive_residual_depth < 0:
        raise ValueError("adaptive residual caps must be nonnegative")
    if args.adaptive_residual_prime_bound < 2:
        raise ValueError("--adaptive-residual-prime-bound must be at least 2")

    from sage.all import QQ, ZZ, NumberField, PolynomialRing, RealField, factor, prime_range

    def integer_csv(text, description):
        try:
            values = [ZZ(item.strip()) for item in text.split(",") if item.strip()]
        except ValueError as exc:
            raise ValueError(f"{description} must be comma-separated integers") from exc
        if not values:
            raise ValueError(f"{description} must not be empty")
        return values

    x = PolynomialRing(QQ, "x").gen()
    if args.field_polynomial_ascending:
        field_coefficients = integer_csv(
            args.field_polynomial_ascending, "field-polynomial-ascending"
        )
        if len(field_coefficients) != 4 or field_coefficients[-1] != 1:
            raise ValueError("field-polynomial-ascending must define a monic cubic")
    else:
        field_coefficients = [
            ZZ(167347710468055045100164888198438918505621536951206),
            ZZ(-5750886029903523759416717668139307),
            ZZ(0),
            ZZ(1),
        ]
    f = sum(
        coefficient * x**index
        for index, coefficient in enumerate(field_coefficients)
    )
    if not f.is_irreducible():
        raise ValueError("the defining cubic must be irreducible over QQ")
    K = NumberField(f, "a")
    th = K.gen()

    if args.selmer_rational_primes:
        bad = sorted(set(integer_csv(
            args.selmer_rational_primes, "selmer-rational-primes"
        )))
        if any(prime < 2 or not prime.is_prime(proof=False) for prime in bad):
            raise ValueError("selmer-rational-primes must all be primes")
    else:
        bad = sorted(
            set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(ZZ(f.discriminant())))])
        )

    # ---- identical fixed factor base to R20FIXEDQ2 -------------------------
    fb, s_cols = [], []
    fb_index_by_hnf = {}

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

    if args.norm_factor_mode == "trial":
        trial_division_bound = args.factor_base_bound
    elif args.trial_prime_bound is None:
        trial_division_bound = min(args.factor_base_bound, 10_000)
    else:
        trial_division_bound = min(args.factor_base_bound, args.trial_prime_bound)
    trial_primes = sorted(q for q in fb_by_q if q <= trial_division_bound)

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

    special_depth = {pair: 0 for pair in specials}

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
    sparse_large_primes = SparseLargePrimeEliminator()
    generators = []
    closed_relations = []
    accepted_generator_coordinates = set()
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

    residual_root_cache = {}

    def residual_degree_one_special(signature):
        """Recover ``(ell,r)`` for a degree-one residual prime-ideal vertex."""

        ell, (degree, hnf) = signature
        if degree != 1 or ell > args.adaptive_residual_prime_bound:
            return None
        key = (ell, hnf)
        if key in residual_root_cache:
            return residual_root_cache[key]
        quotient = ZZ.quotient(ZZ(ell))
        for residue in f.change_ring(quotient).roots(multiplicities=False):
            pair = (ZZ(ell), ZZ(int(residue)))
            prime = special_prime(*pair)
            if prime is not None and str(prime.pari_hnf()) == hnf:
                residual_root_cache[key] = pair
                return pair
        residual_root_cache[key] = None
        return None

    def find_path(src, dst):
        if src == dst:
            return 0, set()
        dq = deque([src])
        parent = {src: None}
        edge = {}
        while dq:
            v = dq.popleft()
            for w, row, generator_index in graph[v]:
                if w in parent:
                    continue
                parent[w] = v
                edge[w] = (row, generator_index)
                if w == dst:
                    acc = 0
                    provenance = set()
                    cur = dst
                    while parent[cur] is not None:
                        edge_row, edge_generator = edge[cur]
                        acc ^= edge_row
                        if edge_generator in provenance:
                            provenance.remove(edge_generator)
                        else:
                            provenance.add(edge_generator)
                        cur = parent[cur]
                    return acc, provenance
                dq.append(w)
        return None, None

    def add_partial(v1, v2, row, generator_index):
        path, provenance = find_path(v1, v2)
        graph[v1].append((v2, row, generator_index))
        graph[v2].append((v1, row, generator_index))
        if path is None:
            return None, None
        if generator_index in provenance:
            provenance.remove(generator_index)
        else:
            provenance.add(generator_index)
        return path ^ row, provenance

    adaptive_added = 0
    for sq_index, (q, r) in enumerate(specials, 1):
        current_depth = special_depth[(q, r)]
        Pq = special_prime(q, r)
        if Pq is None:
            continue

        before_rank = len(pivots)
        before_after = afterS()

        candidates = p1 = p2 = phyper = cycles = indep = after_gain = 0
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
                if key in seen_abc or key in accepted_generator_coordinates:
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

                residual_primes = None
                if args.norm_factor_mode in {"trial", "hybrid"}:
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
                else:
                    co = ZZ(1)
                    used = []
                    residual_primes = []
                    for p, exponent in factor(N):
                        p = ZZ(p)
                        exponent = ZZ(exponent)
                        if p == q:
                            continue
                        if p <= args.factor_base_bound:
                            used.append(p)
                        elif int(exponent) & 1:
                            residual_primes.append(p)

                if args.norm_factor_mode == "hybrid" and co != 1 and co <= args.residual_factor_limit:
                    residual = ZZ(1)
                    for p, exponent in factor(co):
                        p = ZZ(p)
                        exponent = ZZ(exponent)
                        if p <= args.factor_base_bound:
                            used.append(p)
                        else:
                            residual *= p**exponent
                    co = residual

                row = exact_fb_row(alpha, used)
                generator_index = len(generators)
                generators.append(
                    {
                        "power_basis": [str(a), str(b), str(c)],
                        "source_special_q": [int(q), int(r)],
                        "norm": str(N),
                    }
                )
                vertices = []

                if vq & 1:
                    vertices.append(
                        (int(q), (1, str(Pq.pari_hnf())))
                    )

                if residual_primes is not None:
                    if len(residual_primes) > args.max_residual_primes:
                        continue
                    residual_signatures = []
                    for residual_prime in residual_primes:
                        if residual_prime > args.large_prime_bound:
                            residual_signatures = []
                            break
                        signature = prime_sig(alpha, residual_prime)
                        if signature is None:
                            residual_signatures = []
                            break
                        residual_signatures.append(signature)
                    if len(residual_signatures) != len(residual_primes):
                        continue
                    accepted_residual_prime += len(residual_signatures)
                    vertices.extend(residual_signatures)
                    if (
                        args.adaptive_residual_specials
                        and current_depth < args.adaptive_residual_depth
                    ):
                        for signature in residual_signatures:
                            if adaptive_added >= args.adaptive_residual_specials:
                                break
                            pair = residual_degree_one_special(signature)
                            if pair is None or pair in seen:
                                continue
                            specials.append(pair)
                            seen.add(pair)
                            special_depth[pair] = current_depth + 1
                            adaptive_added += 1
                elif co != 1:
                    # Trial/hybrid mode retains the historical one-large-prime
                    # behavior because the unfactored cofactor has no exact
                    # multi-prime support yet.
                    if co > args.large_prime_bound or not co.is_prime(proof=False):
                        continue
                    accepted_residual_prime += 1
                    sig = prime_sig(alpha, co)
                    if sig is None:
                        continue
                    vertices.append(sig)

                cycle = None
                provenance = None
                if len(vertices) == 1:
                    p1 += 1
                elif len(vertices) == 2:
                    p2 += 1
                elif len(vertices) > 2:
                    phyper += 1
                accepted_generator_coordinates.add(key)
                if args.large_prime_merge_engine == "sparse-hypergraph":
                    cycle, provenance = sparse_large_primes.add(
                        vertices, row, generator_index
                    )
                elif len(vertices) == 0:
                    cycle = row
                    provenance = {generator_index}
                elif len(vertices) == 1:
                    cycle, provenance = add_partial(
                        ROOT, vertices[0], row, generator_index
                    )
                elif len(vertices) == 2:
                    cycle, provenance = add_partial(
                        vertices[0], vertices[1], row, generator_index
                    )
                else:
                    continue

                if cycle is not None:
                    cycles += 1
                    closed_relations.append(
                        {
                            "fb_parity_mask_hex": hex(cycle),
                            "generator_indices": sorted(provenance),
                            "kind": "unit_dependency" if cycle == 0 else "fb_relation",
                        }
                    )
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
            "partial_hyperedges": phyper,
            "cycles": cycles,
            "rank_gain": len(pivots) - before_rank,
            "afterS_gain": before_after - afterS(),
            "independent": indep,
            "independent_afterS": after_gain,
            "min_norm_bits": norm_bits_min,
            "avg_norm_bits": (
                norm_bits_sum / candidates if candidates else None
            ),
            "special_depth": current_depth,
        }
        history.append(stat)
        save_checkpoint()

        print(
            f"{PROTOCOL}|stage=special_q_done|index={sq_index}/{len(specials)}"
            f"|q={q}|r={r}|rank_gain={stat['rank_gain']}"
            f"|depth={current_depth}|adaptive_added={adaptive_added}"
            f"|afterS_gain={stat['afterS_gain']}"
            f"|rank={len(pivots)}|afterS={afterS()}"
            f"|candidates={candidates}|p1={p1}|p2={p2}|cycles={cycles}"
            f"|phyper={phyper}"
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
        f"|afterS={afterS()}|checkpoint={args.checkpoint}"
        f"|lp_vertices={len(sparse_large_primes.vertex_columns)}"
        f"|lp_edges={sparse_large_primes.edge_count}"
        f"|lp_rank={len(sparse_large_primes.pivots)}"
        f"|lp_nullity={sparse_large_primes.edge_count-len(sparse_large_primes.pivots)}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=recommended_seeds|top={ranked[:15]}",
        flush=True,
    )
    if args.relation_ledger:
        ledger = {
            "schema": "elliptic-curves.bnf-free-principal-relation-ledger.v1",
            "status": "exact_new_relations_not_class_group_completion",
            "field_polynomial": str(f),
            "defining_polynomial_ascending": [
                str(value) for value in field_coefficients
            ],
            "field_discriminant": str(K.discriminant()),
            "generator_coordinate_order": ["1", "theta", "theta^2"],
            "factor_base_bound": args.factor_base_bound,
            "norm_factor_mode": args.norm_factor_mode,
            "large_prime_merge_engine": args.large_prime_merge_engine,
            "max_residual_primes": args.max_residual_primes,
            "adaptive_residual_specials": args.adaptive_residual_specials,
            "adaptive_residual_prime_bound": args.adaptive_residual_prime_bound,
            "adaptive_residual_depth": args.adaptive_residual_depth,
            "large_prime_elimination": {
                "vertex_count": len(sparse_large_primes.vertex_columns),
                "edge_count": sparse_large_primes.edge_count,
                "rank": len(sparse_large_primes.pivots),
                "dependency_count": sparse_large_primes.dependency_count,
                "nullity": (
                    sparse_large_primes.edge_count - len(sparse_large_primes.pivots)
                ),
            },
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
        args.relation_ledger.write_text(
            json.dumps(ledger, indent=2, sort_keys=True) + "\n"
        )
        print(
            f"{PROTOCOL}|stage=write_relation_ledger"
            f"|path={args.relation_ledger}|generators={len(generators)}"
            f"|closed_relations={len(closed_relations)}",
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
