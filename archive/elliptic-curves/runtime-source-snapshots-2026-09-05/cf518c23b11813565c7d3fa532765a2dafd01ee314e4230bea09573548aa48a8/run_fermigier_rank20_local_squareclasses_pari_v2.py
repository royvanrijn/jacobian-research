#!/usr/bin/env python3
"""PARI-backed odd-local square-class fingerprints for Fermigier rank-20.

Avoids Sage residue fields entirely (important on macOS/OpenBLAS setups where
NumberFieldIdeal.residue_symbol can SIGILL in matrix_modn_dense_float).

Sage is used only to load/transport the pinned 20 points. PARI handles:
  * nfinit
  * idealprimedec
  * idealval
  * idealappr (local uniformizer)
  * nfmodprinit / nfmodpr
  * issquare in the residue field

This is still not a full Selmer computation; p=2 is intentionally deferred.

Run:
  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_local_squareclasses_pari.py \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_local_squareclasses_pari.log
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import sys
import time

CAS_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS_ROOT))

from run_fermigier_rank20_pari_descent import (
    KNOWN_RANK,
    load_descent_basis,
    sage_q,
)

PROTOCOL = "R20LOCALPARI"


def f2_rank(rows):
    """Pure-Python GF(2) rank; avoids any dense finite-field matrix backend."""
    pivots = {}
    rank = 0
    for row in rows:
        v = 0
        for j, bit in enumerate(row):
            if bit:
                v ^= 1 << j
        while v:
            k = v.bit_length() - 1
            if k in pivots:
                v ^= pivots[k]
            else:
                pivots[k] = v
                rank += 1
                break
    return rank


def qstr(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"({q.numerator}/{q.denominator})"


def pari_call(pari, name, *args):
    """Call GP function robustly across cypari2 bound/global styles."""
    fn = getattr(pari, name)
    return fn(*args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/"
            "fermigier_rank20_near_miss_v1.json"
        ),
    )
    ap.add_argument(
        "--candidate-record",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/"
            "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
        ),
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier_rank20_local_squareclasses_pari.txt"
        ),
    )
    args = ap.parse_args()

    from sage.all import EllipticCurve, QQ, ZZ, pari

    basis = load_descent_basis(args.manifest, args.candidate_record)
    E = EllipticCurve(QQ, [sage_q(v) for v in basis.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis.points]
    if len(known) != KNOWN_RANK:
        raise ArithmeticError("expected 20 known points")

    short = E.integral_model().short_weierstrass_model()
    a1, a2, a3, a4, a6 = short.a_invariants()
    if a1 != 0 or a3 != 0:
        raise ArithmeticError("expected long model y^2=x^3+A*x^2+B*x+C")
    A, B, C = map(ZZ, (a2, a4, a6))
    iso = E.isomorphism_to(short)

    # PARI variable y is intentionally used for the cubic generator.
    ftxt = f"y^3+({A})*y^2+({B})*y+({C})"
    f = pari(ftxt)

    print(
        f"{PROTOCOL}|stage=input|known={KNOWN_RANK}|f={f}",
        flush=True,
    )

    t0 = time.monotonic()
    nf = pari.nfinit(f)
    theta = pari(f"Mod(y,{f})")
    print(
        f"{PROTOCOL}|stage=nfinit|status=complete|seconds={time.monotonic()-t0:.6f}",
        flush=True,
    )

    # Use the same rational bad-prime support established by the skeleton.
    # Factor the polynomial discriminant with PARI; no BNF involved.
    disc = int(pari.poldisc(f))
    fac = pari.factor(abs(disc))
    rat_bad = {2}
    for i in range(int(fac.nrows())):
        rat_bad.add(int(fac[i, 0]))
    rat_bad = sorted(rat_bad)

    print(
        f"{PROTOCOL}|stage=bad_primes|count={len(rat_bad)}|primes={rat_bad}",
        flush=True,
    )

    places = []
    for p in rat_bad:
        if p == 2:
            continue
        print(f"{PROTOCOL}|stage=place_setup|p={p}|status=start", flush=True)
        t1 = time.monotonic()
        prs = pari.idealprimedec(nf, p)
        for k in range(len(prs)):
            pr = prs[k]
            # idealappr(pr) gives an element with the prime-ideal valuation
            # pattern needed for a local uniformizer. Verify v_pr(pi)=1.
            pi_col = pari.idealappr(nf, pr)
            pi = pari.nfbasistoalg(nf, pi_col)
            vpi = int(pari.idealval(nf, pi, pr))
            if vpi != 1:
                raise ArithmeticError(
                    f"idealappr did not produce a uniformizer: p={p}, k={k}, v={vpi}"
                )
            modpr = pari.nfmodprinit(nf, pr)
            # pr[3] / pr[4] are e/f in GP's prime-ideal representation
            e = int(pr[2])
            ff = int(pr[3])
            places.append((p, k, pr, pi, modpr, e, ff))
        print(
            f"{PROTOCOL}|stage=place_setup|p={p}|status=complete"
            f"|seconds={time.monotonic()-t1:.6f}|places={len(prs)}",
            flush=True,
        )

    print(
        f"{PROTOCOL}|stage=places|odd_finite={len(places)}",
        flush=True,
    )

    val_rows = []
    local_rows = []

    for i, P0 in enumerate(known):
        print(f"{PROTOCOL}|stage=point|i={i}|status=start", flush=True)
        t1 = time.monotonic()
        Q = iso(P0)
        xq = QQ(Q[0])

        # alpha = x(P)-theta.  Construct it algebraically; do not send a
        # dynamically formatted expression back through the GP parser.
        xn = pari(int(xq.numerator()))
        xd = pari(int(xq.denominator()))
        alpha = xn / xd - theta

        vb = []
        lb = []
        for j, (p, k, pr, pi, modpr, e, ff) in enumerate(places):
            tv = time.monotonic()
            v = int(pari.idealval(nf, alpha, pr))
            vp = v & 1

            # Exact local unit. Negative valuations are fine: pi^v is in K*.
            unit = alpha / (pi ** v)
            vu = int(pari.idealval(nf, unit, pr))
            if vu != 0:
                raise ArithmeticError(
                    f"unit normalization failed point={i}, p={p}, k={k}, v={v}, vu={vu}"
                )

            residue = pari.nfmodpr(nf, unit, modpr)
            sq = bool(pari.issquare(residue))
            ub = 0 if sq else 1

            vb.append(vp)
            lb.extend((vp, ub))

            elapsed = time.monotonic() - tv
            if elapsed >= 0.25:
                print(
                    f"{PROTOCOL}|stage=local_place|point={i}|column={j}"
                    f"|p={p}|prime_index={k}|e={e}|f={ff}"
                    f"|valuation={v}|unit_square={str(sq).lower()}"
                    f"|seconds={elapsed:.6f}",
                    flush=True,
                )

        val_rows.append(vb)
        local_rows.append(lb)

        print(
            f"{PROTOCOL}|stage=point|i={i}|status=complete"
            f"|seconds={time.monotonic()-t1:.6f}"
            f"|valuation_weight={sum(vb)}|local_weight={sum(lb)}",
            flush=True,
        )

    rv = f2_rank(val_rows)
    rl = f2_rank(local_rows)

    incremental = [f2_rank(local_rows[:i+1]) for i in range(len(local_rows))]

    print(
        f"{PROTOCOL}|stage=rank|valuation_only={rv}"
        f"|odd_local_squareclasses={rl}|known_mod2_rank={KNOWN_RANK}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=incremental_rank|ranks={incremental}",
        flush=True,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "\n".join([
            f"valuation_only_rank={rv}",
            f"odd_local_squareclass_rank={rl}",
            f"known_mod2_rank={KNOWN_RANK}",
            f"odd_finite_places={len(places)}",
            f"incremental_ranks={incremental}",
        ]) + "\n"
    )

    print(
        f"{PROTOCOL}|stage=summary|visible_known_rank={rl}"
        f"|missing_known_bits={KNOWN_RANK-rl}|output={args.output}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=not_selmer"
        f"|missing=2_adic_squareclasses,real_signs,global_S_unit_relations,classgroup_2part,local_selmer_intersection",
        flush=True,
    )


if __name__ == "__main__":
    main()
