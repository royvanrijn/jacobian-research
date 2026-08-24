#!/usr/bin/env python3
"""Local square-class fingerprint for the Fermigier rank-20 Kummer images.

No class group, regulator, BNF, S-unit group, or full Selmer computation.

For every known P:
  alpha = x(P) - theta in K = Q(theta).

At each odd bad prime ideal p of K we record the two coordinates of
K_p^*/K_p^{*2}:
  * v_p(alpha) mod 2
  * quadratic character of the unit part alpha / pi^v modulo p

We also record signs at all real embeddings.  The 2-adic local square class is
deliberately omitted and should be added in a later experiment.

Run from jacobian-research root:

  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_local_squareclasses.py \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_local_squareclasses.log
"""

from __future__ import annotations

import argparse
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

PROTOCOL = "R20LOCAL2"


def f2_rank(rows):
    from sage.all import GF, Matrix
    if not rows:
        return 0
    return int(Matrix(GF(2), rows).rank())


def bit_from_residue_symbol(s):
    # For m=2 the symbol is +/-1.
    return 0 if s == 1 else 1


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
            "fermigier_rank20_local_squareclasses.txt"
        ),
    )
    args = ap.parse_args()

    from sage.all import EllipticCurve, QQ, ZZ, PolynomialRing, NumberField, factor

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != KNOWN_RANK or not basis.mod2_certified:
        raise ArithmeticError("missing exact rank-20 mod-2 certificate")

    E = EllipticCurve(QQ, [sage_q(v) for v in basis.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis.points]

    short = E.integral_model().short_weierstrass_model()
    a1, a2, a3, a4, a6 = short.a_invariants()
    if a1 != 0 or a3 != 0:
        raise ArithmeticError("expected a1=a3=0")
    A, B, C = ZZ(a2), ZZ(a4), ZZ(a6)

    R = PolynomialRing(QQ, "x")
    x = R.gen()
    f = x**3 + A*x**2 + B*x + C
    disc = ZZ(f.discriminant())

    K = NumberField(f, "th")
    th = K.gen()
    iso = E.isomorphism_to(short)

    rat_bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(disc))]))

    print(
        f"{PROTOCOL}|stage=input|known={KNOWN_RANK}"
        f"|disc_bits={abs(disc).nbits()}|field_disc_bits={abs(ZZ(K.discriminant())).nbits()}"
        f"|bad_primes={rat_bad}",
        flush=True,
    )

    # Build odd local places.  A local uniformizer is enough; we do not need
    # a principal generator of the global prime ideal.
    places = []
    for p in rat_bad:
        if p == 2:
            continue
        print(f"{PROTOCOL}|stage=place_setup|p={p}|status=start", flush=True)
        t0 = time.monotonic()
        for P in K.primes_above(p):
            pi = K.uniformizer(P)
            if int(pi.valuation(P)) != 1:
                raise ArithmeticError(f"bad uniformizer at {P}")
            places.append((int(p), P, pi))
        print(
            f"{PROTOCOL}|stage=place_setup|p={p}|status=complete"
            f"|seconds={time.monotonic()-t0:.6f}"
            f"|places={len(K.primes_above(p))}",
            flush=True,
        )

    # K is expected to be totally real here.  These are square-class sign
    # characters at infinity.
    real_embeddings = K.real_embeddings(prec=128)
    print(
        f"{PROTOCOL}|stage=places|odd_finite={len(places)}"
        f"|real={len(real_embeddings)}",
        flush=True,
    )

    val_rows = []
    odd_rows = []
    full_rows = []

    for i, P0 in enumerate(known):
        print(f"{PROTOCOL}|stage=point|i={i}|status=start", flush=True)
        t0 = time.monotonic()
        Q = iso(P0)
        xq = QQ(Q[0])
        alpha = K(xq) - th

        val_bits = []
        local_bits = []

        for j, (p, P, pi) in enumerate(places):
            v = int(alpha.valuation(P))
            vp = v & 1

            # Remove the entire valuation, not just its parity.  The remaining
            # element is a P-adic unit and its residue square/nonsquare bit is
            # the second coordinate of the odd local square class.
            unit = alpha / (pi ** v)
            if int(unit.valuation(P)) != 0:
                raise ArithmeticError(
                    f"normalisation failed for point={i}, place={j}, p={p}"
                )

            try:
                symbol = P.residue_symbol(unit, 2)
            except Exception as exc:
                print(
                    f"{PROTOCOL}|stage=local_symbol|status=error|point={i}"
                    f"|place={j}|p={p}|valuation={v}|error={exc}",
                    flush=True,
                )
                raise

            ub = bit_from_residue_symbol(symbol)
            val_bits.append(vp)
            local_bits.extend((vp, ub))

        sign_bits = []
        for emb in real_embeddings:
            z = emb(alpha)
            if z == 0:
                raise ArithmeticError("Kummer element vanished at real embedding")
            sign_bits.append(1 if z < 0 else 0)

        val_rows.append(val_bits)
        odd_rows.append(local_bits)
        full_rows.append(local_bits + sign_bits)

        print(
            f"{PROTOCOL}|stage=point|i={i}|status=complete"
            f"|seconds={time.monotonic()-t0:.6f}"
            f"|valuation_weight={sum(val_bits)}"
            f"|odd_local_weight={sum(local_bits)}"
            f"|real_signs={sign_bits}",
            flush=True,
        )

    rv = f2_rank(val_rows)
    ro = f2_rank(odd_rows)
    rf = f2_rank(full_rows)

    print(
        f"{PROTOCOL}|stage=rank"
        f"|valuation_only={rv}"
        f"|odd_local_squareclasses={ro}"
        f"|odd_plus_real={rf}"
        f"|known_mod2_rank={KNOWN_RANK}",
        flush=True,
    )

    # Incremental rank identifies which known MW directions become visible
    # only after adding unit-residue / real characters.
    from sage.all import GF, Matrix
    incremental = []
    for i in range(len(full_rows)):
        r = int(Matrix(GF(2), full_rows[:i+1]).rank())
        incremental.append(r)
    print(
        f"{PROTOCOL}|stage=incremental_rank|ranks={incremental}",
        flush=True,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "\n".join([
            f"valuation_only_rank={rv}",
            f"odd_local_squareclass_rank={ro}",
            f"odd_plus_real_rank={rf}",
            f"known_mod2_rank={KNOWN_RANK}",
            f"odd_finite_places={len(places)}",
            f"real_places={len(real_embeddings)}",
            f"incremental_ranks={incremental}",
        ]) + "\n"
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|visible_known_rank={rf}|missing_known_bits={KNOWN_RANK-rf}"
        f"|output={args.output}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=not_selmer"
        f"|missing=2_adic_squareclasses,global_S_unit_relations,classgroup_2part,local_selmer_intersection",
        flush=True,
    )


if __name__ == "__main__":
    main()
