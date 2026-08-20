#!/usr/bin/env python3
"""Cheap Kummer skeleton for the Fermigier rank-20 anchor.

This deliberately avoids bnfinit/class-group/regulator computations.

It:
  * loads the pinned rank-20 basis from the repository,
  * builds the cubic 2-division field K = Q(theta),
  * factors only rational primes dividing 2*disc(f),
  * computes alpha_P = x(P) - theta for the 20 known points,
  * records v_p(alpha_P) mod 2 at primes p of K over the bad rational primes,
  * computes the F2-rank of those 20 parity vectors.

This is NOT a full 2-Selmer computation.  It measures how much of the obvious
S-valuation square-class space is already occupied by the known Mordell-Weil
subgroup.

Run from jacobian-research root:

  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_kummer_skeleton.py \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_kummer_skeleton.log
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

PROTOCOL = "R20KUMMER"


def parity_rank(rows):
    """Rank over F2 of integer 0/1 rows."""
    if not rows:
        return 0
    from sage.all import GF, Matrix
    return int(Matrix(GF(2), rows).rank())


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
            "artifacts/generated-results/"
            "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
        ),
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier_rank20_kummer_skeleton.txt"
        ),
    )
    args = ap.parse_args()

    from sage.all import (
        EllipticCurve, QQ, ZZ, PolynomialRing, NumberField, factor
    )

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != KNOWN_RANK or not basis.mod2_certified:
        raise ArithmeticError("missing exact rank-20 mod-2 certificate")

    E = EllipticCurve(QQ, [sage_q(v) for v in basis.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis.points]
    if len(known) != KNOWN_RANK:
        raise ArithmeticError("expected exactly 20 known points")

    # Move to the integral long Weierstrass model used by Simon:
    # Y^2 = X^3 + A X^2 + B X + C
    short = E.integral_model().short_weierstrass_model()
    a1, a2, a3, a4, a6 = short.a_invariants()
    if a1 != 0 or a3 != 0:
        raise ArithmeticError("expected a model with a1=a3=0")

    A = ZZ(a2)
    B = ZZ(a4)
    C = ZZ(a6)

    R = PolynomialRing(QQ, "x")
    x = R.gen()
    f = x**3 + A*x**2 + B*x + C

    print(
        f"{PROTOCOL}|stage=input|known={KNOWN_RANK}"
        f"|basis_sha256={basis.basis_sha256}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=field_polynomial|f={f}"
        f"|disc_bits={ZZ(f.discriminant()).abs().nbits()}",
        flush=True,
    )

    # Rational bad-prime support for K(S,2): primes dividing 2*disc(f).
    disc = ZZ(f.discriminant())
    rat_bad = sorted(set([ZZ(2)] + [ZZ(p) for p, _ in factor(abs(disc))]))
    print(
        f"{PROTOCOL}|stage=bad_primes|count={len(rat_bad)}"
        f"|primes={rat_bad}",
        flush=True,
    )

    t0 = time.monotonic()
    K = NumberField(f, "th")
    th = K.gen()
    print(
        f"{PROTOCOL}|stage=number_field|status=complete"
        f"|seconds={time.monotonic()-t0:.6f}"
        f"|degree={K.degree()}|disc={K.discriminant()}",
        flush=True,
    )

    prime_ideals = []
    prime_meta = []

    for p in rat_bad:
        t1 = time.monotonic()
        fac = K.ideal(p).factor()
        elapsed = time.monotonic() - t1

        local = []
        for P, e in fac:
            idx = len(prime_ideals)
            prime_ideals.append(P)
            meta = {
                "index": idx,
                "p": int(p),
                "e": int(e),
                "f": int(P.residue_class_degree()),
                "norm": int(P.norm()),
                "gens": str(P.gens_reduced()),
            }
            prime_meta.append(meta)
            local.append((idx, meta["e"], meta["f"], meta["norm"]))

        print(
            f"{PROTOCOL}|stage=factor_prime|p={p}|seconds={elapsed:.6f}"
            f"|factors={local}",
            flush=True,
        )

    print(
        f"{PROTOCOL}|stage=S_prime_ideals|count={len(prime_ideals)}",
        flush=True,
    )

    # Transport points onto the chosen short model exactly.
    # Sage provides an isomorphism from E to short.
    iso = E.isomorphism_to(short)

    rows = []
    row_weights = []

    for i, P in enumerate(known):
        Q = iso(P)
        if Q.is_zero():
            raise ArithmeticError("known infinite-order point mapped to infinity")

        xq = QQ(Q[0])
        alpha = K(xq) - th

        vals = []
        nz = []
        for j, Pideal in enumerate(prime_ideals):
            v = int(alpha.valuation(Pideal))
            bit = v & 1
            vals.append(bit)
            if bit:
                nz.append(j)

        rows.append(vals)
        row_weights.append(len(nz))
        print(
            f"{PROTOCOL}|stage=kummer_point|i={i}"
            f"|x={xq}|weight={len(nz)}|support={nz}",
            flush=True,
        )

    rank_known = parity_rank(rows)
    ambient = len(prime_ideals)
    naive_residual = ambient - rank_known

    print(
        f"{PROTOCOL}|stage=valuation_matrix"
        f"|rows={len(rows)}|cols={ambient}"
        f"|rank={rank_known}"
        f"|naive_residual_dim={naive_residual}"
        f"|weights={row_weights}",
        flush=True,
    )

    # Also show per-column occupancy: useful for spotting dead/local-only directions.
    col_weights = [
        sum(rows[i][j] for i in range(len(rows)))
        for j in range(ambient)
    ]
    print(
        f"{PROTOCOL}|stage=column_weights|weights={col_weights}",
        flush=True,
    )

    # Human-readable mapping from column index to prime ideal metadata.
    for meta in prime_meta:
        print(
            f"{PROTOCOL}|stage=S_column"
            f"|index={meta['index']}|p={meta['p']}|e={meta['e']}"
            f"|f={meta['f']}|norm={meta['norm']}|gens={meta['gens']}",
            flush=True,
        )

    summary = (
        f"known_rank={rank_known}\n"
        f"S_prime_ideal_count={ambient}\n"
        f"naive_residual_dim={naive_residual}\n"
        f"rational_bad_primes={rat_bad}\n"
        f"column_weights={col_weights}\n"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(summary)

    print(
        f"{PROTOCOL}|stage=summary|known_valuation_rank={rank_known}"
        f"|S_prime_ideal_count={ambient}"
        f"|naive_residual_dim={naive_residual}"
        f"|output={args.output}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=not_full_selmer"
        f"|missing=units_mod_squares,classgroup_2part,norm_condition,local_kummer_conditions",
        flush=True,
    )


if __name__ == "__main__":
    main()
