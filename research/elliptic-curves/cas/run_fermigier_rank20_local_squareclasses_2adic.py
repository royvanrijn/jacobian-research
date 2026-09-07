#!/usr/bin/env python3
"""Add exact 2-adic local square-class coordinates to the Fermigier rank-20 Kummer fingerprint.

This avoids BNF/class groups entirely.

For alpha_i = x(P_i)-theta in the cubic field K, we:
  * compute the odd-local square-class coordinates as before;
  * compute the class of alpha_i in the product over primes P|2 of K_P^*/K_P^{*2}
    using PARI nfislocalpower;
  * derive explicit F2 coordinates incrementally by subgroup membership;
  * add real-sign characters;
  * report ranks of odd-only, 2-adic-only, odd+2-adic, and odd+2-adic+real fingerprints.

Run:
  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_local_squareclasses_2adic.py \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_local_squareclasses_2adic.log
"""

from __future__ import annotations

from research_runtime.pari_context import prepared_prime_ideals, prepared_factor

import argparse
from fractions import Fraction
from itertools import product
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

PROTOCOL = "R20LOCALFULL"


def f2_rank(rows):
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


def qpari(pari, q):
    try:
        n = q.numerator()
        d = q.denominator()
    except TypeError:
        n = q.numerator
        d = q.denominator
    return pari(int(n)) / pari(int(d))


def local_square_everywhere_2(pari, nf, two_primes, a):
    for pr in two_primes:
        if not bool(pari.nfislocalpower(nf, pr, a, 2)):
            return False
    return True


def class_coordinates_2adic(pari, nf, two_primes, alphas):
    """Incremental basis for product_{P|2} K_P^*/K_P^{*2}.

    basis_alphas contains actual K* representatives.  For each alpha we test
    all products of the current basis.  The local space is tiny (expected <=7
    dimensions here), so exhaustive subset membership is cheap and robust.
    """
    basis = []
    coords = []

    for i, alpha in enumerate(alphas):
        r = len(basis)
        found = None

        # alpha has coordinates mask iff alpha / product(b_j^mask_j) is a
        # square at every P|2.  Mod squares, division == multiplication, but
        # division keeps representative growth milder.
        for mask in range(1 << r):
            candidate = alpha
            bits = [0] * r
            for j in range(r):
                if (mask >> j) & 1:
                    candidate = candidate / basis[j]
                    bits[j] = 1
            if local_square_everywhere_2(pari, nf, two_primes, candidate):
                found = bits
                break

        if found is None:
            # New independent 2-adic direction.
            basis.append(alpha)
            # Existing coordinate rows get a zero in the new column.
            for row in coords:
                row.append(0)
            row = [0] * len(basis)
            row[-1] = 1
            coords.append(row)
            print(
                f"{PROTOCOL}|stage=2adic_basis|point={i}|action=new_basis"
                f"|local_rank={len(basis)}",
                flush=True,
            )
        else:
            # Current basis dimension unchanged.
            coords.append(found)
            print(
                f"{PROTOCOL}|stage=2adic_basis|point={i}|action=dependent"
                f"|coords={found}|local_rank={len(basis)}",
                flush=True,
            )

    return basis, coords


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--manifest", type=Path,
        default=Path("artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json")
    )
    ap.add_argument(
        "--candidate-record", type=Path,
        default=Path("artifacts/generated-results/elliptic-curves/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json")
    )
    ap.add_argument(
        "--output", type=Path,
        default=Path("artifacts/local/elliptic-curves/fermigier_rank20_local_squareclasses_2adic.txt")
    )
    args = ap.parse_args()

    from sage.all import EllipticCurve, QQ, ZZ, pari

    basis_data = load_descent_basis(args.manifest, args.candidate_record)
    E = EllipticCurve(QQ, [sage_q(v) for v in basis_data.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in basis_data.points]

    short = E.integral_model().short_weierstrass_model()
    a1, A, a3, B, C = short.a_invariants()
    if a1 != 0 or a3 != 0:
        raise ArithmeticError("expected a1=a3=0")
    A, B, C = map(ZZ, (A, B, C))
    iso = E.isomorphism_to(short)

    f = pari(f"y^3+({A})*y^2+({B})*y+({C})")
    from research_runtime.pari_context import prepared_nf
    nf = prepared_nf(f)
    theta = pari(f"Mod(y,{f})")

    disc = abs(int(pari.poldisc(f)))
    ff = prepared_factor(disc)
    bad = {2}
    for i in range(int(ff.nrows())):
        bad.add(int(ff[i, 0]))
    bad = sorted(bad)

    print(
        f"{PROTOCOL}|stage=input|known={len(known)}|bad_primes={bad}",
        flush=True,
    )

    # Prepare Kummer elements.
    alphas = []
    for P in known:
        Q = iso(P)
        xq = QQ(Q[0])
        alphas.append(qpari(pari, xq) - theta)

    # ----- Odd local coordinates: valuation parity + unit residue square bit.
    odd_places = []
    for p in bad:
        if p == 2:
            continue
        for pr in prepared_prime_ideals(nf, p):
            pi_col = pari.idealappr(nf, pr)
            pi = pari.nfbasistoalg(nf, pi_col)
            if int(pari.idealval(nf, pi, pr)) != 1:
                raise ArithmeticError(f"bad uniformizer at p={p}")
            modpr = pari.nfmodprinit(nf, pr)
            odd_places.append((p, pr, pi, modpr))

    odd_rows = []
    for i, alpha in enumerate(alphas):
        row = []
        for p, pr, pi, modpr in odd_places:
            v = int(pari.idealval(nf, alpha, pr))
            unit = alpha / (pi ** v)
            residue = pari.nfmodpr(nf, unit, modpr)
            row.extend((v & 1, 0 if bool(pari.issquare(residue)) else 1))
        odd_rows.append(row)

    odd_rank = f2_rank(odd_rows)
    print(
        f"{PROTOCOL}|stage=odd_local|places={len(odd_places)}"
        f"|columns={len(odd_rows[0])}|rank={odd_rank}",
        flush=True,
    )

    # ----- Exact product of 2-adic square-class groups.
    two_primes = list(prepared_prime_ideals(nf, 2))
    two_meta = []
    for k, pr in enumerate(two_primes):
        # PARI prime ideal structure has e,f in components 3,4 in GP indexing.
        e = int(pr[2])
        residue_degree = int(pr[3])
        two_meta.append((k, e, residue_degree))
    print(
        f"{PROTOCOL}|stage=2adic_setup|primes={len(two_primes)}|meta={two_meta}",
        flush=True,
    )

    t2 = time.monotonic()
    two_basis, two_rows = class_coordinates_2adic(pari, nf, two_primes, alphas)
    two_seconds = time.monotonic() - t2
    two_rank = f2_rank(two_rows)

    print(
        f"{PROTOCOL}|stage=2adic|rank={two_rank}|basis_size={len(two_basis)}"
        f"|seconds={two_seconds:.6f}",
        flush=True,
    )

    # ----- Combine odd and 2-adic.
    odd2_rows = [odd_rows[i] + two_rows[i] for i in range(KNOWN_RANK)]
    odd2_rank = f2_rank(odd2_rows)

    # ----- Real sign characters, obtained from the three real roots theta_j.
    # f is totally real for this field.
    roots = list(pari.polrootsreal(f))
    real_rows = []
    for P in known:
        Q = iso(P)
        xq = float(Q[0])
        signs = []
        for root in roots:
            # Use PARI high-precision comparison rather than converting alpha
            # through a Sage number-field embedding.
            val = pari(qpari(pari, QQ(Q[0])) - root)
            signs.append(1 if val < 0 else 0)
        real_rows.append(signs)

    all_rows = [odd2_rows[i] + real_rows[i] for i in range(KNOWN_RANK)]
    all_rank = f2_rank(all_rows)

    print(
        f"{PROTOCOL}|stage=ranks|odd={odd_rank}|two_adic={two_rank}"
        f"|odd_plus_2={odd2_rank}|real_places={len(roots)}"
        f"|odd_plus_2_plus_real={all_rank}|known_mod2_rank={KNOWN_RANK}",
        flush=True,
    )

    incremental = [f2_rank(all_rows[:i+1]) for i in range(KNOWN_RANK)]
    print(
        f"{PROTOCOL}|stage=incremental_rank|ranks={incremental}",
        flush=True,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "\n".join([
            f"odd_rank={odd_rank}",
            f"two_adic_rank={two_rank}",
            f"odd_plus_2_rank={odd2_rank}",
            f"real_places={len(roots)}",
            f"odd_plus_2_plus_real_rank={all_rank}",
            f"known_mod2_rank={KNOWN_RANK}",
            f"two_prime_meta={two_meta}",
            f"incremental_ranks={incremental}",
        ]) + "\n"
    )

    print(
        f"{PROTOCOL}|stage=summary|visible_known_rank={all_rank}"
        f"|missing_known_bits={KNOWN_RANK-all_rank}|output={args.output}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|warning=not_full_selmer"
        f"|missing=global_S_unit_relations,classgroup_2part,actual_local_selmer_intersection",
        flush=True,
    )


if __name__ == "__main__":
    main()
