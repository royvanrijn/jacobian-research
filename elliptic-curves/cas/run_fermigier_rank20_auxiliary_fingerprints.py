#!/usr/bin/env python3
"""Find cheap auxiliary good-prime fingerprints separating the 20 known Kummer classes.

Start from the same Kummer elements alpha_i = x(P_i)-theta.  At small rational
good primes q (q not dividing 2*disc(f)), compute local square-class characters
at all primes P|q in K:
    v_P(alpha_i) mod 2
    square/nonsquare of the unit part mod P

Greedily keep only rational primes whose local characters increase the F2-rank.
This is NOT adding Selmer conditions: these are witness coordinates only.

Goal: see whether the known rank-20 Kummer image can be represented faithfully
without bnfinit/global unit computations.

Run:
  PYTHONUNBUFFERED=1 caffeinate -i \
    sage -python elliptic-curves/cas/run_fermigier_rank20_auxiliary_fingerprints.py \
    --prime-bound 5000 \
    2>&1 | tee artifacts/local/elliptic-curves/fermigier_rank20_auxiliary_fingerprints.log
"""

from __future__ import annotations

from research_runtime.pari_context import prepared_prime_ideals, prepared_factor

import argparse
from fractions import Fraction
import json
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

PROTOCOL = "R20AUXFP"


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


def f2_mask(row):
    """Pack an explicitly ordered GF(2) row for the residual quotient ledger."""

    return sum(int(bit) << index for index, bit in enumerate(row))


def qpari(pari, q):
    try:
        n = q.numerator()
        d = q.denominator()
    except TypeError:
        n = q.numerator
        d = q.denominator
    return pari(int(n)) / pari(int(d))


def local_square_everywhere_2(pari, nf, two_primes, a):
    return all(bool(pari.nfislocalpower(nf, pr, a, 2)) for pr in two_primes)


def two_adic_coords(pari, nf, two_primes, alphas):
    basis = []
    basis_origins = []
    coords = []
    for i, alpha in enumerate(alphas):
        r = len(basis)
        found = None
        for mask in range(1 << r):
            c = alpha
            bits = [0] * r
            for j in range(r):
                if (mask >> j) & 1:
                    c /= basis[j]
                    bits[j] = 1
            if local_square_everywhere_2(pari, nf, two_primes, c):
                found = bits
                break
        if found is None:
            basis.append(alpha)
            basis_origins.append(i)
            for row in coords:
                row.append(0)
            row = [0] * len(basis)
            row[-1] = 1
            coords.append(row)
        else:
            coords.append(found)
    return basis, basis_origins, coords


def prime_local_rows(pari, nf, alphas, q):
    """Return 2 bits/place for rational prime q, using PARI only."""
    places = []
    for pr in prepared_prime_ideals(nf, q):
        pi_col = pari.idealappr(nf, pr)
        pi = pari.nfbasistoalg(nf, pi_col)
        if int(pari.idealval(nf, pi, pr)) != 1:
            raise ArithmeticError(f"bad uniformizer q={q}")
        modpr = pari.nfmodprinit(nf, pr)
        places.append((pr, pi, modpr))

    rows = []
    for alpha in alphas:
        row = []
        for pr, pi, modpr in places:
            v = int(pari.idealval(nf, alpha, pr))
            unit = alpha / (pi ** v)
            residue = pari.nfmodpr(nf, unit, modpr)
            row.extend((v & 1, 0 if bool(pari.issquare(residue)) else 1))
        rows.append(row)
    return rows, tuple(str(pr) for pr, _, _ in places)


def append_columns(rows, extra):
    return [rows[i] + extra[i] for i in range(len(rows))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=Path(
        "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"))
    ap.add_argument("--candidate-record", type=Path, default=Path(
        "artifacts/generated-results/elliptic-curves/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"))
    ap.add_argument("--prime-bound", type=int, default=5000)
    ap.add_argument("--target-rank", type=int, default=20)
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier_rank20_signature_map.json"
        ),
        help=(
            "write the known Kummer image in the BNF-free local/fingerprint "
            "quotient format"
        ),
    )
    args = ap.parse_args()

    from sage.all import AA, EllipticCurve, PolynomialRing, QQ, ZZ, pari, prime_range

    data = load_descent_basis(args.manifest, args.candidate_record)
    E = EllipticCurve(QQ, [sage_q(v) for v in data.model])
    known = [E(sage_q(x), sage_q(y)) for x, y in data.points]

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

    alphas = []
    xqs = []
    for P in known:
        Q = iso(P)
        xq = QQ(Q[0])
        xqs.append(xq)
        alphas.append(qpari(pari, xq) - theta)

    # Baseline: odd bad-prime local square classes.
    rows = [[] for _ in alphas]
    local_coordinates = []
    for p in sorted(bad):
        if p == 2:
            continue
        extra, places = prime_local_rows(pari, nf, alphas, p)
        rows = append_columns(rows, extra)
        for place_index, prime_ideal in enumerate(places):
            local_coordinates.extend(
                [
                    {
                        "kind": "odd_valuation_parity",
                        "rational_prime": p,
                        "prime_ideal": prime_ideal,
                        "place_index": place_index,
                    },
                    {
                        "kind": "odd_unit_squareclass",
                        "rational_prime": p,
                        "prime_ideal": prime_ideal,
                        "place_index": place_index,
                    },
                ]
            )

    rank_odd = f2_rank(rows)

    # Add exact 2-adic coordinates.
    two_primes = list(prepared_prime_ideals(nf, 2))
    two_basis, two_basis_origins, tw = two_adic_coords(
        pari, nf, two_primes, alphas
    )
    rows = append_columns(rows, tw)
    local_coordinates.extend(
        {
            "kind": "two_adic_product_basis",
            "basis_index": index,
            "generator_power_basis": f"({xqs[two_basis_origins[index]]},-1,0)",
            "two_adic_primes": [str(pr) for pr in two_primes],
        }
        for index in range(len(two_basis))
    )
    rank_odd2 = f2_rank(rows)

    # Add real signs.
    polynomial_ring = PolynomialRing(QQ, "z")
    z = polynomial_ring.gen()
    exact_roots = list((z**3 + QQ(A) * z**2 + QQ(B) * z + QQ(C)).roots(
        AA, multiplicities=False
    ))
    real_rows = []
    for xq in xqs:
        real_rows.append([1 if QQ(xq) - root < 0 else 0 for root in exact_roots])
    rows = append_columns(rows, real_rows)
    local_coordinates.extend(
        {
            "kind": "real_sign",
            "embedding_index": index,
            "root_order": "increasing_real_root",
        }
        for index in range(len(exact_roots))
    )
    baseline_rank = f2_rank(rows)
    local_rows = [list(row) for row in rows]

    print(
        f"{PROTOCOL}|stage=baseline|odd_rank={rank_odd}|odd_plus_2={rank_odd2}"
        f"|odd_plus_2_plus_real={baseline_rank}|target={args.target_rank}",
        flush=True,
    )

    selected = []
    fingerprint_coordinates = []
    current = rows

    for q in prime_range(3, args.prime_bound + 1):
        q = int(q)
        if q in bad:
            continue

        t0 = time.monotonic()
        try:
            extra, places = prime_local_rows(pari, nf, alphas, q)
        except Exception as exc:
            print(
                f"{PROTOCOL}|stage=prime|q={q}|status=error|error={exc}",
                flush=True,
            )
            continue

        trial = append_columns(current, extra)
        old_rank = f2_rank(current)
        new_rank = f2_rank(trial)
        elapsed = time.monotonic() - t0

        if new_rank > old_rank:
            selected.append((q, len(places), new_rank - old_rank))
            current = trial
            for place_index, prime_ideal in enumerate(places):
                fingerprint_coordinates.extend(
                    [
                        {
                            "kind": "auxiliary_valuation_parity",
                            "rational_prime": q,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                        {
                            "kind": "auxiliary_unit_squareclass",
                            "rational_prime": q,
                            "prime_ideal": prime_ideal,
                            "place_index": place_index,
                        },
                    ]
                )
            print(
                f"{PROTOCOL}|stage=prime|q={q}|status=selected"
                f"|places={len(places)}|gain={new_rank-old_rank}|rank={new_rank}"
                f"|seconds={elapsed:.6f}",
                flush=True,
            )
        elif elapsed >= 0.1:
            print(
                f"{PROTOCOL}|stage=prime|q={q}|status=no_gain"
                f"|places={len(places)}|rank={old_rank}|seconds={elapsed:.6f}",
                flush=True,
            )

        if f2_rank(current) >= args.target_rank:
            break

    final_rank = f2_rank(current)
    print(
        f"{PROTOCOL}|stage=summary|baseline_rank={baseline_rank}|final_rank={final_rank}"
        f"|selected={selected}|prime_bound={args.prime_bound}",
        flush=True,
    )

    if final_rank == KNOWN_RANK:
        print(
            f"{PROTOCOL}|result=faithful_known_kummer_fingerprint"
            f"|dimensions={KNOWN_RANK}|auxiliary_primes={[q for q,_,_ in selected]}",
            flush=True,
        )
    else:
        print(
            f"{PROTOCOL}|result=incomplete_fingerprint"
            f"|missing={KNOWN_RANK-final_rank}|action=raise_prime_bound",
            flush=True,
        )

    # This is a faithful coordinate representation of the *known* Kummer
    # image, not a Selmer computation.  A relation collector can append an
    # unexplained global squareclass in these exact coordinates and immediately
    # reduce it modulo the displayed known-MW target image.
    fingerprint_width = len(fingerprint_coordinates)
    if any(len(row) != len(local_coordinates) + fingerprint_width for row in current):
        raise ArithmeticError("signature-coordinate bookkeeping lost alignment")
    signature_map = {
        "schema": "elliptic-curves.bnf-free-signature-map.v1",
        "status": "known_kummer_image_only_not_a_selmer_bound",
        "field_generator": "theta",
        "generator_coordinate_order": ["1", "theta", "theta^2"],
        "defining_polynomial_ascending": [str(C), str(B), str(A), "1"],
        "local_dimension": len(local_coordinates),
        "fingerprint_dimension": fingerprint_width,
        "local_coordinates": local_coordinates,
        "fingerprint_coordinates": fingerprint_coordinates,
        "known_mw_images": [
            {
                "label": f"P{index + 1}",
                "generator": f"({xqs[index]},-1,0)",
                "generator_coefficients": [str(xqs[index]), "-1", "0"],
                "local": f"0x{f2_mask(local_rows[index]):x}",
                "fingerprint": f"0x{f2_mask(current[index][len(local_coordinates):]):x}",
            }
            for index in range(KNOWN_RANK)
        ],
        "selected_auxiliary_primes": [q for q, _, _ in selected],
        "known_mw_target_rank": final_rank,
        "class_quotient_certification": {
            "method": "none",
            "remaining_dimension_upper_bound": None,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(signature_map, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=write_signature_map|path={args.output}"
        f"|local_dimensions={len(local_coordinates)}"
        f"|fingerprint_dimensions={fingerprint_width}"
        f"|known_mw_target_rank={final_rank}",
        flush=True,
    )


if __name__ == "__main__":
    main()
