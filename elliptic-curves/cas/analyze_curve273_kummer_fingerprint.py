#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

from sage.all import QQ, ZZ, pari, prime_range

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import (
    POINTS,
    SHORT_POINTS,
    short_coefficients,
)

PROTOCOL = "R30KUMMER"
KNOWN_RANK = 30


def sage_q(value):
    """Convert Python Fraction / Sage rational / integer safely to Sage QQ."""
    n = value.numerator
    d = value.denominator

    if callable(n):
        n = n()
    if callable(d):
        d = d()

    return QQ(ZZ(n)) / QQ(ZZ(d))


def f2_rank(rows):
    pivots = {}

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
                break

    return len(pivots)


def f2_mask(row):
    return sum(int(bit) << index for index, bit in enumerate(row))


def qpari(pari_obj, q):
    q = sage_q(q)
    return pari_obj(int(q.numerator())) / pari_obj(int(q.denominator()))


def append_columns(rows, extra):
    return [
        rows[i] + extra[i]
        for i in range(len(rows))
    ]


def local_square_everywhere_2(pari_obj, nf, two_primes, alpha):
    return all(
        bool(
            pari_obj.nfislocalpower(
                nf,
                prime,
                alpha,
                2,
            )
        )
        for prime in two_primes
    )


def two_adic_coords(pari_obj, nf, two_primes, alphas):
    """
    Construct coordinates in the product of local K*/K*^2 at primes over 2.
    """

    basis = []
    basis_origins = []
    coords = []

    for index, alpha in enumerate(alphas):
        r = len(basis)
        found = None

        for mask in range(1 << r):
            c = alpha
            bits = [0] * r

            for j in range(r):
                if (mask >> j) & 1:
                    c /= basis[j]
                    bits[j] = 1

            if local_square_everywhere_2(
                pari_obj,
                nf,
                two_primes,
                c,
            ):
                found = bits
                break

        if found is None:
            basis.append(alpha)
            basis_origins.append(index)

            for row in coords:
                row.append(0)

            row = [0] * len(basis)
            row[-1] = 1
            coords.append(row)
        else:
            coords.append(found)

    return basis, basis_origins, coords


def prime_local_rows(pari_obj, nf, alphas, q):
    """
    For every prime ideal P|q, append:

        v_P(alpha) mod 2
        square/nonsquare(unit part mod P)
    """

    places = []

    for pr in pari_obj.idealprimedec(nf, q):
        pi_col = pari_obj.idealappr(nf, pr)
        pi = pari_obj.nfbasistoalg(nf, pi_col)

        if int(pari_obj.idealval(nf, pi, pr)) != 1:
            raise ArithmeticError(
                f"bad uniformizer at q={q}"
            )

        modpr = pari_obj.nfmodprinit(nf, pr)

        places.append(
            (pr, pi, modpr)
        )

    rows = []

    for alpha in alphas:
        row = []

        for pr, pi, modpr in places:
            v = int(
                pari_obj.idealval(
                    nf,
                    alpha,
                    pr,
                )
            )

            unit = alpha / (pi ** v)

            residue = pari_obj.nfmodpr(
                nf,
                unit,
                modpr,
            )

            row.extend(
                (
                    v & 1,
                    0
                    if bool(pari_obj.issquare(residue))
                    else 1,
                )
            )

        rows.append(row)

    return rows, tuple(str(pr) for pr, _, _ in places)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--prime-bound",
        type=int,
        default=5000,
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/curve273_signature_map.json"
        ),
        help=(
            "write known Kummer images in the BNF-free local/fingerprint "
            "quotient format"
        ),
    )

    args = ap.parse_args()

    if len(POINTS) != KNOWN_RANK:
        raise AssertionError(
            f"expected {KNOWN_RANK} points"
        )

    coeffs = short_coefficients()

    if any(coeffs[:3]):
        raise AssertionError(
            "expected integral short model"
        )

    A = ZZ(sage_q(coeffs[3]))
    B = ZZ(sage_q(coeffs[4]))

    # --------------------------------------------------------
    # 2-division field
    #
    # short model:
    #
    #     Y^2 = X^3 + A X + B
    #
    # so theta is a root of:
    #
    #     f(T)=T^3+A*T+B.
    # --------------------------------------------------------

    f = pari(
        f"z^3+({A})*z+({B})"
    )

    print(
        f"{PROTOCOL}|stage=input"
        f"|points={len(SHORT_POINTS)}"
        f"|A_bits={abs(int(A)).bit_length()}"
        f"|B_bits={abs(int(B)).bit_length()}",
        flush=True,
    )

    t0 = time.monotonic()

    nf = pari.nfinit(f)

    print(
        f"{PROTOCOL}|stage=nfinit"
        f"|seconds={time.monotonic()-t0:.6f}",
        flush=True,
    )

    theta = pari(
        f"Mod(z,{f})"
    )

    disc = abs(
        int(
            pari.poldisc(f)
        )
    )

    ff = pari.factor(disc)

    bad = {2}

    factors = []

    for i in range(int(ff.nrows())):
        p = int(ff[i, 0])
        e = int(ff[i, 1])

        bad.add(p)
        factors.append((p, e))

    print(
        f"{PROTOCOL}|stage=division_field"
        f"|disc_bits={disc.bit_length()}"
        f"|bad_rational_primes={len(bad)}"
        f"|bad={','.join(map(str, sorted(bad)))}",
        flush=True,
    )

    print(
        f"{PROTOCOL}|stage=division_disc"
        f"|factors="
        + ",".join(
            f"{p}^{e}"
            for p, e in factors
        ),
        flush=True,
    )

    # --------------------------------------------------------
    # Known Kummer image
    # --------------------------------------------------------

    alphas = []
    xqs = []

    for index, P in enumerate(SHORT_POINTS, 1):
        xq = sage_q(P[0])

        alpha = (
            qpari(pari, xq)
            - theta
        )

        alphas.append(alpha)
        xqs.append(xq)

    print(
        f"{PROTOCOL}|stage=kummer"
        f"|classes={len(alphas)}",
        flush=True,
    )

    # --------------------------------------------------------
    # Odd bad-prime local square classes
    # --------------------------------------------------------

    rows = [
        []
        for _ in alphas
    ]

    odd_places = 0
    local_coordinates = []

    for p in sorted(bad):
        if p == 2:
            continue

        start = time.monotonic()

        extra, places = prime_local_rows(
            pari,
            nf,
            alphas,
            p,
        )

        rows = append_columns(
            rows,
            extra,
        )

        odd_places += len(places)
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

        print(
            f"{PROTOCOL}|stage=bad_prime"
            f"|p={p}"
            f"|places={len(places)}"
            f"|rank={f2_rank(rows)}"
            f"|seconds={time.monotonic()-start:.6f}",
            flush=True,
        )

    rank_odd = f2_rank(rows)

    # --------------------------------------------------------
    # 2-adic local square classes
    # --------------------------------------------------------

    two_primes = list(
        pari.idealprimedec(
            nf,
            2,
        )
    )

    start = time.monotonic()

    two_basis, two_basis_origins, tw = two_adic_coords(
        pari,
        nf,
        two_primes,
        alphas,
    )

    rows = append_columns(
        rows,
        tw,
    )
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

    print(
        f"{PROTOCOL}|stage=two_adic"
        f"|places={len(two_primes)}"
        f"|gain={rank_odd2-rank_odd}"
        f"|rank={rank_odd2}"
        f"|seconds={time.monotonic()-start:.6f}",
        flush=True,
    )

    # --------------------------------------------------------
    # Archimedean signs
    # --------------------------------------------------------

    roots = list(
        pari.polrootsreal(f)
    )

    real_rows = []

    for xq in xqs:
        xx = qpari(
            pari,
            xq,
        )

        real_rows.append(
            [
                1 if xx - root < 0 else 0
                for root in roots
            ]
        )

    rows = append_columns(
        rows,
        real_rows,
    )
    local_coordinates.extend(
        {
            "kind": "real_sign",
            "embedding_index": index,
        }
        for index in range(len(roots))
    )

    baseline_rank = f2_rank(rows)
    local_rows = [list(row) for row in rows]

    print(
        f"{PROTOCOL}|stage=baseline"
        f"|odd_rank={rank_odd}"
        f"|odd_plus_2={rank_odd2}"
        f"|odd_plus_2_plus_real={baseline_rank}"
        f"|real_places={len(roots)}"
        f"|target={KNOWN_RANK}",
        flush=True,
    )

    # --------------------------------------------------------
    # Good-prime witness characters
    #
    # These do NOT add Selmer conditions.
    #
    # They only build a cheap faithful coordinate system for
    # the already-known global Kummer classes.
    # --------------------------------------------------------

    selected = []
    fingerprint_coordinates = []
    current = rows

    for q in prime_range(
        3,
        args.prime_bound + 1,
    ):
        q = int(q)

        if q in bad:
            continue

        old_rank = f2_rank(
            current
        )

        if old_rank >= KNOWN_RANK:
            break

        start = time.monotonic()

        try:
            extra, places = prime_local_rows(
                pari,
                nf,
                alphas,
                q,
            )
        except Exception as exc:
            print(
                f"{PROTOCOL}|stage=aux"
                f"|q={q}"
                f"|status=ERROR"
                f"|error={type(exc).__name__}:{exc}",
                flush=True,
            )
            continue

        trial = append_columns(
            current,
            extra,
        )

        new_rank = f2_rank(
            trial
        )

        if new_rank > old_rank:
            gain = new_rank - old_rank

            selected.append(
                (
                    q,
                    len(places),
                    gain,
                )
            )

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
                f"{PROTOCOL}|stage=aux"
                f"|q={q}"
                f"|status=SELECTED"
                f"|places={len(places)}"
                f"|gain={gain}"
                f"|rank={new_rank}"
                f"|seconds={time.monotonic()-start:.6f}",
                flush=True,
            )

    final_rank = f2_rank(
        current
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|baseline_rank={baseline_rank}"
        f"|final_rank={final_rank}"
        f"|target={KNOWN_RANK}"
        f"|selected={selected}"
        f"|prime_bound={args.prime_bound}",
        flush=True,
    )

    if final_rank == KNOWN_RANK:
        print(
            f"{PROTOCOL}|result="
            "FAITHFUL_KNOWN_KUMMER_FINGERPRINT"
            f"|dimensions={KNOWN_RANK}"
            f"|auxiliary_primes="
            f"{[q for q,_,_ in selected]}",
            flush=True,
        )
    else:
        print(
            f"{PROTOCOL}|result="
            "INCOMPLETE_FINGERPRINT"
            f"|missing={KNOWN_RANK-final_rank}"
            f"|action=raise_prime_bound",
            flush=True,
        )

    fingerprint_width = len(fingerprint_coordinates)
    if any(len(row) != len(local_coordinates) + fingerprint_width for row in current):
        raise ArithmeticError("signature-coordinate bookkeeping lost alignment")
    signature_map = {
        "schema": "elliptic-curves.bnf-free-signature-map.v1",
        "status": "known_kummer_image_only_not_a_selmer_bound",
        "field_generator": "theta",
        "generator_coordinate_order": ["1", "theta", "theta^2"],
        "defining_polynomial_ascending": [str(B), str(A), "0", "1"],
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
