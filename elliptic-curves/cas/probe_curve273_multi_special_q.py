#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients


PROTOCOL = "R30MULTIQ"

BAD = (
    2, 3, 5, 7, 13, 31, 41, 47, 53, 67,
    379, 4349, 25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
)


def sage_q(QQ, ZZ, value):
    n = value.numerator
    d = value.denominator

    if callable(n):
        n = n()
    if callable(d):
        d = d()

    return QQ(ZZ(n)) / QQ(ZZ(d))


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--special-q",
        type=int,
        action="append",
        required=True,
        help="repeat for every forced rational prime",
    )

    ap.add_argument(
        "--seed-m",
        type=int,
        required=True,
    )

    ap.add_argument(
        "--half-width",
        type=int,
        default=1_000_000,
    )

    ap.add_argument(
        "--factor-base-bound",
        type=int,
        default=1_000_000,
    )

    ap.add_argument(
        "--top",
        type=int,
        default=500,
    )

    ap.add_argument(
        "--print-best",
        type=int,
        default=25,
    )

    ap.add_argument(
        "--factor-max-bits",
        type=int,
        default=140,
    )

    args = ap.parse_args()

    from sage.all import (
        GF,
        QQ,
        ZZ,
        PolynomialRing,
        RealField,
        gcd,
        prime_range,
    )

    coeffs = short_coefficients()

    A = ZZ(sage_q(QQ, ZZ, coeffs[3]))
    B = ZZ(sage_q(QQ, ZZ, coeffs[4]))

    R = PolynomialRing(ZZ, "x")
    x = R.gen()

    f = x**3 + A*x + B

    seed_m = ZZ(args.seed_m)

    special = tuple(
        ZZ(q)
        for q in args.special_q
    )

    if len(set(special)) != len(special):
        raise SystemExit("special q values must be distinct")

    for q in special:
        if q <= args.factor_base_bound:
            raise SystemExit(
                f"special q={q} lies inside ordinary factor base"
            )

        if not q.is_prime(proof=True):
            raise SystemExit(
                f"special q={q} is not prime"
            )

        if f(seed_m) % q:
            raise SystemExit(
                f"special q={q} does not divide f(seed_m)"
            )

    modulus = ZZ(1)

    for q in special:
        modulus *= q

    # Because seed_m satisfies every congruence, the simultaneous
    # CRT class is simply seed_m mod product(q_i).
    residue = seed_m % modulus

    RF = RealField(256)

    roots = f.change_ring(RF).roots()

    if len(roots) != 1:
        raise RuntimeError(
            f"expected one real root, got {len(roots)}"
        )

    real_root = roots[0][0]

    center_k = ZZ(
        (
            (real_root - RF(residue))
            / RF(modulus)
        ).round()
    )

    center_m = (
        residue
        + modulus*center_k
    )

    seed_k = (
        seed_m - residue
    ) // modulus

    if f(center_m) % modulus:
        raise AssertionError(
            "combined special-q congruence failed"
        )

    quotient0 = abs(
        ZZ(f(center_m) // modulus)
    )

    print(
        f"{PROTOCOL}|stage=input"
        f"|special={','.join(map(str,special))}"
        f"|count={len(special)}"
        f"|modulus_bits={modulus.nbits()}"
        f"|seed_m={seed_m}"
        f"|center_m={center_m}"
        f"|seed_delta_k={seed_k-center_k}"
        f"|center_quotient_bits={quotient0.nbits()}"
        f"|half_width={args.half_width}"
        f"|fb_bound={args.factor_base_bound}",
        flush=True,
    )

    M = args.half_width
    length = 2*M + 1

    base_k = int(
        center_k - M
    )

    scores = np.zeros(
        length,
        dtype=np.float32,
    )

    primes = [
        int(p)
        for p in prime_range(
            2,
            args.factor_base_bound + 1,
        )
    ]

    special_set = {
        int(q)
        for q in special
    }

    primes = [
        p
        for p in primes
        if p not in special_set
    ]

    print(
        f"{PROTOCOL}|stage=sieve"
        f"|points={length}"
        f"|primes={len(primes)}"
        f"|status=START",
        flush=True,
    )

    t0 = time.monotonic()

    residue_int = int(residue)
    modulus_int = int(modulus)

    roots_seen = 0

    for pi, p in enumerate(
        primes,
        1,
    ):

        Fp = GF(p)

        roots_mod_p = [
            int(root)
            for root, multiplicity
            in f.change_ring(Fp).roots()
        ]

        if not roots_mod_p:
            continue

        Minv = pow(
            modulus_int % p,
            -1,
            p,
        )

        base_mod = (
            base_k % p
        )

        weight = np.float32(
            math.log2(p)
        )

        for root in roots_mod_p:

            kres = (
                (root - residue_int)
                * Minv
            ) % p

            start = (
                kres - base_mod
            ) % p

            scores[start::p] += weight

            roots_seen += 1

        if pi % 10000 == 0:

            print(
                f"{PROTOCOL}|stage=sieve_progress"
                f"|prime_index={pi}/{len(primes)}"
                f"|p={p}"
                f"|roots_seen={roots_seen}"
                f"|seconds={time.monotonic()-t0:.2f}",
                flush=True,
            )

    # --------------------------------------------------------
    # Approximate log2(|f(m)| / product(q)).
    #
    # m = center_m + modulus*d
    # --------------------------------------------------------

    d = np.arange(
        -M,
        M + 1,
        dtype=np.float64,
    )

    m0 = ZZ(center_m)

    c0 = ZZ(
        f(m0) // modulus
    )

    c1 = ZZ(
        3*m0*m0 + A
    )

    c2 = ZZ(
        3*m0*modulus
    )

    c3 = ZZ(
        modulus*modulus
    )

    approximate = np.abs(
        float(c0)
        + float(c1)*d
        + float(c2)*d*d
        + float(c3)*d*d*d
    )

    residual_estimate = (
        np.log2(approximate)
        - scores
    )

    top = min(
        args.top,
        length,
    )

    indices = np.argpartition(
        residual_estimate,
        top - 1,
    )[:top]

    indices = sorted(
        indices,
        key=lambda i: (
            float(residual_estimate[i]),
            int(i),
        ),
    )

    print(
        f"{PROTOCOL}|stage=exact"
        f"|candidates={len(indices)}",
        flush=True,
    )

    records = []

    full_smooth = 0
    le32 = 0
    le48 = 0
    le64 = 0
    le80 = 0
    le96 = 0

    all_factors_le24 = 0
    all_factors_le28 = 0
    all_factors_le32 = 0
    all_factors_le36 = 0
    all_factors_le40 = 0

    for candidate_index in indices:

        delta = (
            int(candidate_index)
            - M
        )

        k = (
            center_k + delta
        )

        m = (
            residue
            + modulus*k
        )

        N = abs(
            ZZ(f(m))
        )

        if N % modulus:
            raise AssertionError(
                "candidate escaped combined special-q class"
            )

        quotient = (
            N // modulus
        )

        co = quotient

        small_factors = []

        for p in primes:

            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            small_factors.append(
                (p, exponent)
            )

            if co == 1:
                break

        # S-primes above the normal small FB range are irrelevant
        # in the eventual S-class quotient, so strip them too.
        for p0 in BAD:

            p = ZZ(p0)

            if p <= args.factor_base_bound:
                continue

            if p in special:
                continue

            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            small_factors.append(
                (int(p), exponent)
            )

        co_bits = (
            0
            if co == 1
            else int(co.nbits())
        )

        if co == 1:
            full_smooth += 1

        if co_bits <= 32:
            le32 += 1

        if co_bits <= 48:
            le48 += 1

        if co_bits <= 64:
            le64 += 1

        if co_bits <= 80:
            le80 += 1

        if co_bits <= 96:
            le96 += 1

        factors = None
        max_factor_bits = None

        if (
            co not in (0, 1)
            and
            co_bits <= args.factor_max_bits
        ):

            try:
                ff = list(
                    co.factor(
                        proof=False
                    )
                )

                # Upgrade the heuristic factorization to an exact
                # one by proving every returned factor prime.
                certified = True

                for p, exponent in ff:

                    if not ZZ(p).is_prime(
                        proof=True
                    ):
                        certified = False
                        break

                if certified:

                    product_check = ZZ(1)

                    for p, exponent in ff:
                        product_check *= (
                            ZZ(p) ** int(exponent)
                        )

                    if product_check != co:
                        raise AssertionError(
                            "factorization product mismatch"
                        )

                    factors = [
                        (
                            int(p),
                            int(exponent),
                        )
                        for p, exponent in ff
                    ]

                    max_factor_bits = max(
                        int(ZZ(p).nbits())
                        for p, exponent
                        in ff
                    )

                    if max_factor_bits <= 24:
                        all_factors_le24 += 1

                    if max_factor_bits <= 28:
                        all_factors_le28 += 1

                    if max_factor_bits <= 32:
                        all_factors_le32 += 1

                    if max_factor_bits <= 36:
                        all_factors_le36 += 1

                    if max_factor_bits <= 40:
                        all_factors_le40 += 1

            except Exception:
                pass

        records.append(
            (
                co_bits,
                (
                    999
                    if max_factor_bits is None
                    else max_factor_bits
                ),
                (
                    999
                    if factors is None
                    else len(factors)
                ),
                delta,
                int(m),
                int(quotient.nbits()),
                int(co),
                factors,
                float(
                    residual_estimate[
                        candidate_index
                    ]
                ),
            )
        )

    factor_records = sorted(
        (
            row
            for row in records
            if row[7] is not None
        ),
        key=lambda row: (
            row[1],  # max residual prime bits
            row[0],  # total residual bits
            row[2],  # residual factor count
            abs(row[3]),
        ),
    )

    records.sort(
        key=lambda row: (
            row[0],
            row[1],
            abs(row[3]),
        )
    )

    best_factor = (
        factor_records[0]
        if factor_records
        else None
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|special={','.join(map(str,special))}"
        f"|full_smooth={full_smooth}"
        f"|cofactor_le32={le32}"
        f"|cofactor_le48={le48}"
        f"|cofactor_le64={le64}"
        f"|cofactor_le80={le80}"
        f"|cofactor_le96={le96}"
        f"|all_factors_le24={all_factors_le24}"
        f"|all_factors_le28={all_factors_le28}"
        f"|all_factors_le32={all_factors_le32}"
        f"|all_factors_le36={all_factors_le36}"
        f"|all_factors_le40={all_factors_le40}"
        f"|best_cofactor_bits={records[0][0]}"
        f"|best_largest_factor_bits="
        f"{'NA' if best_factor is None else best_factor[1]}"
        f"|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )

    for rank, row in enumerate(
        factor_records[
            :args.print_best
        ],
        1,
    ):

        (
            co_bits,
            max_bits,
            factor_count,
            delta,
            m,
            quotient_bits,
            co,
            factors,
            estimate,
        ) = row

        factor_text = "*".join(
            (
                str(p)
                if exponent == 1
                else f"{p}^{exponent}"
            )
            for p, exponent
            in factors
        )

        print(
            f"{PROTOCOL}|factorbest={rank}"
            f"|delta_k={delta}"
            f"|m={m}"
            f"|quotient_bits={quotient_bits}"
            f"|cofactor_bits={co_bits}"
            f"|largest_factor_bits={max_bits}"
            f"|factor_count={factor_count}"
            f"|factorization={factor_text}"
            f"|estimate={estimate:.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
