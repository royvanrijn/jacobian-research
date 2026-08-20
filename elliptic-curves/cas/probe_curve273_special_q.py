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


PROTOCOL = "R30SPECIALQ"

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

    ap.add_argument("--special-q", type=int, required=True)
    ap.add_argument("--seed-m", type=int, required=True)

    ap.add_argument("--half-width", type=int, default=1_000_000)
    ap.add_argument("--factor-base-bound", type=int, default=1_000_000)

    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("--print-best", type=int, default=20)
    ap.add_argument("--factor-max-bits", type=int, default=128)

    args = ap.parse_args()

    from sage.all import (
        GF,
        QQ,
        ZZ,
        PolynomialRing,
        RealField,
        prime_range,
    )

    coeffs = short_coefficients()

    A = ZZ(sage_q(QQ, ZZ, coeffs[3]))
    B = ZZ(sage_q(QQ, ZZ, coeffs[4]))

    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    f = x**3 + A*x + B

    q = ZZ(args.special_q)
    seed_m = ZZ(args.seed_m)

    if q <= args.factor_base_bound:
        raise SystemExit(
            "special q should be outside the ordinary factor base"
        )

    seed_N = ZZ(f(seed_m))

    if seed_N % q != 0:
        raise SystemExit(
            f"q={q} does not divide f(seed_m)"
        )

    # m = residue + q*k
    residue = seed_m % q
    seed_k = (seed_m - residue) // q

    RF = RealField(256)
    roots = f.change_ring(RF).roots()

    if len(roots) != 1:
        raise RuntimeError("expected one real root")

    real_root = roots[0][0]

    # Closest k such that residue + q*k is near the real root.
    center_k = ZZ(
        ((real_root - RF(residue)) / RF(q)).round()
    )

    center_m = residue + q*center_k

    if f(center_m) % q != 0:
        raise AssertionError("special-q congruence failed")

    center_quotient = abs(
        ZZ(f(center_m) // q)
    )

    print(
        f"{PROTOCOL}|stage=input"
        f"|q={q}"
        f"|q_bits={q.nbits()}"
        f"|seed_m={seed_m}"
        f"|seed_k={seed_k}"
        f"|center_k={center_k}"
        f"|center_m={center_m}"
        f"|seed_delta_k={seed_k-center_k}"
        f"|quotient_bits={center_quotient.nbits()}"
        f"|half_width={args.half_width}"
        f"|fb_bound={args.factor_base_bound}",
        flush=True,
    )

    M = args.half_width
    length = 2*M + 1

    base_k = int(center_k - M)

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
        if int(p) != int(q)
    ]

    print(
        f"{PROTOCOL}|stage=sieve"
        f"|points={length}"
        f"|primes={len(primes)}"
        f"|status=START",
        flush=True,
    )

    t0 = time.monotonic()

    roots_seen = 0

    residue_int = int(residue)
    q_int = int(q)

    for pi, p in enumerate(primes, 1):

        Fp = GF(p)

        roots_mod_p = [
            int(r)
            for r, mult
            in f.change_ring(Fp).roots()
        ]

        if not roots_mod_p:
            continue

        qmod = q_int % p

        # q is prime and q > FB, so p != q.
        qinv = pow(qmod, -1, p)

        weight = np.float32(
            math.log2(p)
        )

        base_mod = base_k % p

        for root in roots_mod_p:

            # residue + q*k == root mod p
            kres = (
                (root - residue_int)
                * qinv
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
    # Approximate log2 |f(residue+q*k)/q|
    # near center_k.
    # --------------------------------------------------------

    d = np.arange(
        -M,
        M + 1,
        dtype=np.float64,
    )

    m0 = ZZ(center_m)

    # Exact Taylor coefficients for:
    #
    # f(m0 + q*d)/q
    #
    c0 = ZZ(f(m0) // q)
    c1 = ZZ(3*m0*m0 + A)
    c2 = ZZ(3*m0*q)
    c3 = ZZ(q*q)

    approximate = np.abs(
        float(c0)
        + float(c1)*d
        + float(c2)*d*d
        + float(c3)*d*d*d
    )

    norm_logs = np.log2(
        approximate
    )

    residual_estimate = (
        norm_logs - scores
    )

    top = min(
        args.top,
        length,
    )

    candidate_indices = np.argpartition(
        residual_estimate,
        top - 1,
    )[:top]

    candidate_indices = sorted(
        candidate_indices,
        key=lambda i: (
            float(residual_estimate[i]),
            int(i),
        ),
    )

    print(
        f"{PROTOCOL}|stage=exact"
        f"|candidates={len(candidate_indices)}",
        flush=True,
    )

    full = 0

    le64 = 0
    le80 = 0
    le96 = 0
    le128 = 0

    factored = 0

    largest_le32 = 0
    largest_le36 = 0
    largest_le40 = 0
    largest_le48 = 0

    records = []

    # Remove ordinary FB factors first.
    for candidate_index in candidate_indices:

        delta = int(candidate_index) - M

        k = center_k + delta
        m = residue + q*k

        N = abs(
            ZZ(f(m))
        )

        if N % q:
            raise AssertionError(
                "candidate left special-q progression"
            )

        quotient = N // q

        co = quotient
        used = []

        for p in primes:

            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            used.append(
                (p, exponent)
            )

            if co == 1:
                break

        # Also remove huge bad/S primes not in the ordinary FB.
        for p in BAD:

            if p <= args.factor_base_bound:
                continue

            if p == q:
                continue

            if co % p:
                continue

            exponent = 0

            while co % p == 0:
                co //= p
                exponent += 1

            used.append(
                (p, exponent)
            )

        co_bits = (
            0
            if co == 1
            else int(co.nbits())
        )

        if co == 1:
            full += 1

        if co_bits <= 64:
            le64 += 1

        if co_bits <= 80:
            le80 += 1

        if co_bits <= 96:
            le96 += 1

        if co_bits <= 128:
            le128 += 1

        factors = None
        max_factor_bits = None

        if (
            co not in (0, 1)
            and co_bits <= args.factor_max_bits
        ):
            try:
                ff = list(
                    co.factor(
                        proof=False
                    )
                )

                factors = [
                    (int(p), int(e))
                    for p, e in ff
                ]

                max_factor_bits = max(
                    int(ZZ(p).nbits())
                    for p, e in ff
                )

                factored += 1

                if max_factor_bits <= 32:
                    largest_le32 += 1

                if max_factor_bits <= 36:
                    largest_le36 += 1

                if max_factor_bits <= 40:
                    largest_le40 += 1

                if max_factor_bits <= 48:
                    largest_le48 += 1

            except Exception:
                pass

        records.append(
            (
                co_bits,
                999
                if max_factor_bits is None
                else max_factor_bits,
                int(abs(ZZ(quotient)).nbits()),
                delta,
                int(k),
                int(m),
                int(co),
                factors,
                used,
                float(
                    residual_estimate[
                        candidate_index
                    ]
                ),
            )
        )

    records.sort(
        key=lambda row: (
            row[0],
            row[1],
            abs(row[3]),
        )
    )

    # Also rank FACTORED NON-SEED candidates by their largest
    # remaining prime factor.  For recursive special-q descent,
    # this is often much more important than total cofactor size.
    factor_records = sorted(
        (
            row
            for row in records
            if row[3] != 0
            and row[7] is not None
            and row[1] != 999
        ),
        key=lambda row: (
            row[1],       # largest factor bits
            row[0],       # total cofactor bits
            len(row[7]),  # number of residual rational primes
            abs(row[3]),
        ),
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|q={q}"
        f"|full_smooth={full}"
        f"|cofactor_le64={le64}"
        f"|cofactor_le80={le80}"
        f"|cofactor_le96={le96}"
        f"|cofactor_le128={le128}"
        f"|factored={factored}"
        f"|all_factors_le32={largest_le32}"
        f"|all_factors_le36={largest_le36}"
        f"|all_factors_le40={largest_le40}"
        f"|all_factors_le48={largest_le48}"
        f"|best_cofactor_bits={records[0][0]}"
        f"|best_largest_factor_bits={records[0][1]}"
        f"|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )

    for rank, row in enumerate(
        factor_records[:args.print_best],
        1,
    ):
        (
            co_bits,
            max_factor_bits,
            quotient_bits,
            delta,
            k,
            m,
            co,
            factors,
            used,
            estimate,
        ) = row

        factor_text = "*".join(
            (
                str(pp)
                if ee == 1
                else f"{pp}^{ee}"
            )
            for pp, ee in factors
        )

        print(
            f"{PROTOCOL}|factorbest={rank}"
            f"|delta_k={delta}"
            f"|m={m}"
            f"|quotient_bits={quotient_bits}"
            f"|cofactor_bits={co_bits}"
            f"|largest_factor_bits={max_factor_bits}"
            f"|factor_count={len(factors)}"
            f"|factorization={factor_text}"
            f"|estimate={estimate:.3f}",
            flush=True,
        )

    for rank, row in enumerate(
        records[:args.print_best],
        1,
    ):

        (
            co_bits,
            max_factor_bits,
            quotient_bits,
            delta,
            k,
            m,
            co,
            factors,
            used,
            estimate,
        ) = row

        factor_text = (
            "?"
            if factors is None
            else "*".join(
                (
                    str(p)
                    if e == 1
                    else f"{p}^{e}"
                )
                for p, e in factors
            )
        )

        print(
            f"{PROTOCOL}|best={rank}"
            f"|delta_k={delta}"
            f"|m={m}"
            f"|quotient_bits={quotient_bits}"
            f"|cofactor_bits={co_bits}"
            f"|largest_factor_bits={max_factor_bits}"
            f"|factorization={factor_text}"
            f"|estimate={estimate:.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
