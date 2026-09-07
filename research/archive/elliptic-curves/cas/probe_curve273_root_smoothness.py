#!/usr/bin/env python3

from pathlib import Path
import argparse
import math
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients


PROTOCOL = "R30ROOTPROBE"


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

    ap.add_argument("--half-width", type=int, default=500_000)
    ap.add_argument("--factor-base-bound", type=int, default=100_000)
    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("--print-best", type=int, default=30)

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

    RF = RealField(256)
    roots = f.change_ring(RF).roots()

    real_roots = [
        root
        for root, multiplicity in roots
    ]

    if len(real_roots) != 1:
        raise RuntimeError(
            f"expected one real root, found {len(real_roots)}"
        )

    root = real_roots[0]
    center = ZZ(root.round())

    f0 = ZZ(f(center))
    derivative = ZZ(3*center*center + A)

    print(
        f"{PROTOCOL}|stage=input"
        f"|center={center}"
        f"|center_norm_bits={abs(f0).nbits()}"
        f"|derivative_bits={abs(derivative).nbits()}"
        f"|half_width={args.half_width}"
        f"|fb_bound={args.factor_base_bound}",
        flush=True,
    )

    M = args.half_width
    length = 2*M + 1

    # index i corresponds to m = center-M+i
    base = int(center - M)

    scores = np.zeros(length, dtype=np.float32)

    primes = [
        int(p)
        for p in prime_range(
            2,
            args.factor_base_bound + 1,
        )
    ]

    print(
        f"{PROTOCOL}|stage=sieve"
        f"|points={length}"
        f"|primes={len(primes)}"
        f"|status=START",
        flush=True,
    )

    t0 = time.monotonic()

    root_count = 0

    for index, p in enumerate(primes, 1):

        Fp = GF(p)
        fp = f.change_ring(Fp)

        roots_mod_p = [
            int(r)
            for r, multiplicity in fp.roots()
        ]

        if not roots_mod_p:
            continue

        root_count += len(roots_mod_p)

        weight = np.float32(math.log2(p))

        for r in roots_mod_p:
            start = (r - (base % p)) % p
            scores[start::p] += weight

        if index % 2000 == 0:
            print(
                f"{PROTOCOL}|stage=sieve_progress"
                f"|prime_index={index}/{len(primes)}"
                f"|p={p}"
                f"|roots_seen={root_count}"
                f"|seconds={time.monotonic()-t0:.2f}",
                flush=True,
            )

    # Approximate log2 |f(center+n)| using double precision.
    # Offsets are tiny compared to center, so this is stable enough
    # for ranking; shortlisted candidates are checked exactly below.
    offsets = np.arange(-M, M + 1, dtype=np.float64)

    c = float(center)
    f0_float = float(f0)
    d_float = float(derivative)

    approximate_values = np.abs(
        f0_float
        + d_float*offsets
        + 3.0*c*offsets*offsets
        + offsets*offsets*offsets
    )

    norm_logs = np.log2(approximate_values)

    residual_estimate = norm_logs - scores

    top = min(args.top, length)

    candidates = np.argpartition(
        residual_estimate,
        top - 1,
    )[:top]

    candidates = sorted(
        candidates,
        key=lambda i: (
            float(residual_estimate[i]),
            int(i),
        ),
    )

    print(
        f"{PROTOCOL}|stage=exact"
        f"|candidates={len(candidates)}",
        flush=True,
    )

    best = []

    full_smooth = 0
    prime_cofactors = 0
    cofactor_le_64 = 0
    cofactor_le_96 = 0
    cofactor_le_128 = 0

    for candidate_index in candidates:

        offset = int(candidate_index) - M
        m = center + offset

        N = abs(ZZ(f(m)))

        cofactor = N
        used = []

        for p in primes:

            if cofactor % p:
                continue

            exponent = 0

            while cofactor % p == 0:
                cofactor //= p
                exponent += 1

            used.append(
                (p, exponent)
            )

            if cofactor == 1:
                break

        bits = (
            0
            if cofactor == 1
            else int(cofactor.nbits())
        )

        is_prime = (
            False
            if cofactor in (0, 1)
            else bool(cofactor.is_prime(proof=False))
        )

        if cofactor == 1:
            full_smooth += 1

        if is_prime:
            prime_cofactors += 1

        if bits <= 64:
            cofactor_le_64 += 1

        if bits <= 96:
            cofactor_le_96 += 1

        if bits <= 128:
            cofactor_le_128 += 1

        best.append(
            (
                bits,
                int(N.nbits()),
                offset,
                int(m),
                int(cofactor),
                is_prime,
                used,
                float(residual_estimate[candidate_index]),
            )
        )

    best.sort(
        key=lambda row: (
            row[0],
            row[1],
            abs(row[2]),
        )
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|full_smooth={full_smooth}"
        f"|prime_cofactor={prime_cofactors}"
        f"|cofactor_le64={cofactor_le_64}"
        f"|cofactor_le96={cofactor_le_96}"
        f"|cofactor_le128={cofactor_le_128}"
        f"|best_cofactor_bits={best[0][0]}"
        f"|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )

    for rank, row in enumerate(
        best[:args.print_best],
        1,
    ):
        (
            cofactor_bits,
            norm_bits,
            offset,
            m,
            cofactor,
            is_prime,
            used,
            estimated,
        ) = row

        factors = ",".join(
            f"{p}^{e}"
            for p, e in used
        )

        print(
            f"{PROTOCOL}|best={rank}"
            f"|offset={offset}"
            f"|m={m}"
            f"|norm_bits={norm_bits}"
            f"|cofactor_bits={cofactor_bits}"
            f"|cofactor_prime={int(is_prime)}"
            f"|estimated_residual={estimated:.3f}"
            f"|small_factors={factors}"
            f"|cofactor={cofactor}",
            flush=True,
        )


if __name__ == "__main__":
    main()
