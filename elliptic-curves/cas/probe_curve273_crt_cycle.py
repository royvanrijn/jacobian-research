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


PROTOCOL = "R30CRT"

S_RATIONAL = {
    2, 3, 5, 7, 13, 31, 41, 47, 53, 67,
    379, 4349, 25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
}


def sage_q(QQ, ZZ, value):
    n = value.numerator
    d = value.denominator

    if callable(n):
        n = n()
    if callable(d):
        d = d()

    return QQ(ZZ(n)) / QQ(ZZ(d))


def parse_ideal(text):
    q, residue = text.split(":", 1)
    return int(q), int(residue)


def crt_pairs(pairs):
    """
    CRT for pairwise-coprime q values.

    Return r,M with r == residue mod q for every pair.
    """
    r = 0
    M = 1

    for q, residue in pairs:
        residue %= q

        if math.gcd(M, q) != 1:
            raise ValueError(
                f"CRT moduli not coprime: M={M}, q={q}"
            )

        t = (
            (residue - r)
            * pow(M % q, -1, q)
        ) % q

        r += M*t
        M *= q
        r %= M

    return r, M


def support_text(support):
    return ",".join(
        f"{q}:{r}"
        for q, r in sorted(support)
    )


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--force",
        action="append",
        required=True,
        help="prime-ideal coordinate q:root",
    )

    ap.add_argument(
        "--target",
        action="append",
        required=True,
        help="LP support q:root to close",
    )

    ap.add_argument(
        "--factor-base-bound",
        type=int,
        default=1_000_000,
    )

    ap.add_argument(
        "--half-width",
        type=int,
        default=1_000_000,
    )

    ap.add_argument(
        "--top",
        type=int,
        default=1000,
    )

    ap.add_argument(
        "--print-best",
        type=int,
        default=25,
    )

    ap.add_argument(
        "--factor-max-bits",
        type=int,
        default=160,
    )

    args = ap.parse_args()

    from sage.all import (
        GF,
        QQ,
        ZZ,
        PolynomialRing,
        RealField,
        prime_range,
    )

    forced = tuple(
        parse_ideal(x)
        for x in args.force
    )

    target = frozenset(
        parse_ideal(x)
        for x in args.target
    )

    if len(set(q for q, r in forced)) != len(forced):
        raise SystemExit("forced rational primes must be distinct")

    coeffs = short_coefficients()

    A = ZZ(
        sage_q(QQ, ZZ, coeffs[3])
    )

    B = ZZ(
        sage_q(QQ, ZZ, coeffs[4])
    )

    R = PolynomialRing(ZZ, "x")
    x = R.gen()

    f = x**3 + A*x + B

    # --------------------------------------------------------
    # Validate actual degree-one prime ideals.
    # --------------------------------------------------------

    for q, residue in forced:

        if q <= args.factor_base_bound:
            raise SystemExit(
                f"forced q={q} lies in ordinary factor base"
            )

        if not ZZ(q).is_prime(proof=True):
            raise SystemExit(
                f"forced q={q} is not prime"
            )

        if f(ZZ(residue)) % q:
            raise SystemExit(
                f"{residue} is not a root mod {q}"
            )

    for q, residue in target:

        if not ZZ(q).is_prime(proof=True):
            raise SystemExit(
                f"target q={q} is not prime"
            )

        if f(ZZ(residue)) % q:
            raise SystemExit(
                f"target {q}:{residue} is not a root"
            )

    residue, modulus = crt_pairs(
        forced
    )

    modulus = ZZ(modulus)
    residue = ZZ(residue)

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

    for q, r in forced:

        if center_m % q != r % q:
            raise AssertionError(
                "CRT center mismatch"
            )

    quotient0 = abs(
        ZZ(f(center_m) // modulus)
    )

    print(
        f"{PROTOCOL}|stage=input"
        f"|forced={support_text(forced)}"
        f"|target={support_text(target)}"
        f"|forced_count={len(forced)}"
        f"|target_count={len(target)}"
        f"|modulus_bits={modulus.nbits()}"
        f"|center_m={center_m}"
        f"|center_quotient_bits={quotient0.nbits()}"
        f"|half_width={args.half_width}"
        f"|fb_bound={args.factor_base_bound}",
        flush=True,
    )

    # --------------------------------------------------------
    # Sieve arithmetic progression:
    #
    #   m = residue + modulus*k
    # --------------------------------------------------------

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

    forced_q = {
        q
        for q, r in forced
    }

    primes = [
        p
        for p in primes
        if p not in forced_q
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

        roots_mod_p = [
            int(root)
            for root, multiplicity
            in f.change_ring(GF(p)).roots()
        ]

        if not roots_mod_p:
            continue

        inv = pow(
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
                * inv
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
    # Approximate log2 |f(m)/modulus|
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

    estimate = (
        np.log2(approximate)
        - scores
    )

    top = min(
        args.top,
        length,
    )

    indices = np.argpartition(
        estimate,
        top - 1,
    )[:top]

    indices = sorted(
        indices,
        key=lambda i: (
            float(estimate[i]),
            int(i),
        ),
    )

    print(
        f"{PROTOCOL}|stage=exact"
        f"|candidates={len(indices)}",
        flush=True,
    )

    records = []

    fully_closed = 0
    improved = 0

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

        # Every candidate is constructed with
        #
        #     m == r_i (mod q_i)
        #
        # for each forced degree-one prime ideal, hence one copy
        # of every forced q_i is known to divide f(m).
        #
        # Remove that declared special-q modulus BEFORE testing
        # residual smoothness.  Extra powers of a forced q remain
        # in the quotient and are handled normally below.
        if N % modulus:
            raise AssertionError(
                "candidate escaped forced CRT divisibility"
            )

        # ----------------------------------------------------
        # Strip ordinary factor base and S support from the
        # quotient after the declared forced factors.
        # ----------------------------------------------------

        co = N // modulus

        for p in primes:

            if co % p:
                continue

            while co % p == 0:
                co //= p

            if co == 1:
                break

        # S primes > ordinary FB.
        for p0 in S_RATIONAL:

            p = ZZ(p0)

            if p <= args.factor_base_bound:
                continue

            if co % p:
                continue

            while co % p == 0:
                co //= p

        # ----------------------------------------------------
        # Factor the remaining LP part.
        # ----------------------------------------------------

        if co == 1:
            ff = []
        elif co.nbits() <= args.factor_max_bits:

            try:
                raw = list(
                    co.factor(
                        proof=False
                    )
                )
            except Exception:
                continue

            # Certify result exactly.
            product_check = ZZ(1)

            good = True

            for p, exponent in raw:

                if not ZZ(p).is_prime(
                    proof=True
                ):
                    good = False
                    break

                product_check *= (
                    ZZ(p) ** int(exponent)
                )

            if not good or product_check != co:
                continue

            ff = [
                (
                    int(p),
                    int(exponent),
                )
                for p, exponent
                in raw
            ]

        else:
            continue

        # ----------------------------------------------------
        # Build candidate LP PRIME-IDEAL support.
        #
        # For alpha=m-theta, a degree-one prime ideal over p is
        # identified by the root theta=m mod p.
        #
        # IMPORTANT: use parity of the COMPLETE valuation in N,
        # including forced q values.
        # ----------------------------------------------------

        candidate_lp = set()

        all_large_q = {
            p
            for p, exponent in ff
            if p > args.factor_base_bound
            and p not in S_RATIONAL
        }

        all_large_q.update(
            q
            for q, r in forced
            if q > args.factor_base_bound
            and q not in S_RATIONAL
        )

        for q in sorted(all_large_q):

            qz = ZZ(q)

            tmp = N
            valuation = 0

            while tmp % qz == 0:
                tmp //= qz
                valuation += 1

            if valuation & 1:
                candidate_lp.add(
                    (
                        int(q),
                        int(m % qz),
                    )
                )

        candidate_lp = frozenset(
            candidate_lp
        )

        cycle = (
            target
            ^ candidate_lp
        )

        cycle_count = len(
            cycle
        )

        max_bits = (
            max(
                q.bit_length()
                for q, r in cycle
            )
            if cycle
            else 0
        )

        total_bits = sum(
            q.bit_length()
            for q, r in cycle
        )

        overlap = len(
            target & candidate_lp
        )

        if cycle_count < len(target):
            improved += 1

        if cycle_count == 0:
            fully_closed += 1

        records.append(
            {
                "cycle_count": cycle_count,
                "max_bits": max_bits,
                "total_bits": total_bits,
                "overlap": overlap,
                "candidate_count": len(
                    candidate_lp
                ),
                "m": int(m),
                "delta": delta,
                "candidate": candidate_lp,
                "cycle": cycle,
                "cofactor_bits": (
                    0
                    if co == 1
                    else int(co.nbits())
                ),
                "estimate": float(
                    estimate[candidate_index]
                ),
            }
        )

    records.sort(
        key=lambda r: (
            r["cycle_count"],
            r["max_bits"],
            r["total_bits"],
            -r["overlap"],
            abs(r["delta"]),
        )
    )

    print(
        f"{PROTOCOL}|stage=summary"
        f"|forced={support_text(forced)}"
        f"|tested={len(records)}"
        f"|improved={improved}"
        f"|fully_closed={fully_closed}"
        f"|best_cycle_remaining="
        f"{records[0]['cycle_count'] if records else 'NA'}"
        f"|best_max_bits="
        f"{records[0]['max_bits'] if records else 'NA'}"
        f"|seconds={time.monotonic()-t0:.3f}",
        flush=True,
    )

    for rank, record in enumerate(
        records[:args.print_best],
        1,
    ):

        print(
            f"{PROTOCOL}|cyclebest={rank}"
            f"|m={record['m']}"
            f"|delta_k={record['delta']}"
            f"|overlap={record['overlap']}"
            f"|candidate_lp={record['candidate_count']}"
            f"|cycle_remaining={record['cycle_count']}"
            f"|max_bits={record['max_bits']}"
            f"|cofactor_bits={record['cofactor_bits']}"
            f"|candidate={support_text(record['candidate'])}"
            f"|remaining={support_text(record['cycle'])}"
            f"|estimate={record['estimate']:.3f}",
            flush=True,
        )

        if record["cycle_count"] == 0:

            print(
                f"{PROTOCOL}|CYCLE"
                f"|status=EXACT_LP_CLOSURE"
                f"|m={record['m']}",
                flush=True,
            )


if __name__ == "__main__":
    main()
