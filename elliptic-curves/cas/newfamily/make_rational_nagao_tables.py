#!/usr/bin/env sage -python
"""Build projective local Nagao tables for rational T=a/b.

This consumes the recovered finite-minimal newfamily builder from the exploratory
archive.  For a short model

    y^2 = x^3 + A(T) x + B(T),  deg A <= 8, deg B <= 12,

a rational parameter T=a/b is represented modulo p by the projective point
(a:b).  If b != 0, its local trace equals the trace at t=a/b mod p; b=0 is the
point at infinity and uses the leading homogeneous coefficients.

The output is intentionally family-independent and is consumed by
scan_rational_nagao_tables.cpp.
"""

from __future__ import annotations

import argparse
from math import log
from pathlib import Path
import sys

from sage.all import GF, ZZ

ROOTS = (-47, -43, -31, 30, 45, 46)

DISCOVERY_PRIMES = (
    401, 409, 419, 421, 431, 433, 439, 443, 449,
    457, 461, 463, 467, 479, 487, 491, 499,
)
HELD_PRIMES = (
    503, 509, 521, 523, 541, 547, 557,
    563, 569, 571, 577, 587, 593, 599,
)


def load_builder():
    """Import the recovered builder from either canonical or archive location."""
    candidates = [
        Path("elliptic-curves/cas/newfamily/newfamily_rank11_minimal_common.py"),
        Path("elliptic-curves/cas/newfamily/archive/newfamily_rank11_minimal_common.py"),
        Path("/tmp/newfamily_rank11_minimal_common.py"),
    ]
    for path in candidates:
        if path.exists():
            sys.path.insert(0, str(path.parent.resolve()))
            from newfamily_rank11_minimal_common import build_finite_minimal_family
            return build_finite_minimal_family
    raise SystemExit(
        "missing newfamily_rank11_minimal_common.py; preserve the recovered /tmp "
        "scripts under elliptic-curves/cas/newfamily/archive first"
    )


def trace_of_frobenius(a: int, b: int, p: int) -> int:
    total = 0
    for x in range(p):
        rhs = (x * x % p * x + a * x + b) % p
        if rhs:
            ls = pow(rhs, (p - 1) // 2, p)
            total += 1 if ls == 1 else -1
    return -total


def polynomial_mod_value(poly, residue: int, p: int) -> int:
    answer = 0
    for coefficient in reversed(poly.list()):
        answer = (answer * residue + int(ZZ(coefficient) % p)) % p
    return answer


def local_symbol(A, B, residue: int | None, p: int) -> tuple[int, int, int]:
    if residue is None:
        a = int(ZZ(A[A.degree()]) % p) if A.degree() == 8 else 0
        b = int(ZZ(B[B.degree()]) % p) if B.degree() == 12 else 0
    else:
        a = polynomial_mod_value(A, residue, p)
        b = polynomial_mod_value(B, residue, p)
    if (4 * a * a % p * a + 27 * b * b) % p == 0:
        return 0, 0, 0
    trace = trace_of_frobenius(a, b, p)
    order = p + 1 - trace
    score = (2.0 - trace) / order * log(float(p))
    units = int(round(score * 1_000_000_000_000.0))
    return 1, trace, units


def write_band(handle, label: str, primes, A, B) -> None:
    handle.write(f"B {label} {len(primes)}\n")
    for p in primes:
        handle.write(f"P {p}\n")
        for residue in range(p):
            good, trace, units = local_symbol(A, B, residue, p)
            handle.write(f"{good} {trace} {units}\n")
        good, trace, units = local_symbol(A, B, None, p)
        handle.write(f"{good} {trace} {units}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    build = load_builder()
    family = build(ROOTS)
    A = family["Amin"]
    B = family["Bmin"]
    if A.degree() > 8 or B.degree() > 12:
        raise SystemExit(f"unexpected finite-minimal degrees {A.degree()}, {B.degree()}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write("RATIONAL_NAGAO_LOCAL_TABLE_V1\n")
        handle.write(f"F NEWFAMILY_RANK11 {A.degree()} {B.degree()}\n")
        write_band(handle, "D", DISCOVERY_PRIMES, A, B)
        write_band(handle, "H", HELD_PRIMES, A, B)
        handle.write("END\n")

    print(f"A_degree={A.degree()} B_degree={B.degree()}")
    print(f"discovery_primes={len(DISCOVERY_PRIMES)} held_primes={len(HELD_PRIMES)}")
    print(f"saved={out}")
    print("DONE")


if __name__ == "__main__":
    main()
