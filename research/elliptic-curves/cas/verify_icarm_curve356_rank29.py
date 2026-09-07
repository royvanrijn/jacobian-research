#!/usr/bin/env python3
"""Independently replay the rank-at-least-29 certificate for ICARM curve 356.

The rank proof uses exact rational arithmetic and exhaustive finite quotients
``E(F_p)/2E(F_p)``.  PARI/GP is used separately for the local-reduction and
root-number replay; it is not used to decide point independence.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve356 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
    SHORT_POINTS,
    on_curve,
    on_short_curve,
    short_coefficients,
    weierstrass_invariants,
)
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    finite_curve_points,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


PROTOCOL = "R29ICARM356"
PUBLIC_SOURCE = "https://elliptic-rank.icarm.cloud/curve/356"

PUBLIC_CONDUCTOR = int(
    "1138311211986944120854212482517370029218279069247307312958252982216018500197158165165915709784372000548257520927912312147452126835606642920"
)
PUBLIC_DISCRIMINANT = -int(
    "23225588805821255971980452316400390985539123780447843034585631625599190483310342757460742687943805554853363674964538200142891911461727552231886578226014220310540000000"
)

C4_FACTORIZATION = (
    (2, 4),
    (31, 1),
    (223, 1),
    (8009, 1),
    (8753, 1),
    (938083, 1),
    (160962024890210860516928562131837019247, 1),
)

DISCRIMINANT_FACTORIZATION = (
    (2, 8),
    (3, 10),
    (5, 7),
    (13, 6),
    (23, 2),
    (29, 3),
    (37, 4),
    (41, 2),
    (139, 2),
    (751, 1),
    (28960331, 1),
    (1204882855601765528877267647500895974865482613, 1),
    (197980272243427555346397293722916980361535459279712115031762027678304939, 1),
)

CONDUCTOR_FACTORIZATION = tuple(
    (prime, 3 if prime == 2 else 1)
    for prime, _exponent in DISCRIMINANT_FACTORIZATION
)

EXPECTED_CERTIFICATE_PRIMES = (
    11,
    53,
    73,
    79,
    83,
    89,
    101,
    109,
    127,
    137,
    157,
    163,
    173,
    179,
    191,
    193,
    229,
    251,
    263,
    281,
)

# prime -> (conductor exponent, PARI Kodaira code, Tamagawa number, root sign)
EXPECTED_LOCAL = {
    2: (3, -5, 4, -1),
    3: (1, 14, 10, -1),
    5: (1, 11, 7, -1),
    13: (1, 10, 6, -1),
    23: (1, 6, 2, 1),
    29: (1, 7, 3, -1),
    37: (1, 8, 4, -1),
    41: (1, 6, 2, 1),
    139: (1, 6, 2, -1),
    751: (1, 5, 1, -1),
    28960331: (1, 5, 1, -1),
    1204882855601765528877267647500895974865482613: (1, 5, 1, 1),
    197980272243427555346397293722916980361535459279712115031762027678304939: (1, 5, 1, -1),
}

LOCAL_PATTERN = re.compile(
    r"^LOCAL\|(\d+)\|(\d+)\|(-?\d+)\|\[([^]]+)\]\|(\d+)\|(-?\d+)$",
    re.MULTILINE,
)


def factor_product(factors: tuple[tuple[int, int], ...]) -> int:
    answer = 1
    for prime, exponent in factors:
        answer *= prime**exponent
    return answer


def valuation(value: int, prime: int) -> int:
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def gp_rational(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def finite_group_order(prime: int) -> int:
    coefficients = short_coefficients()
    coefficient_a = int(coefficients[3]) % prime
    coefficient_b = int(coefficients[4]) % prime
    return len(finite_curve_points(coefficient_a, coefficient_b, prime))


def local_reduction_replay() -> tuple[str, int]:
    ainvs = ",".join(gp_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS)
    commands = ["default(realprecision,80);", f"E=ellinit([{ainvs}]);"]
    for prime in EXPECTED_LOCAL:
        commands.extend(
            [
                f"L=elllocalred(E,{prime});",
                (
                    f'print("LOCAL|{prime}|",L[1],"|",L[2],"|",L[3],"|",'
                    f'L[4],"|",ellrootno(E,{prime}));'
                ),
            ]
        )
    commands.extend(
        [
            'print("MINIMAL|",Vec(ellminimalmodel(E))[1..5]);',
            'print("PARI|",version());',
            "quit",
        ]
    )
    completed = subprocess.run(
        ["gp", "-q"],
        input="\n".join(commands) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"PARI/GP stderr: {completed.stderr.strip()}")

    rows = {}
    for match in LOCAL_PATTERN.finditer(completed.stdout):
        prime = int(match.group(1))
        change = tuple(item.strip() for item in match.group(4).split(","))
        if change != ("1", "0", "0", "0"):
            raise AssertionError(f"public model is not minimal at {prime}: {change}")
        rows[prime] = (
            int(match.group(2)),
            int(match.group(3)),
            int(match.group(5)),
            int(match.group(6)),
        )
    if rows != EXPECTED_LOCAL:
        raise AssertionError(f"local-reduction fingerprint changed: {rows}")

    minimal_match = re.search(r"^MINIMAL\|\[(.*?)\]$", completed.stdout, re.MULTILINE)
    pari_match = re.search(r"^PARI\|\[(.*?)\]$", completed.stdout, re.MULTILINE)
    if minimal_match is None or pari_match is None:
        raise AssertionError("PARI/GP omitted the minimal-model diagnostics")
    minimal = tuple(item.strip() for item in minimal_match.group(1).split(","))
    expected = tuple(str(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS)
    if minimal != expected:
        raise AssertionError("PARI/GP changed the public global minimal model")

    conductor = factor_product(CONDUCTOR_FACTORIZATION)
    if conductor != PUBLIC_CONDUCTOR:
        raise AssertionError("local conductor exponents do not recover the public conductor")
    # The global sign is the product of the finite local signs and -1 at infinity.
    root_number = -1
    for _prime, (_f, _kod, _c, sign) in rows.items():
        root_number *= sign
    return pari_match.group(1), root_number


def verify_primality() -> None:
    try:
        from sympy import isprime
    except ImportError as error:  # pragma: no cover - optional diagnostic
        raise RuntimeError("--verify-primality requires SymPy") from error
    values = {prime for prime, _ in (*C4_FACTORIZATION, *DISCRIMINANT_FACTORIZATION)}
    composites = [value for value in sorted(values) if not isprime(value)]
    if composites:
        raise AssertionError(f"reported prime factors failed primality: {composites}")
    print(f"{PROTOCOL}|stage=primality|factors={len(values)}|status=PASS", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-primality", action="store_true")
    args = parser.parse_args()

    if len(POINTS) != 29:
        raise AssertionError("expected exactly 29 public points")
    if not all(on_curve(point) for point in POINTS):
        raise AssertionError("a public point is off curve 356")
    if not all(on_short_curve(point) for point in SHORT_POINTS):
        raise AssertionError("the integral short-model transport failed")
    print(f"{PROTOCOL}|stage=membership|checked=29|status=PASS", flush=True)

    _b2, _b4, _b6, _b8, c4, _c6, discriminant = weierstrass_invariants()
    if discriminant != PUBLIC_DISCRIMINANT:
        raise AssertionError("the public discriminant changed")
    if factor_product(C4_FACTORIZATION) != c4:
        raise AssertionError("the c4 factorization changed")
    if factor_product(DISCRIMINANT_FACTORIZATION) != abs(discriminant):
        raise AssertionError("the discriminant factorization changed")
    if factor_product(CONDUCTOR_FACTORIZATION) != PUBLIC_CONDUCTOR:
        raise AssertionError("the conductor factorization changed")
    if max(valuation(discriminant, prime) for prime, _ in DISCRIMINANT_FACTORIZATION) >= 12:
        raise AssertionError("the elementary global-minimality check no longer applies")
    if args.verify_primality:
        verify_primality()
    print(
        f"{PROTOCOL}|stage=invariants|minimal=true|bad_primes={len(DISCRIMINANT_FACTORIZATION)}|status=PASS",
        flush=True,
    )

    short = short_coefficients()
    two_torsion_prime = find_two_torsion_certificate_prime(short, prime_bound=100)
    if two_torsion_prime != 17:
        raise AssertionError(f"least no-2-torsion witness changed to {two_torsion_prime}")
    torsion_counts = {7: finite_group_order(7), 17: finite_group_order(17)}
    if torsion_counts != {7: 12, 17: 25} or gcd(*torsion_counts.values()) != 1:
        raise AssertionError(f"finite torsion witnesses changed: {torsion_counts}")
    print(
        f"{PROTOCOL}|stage=torsion|orders=7:12,17:25|status=PASS",
        flush=True,
    )

    signatures = find_mod2_reduction_certificate(short, SHORT_POINTS, prime_bound=400)
    certificate_primes = tuple(signature.prime for signature in signatures)
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES:
        raise AssertionError(f"finite-reduction certificate changed: {certificate_primes}")
    rank = combined_mod2_rank(signatures, len(POINTS))
    row_count = sum(len(signature.rows) for signature in signatures)
    if rank != 29 or row_count != 30:
        raise AssertionError(f"finite reductions reached rank={rank}, rows={row_count}")
    print(
        f"{PROTOCOL}|stage=mod2|rank=29|rows=30|primes=20|status=PASS",
        flush=True,
    )

    pari_version, root_number = local_reduction_replay()
    if root_number != -1:
        raise AssertionError(f"global root number changed to {root_number}")
    print(
        f"{PROTOCOL}|stage=local|pari={pari_version}|root_number=-1|status=PASS",
        flush=True,
    )
    print(
        f"{PROTOCOL}|stage=done|rank_lower_bound=29|exact_rank=false|source={PUBLIC_SOURCE}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
