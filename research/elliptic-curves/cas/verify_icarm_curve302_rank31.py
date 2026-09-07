#!/usr/bin/env python3
"""Build an exact unconditional rank-at-least-31 certificate for ICARM curve 302.

The public leaderboard certifies the submitted points with its Brumer--Cremona
quadratic-character implementation.  This checker is deliberately independent:
it verifies every coordinate over Q and proves independence from exact images in
products of finite quotients E(F_p)/2E(F_p).  Trivial rational 2-torsion then
turns parity of every putative relation into an infinite-descent contradiction.

The conductor/fibre block is ancillary to the rank proof.  It verifies the
published integer factorizations exactly and derives the Kodaira symbols from
minimal-model valuations; ``--verify-primality`` additionally rechecks that the
listed factors are prime with SymPy.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import sys
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve302 import (  # noqa: E402
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


PROTOCOL = "R31ICARM"
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve302_rank31_v1.json"
)
REPRODUCING_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/verify_icarm_curve302_rank31.py"
)
PUBLIC_SOURCE = "https://elliptic-rank.icarm.cloud/curve/302"
PUBLIC_POST = "https://x.com/dwrensha/status/2091619906701689228"

PUBLIC_CONDUCTOR = int(
    "9058040138868690262858791849730511364320517732812998655982537144987090191789878126630171056773171491466723174687496041392810438595512770829603147288102844437462050"
)
PUBLIC_DISCRIMINANT = int(
    "56986667162894850943626331069759320443067423220527305082244725946223212417010574363800489149279120357391198820328667683955054948591887684940751932057558936585975074155693048247805850350226862080000"
)
PUBLIC_NAIVE_HEIGHT = "468.2771"
PUBLIC_FALTINGS_HEIGHT = "36.7425"
PUBLIC_REGULATOR = (
    "5520367374821893536678475926502746956624.072603230810343324240227195436655202216427344786535970988197"
)

C4_FACTORIZATION = (
    (5, 2),
    (31, 1),
    (109, 2),
    (223, 1),
    (5351, 1),
    (10745106459799, 1),
    (522332575866045267821375338901621944927313, 1),
)

DISCRIMINANT_FACTORIZATION = (
    (2, 15),
    (3, 4),
    (5, 4),
    (7, 6),
    (11, 4),
    (13, 5),
    (19, 2),
    (23, 2),
    (29, 3),
    (37, 2),
    (41, 2),
    (73, 2),
    (131, 2),
    (167, 2),
    (7547, 1),
    (632881, 1),
    (966509, 1),
    (18145679437533309132469, 1),
    (767028866604834801397681553, 1),
    (
        30580600452196904409276223329355584892025407195996968868775951126238056443210297,
        1,
    ),
)

CONDUCTOR_FACTORIZATION = tuple(
    (prime, 2 if prime == 5 else 1) for prime, _ in DISCRIMINANT_FACTORIZATION
)

EXPECTED_CERTIFICATE_PRIMES = (
    17,
    47,
    53,
    61,
    67,
    71,
    79,
    83,
    89,
    101,
    107,
    113,
    127,
    137,
    149,
    179,
    191,
    197,
    211,
    233,
    241,
    263,
    269,
    281,
    283,
    293,
    311,
)

EXPECTED_CHARACTER_PRIMES = (
    17,
    47,
    53,
    61,
    67,
    71,
    79,
    83,
    89,
    101,
    107,
    113,
    127,
    137,
    149,
    179,
    191,
    197,
    211,
    229,
    233,
    241,
    263,
    269,
    281,
    283,
)


def factor_product(factors: Iterable[tuple[int, int]]) -> int:
    answer = 1
    for prime, exponent in factors:
        answer *= prime**exponent
    return answer


def valuation(value: int, prime: int) -> int:
    if not value:
        raise ValueError("valuation of zero is not used here")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def rational_text(value: Fraction | int) -> str:
    return str(Fraction(value))


def point_digest() -> str:
    payload = json.dumps(
        [[rational_text(x), rational_text(y)] for x, y in POINTS],
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def finite_group_order(prime: int) -> int:
    coefficients = short_coefficients()
    coefficient_a = int(coefficients[3]) % prime
    coefficient_b = int(coefficients[4]) % prime
    return len(finite_curve_points(coefficient_a, coefficient_b, prime))


def local_reduction_records(c4: int, c6: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for prime, discriminant_exponent in DISCRIMINANT_FACTORIZATION:
        c4_exponent = valuation(c4, prime)
        c6_exponent = valuation(c6, prime)
        if prime == 5:
            if (c4_exponent, c6_exponent, discriminant_exponent) != (2, 2, 4):
                raise AssertionError("the type-IV valuation triple at 5 changed")
            kodaira = "IV"
            conductor_exponent = 2
        else:
            if c4_exponent != 0:
                raise AssertionError(
                    f"expected multiplicative reduction at {prime}, got v(c4)={c4_exponent}"
                )
            kodaira = f"I{discriminant_exponent}"
            conductor_exponent = 1
        records.append(
            {
                "prime": str(prime),
                "v_c4": c4_exponent,
                "v_c6": c6_exponent,
                "v_discriminant": discriminant_exponent,
                "kodaira_symbol": kodaira,
                "conductor_exponent": conductor_exponent,
            }
        )
    return records


def gf2_rank_rows(rows: Iterable[Iterable[int]], column_count: int) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        row = tuple(int(value) for value in row)
        if len(row) != column_count or any(value not in (0, 1) for value in row):
            raise ValueError("quadratic-character rows must be binary and have fixed width")
        packed = sum(value << index for index, value in enumerate(row))
        while packed:
            pivot = packed.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = packed
                break
            packed ^= pivots[pivot]
    return len(pivots)


def small_primes_up_to(bound: int) -> tuple[int, ...]:
    answer: list[int] = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in answer if prime * prime <= candidate):
            answer.append(candidate)
    return tuple(answer)


def quadratic_character_block(prime: int) -> tuple[tuple[int, ...], ...]:
    """Return Brumer--Cremona character rows for one good prime.

    For a root theta of the short 2-division cubic and x(P)=u/w^2 in
    lowest terms, the character is the Legendre symbol of
    ``u-theta*w^2``; at zero use the derivative.  At most two roots are
    needed because the product of all three characters is trivial.
    """

    _, _, _, coefficient_a_q, coefficient_b_q = short_coefficients()
    coefficient_a = int(coefficient_a_q)
    coefficient_b = int(coefficient_b_q)
    short_discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if prime <= 3 or short_discriminant % prime == 0:
        return ()
    roots = tuple(
        value
        for value in range(prime)
        if (value**3 + coefficient_a * value + coefficient_b) % prime == 0
    )
    rows: list[tuple[int, ...]] = []
    for root in roots[:2]:
        entries: list[int] = []
        for x_value, _ in SHORT_POINTS:
            numerator = x_value.numerator
            square_denominator = x_value.denominator
            denominator_root = isqrt(square_denominator)
            if denominator_root * denominator_root != square_denominator:
                raise AssertionError("a short-model x denominator is not a square")
            character_argument = (numerator - root * square_denominator) % prime
            if character_argument == 0:
                character_argument = (3 * root * root + coefficient_a) % prime
            legendre = pow(character_argument, (prime - 1) // 2, prime)
            if legendre == 1:
                entries.append(0)
            elif legendre == prime - 1:
                entries.append(1)
            else:
                raise AssertionError(
                    f"zero quadratic-character argument at good prime {prime}"
                )
        rows.append(tuple(entries))
    return tuple(rows)


def quadratic_character_certificate(
    *, prime_bound: int = 400
) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...], int]:
    used_primes: list[int] = []
    rows: list[tuple[int, ...]] = []
    for prime in small_primes_up_to(prime_bound):
        block = quadratic_character_block(prime)
        if not block:
            continue
        used_primes.append(prime)
        rows.extend(block)
        rank = gf2_rank_rows(rows, len(POINTS))
        if rank == len(POINTS):
            return tuple(used_primes), tuple(rows), rank
    return tuple(used_primes), tuple(rows), gf2_rank_rows(rows, len(POINTS))


def verify_reported_primes() -> None:
    try:
        from sympy import isprime
    except ImportError as error:  # pragma: no cover - optional diagnostic
        raise RuntimeError("--verify-primality requires SymPy") from error
    factors = {prime for prime, _ in (*C4_FACTORIZATION, *DISCRIMINANT_FACTORIZATION)}
    composites = [value for value in sorted(factors) if not isprime(value)]
    if composites:
        raise AssertionError(f"reported prime factors failed primality: {composites}")
    print(f"{PROTOCOL}|stage=primality|factors={len(factors)}|status=PASS", flush=True)


def build_certificate() -> dict[str, object]:
    print(f"{PROTOCOL}|stage=input|points={len(POINTS)}", flush=True)
    if len(POINTS) != 31:
        raise AssertionError("expected exactly 31 public points")
    for index, point in enumerate(POINTS, 1):
        if not on_curve(point):
            raise AssertionError(f"public point {index} is off the announced curve")
    if not all(on_short_curve(point) for point in SHORT_POINTS):
        raise AssertionError("the integral short-model transport failed")
    print(f"{PROTOCOL}|stage=membership|checked=31|status=PASS", flush=True)

    b2, b4, b6, b8, c4, c6, discriminant = weierstrass_invariants()
    if discriminant != PUBLIC_DISCRIMINANT:
        raise AssertionError("the public discriminant changed")
    if factor_product(C4_FACTORIZATION) != c4:
        raise AssertionError("the c4 factorization changed")
    if factor_product(DISCRIMINANT_FACTORIZATION) != discriminant:
        raise AssertionError("the discriminant factorization changed")
    if factor_product(CONDUCTOR_FACTORIZATION) != PUBLIC_CONDUCTOR:
        raise AssertionError("the conductor factorization changed")

    # An integral model can be nonminimal at p only if p^4 divides c4 and
    # p^12 divides Delta.  Every displayed v_p(Delta) is below 12 except at 2,
    # while c4 is odd.  Hence the public equation is globally minimal.
    possible_nonminimal = [
        prime
        for prime, delta_exponent in DISCRIMINANT_FACTORIZATION
        if delta_exponent >= 12 and valuation(c4, prime) >= 4
    ]
    if possible_nonminimal:
        raise AssertionError(f"nonminimal primes survived: {possible_nonminimal}")
    reductions = local_reduction_records(c4, c6)
    derived_conductor = factor_product(
        (int(record["prime"]), int(record["conductor_exponent"]))
        for record in reductions
    )
    if derived_conductor != PUBLIC_CONDUCTOR:
        raise AssertionError("the local fibre data did not recover the conductor")
    print(
        f"{PROTOCOL}|stage=invariants|minimal=true|bad_primes={len(reductions)}|status=PASS",
        flush=True,
    )

    short = short_coefficients()
    two_torsion_prime = find_two_torsion_certificate_prime(short, prime_bound=100)
    if two_torsion_prime != 31:
        raise AssertionError(f"least no-2-torsion witness changed to {two_torsion_prime}")
    torsion_counts = {17: finite_group_order(17), 31: finite_group_order(31)}
    if torsion_counts != {17: 26, 31: 43}:
        raise AssertionError(f"finite torsion witnesses changed: {torsion_counts}")
    if gcd(*torsion_counts.values()) != 1:
        raise AssertionError("the two good reductions no longer force trivial torsion")
    print(
        f"{PROTOCOL}|stage=torsion|orders=17:{torsion_counts[17]},31:{torsion_counts[31]}|status=PASS",
        flush=True,
    )

    signatures = find_mod2_reduction_certificate(short, SHORT_POINTS, prime_bound=400)
    certificate_primes = tuple(signature.prime for signature in signatures)
    if certificate_primes != EXPECTED_CERTIFICATE_PRIMES:
        raise AssertionError(
            "the deterministic finite-reduction certificate changed: "
            f"{certificate_primes}"
        )
    rank = combined_mod2_rank(signatures, len(POINTS))
    if rank != 31:
        raise AssertionError(f"finite reductions reached rank {rank}, not 31")
    row_count = sum(len(signature.rows) for signature in signatures)
    print(
        f"{PROTOCOL}|stage=mod2|rank={rank}|rows={row_count}|primes={len(signatures)}|status=PASS",
        flush=True,
    )

    character_primes, character_rows, character_rank = quadratic_character_certificate(
        prime_bound=400
    )
    if character_primes != EXPECTED_CHARACTER_PRIMES:
        raise AssertionError(
            "the deterministic quadratic-character certificate changed: "
            f"{character_primes}"
        )
    if character_rank != 31:
        raise AssertionError(
            f"quadratic characters reached rank {character_rank}, not 31"
        )
    print(
        f"{PROTOCOL}|stage=characters|rank={character_rank}"
        f"|rows={len(character_rows)}|primes={len(character_primes)}|status=PASS",
        flush=True,
    )

    j_value = Fraction(c4**3, discriminant)
    script_path = Path(__file__).resolve()
    model_path = script_path.with_name("icarm_curve302.py")
    return {
        "schema_version": 1,
        "artifact_kind": "exact_elliptic_curve_rank_lower_bound",
        "curve_id": "icarm_curve_302",
        "claim": "rank E(Q) >= 31",
        "claim_status": "exact unconditional lower bound; no unconditional exact-rank claim",
        "public_sources": [PUBLIC_SOURCE, PUBLIC_POST],
        "public_attribution": "Claude, Levent Alpoge, and Ava Howell",
        "public_conditional_statement": {
            "statement": "rank 31 conditional on BSD and GRH",
            "independently_replayed_here": False,
            "used_for_rank_lower_bound": False,
        },
        "curve": {
            "ainvs": [rational_text(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS],
            "b2": str(b2),
            "b4": str(b4),
            "b6": str(b6),
            "b8": str(b8),
            "c4": str(c4),
            "c6": str(c6),
            "j_invariant": rational_text(j_value),
            "j_invariant_sha256": hashlib.sha256(str(j_value).encode()).hexdigest(),
            "discriminant": str(discriminant),
            "conductor": str(PUBLIC_CONDUCTOR),
            "global_minimal_model": True,
            "minimality_proof": (
                "Nonminimality at p requires v_p(c4)>=4 and v_p(Delta)>=12; "
                "the only listed Delta exponent >=12 is p=2 and c4 is odd."
            ),
            "torsion_subgroup": "trivial",
            "naive_log_height_public": PUBLIC_NAIVE_HEIGHT,
            "faltings_height_public": PUBLIC_FALTINGS_HEIGHT,
        },
        "short_model": {
            "change_of_variables": {
                "X": "36*x+15",
                "Y": "108*(2*y+x+1)",
            },
            "ainvs": [rational_text(value) for value in short],
            "isomorphism_over_Q": True,
            "all_transported_points_checked": True,
        },
        "points": [[rational_text(x), rational_text(y)] for x, y in POINTS],
        "point_count": len(POINTS),
        "points_sha256": point_digest(),
        "exact_membership_checks_passed": len(POINTS),
        "torsion_certificate": {
            "method": "coprime finite-group orders at two good primes",
            "finite_group_orders": {str(prime): order for prime, order in torsion_counts.items()},
            "gcd": gcd(*torsion_counts.values()),
            "no_rational_2_torsion_witness_prime": two_torsion_prime,
            "short_2_division_cubic_has_no_root_mod_witness_prime": True,
        },
        "independence_certificate": {
            "method": "finite good-reduction quotients E(F_p)/2E(F_p)",
            "proof": (
                "Full F2 column rank forces every coefficient of a rational-point "
                "relation to be even. Trivial rational 2-torsion permits division "
                "by two; infinite descent forces every coefficient to vanish."
            ),
            "relation_prime": 2,
            "certificate_prime_bound": 400,
            "certificate_primes": list(certificate_primes),
            "matrix_row_count": row_count,
            "combined_binary_rank": rank,
            "signatures": [
                {
                    "prime": signature.prime,
                    "group_order": signature.group_order,
                    "doubled_subgroup_order": signature.doubled_subgroup_order,
                    "quotient_dimension": signature.quotient_dimension,
                    "matrix_rows": [list(row) for row in signature.rows],
                }
                for signature in signatures
            ],
        },
        "quadratic_character_cross_check": {
            "method": "Brumer--Cremona quadratic characters at roots of the 2-division cubic",
            "role": "second exact independence certificate; independently implemented",
            "certificate_prime_bound": 400,
            "certificate_primes": list(character_primes),
            "matrix_row_count": len(character_rows),
            "combined_binary_rank": character_rank,
            "matrix_rows": [list(row) for row in character_rows],
        },
        "factorizations": {
            "c4": [
                {"factor": str(prime), "exponent": exponent}
                for prime, exponent in C4_FACTORIZATION
            ],
            "discriminant": [
                {"factor": str(prime), "exponent": exponent}
                for prime, exponent in DISCRIMINANT_FACTORIZATION
            ],
            "conductor": [
                {"factor": str(prime), "exponent": exponent}
                for prime, exponent in CONDUCTOR_FACTORIZATION
            ],
            "all_integer_products_checked_exactly": True,
            "primality_required_for_rank_claim": False,
            "optional_primality_command": f"{REPRODUCING_COMMAND} --check --verify-primality",
        },
        "local_reduction": {
            "records": reductions,
            "all_bad_fibres_multiplicative_except_prime_5": True,
            "prime_5_kodaira_symbol": "IV",
            "derived_conductor_equals_public_conductor": True,
        },
        "height_diagnostic": {
            "public_regulator_numeric": PUBLIC_REGULATOR,
            "used_for_rank_claim": False,
        },
        "construction_status": {
            "identified_as_H3_or_R17_specialization": False,
            "reason": (
                "No public family equation or specialization parameter is supplied, "
                "and the reconstructed generic rootless MW17 equation is not yet available."
            ),
            "conditional_extra_specialization_directions_if_generic_rank_17": 14,
        },
        "generation": {
            "command": REPRODUCING_COMMAND,
            "arithmetic": "exact rational and exhaustive finite-field group operations",
            "python": sys.version.split()[0],
            "checker_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
            "model_data_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify that --output already equals the deterministic certificate",
    )
    parser.add_argument(
        "--verify-primality",
        action="store_true",
        help="also recheck all reported factors with SymPy isprime",
    )
    args = parser.parse_args()
    if args.verify_primality:
        verify_reported_primes()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise SystemExit(f"missing pinned certificate: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale pinned certificate: rerun {REPRODUCING_COMMAND}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(
        f"{PROTOCOL}|stage=done|status=PASS|rank_lower_bound=31"
        f"|mode={'check' if args.check else 'write'}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
