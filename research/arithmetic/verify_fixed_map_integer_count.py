#!/usr/bin/env python3
"""Exact bounded controls for FIXED_MAP_INTEGER_COUNT.md; no analytic proof.

Run from the repository root. To refresh the compact receipt explicitly:
    python3 research/arithmetic/verify_fixed_map_integer_count.py --write-receipt
Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
import platform


def sieve(limit: int) -> list[int]:
    smallest = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if smallest[p] == p:
            for n in range(p * p, limit + 1, p):
                if smallest[n] == n:
                    smallest[n] = p
    return smallest


def factor(n: int, smallest: list[int]) -> dict[int, int]:
    assert 1 <= n < len(smallest)
    result: dict[int, int] = {}
    while n > 1:
        p = smallest[n]
        result[p] = result.get(p, 0) + 1
        n //= p
    return result


def strip_three(n: int) -> tuple[int, int]:
    assert n > 0
    exponent = 0
    while n % 3 == 0:
        exponent += 1
        n //= 3
    return exponent, n


def integer_local(a: int, smallest: list[int]) -> bool:
    """Original integer local criterion, before excluding cubes."""
    assert a != 0
    e, unit = strip_three(abs(a))
    return (e % 3 == 0 and unit % 9 in (1, 8)
            and all(p % 3 != 2 or v % 3 == 0
                    for p, v in factor(unit, smallest).items()))


def index_local(h: int, smallest: list[int]) -> bool:
    """Multiplicative height-index restrictions and character projector."""
    e, unit = strip_three(h)
    return ((e == 0 or e % 3 == 2) and unit % 9 in (1, 8)
            and all(p % 3 != 2 or v % 3 == 0
                    for p, v in factor(unit, smallest).items()))


def is_cube(n: int, smallest: list[int]) -> bool:
    return all(v % 3 == 0 for v in factor(abs(n), smallest).values())


def direct_height(a: int) -> int:
    coords = (9, -9, 32 * a, 24 * a + 3)
    return max(map(abs, coords)) // gcd(*coords)


def char_exponent(unit: int) -> int:
    """chi(2)=omega, encoded exactly by its exponent modulo three."""
    return {1: 0, 2: 1, 4: 2, 8: 0, 7: 1, 5: 2}[unit % 9]


def run(bound: int, maximum_height: int) -> dict:
    assert bound >= 1 and maximum_height >= 1000
    limit = max(513, 3 * bound, 3 * maximum_height // 32 + 1)
    smallest = sieve(limit)
    counts = {"height_index_bijection": 0, "signed_primitive_heights": 0,
              "character_extension": 0, "character_projector": 0,
              "euler_coefficient_controls": 0, "rejection_controls": 0}

    for h in range(1, bound + 1):
        a = h if h % 3 else 3 * h
        e, _ = strip_three(h)
        # The inverse is asserted only on the permissible 3-adic index set.
        if e == 0 or e % 3 == 2:
            assert 32 * a // gcd(a, 3) == 32 * h
            assert index_local(h, smallest) == integer_local(a, smallest)
            counts["height_index_bijection"] += 2
        else:
            assert not index_local(h, smallest)
            counts["height_index_bijection"] += 1
        _, unit = strip_three(h)
        via_primes = sum(char_exponent(p) * v for p, v in
                         factor(h, smallest).items() if p != 3) % 3
        assert via_primes == char_exponent(unit)
        counts["character_extension"] += 1

    for a in range(1, bound + 1):
        for signed in (a, -a):
            assert direct_height(signed) == 32 * a // gcd(a, 3)
            assert integer_local(signed, smallest) == integer_local(a, smallest)
            counts["signed_primitive_heights"] += 2

    omega = ((1, 0), (0, 1), (-1, -1))
    for residue in (1, 2, 4, 5, 7, 8):
        k = char_exponent(residue)
        summed = tuple(sum(omega[j * k % 3][i] for j in range(3))
                       for i in range(2))
        assert summed == ((3, 0) if residue in (1, 8) else (0, 0))
        counts["character_projector"] += 1

    # A3(q)*(1-q^3)=1+q^2-q^3, checked coefficientwise through degree 30.
    coefficients = [int(e == 0 or e % 3 == 2) for e in range(31)]
    for e, value in enumerate(coefficients):
        product_coefficient = value - (coefficients[e - 3] if e >= 3 else 0)
        assert product_coefficient == {0: 1, 2: 1, 3: -1}.get(e, 0)
        counts["euler_coefficient_controls"] += 1
    a3 = 1 + Q(1, 9) / (1 - Q(1, 27))
    assert a3 == Q(29, 26)
    assert Q(2, 3 * 32) * a3 == Q(29, 1248)
    assert Q(29, 1248) / Q(1, 96) == Q(29, 13)
    x = Q(1, 3)
    a3_derivative = (2 * x + x**4) / (1 - x**3)**2
    assert x * a3_derivative / a3 == Q(165, 754)
    # Gamma(1/2)/Gamma(-1/2)=-1/2 fixes the correction sign.
    for gamma, logarithmic_derivative in ((Q(1, 7), Q(2, 5)),
                                         (Q(3, 2), Q(-4, 3))):
        coefficient = -Q(1, 2) * (logarithmic_derivative + gamma / 2 - 1)
        assert coefficient == Q(1, 2) - gamma / 4 - logarithmic_derivative / 2
    counts["euler_coefficient_controls"] += 6

    examples = []
    for a, expected_local, expected_hasse in (
            (19, True, True), (-19, True, True), (152, True, True),
            (513, True, True), (7, False, False), (8, True, False),
            (27, True, False), (54, False, False)):
        local = integer_local(a, smallest)
        hasse = local and not is_cube(a, smallest)
        assert (local, hasse) == (expected_local, expected_hasse)
        examples.append({"a": a, "height": direct_height(a),
                         "local": local, "hasse_failure": hasse})
    # Reject deleting inert cubes, discarding the 3-adic content, and using
    # the ordinary Dirichlet character (zero at three) in the projector.
    assert integer_local(152, smallest) and 152 % 2 == 0
    assert direct_height(513) != 32 * 513
    assert index_local(171, smallest) and 171 % 3 == 0
    assert Q(29, 26) != Q(27, 26) and Q(29, 13) != Q(29, 26)
    counts["rejection_controls"] = 4

    enumeration = []
    for height in sorted({1000, 10000, 100000, maximum_height}):
        if height > maximum_height:
            continue
        direct = sum(integer_local(a, smallest) and not is_cube(a, smallest)
                     and direct_height(a) <= height
                     for a in range(1, 3 * height // 32 + 1))
        indexed = sum(index_local(h, smallest)
                      and not is_cube(h if h % 3 else 3 * h, smallest)
                      for h in range(1, height // 32 + 1))
        assert direct == indexed
        enumeration.append({"height": height, "positive_hasse_count": direct,
                            "signed_hasse_count": 2 * direct})

    here = Path(__file__).resolve().parent
    proof = here / "FIXED_MAP_INTEGER_COUNT.md"
    return {
        "status": "pass", "python": platform.python_version(),
        "command": ("python3 research/arithmetic/verify_fixed_map_integer_count.py "
                    f"--bound {bound} --maximum-height {maximum_height}"),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
        "bounds": {"integer_and_index_bound": bound,
                   "maximum_height": maximum_height, "sieve_limit": limit},
        "checks": counts, "examples": examples, "enumeration": enumeration,
        "exact_constants": {"A3_at_1": "29/26", "C_int_rational_factor": "29/1248",
                            "signed_ratio_to_sufficient_without_J": "29/13",
                            "A3_log_derivative_over_minus_log_3": "165/754"},
        "scope": ("Exact bounded arithmetic and coefficient controls only; "
                  "the written proof and imported Selberg-Delange theorem establish "
                  "the asymptotic, not finite extrapolation."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=50000)
    parser.add_argument("--maximum-height", type=int, default=1000000)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run(args.bound, args.maximum_height)
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write_receipt:
        path = Path(__file__).resolve().with_name("fixed_map_integer_count_verification.json")
        path.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
