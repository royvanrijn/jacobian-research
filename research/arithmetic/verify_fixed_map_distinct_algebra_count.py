#!/usr/bin/env python3
"""Bounded exact controls for FIXED_MAP_DISTINCT_ALGEBRA_COUNT.md.

Uses only Python's standard library. Finite counts do not prove asymptotics
or establish the pure-cubic field-isomorphism theorem.
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


def factors(n: int, smallest: list[int]) -> dict[int, int]:
    assert 1 <= n < len(smallest)
    result: dict[int, int] = {}
    while n > 1:
        p = smallest[n]
        result[p] = result.get(p, 0) + 1
        n //= p
    return result


def squarefree_split(n: int, smallest: list[int]) -> bool:
    return all(p % 3 == 1 and e == 1 for p, e in factors(n, smallest).items())


def cube_class(m: int, n: int, smallest: list[int]) -> tuple[tuple[int, int], ...]:
    exponents = factors(abs(m), smallest)
    for p, e in factors(n, smallest).items():
        exponents[p] = exponents.get(p, 0) - e
    return tuple(sorted((p, e % 3) for p, e in exponents.items() if e % 3))


def inverse_class(key: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    return tuple((p, 3 - e) for p, e in key)


def algebra_class_key(m: int, n: int, smallest: list[int]) -> tuple:
    key = cube_class(m, n, smallest)
    return min(key, inverse_class(key))


def height(m: int, n: int) -> int:
    coords = (9 * n, -9 * n, 32 * m, 24 * m + 3 * n)
    return max(map(abs, coords)) // gcd(*coords)


def omega(exponent: int) -> tuple[int, int]:
    return ((1, 0), (0, 1), (-1, -1))[exponent % 3]


def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return x[0] + y[0], x[1] + y[1]


def mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] - x[1] * y[1]


def local_h(u: int, v: int, twist: int) -> tuple[int, int]:
    if u == v == 0:
        return 1, 0
    if u == 0 or v == 0:
        return 0, 0
    z = omega(twist * (u - v))
    sign = (-1) ** (u + v - 1)
    return sign * z[0], sign * z[1]


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def correction(d: int, e: int, smallest: list[int]) -> int:
    fd, fe = factors(d, smallest), factors(e, smallest)
    if fd.keys() != fe.keys() or any(p % 3 != 1 for p in fd):
        return 0
    return (-1) ** sum(fd[p] + fe[p] - 1 for p in fd)


def cubefree_split(n: int, smallest: list[int]) -> bool:
    return all(p % 3 == 1 and e <= 2 for p, e in factors(n, smallest).items())


def canonical_pair(key: tuple) -> tuple[int, int]:
    m = n = 1
    for p, e in key:
        if e == 1:
            m *= p
        else:
            n *= p
    return m, n


def allocation_pairs(key: tuple, free_primes: set[int] | None = None) -> list[tuple[int, int]]:
    pairs = [(1, 1)]
    for p, e in key:
        if free_primes is not None and p not in free_primes:
            choices = [(p, 1)] if e == 1 else [(1, p)]
        else:
            choices = [(p, 1), (1, p*p)] if e == 1 else [(p*p, 1), (1, p)]
        pairs = [(r*rp, s*sp) for r, s in pairs for rp, sp in choices]
    return pairs


def optimized_height(r: int, s: int) -> int:
    return min(max(4*r, 9*s), max(36*r, 3*s))


def minimum_height(key: tuple, free_primes: set[int] | None = None) -> int:
    return min(min(optimized_height(r, s), optimized_height(s, r))
               for r, s in allocation_pairs(key, free_primes))


OPTIMIZED_RECTANGLES = [(Q(1, 4), Q(1, 9)), (Q(1, 36), Q(1, 3)),
                        (Q(1, 9), Q(1, 4)), (Q(1, 3), Q(1, 36))]


def rectangle_union_area(rectangles: list[tuple[Q, Q]]) -> Q:
    previous = answer = Q(0)
    for width in sorted({w for w, _ in rectangles}):
        top = max(h for w, h in rectangles if w >= width)
        answer += (width - previous) * top
        previous = width
    return answer


def constant_coefficient(primes: tuple[int, ...]) -> Q:
    from itertools import product
    total = Q(0)
    for pattern in product((0, 1, 2), repeat=len(primes)):
        key = tuple((p, e) for p, e in zip(primes, pattern) if e)
        rectangles = [(w/r, h/s) for r, s in allocation_pairs(key)
                      for w, h in OPTIMIZED_RECTANGLES]
        total += rectangle_union_area(rectangles)
    denominator = Q(6)
    for p in primes:
        denominator *= 1 + Q(2, p)
    return total / denominator


def local_info(n: int, smallest: list[int]) -> tuple[bool, int, int]:
    exponents = factors(n, smallest)
    inert_ok = all(e % 3 == 0 for p, e in exponents.items() if p % 3 == 2)
    e3 = exponents.get(3, 0)
    return inert_ok, e3, n // 3**e3


def new_class_controls(smallest: list[int]) -> tuple[dict, list, list]:
    counts = {"optimal_height_cube_moves": 0, "constructive_optimal_heights": 0,
              "exact_finite_prime_area_sums": 0, "complete_class_enumerations": 0}
    assert rectangle_union_area(OPTIMIZED_RECTANGLES) == Q(31, 648)
    coefficients = []
    expected = [((), Q(31, 3888)), ((7,), Q(1747, 214326)),
                ((7, 13), Q(26155, 3184272)),
                ((7, 13, 19), Q(18113057423, 2196736908912)),
                ((7, 13, 19, 31), Q(2157454787546027, 261244190971223460))]
    previous = Q(0)
    for primes, value in expected:
        coefficient = constant_coefficient(primes)
        assert coefficient == value and coefficient >= previous
        previous = coefficient
        coefficients.append({"primes": list(primes), "C_P_over_kappa_sf": str(coefficient)})
        counts["exact_finite_prime_area_sums"] += 1

    for m, n in [(7, 1), (13, 1), (19, 1), (7, 13), (1729, 1)]:
        key = cube_class(m, n, smallest)
        assert canonical_pair(key) == (m, n)
        claimed = minimum_height(key)
        brute = None
        for base_m, base_n in [(m, n), (n, m)]:
            for x in range(1, 43):
                for y in range(1, 43):
                    parameter = Q(base_m*x**3, base_n*y**3)
                    actual = height(-parameter.numerator, parameter.denominator)
                    assert actual >= claimed
                    brute = actual if brute is None else min(brute, actual)
                    counts["optimal_height_cube_moves"] += 1
        assert brute == claimed
        for r, s in allocation_pairs(key):
            a0, a1 = Q(-r, 8*s), Q(-27*r, 8*s)
            assert height(a0.numerator, a0.denominator) == max(4*r, 9*s)
            assert height(a1.numerator, a1.denominator) == max(36*r, 3*s)
            counts["constructive_optimal_heights"] += 2
    assert minimum_height(cube_class(19, 1, smallest)) == 57
    assert minimum_height(cube_class(1729, 1, smallest)) == 988

    full_counts = []
    for bound in (100, 300, 1000):
        info = {n: local_info(n, smallest) for n in range(1, 8*bound // 3 + 1)}
        direct_classes = set()
        for m in range(1, 3*bound // 4 + 1):
            gm, em, um = info[m]
            if not gm:
                continue
            for n in range(1, 8*bound // 3 + 1):
                gn, en, un = info[n]
                if not gn or (em-en) % 3 or gcd(m, n) != 1:
                    continue
                if um * pow(un, -1, 9) % 9 not in (1, 8) or height(-m, n) > bound:
                    continue
                key = algebra_class_key(m, n, smallest)
                if key:
                    direct_classes.add(key)
        cf = [n for n in range(1, bound // 3 + 1) if cubefree_split(n, smallest)]
        allocation_classes = set()
        for r in cf:
            for s in cf:
                if gcd(r, s) != 1 or min(optimized_height(r, s), optimized_height(s, r)) > bound:
                    continue
                key = algebra_class_key(r, s, smallest)
                m, n = canonical_pair(key)
                if key and (m-n) % 9 == 0:
                    allocation_classes.add(key)
        assert direct_classes == allocation_classes
        previous_set = set()
        finite_counts = []
        for primes in (set(), {7}, {7, 13}, {7, 13, 19, 31}):
            captured = {key for key in direct_classes if minimum_height(key, primes) <= bound}
            assert previous_set <= captured <= direct_classes
            previous_set = captured
            finite_counts.append({"primes": sorted(primes), "classes": len(captured)})
        assert all(minimum_height(key) <= bound for key in direct_classes)
        full_counts.append({"height": bound, "direct_and_normal_form_classes": len(direct_classes),
                            "finite_prime_classes": finite_counts})
        counts["complete_class_enumerations"] += 1
    return counts, coefficients, full_counts

def run(pair_bound: int, coefficient_bound: int, maximum_height: int) -> dict:
    assert pair_bound >= 1 and coefficient_bound >= 1 and maximum_height >= 1000
    limit = max(2667, pair_bound, coefficient_bound, maximum_height // 9 + 1)
    smallest = sieve(limit)
    allowed = [n for n in range(1, limit + 1) if squarefree_split(n, smallest)]
    controls = {"twisted_euler_coefficients": 0, "global_coprime_convolution": 0,
                "cube_class_representatives": 0, "signed_heights": 0,
                "constant_algebra": 0, "rejection_controls": 0}

    for twist in range(3):
        factor = {(0, 0): (1, 0), (1, 0): omega(twist),
                  (0, 1): omega(-twist), (1, 1): (1, 0)}
        numerator = {(0, 0): (1, 0), (1, 0): omega(twist),
                     (0, 1): omega(-twist)}
        for u in range(7):
            for v in range(7):
                product = (0, 0)
                for (i, j), value in factor.items():
                    if i <= u and j <= v:
                        product = add(product, mul(value, local_h(u - i, v - j, twist)))
                assert product == numerator.get((u, v), (0, 0))
                controls["twisted_euler_coefficients"] += 1

    divisor_table = {n: divisors(n) for n in range(1, coefficient_bound + 1)}
    for m in range(1, coefficient_bound + 1):
        for n in range(1, coefficient_bound + 1):
            convolution = sum(correction(d, e, smallest)
                              * squarefree_split(m // d, smallest)
                              * squarefree_split(n // e, smallest)
                              for d in divisor_table[m] for e in divisor_table[n])
            expected = int(squarefree_split(m, smallest)
                           and squarefree_split(n, smallest) and gcd(m, n) == 1)
            assert convolution == expected
            controls["global_coprime_convolution"] += 1

    representatives: dict[tuple, tuple[int, int]] = {}
    small_allowed = [n for n in allowed if n <= pair_bound]
    for m in small_allowed:
        for n in small_allowed:
            if gcd(m, n) != 1:
                continue
            key = cube_class(m, n, smallest)
            reconstructed_m = reconstructed_n = 1
            for p, e in key:
                if e == 1:
                    reconstructed_m *= p
                else:
                    reconstructed_n *= p
            assert (reconstructed_m, reconstructed_n) == (m, n)
            assert inverse_class(key) == cube_class(n, m, smallest)
            class_key = algebra_class_key(m, n, smallest)
            sorted_pair = tuple(sorted((m, n)))
            assert representatives.setdefault(class_key, sorted_pair) == sorted_pair
            controls["cube_class_representatives"] += 3
            assert height(-m, n) == max(32 * m, 9 * n)
            assert height(m, n) >= height(-m, n)
            controls["signed_heights"] += 2

    assert algebra_class_key(19, 1, smallest) == algebra_class_key(152, 1, smallest)
    assert algebra_class_key(19, 1, smallest) == algebra_class_key(1, 19, smallest)
    assert algebra_class_key(19, 1, smallest) != algebra_class_key(37, 1, smallest)
    controls["rejection_controls"] = 3

    for p in (7, 13, 19, 31):
        q = Q(1, p)
        assert (1 - q**2)**2 * (1 + 2 * q) / (1 + q)**2 == 1 - 3 * q**2 + 2 * q**3
        assert 0 < 1 - 3 * q**2 + 2 * q**3 < 1
        controls["constant_algebra"] += 2
    area = 2 * Q(1, 32 * 9) - Q(1, 32**2)
    assert area == Q(55, 9216)
    assert area / 6 == Q(55, 55296)
    controls["constant_algebra"] += 2

    enumeration = []
    for bound in sorted({1000, 10000, maximum_height}):
        if bound > maximum_height:
            continue
        rectangle = {(m, n) for m in allowed if m <= bound // 32
                     for n in allowed if n <= bound // 9
                     if gcd(m, n) == 1 and (m - n) % 9 == 0}
        square = {(m, n) for m in allowed if m <= bound // 32
                  for n in allowed if n <= bound // 32
                  if gcd(m, n) == 1 and (m - n) % 9 == 0}
        union = rectangle | {(n, m) for m, n in rectangle}
        assert len(union) == 2 * len(rectangle) - len(square)
        classes = {algebra_class_key(m, n, smallest) for m, n in rectangle
                   if (m, n) != (1, 1)}
        assert 2 * len(classes) == len(union) - 1
        assert len(classes) == len(rectangle) - (len(square) + 1) // 2
        assert all(height(-m, n) <= bound for m, n in rectangle)
        enumeration.append({"height": bound, "rectangle_pairs_including_cube": len(rectangle),
                            "square_pairs_including_cube": len(square),
                            "distinct_restricted_classes": len(classes)})

    extra_controls, finite_coefficients, full_counts = new_class_controls(smallest)
    controls.update(extra_controls)
    here = Path(__file__).resolve().parent
    proof = here / "FIXED_MAP_DISTINCT_ALGEBRA_COUNT.md"
    return {
        "status": "pass", "python": platform.python_version(),
        "command": ("python3 research/arithmetic/verify_fixed_map_distinct_algebra_count.py "
                    f"--pair-bound {pair_bound} --coefficient-bound {coefficient_bound} "
                    f"--maximum-height {maximum_height}"),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
        "bounds": {"squarefree_pair_bound": pair_bound,
                   "global_coefficient_bound": coefficient_bound,
                   "maximum_height": maximum_height, "sieve_limit": limit,
                   "full_class_maximum_height": 1000,
                   "cube_multiplier_maximum": 42, "finite_prime_area_cutoff": 31},
        "checks": controls, "enumeration": enumeration,
        "finite_prime_coefficients": finite_coefficients,
        "full_class_enumeration": full_counts,
        "exact_constants": {"union_area": "55/9216", "restricted_coefficient": "55/55296",
                            "split_local_product": "1-3/p^2+2/p^3"},
        "scope": ("Exact bounded Euler, coprimality convolution, cube-class and height "
                  "controls. Cube-class computations do not independently prove the "
                  "pure-cubic field-isomorphism theorem; finite counts do not prove "
                  "the asymptotic or either uniform tail bound. The finite-prime constants are exact area sums; bounded cube searches do not prove the global minimum-height theorem."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-bound", type=int, default=500)
    parser.add_argument("--coefficient-bound", type=int, default=80)
    parser.add_argument("--maximum-height", type=int, default=32000)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run(args.pair_bound, args.coefficient_bound, args.maximum_height)
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("fixed_map_distinct_algebra_count_verification.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
