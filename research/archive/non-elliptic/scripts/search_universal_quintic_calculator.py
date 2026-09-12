#!/usr/bin/env python3
"""Bounded projective-height search for the five quintic Galois groups.

Enumerate primitive targets ``[W:P:B:C]`` for

    E = (P/W)^5*S^5 - 5(P/W)*S^3 - 2(B/W)*S^2 + 4*S - 2(C/W).

Cheap exact discriminant and Frobenius-pattern filters leave a small list
for PARI/GP's ``polgalois``.  The search is discovery/minimality evidence;
``verify_universal_quintic_calculator.py`` checks the selected rows without
using a Galois-group oracle.
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import subprocess
import warnings
from dataclasses import dataclass
from functools import reduce

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


SCREEN_PRIMES = (
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
    101,
    103,
    107,
    109,
    113,
    127,
    131,
    137,
    139,
    149,
    151,
    157,
    163,
    167,
    173,
)


@dataclass(frozen=True)
class Target:
    w: int
    pi: int
    b: int
    c: int

    @property
    def height(self) -> int:
        return max(self.w, abs(self.pi), abs(self.b), abs(self.c))

    @property
    def primitive_inverse_coefficients(self) -> tuple[int, ...]:
        coefficients = (
            self.pi**5,
            0,
            -5 * self.pi * self.w**4,
            -2 * self.b * self.w**4,
            4 * self.w**5,
            -2 * self.c * self.w**4,
        )
        content = reduce(math.gcd, (abs(value) for value in coefficients))
        coefficients = tuple(value // content for value in coefficients)
        if coefficients[0] < 0:
            coefficients = tuple(-value for value in coefficients)
        return coefficients

    @property
    def coefficient_height(self) -> int:
        return max(abs(value) for value in self.primitive_inverse_coefficients)


def discriminant_square_class(target: Target) -> int:
    """Square class of the discriminant, after removing ``(4 P^4 W^8)^2``."""

    w, p, b, c = target.w, target.pi, target.b, target.c
    return (
        432 * b**5 * c * p**2 * w**8
        - 432 * b**4 * p**2 * w**10
        + 12600 * b**3 * c * p**3 * w**9
        - 2000 * b**3 * c * w**12
        + 9000 * b**2 * c**2 * p**7 * w**5
        + 20625 * b**2 * c**2 * p**4 * w**8
        - 11520 * b**2 * p**3 * w**11
        + 2000 * b**2 * w**14
        + 18750 * b * c**3 * p**8 * w**4
        - 25600 * b * c * p**7 * w**7
        + 56000 * b * c * p**4 * w**10
        - 45000 * b * c * p * w**13
        + 3125 * c**4 * p**12
        - 40000 * c**2 * p**8 * w**6
        + 112500 * c**2 * p**5 * w**9
        - 84375 * c**2 * p**2 * w**12
        + 16384 * p**7 * w**9
        - 51200 * p**4 * w**12
        + 40000 * p * w**15
    )


def is_square(value: int) -> bool:
    if value < 0:
        return False
    root = math.isqrt(value)
    return root * root == value


def trim(poly: list[int]) -> list[int]:
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_remainder(dividend: list[int], divisor: list[int], prime: int) -> list[int]:
    dividend = dividend[:]
    trim(dividend)
    divisor = divisor[:]
    trim(divisor)
    inverse_lead = pow(divisor[-1], -1, prime)
    while len(dividend) >= len(divisor):
        quotient = dividend[-1] * inverse_lead % prime
        shift = len(dividend) - len(divisor)
        for index, coefficient in enumerate(divisor):
            dividend[shift + index] = (
                dividend[shift + index] - quotient * coefficient
            ) % prime
        trim(dividend)
    return dividend


def poly_multiply_mod(
    left: list[int], right: list[int], modulus: list[int], prime: int
) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            product[i + j] = (
                product[i + j] + left_coefficient * right_coefficient
            ) % prime
    return poly_remainder(product, modulus, prime)


def x_power_mod(exponent: int, modulus: list[int], prime: int) -> list[int]:
    result = [1]
    base = [0, 1]
    while exponent:
        if exponent & 1:
            result = poly_multiply_mod(result, base, modulus, prime)
        base = poly_multiply_mod(base, base, modulus, prime)
        exponent >>= 1
    return result


def gcd_degree(left: list[int], right: list[int], prime: int) -> int:
    trim(left)
    trim(right)
    while right:
        left, right = right, poly_remainder(left, right, prime)
    return len(left) - 1


def local_f20_allowed(coefficients: tuple[int, ...], prime: int) -> bool:
    """Test whether a squarefree factor pattern is possible in F_20."""

    polynomial = [coefficient % prime for coefficient in reversed(coefficients)]
    trim(polynomial)
    if len(polynomial) != 6:
        return True
    derivative = [
        index * polynomial[index] % prime for index in range(1, len(polynomial))
    ]
    if gcd_degree(polynomial[:], derivative, prime) > 0:
        return True

    x_prime = x_power_mod(prime, polynomial, prime)
    if len(x_prime) < 2:
        x_prime.extend([0] * (2 - len(x_prime)))
    x_prime[1] = (x_prime[1] - 1) % prime
    linear_degree = gcd_degree(polynomial[:], x_prime, prime)
    if linear_degree in (2, 3):
        return False  # (3,1,1) or (2,1,1,1)

    if linear_degree == 0:
        x_prime_squared = x_power_mod(prime**2, polynomial, prime)
        if len(x_prime_squared) < 2:
            x_prime_squared.extend([0] * (2 - len(x_prime_squared)))
        x_prime_squared[1] = (x_prime_squared[1] - 1) % prime
        if gcd_degree(polynomial[:], x_prime_squared, prime) > 0:
            return False  # (3,2)
    return True


def projective_targets(bound: int):
    """Enumerate targets modulo ``(B,C) -> (-B,-C)``."""

    for w in range(1, bound + 1):
        for pi, b, c in itertools.product(
            range(-bound, bound + 1),
            range(-bound, bound + 1),
            range(-bound, bound + 1),
        ):
            if not pi or b < 0 or (b == 0 and c < 0):
                continue
            target = Target(w, pi, b, c)
            if target.height > bound:
                continue
            if reduce(math.gcd, (w, abs(pi), abs(b), abs(c))) != 1:
                continue
            yield target


def pari_classify(candidates: set[Target], gp: str) -> dict[Target, str]:
    def target_key(item: Target) -> tuple[int, int, int, int]:
        return item.w, item.pi, item.b, item.c

    lines = []
    ordered = sorted(candidates, key=target_key)
    for index, target in enumerate(ordered):
        coefficients = target.primitive_inverse_coefficients
        polynomial = " + ".join(
            f"({coefficient})*x^{5 - position}"
            for position, coefficient in enumerate(coefficients)
            if coefficient
        )
        lines.append(
            f'f={polynomial};if(polisirreducible(f),'
            f'print("{index}|",polgalois(f)[4]))'
        )
    completed = subprocess.run(
        [gp, "-q"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    result: dict[Target, str] = {}
    for line in completed.stdout.splitlines():
        index_text, pari_name = line.split("|", 1)
        target = ordered[int(index_text)]
        if pari_name == "S5":
            result[target] = "S5"
        elif pari_name == "A5":
            result[target] = "A5"
        elif pari_name.startswith("C(5)"):
            result[target] = "C5"
        elif pari_name.startswith("D(5)"):
            result[target] = "D5"
        elif pari_name.startswith("F(5)"):
            result[target] = "F20"
    return result


def s5_witness_score(target: Target) -> tuple[int, int, int]:
    variable = sp.symbols("x")
    coefficients = target.primitive_inverse_coefficients
    expression = sum(
        coefficient * variable ** (5 - index)
        for index, coefficient in enumerate(coefficients)
    )
    polynomial = sp.Poly(expression, variable, domain=sp.ZZ)
    discriminant = int(sp.discriminant(expression, variable))

    def pattern(prime: int) -> tuple[int, ...] | None:
        if int(polynomial.LC()) % prime == 0 or discriminant % prime == 0:
            return None
        factors = sp.factor_list(expression, modulus=prime)[1]
        return tuple(
            sorted(
                (
                    int(sp.degree(factor, variable))
                    for factor, exponent in factors
                    for _ in range(exponent)
                ),
                reverse=True,
            )
        )

    primes = list(sp.primerange(2, 500))
    irreducible = next(prime for prime in primes if pattern(prime) == (5,))
    transposition = next(
        prime for prime in primes if pattern(prime) == (2, 1, 1, 1)
    )
    return max(irreducible, transposition), irreducible + transposition, irreducible


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=21)
    parser.add_argument("--gp", default="gp")
    args = parser.parse_args()
    assert args.bound >= 1

    square_candidates: set[Target] = set()
    f20_candidates: set[Target] = set()
    s5_candidates: set[Target] = set()
    examined = 0

    for target in projective_targets(args.bound):
        examined += 1
        if target.height == 1:
            s5_candidates.add(target)
        square_class = discriminant_square_class(target)
        if not square_class:
            continue
        if square_class > 0 and is_square(square_class):
            square_candidates.add(target)
            continue
        if square_class < 0:
            continue
        coefficients = target.primitive_inverse_coefficients
        if all(
            local_f20_allowed(coefficients, prime) for prime in SCREEN_PRIMES
        ):
            f20_candidates.add(target)

    candidates = square_candidates | f20_candidates | s5_candidates
    classifications = pari_classify(candidates, args.gp)
    by_group: dict[str, list[Target]] = {
        name: [] for name in ("S5", "A5", "C5", "D5", "F20")
    }
    for target, group in classifications.items():
        by_group[group].append(target)

    print(f"examined {examined} primitive targets modulo (B,C)~(-B,-C)")
    print(
        f"PARI candidates: square={len(square_candidates)}, "
        f"F20-screen={len(f20_candidates)}, height-one={len(s5_candidates)}"
    )
    for group, entries in by_group.items():
        if not entries:
            print(f"{group}: not found through height {args.bound}")
            continue
        if group == "S5":
            best = min(
                entries,
                key=lambda target: (
                    target.height,
                    target.coefficient_height,
                    s5_witness_score(target),
                ),
            )
        else:
            best = min(
                entries,
                key=lambda target: (target.height, target.coefficient_height),
            )
        print(
            f"{group}: [{best.w}:{best.pi}:{best.b}:{best.c}], "
            f"height={best.height}, coefficient_height={best.coefficient_height}"
        )


if __name__ == "__main__":
    main()
