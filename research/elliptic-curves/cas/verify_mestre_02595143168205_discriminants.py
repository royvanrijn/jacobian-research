#!/usr/bin/env python3
"""Pin the direct and split-infinity discriminant cores for family 2.

The direct polynomial is the fixed-square-content-free binary-quartic
discriminant in ``T``, divided by its integer coefficient content.  The
pullback is

    (2u)^20 * D((39146-u^2)/(2u)),

again divided by integer coefficient content.  Coefficients are stored only
transiently and hashed in ascending order as compact JSON integer lists.
PARI/GP supplies an exact Q-factorization/irreducibility replay; no numerical
root calculation is involved.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
import hashlib
import json
from math import gcd
import shutil
import subprocess
from typing import Sequence

from search_mestre_dsquare_four import FAMILIES


BASE_CHANGE_CONSTANT = 39_146
EXPECTED_DIRECT_CONTENT = 16
EXPECTED_PULLBACK_CONTENT = 4_096
EXPECTED_DIRECT_SHA256 = (
    "fc36f00ad71a6b30126402aae310cdd2c9d35553e9f22910334c2ba4b9a05590"
)
EXPECTED_PULLBACK_SHA256 = (
    "876a5e46a21c20cf531eb63469b55fe2cecf58d4fd5fdedfedacb3950a0e3a41"
)


def multiply(left: Sequence[int], right: Sequence[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return answer


def integer_content(coefficients: Sequence[int]) -> int:
    return reduce(gcd, (abs(value) for value in coefficients))


def coefficient_digest(coefficients: Sequence[int]) -> str:
    encoded = json.dumps(list(coefficients), separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def trim(coefficients: Sequence[Fraction]) -> list[Fraction]:
    answer = list(coefficients)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def polynomial_remainder(
    dividend: Sequence[Fraction], divisor: Sequence[Fraction]
) -> list[Fraction]:
    answer = trim(dividend)
    divisor = trim(divisor)
    if divisor == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    while len(answer) >= len(divisor) and answer != [0]:
        shift = len(answer) - len(divisor)
        quotient = answer[-1] / divisor[-1]
        for index, coefficient in enumerate(divisor):
            answer[index + shift] -= quotient * coefficient
        answer = trim(answer)
    return answer


def squarefree_over_q(coefficients: Sequence[int]) -> bool:
    polynomial = [Fraction(value) for value in coefficients]
    derivative = [Fraction(index * value) for index, value in enumerate(coefficients)][1:]
    left = trim(polynomial)
    right = trim(derivative)
    while right != [0]:
        left, right = right, polynomial_remainder(left, right)
    return len(trim(left)) == 1


def gp_polynomial(coefficients: Sequence[int]) -> str:
    return "+".join(
        f"({coefficient})*x^{degree}"
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    )


def exact_irreducible_over_q(coefficients: Sequence[int]) -> bool:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP is required for exact irreducibility replay")
    command = f"p={gp_polynomial(coefficients)};print(polisirreducible(p));\n"
    completed = subprocess.run(
        [executable, "-fq"],
        input=command,
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"PARI/GP wrote stderr: {completed.stderr.strip()}")
    return completed.stdout.strip() == "1"


def direct_core() -> tuple[list[int], int]:
    construction = FAMILIES[2].construction
    raw = construction.primitive_discriminant_polynomial
    if len(raw) != 21 or any(value.denominator != 1 for value in raw):
        raise AssertionError("the direct discriminant ceased to be integral degree 20")
    integers = [int(value) for value in raw]
    content = integer_content(integers)
    return [value // content for value in integers], content


def pulled_back_core(direct: Sequence[int]) -> tuple[list[int], int]:
    numerator = [BASE_CHANGE_CONSTANT, 0, -1]
    denominator = [0, 2]
    numerator_powers = [[1]]
    denominator_powers = [[1]]
    for _ in range(20):
        numerator_powers.append(multiply(numerator_powers[-1], numerator))
        denominator_powers.append(multiply(denominator_powers[-1], denominator))
    answer = [0] * 41
    for power, coefficient in enumerate(direct):
        term = multiply(
            numerator_powers[power],
            denominator_powers[20 - power],
        )
        for index, value in enumerate(term):
            answer[index] += coefficient * value
    content = integer_content(answer)
    return [value // content for value in answer], content


def replay() -> dict[str, object]:
    family = FAMILIES[2]
    if family.roots != (0, 25, 95, 143, 168, 205):
        raise AssertionError("family 2 roots changed")
    direct, direct_content = direct_core()
    pullback, pullback_content = pulled_back_core(direct)
    if direct_content != EXPECTED_DIRECT_CONTENT:
        raise AssertionError("direct discriminant content changed")
    if pullback_content != EXPECTED_PULLBACK_CONTENT:
        raise AssertionError("pullback discriminant content changed")
    if integer_content(direct) != 1 or integer_content(pullback) != 1:
        raise AssertionError("a declared primitive discriminant has nontrivial content")
    if len(direct) != 21 or len(pullback) != 41:
        raise AssertionError("a discriminant degree changed")
    if any(direct[1::2]) or any(pullback[1::2]):
        raise AssertionError("an even discriminant acquired an odd coefficient")
    if coefficient_digest(direct) != EXPECTED_DIRECT_SHA256:
        raise AssertionError("direct coefficient digest changed")
    if coefficient_digest(pullback) != EXPECTED_PULLBACK_SHA256:
        raise AssertionError("pullback coefficient digest changed")
    if not squarefree_over_q(direct) or not squarefree_over_q(pullback):
        raise AssertionError("a discriminant core ceased to be squarefree")
    direct_irreducible = exact_irreducible_over_q(direct)
    pullback_irreducible = exact_irreducible_over_q(pullback)
    if not direct_irreducible or not pullback_irreducible:
        raise AssertionError("a discriminant core ceased to be irreducible over Q")
    return {
        "status": "verified exact primitive discriminant cores",
        "family_roots": list(family.roots),
        "base_change": "T=(39146-u^2)/(2u)",
        "coefficient_hash_encoding": "SHA256(compact JSON ascending integer list)",
        "direct": {
            "degree": 20,
            "raw_integer_content": direct_content,
            "primitive_content": integer_content(direct),
            "squarefree_over_Q": True,
            "irreducible_over_Q": direct_irreducible,
            "factor_degrees_over_Q": [20],
            "coefficient_sha256": coefficient_digest(direct),
        },
        "pullback": {
            "degree": 40,
            "cleared_integer_content": pullback_content,
            "primitive_content": integer_content(pullback),
            "squarefree_over_Q": True,
            "irreducible_over_Q": pullback_irreducible,
            "factor_degrees_over_Q": [40],
            "coefficient_sha256": coefficient_digest(pullback),
        },
    }


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
