#!/usr/bin/env python3
"""Exact factorial-moment translation and sharp finite witness audit.

For the standard factorial functional

    L(U_1^a_1 ... U_n^a_n) = a_1! ... a_n!,

the script checks three related calculations.

1. Torus-diagonal extraction of the Dvorsky--Long GVC(5) polynomial gives
   an all-order sequence h_m with L(h_m)=0.  It also verifies that diagonal
   extraction is not multiplicative, so h_m is not the m-th power of h_1.
2. Cyclotomic homogeneous linear forms have exactly r-1 initial zero
   factorial moments.
3. Tensoring an odd cyclotomic filter with the Dvorsky diagonal shadow gives
   a 2r-term quartic whose first 2r-1 factorial moments vanish.

The bounded sparse and cyclotomic calculations below replay the displayed
all-order identities in the accompanying note.  They do not claim a
counterexample to the Factorial Conjecture.
"""

from __future__ import annotations

import json
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "factorial_moment_witnesses.json"

Exponent = tuple[int, ...]
BiExponent = tuple[Exponent, Exponent]
BiPolynomial = dict[BiExponent, int]
Polynomial = dict[Exponent, int]
CyclotomicElement = tuple[int, ...]


def multiply_bipolynomials(left: BiPolynomial, right: BiPolynomial) -> BiPolynomial:
    result: BiPolynomial = {}
    for (left_w, left_z), left_coefficient in left.items():
        for (right_w, right_z), right_coefficient in right.items():
            w_exponent = tuple(a + b for a, b in zip(left_w, right_w))
            z_exponent = tuple(a + b for a, b in zip(left_z, right_z))
            key = (w_exponent, z_exponent)
            result[key] = (
                result.get(key, 0) + left_coefficient * right_coefficient
            )
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def power_bipolynomial(polynomial: BiPolynomial, exponent: int) -> BiPolynomial:
    variable_count = len(next(iter(polynomial))[0])
    zero = (0,) * variable_count
    result: BiPolynomial = {(zero, zero): 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_bipolynomials(result, base)
        base = multiply_bipolynomials(base, base)
        remaining //= 2
    return result


def diagonal(polynomial: BiPolynomial) -> Polynomial:
    return {
        w_exponent: coefficient
        for (w_exponent, z_exponent), coefficient in polynomial.items()
        if w_exponent == z_exponent
    }


def multiply_polynomials(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(a + b for a, b in zip(left_exponent, right_exponent))
            result[exponent] = (
                result.get(exponent, 0) + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def factorial_functional(polynomial: Polynomial) -> int:
    return sum(
        coefficient
        * product(factorial(component) for component in exponent)
        for exponent, coefficient in polynomial.items()
    )


def product(values) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def dvorsky_bipolynomial() -> BiPolynomial:
    """Return p(w,z)=lambda(w)P(z), ordered as (t,a,b,c,d)."""
    zero = (0, 0, 0, 0, 0)
    lambda_terms = {
        (1, 1, 0, 0, 1): 1,
        (1, 0, 1, 1, 0): -1,
    }
    p_terms = {
        (1, 1, 0, 0, 1): 1,  # tad
        (0, 1, 0, 1, 1): 1,  # cad
        (2, 0, 1, 0, 0): 1,  # bt^2
        (1, 0, 1, 1, 0): 1,  # bct
    }
    return {
        (w_exponent, z_exponent): w_coefficient * z_coefficient
        for w_exponent, w_coefficient in lambda_terms.items()
        for z_exponent, z_coefficient in p_terms.items()
        if w_exponent != zero and z_exponent != zero
    }


def expected_diagonal_power(m: int) -> Polynomial:
    """Closed form for Diag(p^m)."""
    result: Polynomial = {}
    for k in range(m + 1):
        exponent = (m, m - k, k, k, m - k)
        result[exponent] = (-1) ** k * comb(m, k) ** 3
    return result


def dvorsky_shadow_moment(m: int) -> int:
    """L(A^m), A=U_t(U_a U_d-U_b U_c)."""
    total = 0
    for k in range(m + 1):
        total += (
            (-1) ** k
            * comb(m, k)
            * factorial(m)
            * factorial(m - k) ** 2
            * factorial(k) ** 2
        )
    return total


def expected_dvorsky_shadow_moment(m: int) -> int:
    if m % 2:
        return 0
    return factorial(m) ** 3 * 2 * (m + 1) // (m + 2)


def cyclotomic_reduce_prime(coefficients: list[int], prime: int) -> CyclotomicElement:
    """Reduce in Z[x]/(1+x+...+x^(prime-1)) for prime `prime`."""
    circular = [0] * prime
    for exponent, coefficient in enumerate(coefficients):
        circular[exponent % prime] += coefficient
    top = circular[-1]
    return tuple(circular[index] - top for index in range(prime - 1))


def cyclotomic_add(
    left: CyclotomicElement, right: CyclotomicElement
) -> CyclotomicElement:
    return tuple(a + b for a, b in zip(left, right))


def cyclotomic_multiply(
    left: CyclotomicElement, right: CyclotomicElement, prime: int
) -> CyclotomicElement:
    raw = [0] * (2 * prime - 3)
    for left_exponent, left_coefficient in enumerate(left):
        for right_exponent, right_coefficient in enumerate(right):
            raw[left_exponent + right_exponent] += (
                left_coefficient * right_coefficient
            )
    return cyclotomic_reduce_prime(raw, prime)


def cyclotomic_power_of_root(exponent: int, prime: int) -> CyclotomicElement:
    raw = [0] * (exponent % prime + 1)
    raw[exponent % prime] = 1
    return cyclotomic_reduce_prime(raw, prime)


def cyclotomic_linear_normalized_moments(
    prime: int, maximum_order: int
) -> list[CyclotomicElement]:
    """Compute h_m(1,zeta,...,zeta^(prime-1)) by truncated products."""
    zero = (0,) * (prime - 1)
    one = (1,) + (0,) * (prime - 2)
    series = [zero for _ in range(maximum_order + 1)]
    series[0] = one
    for root_exponent in range(prime):
        factor = [
            cyclotomic_power_of_root(root_exponent * degree, prime)
            for degree in range(maximum_order + 1)
        ]
        updated = [zero for _ in range(maximum_order + 1)]
        for left_degree in range(maximum_order + 1):
            for right_degree in range(maximum_order - left_degree + 1):
                term = cyclotomic_multiply(
                    series[left_degree], factor[right_degree], prime
                )
                updated[left_degree + right_degree] = cyclotomic_add(
                    updated[left_degree + right_degree], term
                )
        series = updated
    return series


def cyclotomic_linear_moment(order: int, root_order: int) -> int:
    return factorial(order) if order % root_order == 0 else 0


def filtered_quartic_moment(order: int, odd_root_order: int) -> int:
    return (
        dvorsky_shadow_moment(order)
        * cyclotomic_linear_moment(order, odd_root_order)
    )


def main() -> None:
    p = dvorsky_bipolynomial()
    diagonal_replay = []
    for m in range(1, 9):
        actual = diagonal(power_bipolynomial(p, m))
        expected = expected_diagonal_power(m)
        assert actual == expected
        assert factorial_functional(actual) == 0
        diagonal_replay.append(
            {
                "order": m,
                "diagonal_terms": len(actual),
                "factorial_moment": 0,
            }
        )

    shadow = diagonal(p)
    shadow_square = multiply_polynomials(shadow, shadow)
    diagonal_square = diagonal(power_bipolynomial(p, 2))
    assert shadow_square != diagonal_square
    assert factorial_functional(shadow_square) == 12
    assert factorial_functional(diagonal_square) == 0

    for m in range(1, 21):
        assert dvorsky_shadow_moment(m) == expected_dvorsky_shadow_moment(m)

    cyclotomic_replay = []
    for prime in (2, 3, 5, 7, 11):
        maximum_order = 2 * prime
        normalized = cyclotomic_linear_normalized_moments(prime, maximum_order)
        zero = (0,) * (prime - 1)
        one = (1,) + (0,) * (prime - 2)
        for m, value in enumerate(normalized):
            assert value == (one if m % prime == 0 else zero)
        cyclotomic_replay.append(
            {
                "root_order": prime,
                "term_count": prime,
                "degree": 1,
                "initial_zero_orders": list(range(1, prime)),
                "first_nonzero_order": prime,
                "first_nonzero_moment": factorial(prime),
                "checked_through": maximum_order,
            }
        )

    filtered_replay = []
    for odd_prime in (3, 5, 7, 11):
        for m in range(1, 2 * odd_prime):
            assert filtered_quartic_moment(m, odd_prime) == 0
        first_nonzero = filtered_quartic_moment(2 * odd_prime, odd_prime)
        expected_first_nonzero = (
            factorial(2 * odd_prime) ** 4
            * (2 * odd_prime + 1)
            // (odd_prime + 1)
        )
        assert first_nonzero == expected_first_nonzero
        assert first_nonzero != 0
        filtered_replay.append(
            {
                "root_order": odd_prime,
                "term_count": 2 * odd_prime,
                "degree": 4,
                "initial_zero_range": [1, 2 * odd_prime - 1],
                "first_nonzero_order": 2 * odd_prime,
                "first_nonzero_moment": str(first_nonzero),
            }
        )

    artifact = {
        "format": "factorial-moment-witnesses-v1",
        "field": "complex numbers",
        "factorial_functional": "L(U^alpha)=product_i alpha_i!",
        "source_translation": {
            "source": "Dvorsky--Long GVC(5)/SIC(5) witness",
            "p": (
                "w_t(w_a w_d-w_b w_c)(t+c)(ad+bt), "
                "with variable order (t,a,b,c,d)"
            ),
            "diagonal_formula": (
                "Diag(p^m)=U_t^m sum_(k=0)^m "
                "(-1)^k binom(m,k)^3 "
                "(U_a U_d)^(m-k)(U_b U_c)^k"
            ),
            "all_order_factorial_value": "L(Diag(p^m))=0 for every m>=1",
            "multiplicativity_obstruction": {
                "L(Diag(p)^2)": 12,
                "L(Diag(p^2))": 0,
            },
            "bounded_sparse_replay": diagonal_replay,
        },
        "linear_search": {
            "family": "G_r=sum_(j=0)^(r-1) zeta_r^j X_j",
            "all_order_formula": (
                "L(G_r^m)=m! if r divides m, and 0 otherwise"
            ),
            "minimality_scope": (
                "Among nonzero homogeneous linear forms with s nonzero terms, "
                "at most s-1 initial factorial moments vanish"
            ),
            "bounded_exact_replay": cyclotomic_replay,
        },
        "witness_derived_quartics": {
            "family": (
                "F_r=U_t(U_a U_d-U_b U_c)G_r on disjoint variables, "
                "for odd r"
            ),
            "all_order_formula": (
                "L(F_r^m)=L(A^m)L(G_r^m), "
                "L(A^m)=(m!)^3(1+(-1)^m)(m+1)/(m+2)"
            ),
            "bounded_exact_replay": filtered_replay,
        },
        "scope": {
            "factorial_conjecture": "no counterexample; each family has a nonzero later moment",
            "strong_factorial_threshold": (
                "sharp: the displayed families attain N(f)-1 initial zero moments"
            ),
            "degree_only_cutoff": (
                "fails uniformly in ambient dimension already for degree one"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS factorial translation: Diag(p^m), m=1..8")
    print("PASS factorial translation: nonmultiplicativity at m=2")
    print("PASS cyclotomic linear moments: prime orders 2,3,5,7,11")
    print("PASS witness-derived quartics: exact first nonzero moments")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
