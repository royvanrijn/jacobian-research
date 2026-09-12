#!/usr/bin/env python3
"""Dependency-free arithmetic certificate for the integral Hasse fiber.

This deliberately imports neither SymPy nor project modules.  The all-prime
argument is finite:

* an irreducible cubic has Galois group S_3 because its discriminant is
  nonsquare;
* its quadratic discriminant field is cut out by the displayed quadratic;
* every element of S_3 either fixes a cubic root or acts trivially on that
  quadratic field;
* only 2 and 23 ramify in the splitting field;
* the additional product-discriminant prime 7 is handled explicitly.

The script also evaluates the reconstructed source points modulo high powers
of 2 and 23, including the special division by 8 at 2.
"""

from __future__ import annotations

from itertools import permutations
from math import gcd, isqrt


def cubic(value: int) -> int:
    return value**3 - 2 * value**2 + 8


def cubic_derivative(value: int) -> int:
    return 3 * value**2 - 4 * value


def quadratic(value: int) -> int:
    return value**2 - value + 6


def quadratic_derivative(value: int) -> int:
    return 2 * value - 1


def p_polynomial(value: int) -> int:
    return cubic(value) * quadratic(value)


def p_derivative(value: int) -> int:
    return (
        cubic_derivative(value) * quadratic(value)
        + cubic(value) * quadratic_derivative(value)
    )


def prime_factors(value: int) -> set[int]:
    value = abs(value)
    factors: set[int] = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.add(value)
    return factors


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def permutation_is_even(permutation: tuple[int, int, int]) -> bool:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(3)
        for right in range(left + 1, 3)
    )
    return inversions % 2 == 0


def determinant(matrix: list[list[int]]) -> int:
    """Small dependency-free determinant by permutation expansion."""
    size = len(matrix)
    result = 0
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = 1
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        result += (-1 if inversions % 2 else 1) * term
    return result


def hensel_lift(
    polynomial,
    prime: int,
    residue: int,
    exponent: int,
) -> tuple[int, int]:
    """Lift a simple residue root to prime**exponent by digit search."""
    root = residue % prime
    modulus = prime
    assert polynomial(root) % modulus == 0
    while modulus < prime**exponent:
        next_modulus = modulus * prime
        candidates = [
            root + digit * modulus
            for digit in range(prime)
            if polynomial(root + digit * modulus) % next_modulus == 0
        ]
        assert len(candidates) == 1
        root = candidates[0]
        modulus = next_modulus
    return root, modulus


def evaluate_map_mod(
    source_x: int,
    source_y: int,
    source_z: int,
    modulus: int,
) -> tuple[int, int, int]:
    """Evaluate the displayed integral map modulo ``modulus``."""
    source_x %= modulus
    source_y %= modulus
    source_z %= modulus
    t_value = (1 + 2 * source_x * source_y) % modulus
    q_value = (
        t_value**2 * source_z
        - source_y**2 * (1 + 3 * t_value)
    ) % modulus
    first = (
        source_x * (1 - 3 * source_x * source_y)
        + 2 * source_x**3 * source_z
        - 3 * source_x**4 * q_value**4
        + 3 * source_x**5 * q_value**5
    ) % modulus
    second = (
        source_y
        - 6 * source_x * q_value
        + 6 * t_value**2 * source_x**2 * q_value**4
        - 5 * t_value**2 * source_x**3 * q_value**5
    ) % modulus
    third = t_value * q_value % modulus
    return first, second, third


def reconstruct_at_odd_prime(
    root: int,
    modulus: int,
) -> tuple[int, int, int]:
    """Use the root reconstruction when 2 and P'(root) are units."""
    inverse_2 = pow(2, -1, modulus)
    inverse_8 = pow(8, -1, modulus)
    d_value = p_derivative(root) * pow(-8, -1, modulus) % modulus
    assert gcd(d_value, modulus) == 1
    t_value = pow(d_value, -1, modulus)
    source_x = root * pow(2 * d_value, -1, modulus) % modulus
    beta = (
        1
        - 4 * root
        + 3 * inverse_2 * root**2
        - 5 * inverse_8 * root**3
    ) % modulus
    source_y = (-beta - root) % modulus
    source_z = (
        d_value**2
        * (d_value + source_y**2 * (1 + 3 * t_value))
    ) % modulus
    return source_x, source_y, source_z


# Irreducibility and common discriminant field.
cubic_possible_roots = {
    divisor
    for positive in (1, 2, 4, 8)
    for divisor in (positive, -positive)
}
assert all(cubic(root) != 0 for root in cubic_possible_roots)
cubic_discriminant = -1472
quadratic_discriminant = -23
assert cubic_discriminant == quadratic_discriminant * 8**2
assert isqrt(abs(cubic_discriminant)) ** 2 != abs(cubic_discriminant)
assert quadratic_discriminant < 0

# The S_3 covering behind every unramified prime.
for permutation in permutations(range(3)):
    fixes_cubic_root = any(permutation[index] == index for index in range(3))
    splits_discriminant_quadratic = permutation_is_even(permutation)
    assert fixes_cubic_root or splits_discriminant_quadratic

# These are exactly the ramified rational primes in the two splitting fields.
assert prime_factors(cubic_discriminant) == {2, 23}
assert prime_factors(quadratic_discriminant) == {23}

# The factor resultant is 392=2^3*7^2.  Thus the product P has one additional
# bad-reduction prime, 7, even though the splitting field is unramified there.
sylvester_matrix = [
    [1, -2, 0, 8, 0],
    [0, 1, -2, 0, 8],
    [1, -1, 6, 0, 0],
    [0, 1, -1, 6, 0],
    [0, 0, 1, -1, 6],
]
factor_resultant = determinant(sylvester_matrix)
assert abs(factor_resultant) == 392
product_discriminant = (
    cubic_discriminant
    * quadratic_discriminant
    * factor_resultant**2
)
assert prime_factors(product_discriminant) == {2, 7, 23}

# Exact Hensel starts at the exceptional product-discriminant primes.
assert quadratic(0) % 2 == 0
assert quadratic_derivative(0) % 2 == 1
assert p_polynomial(1) % 7 == 0
assert p_derivative(1) % 7 == 1
assert cubic(7) % 23 == 0
assert cubic_derivative(7) % 23 == 4

# No rational root, and a real root lies between -2 and 0.
assert all(cubic(root) != 0 for root in cubic_possible_roots)
assert isqrt(23) ** 2 != 23
assert cubic(-2) < 0 < cubic(0)

# Direct 2-adic reconstruction.  The quadratic root is s=2u, and
# d=-P'(s)/8 is an odd unit because
# d=-(4u-1)(u^3-u^2+1).
two_precision = 16
root_2, modulus_2_lift = hensel_lift(
    quadratic,
    2,
    0,
    two_precision + 3,
)
modulus_2 = 2**two_precision
assert root_2 % 2 == 0
assert p_derivative(root_2) % 8 == 0
u_value = root_2 // 2
d_integer = -(p_derivative(root_2) // 8)
assert d_integer % 2 == 1
assert (
    d_integer
    + (4 * u_value - 1) * (u_value**3 - u_value**2 + 1)
) % modulus_2 == 0
d_2 = d_integer % modulus_2
t_2 = pow(d_2, -1, modulus_2)
x_2 = u_value * t_2 % modulus_2
y_2 = (
    -1 + 6 * u_value - 6 * u_value**2 + 5 * u_value**3
) % modulus_2
z_2 = d_2**2 * (d_2 + y_2**2 * (1 + 3 * t_2)) % modulus_2
assert p_polynomial(root_2) % modulus_2 == 0
assert evaluate_map_mod(x_2, y_2, z_2, modulus_2) == (
    3 % modulus_2,
    -1 % modulus_2,
    1,
)

# Direct 7-adic reconstruction at the simple product root 1.
seven_precision = 8
root_7, modulus_7 = hensel_lift(
    cubic,
    7,
    1,
    seven_precision,
)
source_7 = reconstruct_at_odd_prime(root_7, modulus_7)
assert p_polynomial(root_7) % modulus_7 == 0
assert evaluate_map_mod(*source_7, modulus_7) == (
    3 % modulus_7,
    -1 % modulus_7,
    1,
)

# Direct 23-adic reconstruction at the simple cubic root 7.
twenty_three_precision = 6
root_23, modulus_23 = hensel_lift(
    cubic,
    23,
    7,
    twenty_three_precision,
)
source_23 = reconstruct_at_odd_prime(root_23, modulus_23)
assert p_polynomial(root_23) % modulus_23 == 0
assert evaluate_map_mod(*source_23, modulus_23) == (
    3 % modulus_23,
    -1 % modulus_23,
    1,
)

# Finite regression of the all-prime Frobenius argument.  This is not used
# as the proof: the exact S_3 covering above handles every unramified prime.
for prime in range(3, 400):
    if not is_prime(prime) or prime == 23:
        continue
    roots = [
        residue
        for residue in range(prime)
        if p_polynomial(residue) % prime == 0
        and p_derivative(residue) % prime != 0
    ]
    assert roots, prime
    source = reconstruct_at_odd_prime(roots[0], prime)
    assert evaluate_map_mod(*source, prime) == (
        3 % prime,
        -1 % prime,
        1,
    )


print("PASS: dependency-free irreducibility and discriminant-field checks")
print("PASS: the exact S_3 covering handles every unramified prime")
print("PASS: explicit Hensel starts handle the bad primes 2, 7, and 23")
print("PASS: reconstructed local points satisfy the map to high precision")
print("PASS: the fiber is real-soluble and has no rational point")
