"""Replay the exact cubic field arithmetic for the fixed 563 example.

The Hecke trace inputs are the three entries retained in
CM_NEWTON_DISCONNECTED_RM_ORDINARITY.md.
"""

from itertools import product


def multiply(left, right):
    coefficients = [0] * 5
    for i, j in product(range(3), repeat=2):
        coefficients[i + j] += left[i] * right[j]
    for i in (4, 3):
        leading = coefficients[i]
        coefficients[i - 1] += leading
        coefficients[i - 2] += 6 * leading
        coefficients[i - 3] -= 7 * leading
    return tuple(coefficients[:3])


def trace(element):
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    return sum(multiply(element, vector)[i] for i, vector in enumerate(basis))


def symmetric_sums(element):
    square = multiply(element, element)
    first = trace(element)
    second_power = trace(square)
    third_power = trace(multiply(square, element))
    second = (first**2 - second_power) // 2
    third = (first**3 - 3 * first * second_power + 2 * third_power) // 6
    return first, second, third


def frobenius_polynomial(element, prime):
    first, second, third = symmetric_sums(element)
    return (1, -first, 3 * prime + second, -2 * prime * first - third,
            3 * prime**2 + prime * second, -prime**2 * first, prime**3)


def evaluate(coefficients, value):
    return sum(coefficient * value**i for i, coefficient in enumerate(coefficients))


def roots(coefficients, prime):
    return tuple(value for value in range(prime) if evaluate(coefficients, value) % prime == 0)


def main():
    coefficient_polynomial = (7, -6, -1, 1)
    base_polynomial = (-8, -14, -1, 1)
    for prime, element, middle, norm, auxiliary, discriminant_root in (
        (2, (0, 1, 0), 3, -31, 31, 16),
        (11, (0, 0, -1), 335, -7 * 449, 7, 5),
    ):
        assert frobenius_polynomial(element, prime)[3] == middle
        assert middle % prime != 0
        discriminant = list(multiply(element, element))
        discriminant[0] -= 4 * prime
        assert symmetric_sums(discriminant)[2] == norm
        coefficient_roots = roots(coefficient_polynomial, auxiliary)
        assert len(coefficient_roots) == 3
        assert tuple(x for x in coefficient_roots if evaluate(discriminant, x) % auxiliary == 0
                     ) == (discriminant_root,)
        print(f"At {prime}, ordinary reduction and relative discriminant norm {norm}.")

    prime = 563
    coefficient_roots = roots(coefficient_polynomial, prime)
    assert coefficient_roots == (72, 113, 379)
    assert roots(base_polynomial, prime) == (145, 157, 262)
    for base_root, prime_generator, element, residues in (
        (262, (23, 35, -9), (-57, 2, 13), (482, 81, 0)),
        (157, (-171, -21, 13), (5, 11, -2), (0, 482, 81)),
        (145, (-9, -14, -4), (52, -13, -11), (81, 0, 482)),
    ):
        assert tuple(x for x in roots(base_polynomial, prime)
                     if evaluate(prime_generator, x) % prime == 0) == (base_root,)
        assert tuple(evaluate(element, x) % prime for x in coefficient_roots) == residues
        assert symmetric_sums(element) == (0, -931, -19 * prime)
        polynomial = frobenius_polynomial(element, prime)
        assert polynomial == (1, 0, 758, 10697, 426754, 0, 178453547)
        assert tuple(coefficient % prime for coefficient in polynomial) == (1, 0, 195, 0, 0, 0, 0)
        print(f"At {prime} with w={base_root}, the trace residues are {residues} and p-rank is two.")


if __name__ == "__main__":
    main()
