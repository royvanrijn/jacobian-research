"""Replay the branch and lattice arithmetic in the CM moment face note."""

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import comb, gcd


@dataclass(frozen=True, order=True)
class Face:
    degree: int
    n: int
    q: int
    rho_numerator: int
    left: int
    interior: tuple[int, ...]
    right: int

    @property
    def angles(self):
        return range(self.left - self.n, self.n - self.right + 1, self.q)

    @property
    def branches(self):
        return (self.left, *self.interior, self.right)

    def height(self, angle):
        return Fraction(self.n - self.rho_numerator * angle, self.q)


@cache
def partitions(total, length, smallest, largest):
    if length == 0:
        return ((),) if total == 0 else ()
    return tuple(
        (part, *tail)
        for part in range(max(smallest, total - (length - 1) * largest),
                          min(largest, total // length) + 1)
        for tail in partitions(total - part, length - 1, part, largest)
    )


def unit_characters(n):
    return tuple(i for i in range(1, n) if gcd(i, n) == 1)


def coherent_dimensions(n, branches, characters):
    if not all(0 < exponent < n for exponent in branches):
        raise ValueError("Reduced nonzero branch exponents are required.")
    if sum(branches) != 2 * n or gcd(n, *branches) != 1:
        raise ValueError("The cover must be connected and have branch sum 2n.")
    dimensions = tuple(
        sum(character * exponent % n for exponent in branches) // n - 1
        for character in characters
    )
    if any(value + dual != len(branches) - 2
           for value, dual in zip(dimensions, reversed(dimensions))):
        raise ValueError("Conjugate coherent dimensions do not sum to the rank.")
    return dimensions


def minimum_face(n, branches):
    values = set()
    for i, j in permutations(range(len(branches)), 2):
        interior = tuple(value for k, value in enumerate(branches) if k not in (i, j))
        values.add((branches[i], interior, branches[j]))
    for left, interior, right in sorted(values):
        b, a = n - left, n - right
        for q in range(1, gcd(*interior) + 1):
            if gcd(q, n * b) != 1 or any(value % q for value in interior):
                continue
            lower = max(q + 1, 2 * q - n // b)
            for rho in range(lower, n // a + 1):
                if (n + rho * b) % q == 0:
                    degree = (n + (rho - q) * b) // q
                    yield Face(degree, n, q, rho, left, interior, right)
                    break


def lattice_check(face):
    heights = tuple(face.height(angle) for angle in face.angles)
    assert face.rho_numerator > face.q
    assert gcd(face.n, face.q) == gcd(face.rho_numerator, face.q) == 1
    assert face.angles.stop - 1 == face.angles[-1]
    assert all(height.denominator == 1 and height >= max(0, -2 * angle)
               for angle, height in zip(face.angles, heights))
    assert face.degree == max(angle + height for angle, height in zip(face.angles, heights))
    assert sum(face.interior) == face.angles[-1] - face.angles[0]


def convolve(first, second):
    output = defaultdict(Fraction)
    for (angle_a, height_a), a in first.items():
        for (angle_b, height_b), b in second.items():
            output[angle_a + angle_b, height_a + height_b] += a * b
    return output


def multiply_by_power(polynomial, right, exponent):
    for _ in range(exponent):
        polynomial = convolve(polynomial, right)
    return polynomial


def check_prime_coefficient(face, root_multiplicities, prime):
    polynomial = {(0, 0): Fraction(1)}
    for root, multiplicity in root_multiplicities:
        polynomial = multiply_by_power(polynomial, {(1, 0): 1, (0, 0): -root},
                                       multiplicity)
    b = face.n - face.left
    laurent = defaultdict(Fraction)
    for (power, _), value in polynomial.items():
        angle = face.q * power - b
        height = int(face.height(angle))
        if angle >= 0:
            laurent[angle, height] += value
        else:
            negative_angle = -angle
            for exponent in range(negative_angle + 1):
                actual_height = height - 2 * negative_angle + 2 * exponent
                lift_coefficient = (-1) ** (negative_angle + exponent) * comb(
                    negative_angle, exponent
                )
                laurent[angle, actual_height] += value * lift_coefficient
    laurent[0, 0] += Fraction(7, 2)
    laurent[-1, 2] -= Fraction(3, 2)
    laurent[-1, 0] += Fraction(3, 2)
    k = (prime - 1) // face.n
    assert prime == face.n * k + 1
    power = multiply_by_power({(0, 0): Fraction(1)}, laurent, face.q * k)
    moment = sum(value / (height + 1)
                 for (angle, height), value in power.items()
                 if angle == 0 and height % 2 == 0)
    coefficient = multiply_by_power({(0, 0): Fraction(1)}, polynomial, face.q * k)[b * k, 0]
    difference = prime * moment - coefficient
    assert difference.denominator % prime != 0 and difference.numerator % prime == 0
    return int(coefficient) % prime


def check_examples():
    degree_six = Face(6, 4, 1, 2, 2, (1, 1, 1, 1), 2)
    degree_eight = Face(8, 10, 3, 5, 3, (3, 3, 3), 8)
    order_nine = Face(13, 9, 1, 2, 5, (2, 3), 8)
    for face, expected in [(degree_six, (1, 3)), (degree_eight, (1, 3, 0, 2)),
                           (order_nine, (1, 1, 1, 1, 1, 1))]:
        lattice_check(face)
        assert coherent_dimensions(face.n, face.branches,
                                   unit_characters(face.n)) == expected
    orbit = tuple(pow(4, index, 9) for index in range(3))
    assert len({min(i, 9 - i) for i in orbit}) == 3
    assert pow(4, 3, 9) == 1
    assert tuple(4 * exponent % 9 for exponent in (5, 2, 8, 3)) == (2, 8, 5, 3)
    k_values = (
        check_prime_coefficient(degree_six, [(1, 1), (2, 1), (3, 1), (4, 1)], 13),
        check_prime_coefficient(degree_eight, [(1, 1), (2, 1), (4, 1)], 11),
        check_prime_coefficient(order_nine, [(1, 2), (8, 3)], 19),
    )
    assert k_values[2] != 0 and (8 * 8 - 8 + 1) % 19 == 0
    print("Prime congruences at p = 13, 11, 19; isolated residues", k_values)


def check_single_bad_pair():
    prime = 19
    k = (prime - 1) // 9
    values = []
    for character in unit_characters(9):
        left, first, second = (k * (character * exponent % 9) for exponent in (5, 2, 3))
        target = prime - 1 - left
        coefficient = sum(
            comb(first, j) * (-1) ** (first - j) *
            comb(second, target - j) * (-2) ** (second - target + j)
            for j in range(max(0, target - second), min(first, target) + 1)
        )
        values.append(coefficient % prime)
    assert values == [0, 7, 16, 7, 5, 0]
    print("Order-nine characters at p = 19, lambda = 2:", tuple(values))


def census():
    eligible = 0
    best = {}
    for n in range(3, 31):
        characters = unit_characters(n)
        e = len(characters) // 2
        for rank in range(3, 30 // e + 1):
            for branches in partitions(2 * n, rank + 2, 1, n - 1):
                if gcd(n, *branches) != 1:
                    continue
                signature = coherent_dimensions(n, branches, characters)
                if sum(0 < value < rank for value in signature[:e]) != 1:
                    continue
                candidate = min(minimum_face(n, branches), default=None)
                if candidate is None:
                    continue
                lattice_check(candidate)
                eligible += 1
                if (n, rank) not in best or candidate < best[n, rank]:
                    best[n, rank] = candidate
    return eligible, best


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    brute = {row for row in product(range(1, 6), repeat=4)
             if tuple(sorted(row)) == row and sum(row) == 12}
    assert set(partitions(12, 4, 1, 5)) == brute
    check_examples()
    check_single_bad_pair()
    eligible, best = census()
    print("Feasible branch multisets of single-mixed-pair type:", eligible)
    print("Distinct (order, relative rank) combinations:", len(best))
    for (n, rank), face in sorted(best.items()):
        print(n, rank, face.degree, str(Fraction(n, face.q)),
              str(Fraction(face.rho_numerator, face.q)), face.branches)
    assert eligible == 69 and len(best) == 20
    assert best[4, 4] == Face(6, 4, 1, 2, 2, (1, 1, 1, 1), 2)
    assert best[10, 3] == Face(8, 10, 3, 5, 3, (3, 3, 3), 8)


if __name__ == "__main__":
    main()
