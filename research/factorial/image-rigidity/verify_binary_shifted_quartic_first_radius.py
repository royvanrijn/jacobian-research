"""First-radius certificates for the two-adjacent-one foreign exponent."""

from itertools import product
from math import comb, factorial

from verify_binary_two_digit_second_radius import (
    TwoDigitMoments, factorial_unit, factorial_value, free_products,
)
from verify_binary_xy_second_radius import (
    C1, CONSTANT, I, ONE, ZERO, Factorials, Ring, gaussian_power,
    gaussian_product, reducible_term,
)


A3 = (3, 0, 0)
B3 = (0, 3, 0)


def g(n, modulus):
    blocks, remainder = divmod(n, modulus)
    return (-1)**blocks * free_products(modulus)[remainder] % modulus


def axis_binomial(height, exponent, m, x, accuracy, count_digits):
    modulus = 3**accuracy
    n = m*exponent
    d = (n // 3**height) % 3**(accuracy + count_digits)
    coefficient = comb(x + d, x) % modulus
    for j in range(max(0, height - accuracy + 1), height):
        axis = x * 3**(height - j)
        foreign = n // 3**j
        unit = g(axis + foreign, modulus)
        denominator = g(axis, modulus) * g(foreign, modulus)
        coefficient = coefficient * unit * pow(denominator, -1, modulus) % modulus
    return coefficient


def foreign_unit(m, exponent, modulus):
    valuation = factorial_value(m*exponent) - m*factorial_value(exponent) - 2*factorial_value(m)
    assert valuation >= 0
    denominator = pow(factorial_unit(exponent, modulus), m, modulus) * factorial_unit(m, modulus)**2
    return 3**valuation * factorial_unit(m*exponent, modulus) * pow(denominator, -1, modulus) % modulus


class FirstRadiusMoments(TwoDigitMoments):
    def __init__(self, x_height, y_height, exponent_shift, accuracy):
        self.start = 5
        self.ring = ring = Ring(3*accuracy, 3)
        self.factorials = Factorials(ring.moduli[0])
        self.previous_offsets = 9*(accuracy - 1)
        largest_index = 8 + self.previous_offsets
        self.b = ring.add(ring.term(coefficient=I), ring.term(1, (0, 1, 0)))
        u = ring.add(ring.term(1, coefficient=(-1, -1)),
                     ring.term(2, (1, 0, 0), (-1, 0)), ring.term(4, C1))
        exponent = 4 * 3**exponent_shift
        self.native, self.d_native = self.build_rows(
            x_height, y_height, exponent, largest_index,
            self.powers(self.b, largest_index), self.powers(u, 8),
        )
        self.kernel_cache = {}
        self.unit_powers = {}

    def unit_power(self, offset):
        ring = self.ring
        modulus = ring.moduli[0]
        if offset not in self.unit_powers:
            terms = []
            for j in range(ring.precision):
                if offset >= 0:
                    if j > offset:
                        break
                    binomial = comb(offset, j)
                else:
                    binomial = (-1)**j * comb(j - offset - 1, j)
                unit = gaussian_power((1, 1), offset - j, modulus)
                terms.append(ring.term(j, (j, 0, 0),
                                      gaussian_product((binomial, 0), unit, modulus)))
            self.unit_powers[offset] = ring.add(*terms)
        return self.unit_powers[offset]

    def scalar_order(self, high_block, offset):
        base = 9*high_block
        if base + offset < 0:
            return None, 0
        value, unit = self.factorials.ratio([(base, 1)], [(base + offset, 1)])
        order = offset + 3*value
        assert order >= 0
        return order, unit

    def build_rows(self, x_height, y_height, exponent, largest, b_powers, u_powers):
        ring = self.ring
        native, d_native = ([{} for _ in range(largest + 1)] for _ in range(2))
        pure = [1]
        for j in range(largest):
            pure.append(pure[-1] * ((3*j + 1)*(3*j + 2)//2) % ring.moduli[0])
        count_digits = 0
        while 3**(count_digits + 1) <= largest:
            count_digits += 1
        for m, u_power in enumerate(u_powers):
            remaining = ring.precision - m - 3*factorial_value(m)
            if remaining <= 0:
                continue
            accuracy = (remaining + 2) // 3
            modulus = 3**accuracy
            factor = foreign_unit(m, exponent, modulus)
            if not factor:
                continue
            x_row = [(pure[x] if x_height else 1)
                     * axis_binomial(x_height, exponent, m, x, accuracy, count_digits) % modulus
                     for x in range(largest - m + 1)]
            y_row = [(pure[y] if y_height else 1) % modulus for y in range(largest - m + 1)]
            x_d = [x_row[0], *((b-a) % modulus for a,b in zip(x_row, x_row[1:]))]
            y_d = [y_row[0], *((b-a) % modulus for a,b in zip(y_row, y_row[1:]))]
            for destination, x_axis, y_axis in ((native, x_row, y_row), (d_native, x_d, y_d)):
                for n in range(m, largest + 1):
                    terms = []
                    for y in range(n-m+1):
                        coefficient = x_axis[n-m-y] * y_axis[y] * factor % modulus
                        if coefficient:
                            terms.append(ring.scale(b_powers[y], (coefficient, 0)))
                    contribution = ring.scale(ring.add(*terms), (factorial(m), 0))
                    destination[n] = ring.add(destination[n], ring.multiply(contribution, u_power))
        return native, d_native


def normalized_cubes(ring, x, y):
    order, first = ring.initial(x)
    other_order, second = ring.initial(y)
    assert order == other_order == 3
    assert set(first) <= {CONSTANT, A3, B3} and set(second) <= {CONSTANT, A3, B3}
    normalized = []
    for target, remove in ((A3, B3), (B3, A3)):
        base, head, other, other_head = (y, second, x, first) if remove in second else (
            x, first, y, second)
        inverse = gaussian_power(head[remove], -1, 3)
        ratio = gaussian_product(other_head.get(remove, ZERO), inverse, 3)
        row = ring.add(other, ring.scale(base, tuple(-part for part in ratio)))
        level, digit = ring.initial(row)
        assert level == order and target in digit and set(digit) <= {target, CONSTANT}
        normalized.append(row)
    return tuple(normalized)


def at_residues(digit, residues):
    value = ZERO
    for exponents, coefficient in digit.items():
        for exponent, residue in zip(exponents, residues):
            coefficient = gaussian_product(coefficient, gaussian_power(residue, exponent, 3), 3)
        value = tuple((a+b) % 3 for a,b in zip(value, coefficient))
    return value


def filter_rows(system, block):
    ring = system.ring
    e5, e6, z7, remainder = system.forms(block)
    x, y = normalized_cubes(ring, e5, e6)
    pivots, residues = [], []
    for exponent, row in ((A3, x), (B3, y), (C1, z7)):
        order, digit = ring.initial(row)
        assert exponent in digit and set(digit) <= {CONSTANT, exponent}
        inverse = gaussian_power(digit[exponent], -1, 3)
        constant = digit.get(CONSTANT, ZERO)
        residue = gaussian_product(tuple(-part for part in constant), inverse, 3)
        residues.append(gaussian_power(residue, 3, 3) if exponent != C1 else residue)
        pivots.append((exponent, order, row, inverse))
    for order in range(ring.precision):
        while True:
            digit = ring.digit(remainder, order)
            multiplier, row = reducible_term(ring, digit, order, pivots)
            if not multiplier:
                break
            remainder = ring.add(remainder, ring.multiply(multiplier, row))
        if digit:
            if at_residues(digit, residues) != ZERO:
                return False
            assert order == 11 and digit.get((0, 1, 0)) == ONE, (block, order, digit)
            assert set(digit) <= {CONSTANT, (0, 1, 0)}
            assert at_residues(ring.digit(y, 4), residues) != ZERO
            return False
    return True


def verify():
    low = [(r+gap, other, r) for r,gap,other in product((0, 1), (3, 4, 5), (0, 1))]
    upper = [(q, other, q+gap) for q,gap,other in product((0, 1), range(2, 8), (0, 1))]
    retained = 0
    for qx,qy,r in [*low, *upper]:
        system = FirstRadiusMoments(qx, qy, r, 4)
        for block in range(27):
            keep = filter_rows(system, block)
            expected = qx == qy == 0 and r >= 3 and block % 9 == (7 if r == 3 else 1)
            assert keep == expected, (qx, qy, r, block, keep)
            retained += keep
    print(f"Precision four: {len(low+upper)*27-retained} excluded, {retained} retained.")

    tested, retained = 0, 0
    for r in range(3, 11):
        system = FirstRadiusMoments(0, 0, r, 7)
        low_residue = 7 if r == 3 else 1
        for block in range(low_residue, 729, 9):
            keep = filter_rows(system, block)
            expected = r >= 6 and block % 243 == (136 if r == 6 else 217)
            assert keep == expected, (r, block, keep)
            tested += 1
            retained += keep
    print(f"Precision seven: {tested-retained} excluded, {retained} retained.")


if __name__ == "__main__":
    verify()
