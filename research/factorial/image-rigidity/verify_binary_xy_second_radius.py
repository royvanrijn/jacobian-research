"""Check the finite polynomial certificates for the second-radius XY exclusion."""

from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

ZERO = (0, 0)
ONE = (1, 0)
I = (0, 1)
CONSTANT = (0, 0, 0)
A9 = (9, 0, 0)
B9 = (0, 9, 0)
C1 = (0, 0, 1)


def gaussian_product(left, right, modulus):
    a, b = left
    c, d = right
    return (a * c - b * d) % modulus, (a * d + b * c) % modulus


def gaussian_power(value, exponent, modulus):
    if exponent < 0:
        real, imag = value
        inverse_norm = pow(real**2 + imag**2, -1, modulus)
        value = real * inverse_norm, -imag * inverse_norm
        exponent = -exponent
    result = ONE
    while exponent:
        if exponent % 2:
            result = gaussian_product(result, value, modulus)
        value = gaussian_product(value, value, modulus)
        exponent //= 2
    return result


class Ring:
    """Sparse polynomials in A,B,C over Z[i,t]/(i^2+1,t^e-3,t^N)."""

    def __init__(self, precision, ramification):
        self.precision = precision
        self.ramification = ramification
        self.moduli = [
            3 ** ((precision - j + ramification - 1) // ramification)
            for j in range(ramification)
        ]

    def term(self, order=0, exponents=CONSTANT, coefficient=ONE):
        if order >= self.precision:
            return {}
        shift, power = divmod(order, self.ramification)
        value = tuple(3**shift * part % self.moduli[power] for part in coefficient)
        return {(power, *exponents): value} if any(value) else {}

    def add(self, *polynomials):
        result = {}
        for polynomial in polynomials:
            for key, (x, y) in polynomial.items():
                a, b = result.get(key, ZERO)
                modulus = self.moduli[key[0]]
                value = (a + x) % modulus, (b + y) % modulus
                if any(value):
                    result[key] = value
                else:
                    result.pop(key, None)
        return result

    def multiply(self, left, right):
        result = {}
        for (p, a, b, c), value in left.items():
            for (q, d, e, f), other in right.items():
                shift, power = divmod(p + q, self.ramification)
                key = power, a + d, b + e, c + f
                modulus = self.moduli[power]
                x, y = gaussian_product(value, other, modulus)
                r, s = result.get(key, ZERO)
                factor = 3**shift
                total = (r + factor * x) % modulus, (s + factor * y) % modulus
                if any(total):
                    result[key] = total
                else:
                    result.pop(key, None)
        return result

    def scale(self, poly, coefficient, order=0):
        return self.multiply(self.term(order, coefficient=coefficient), poly)

    def digit(self, polynomial, order):
        shift, power = divmod(order, self.ramification)
        factor = 3**shift
        result = {}
        for (index, a, b, c), (x, y) in polynomial.items():
            if index == power:
                value = x // factor % 3, y // factor % 3
                if any(value):
                    result[a, b, c] = value
        return result

    def initial(self, polynomial):
        for order in range(self.precision):
            digit = self.digit(polynomial, order)
            if digit:
                return order, digit
        return self.precision, {}


class Factorials:
    """Valuation and unit part of n!, calculated modulo a power of three."""

    def __init__(self, modulus):
        self.modulus = modulus
        self.free_products = [1]

    @lru_cache(None)
    def parts(self, n):
        if n == 0:
            return 0, 1
        quotient = n // 3
        valuation, unit = self.parts(quotient)
        while len(self.free_products) <= n:
            j = len(self.free_products)
            factor = j if j % 3 else 1
            self.free_products.append(self.free_products[-1] * factor % self.modulus)
        return valuation + quotient, unit * self.free_products[n] % self.modulus

    def ratio(self, numerator, denominator):
        modulus = self.modulus
        valuation, unit = 0, 1
        for sign, factors in ((1, numerator), (-1, denominator)):
            for n, multiplicity in factors:
                power, n_unit = self.parts(n)
                exponent = sign * multiplicity
                valuation += exponent * power
                unit = unit * pow(n_unit, exponent, modulus) % modulus
        return valuation, unit


class Moments:
    def __init__(self, start, axes, precision):
        self.start = start
        self.axes = axes
        ramification = 9 if start == 23 else 27
        self.ring = ring = Ring(precision, ramification)
        self.factorials = Factorials(ring.moduli[0])
        self.root_order = 4 * ramification // 9
        self.step = {17: 1, 20: 2, 23: 1}[start]
        coefficient_order = {17: 27, 20: 18, 23: 10}[start]
        self.b = ring.add(ring.term(coefficient=I), ring.term(self.step, (0, 1, 0)))
        k = ring.add(ring.term(self.root_order, coefficient=(-1, -1)),
                     ring.term(self.root_order + self.step, (1, 0, 0)))
        u = ring.add(ring.scale(k, (-1, 0)), ring.term(coefficient_order, C1))
        self.previous_offsets = 27 * ((precision - 1) // ramification)
        largest_index = self.previous_offsets + start + 3
        b_powers = self.powers(self.b, largest_index)
        # For XY the native coefficients with six foreign insertions
        # already have value at least 6*(4/9)+v_3(6!)>4.
        u_powers = self.powers(u, 5)
        self.native = [
            self.native_row(index, b_powers, u_powers)
            for index in range(largest_index + 1)
        ]

    def powers(self, polynomial, last):
        powers = [self.ring.term()]
        for _ in range(last):
            powers.append(self.ring.multiply(powers[-1], polynomial))
        return powers

    def native_row(self, index, b_powers, u_powers):
        qx, qy = self.axes
        ring = self.ring
        modulus = ring.moduli[0]
        terms = []
        for m in range(min(index, len(u_powers) - 1) + 1):
            for y in range(index - m + 1):
                x = index - m - y
                valuation, unit = self.factorials.ratio(
                    [(qx*x + m, 1), (qy*y + m, 1)],
                    [(x, 1), (y, 1), (m, 1), (qx, x), (qy, y)],
                )
                assert valuation >= 0
                coefficient = 3**valuation * unit % modulus
                if coefficient:
                    term = ring.multiply(b_powers[y], u_powers[m])
                    terms.append(ring.scale(term, (coefficient, 0)))
        return ring.add(*terms)

    @lru_cache(None)
    def unit_power(self, offset):
        ring = self.ring
        terms = []
        for j in range((ring.precision - 1) // self.step + 1):
            if offset >= 0:
                if j > offset:
                    break
                coefficient = comb(offset, j)
            else:
                coefficient = (-1)**j * comb(j - offset - 1, j)
            unit = gaussian_power((-1, -1), offset - j, ring.moduli[0])
            value = gaussian_product((coefficient, 0), unit, ring.moduli[0])
            terms.append(ring.term(self.step*j, (j, 0, 0), value))
        return ring.add(*terms)

    def scalar(self, high_block, offset):
        base = 27 * high_block
        if base + offset < 0:
            return {}
        ring = self.ring
        valuation, unit = self.factorials.ratio([(base, 1)], [(base + offset, 1)])
        order = self.root_order*offset + ring.ramification*valuation
        assert order >= 0
        return ring.scale(self.unit_power(offset), (unit, 0), order)

    def forms(self, high_block):
        ring = self.ring
        scalars = {
            offset: self.scalar(high_block, offset)
            for offset in range(-self.previous_offsets, self.start + 4)
        }
        rows = [
            ring.add(*(ring.multiply(scalars[offset], self.native[index - offset])
                       for offset in range(-self.previous_offsets, index + 1)))
            for index in range(self.start, self.start + 4)
        ]
        minus_one_minus_b = ring.scale(ring.add(ring.term(), self.b), (-1, 0))

        def combine(j):
            return ring.add(rows[j], ring.multiply(minus_one_minus_b, rows[j - 1]),
                            ring.multiply(self.b, rows[j - 2]))

        return rows[0], rows[1], combine(2), combine(3)


def substitute(source, target, polynomial, centers, variable_order):
    """Send the uniformizer to t^3, A,B to center+t^order times a variable."""
    assert target.ramification == 3 * source.ramification

    @lru_cache(None)
    def power(variable, exponent):
        return target.add(*(
            target.term(
                variable_order * j,
                tuple(j if d == variable else 0 for d in range(3)),
                gaussian_product(
                    (comb(exponent, j), 0),
                    gaussian_power(centers[variable], exponent - j, target.moduli[0]),
                    target.moduli[0],
                ),
            )
            for j in range(exponent + 1)
        ))

    return target.add(*(
        target.scale(target.multiply(target.multiply(power(0, a), power(1, b)),
                                     target.term(exponents=(0, 0, c))),
                     coefficient, 3*order)
        for (order, a, b, c), coefficient in polynomial.items()
    ))


def normalize_ninth_rows(ring, x, y):
    order, first = ring.initial(x)
    other_order, second = ring.initial(y)
    assert order == other_order
    assert set(first) <= {CONSTANT, A9, B9} and set(second) <= {CONSTANT, A9, B9}
    rows = []
    for pivot, eliminate in ((A9, B9), (B9, A9)):
        left, right = x, y
        p, q = first, second
        if eliminate not in q:
            left, right, p, q = right, left, q, p
        assert eliminate in q
        factor = gaussian_product(p.get(eliminate, ZERO),
                                  gaussian_power(q[eliminate], -1, 3), 3)
        row = ring.add(left, ring.scale(right, tuple(-part for part in factor)))
        level, digit = ring.initial(row)
        assert level == order and pivot in digit and set(digit) <= {pivot, CONSTANT}
        rows.append(row)
    return rows


def reduce_last_row(ring, rows):
    """Subtract multiples of the first three rows in increasing valuation order."""
    x, y, z, remainder = rows
    x, y = normalize_ninth_rows(ring, x, y)
    pivots = []
    for exponent, row in ((A9, x), (B9, y), (C1, z)):
        order, digit = ring.initial(row)
        assert exponent in digit and set(digit) <= {CONSTANT, exponent}
        inverse = gaussian_power(digit[exponent], -1, 3)
        pivots.append((exponent, order, row, inverse))

    for level in range(ring.precision):
        while True:
            digit = ring.digit(remainder, level)
            multiplier, row = reducible_term(ring, digit, level, pivots)
            if multiplier:
                remainder = ring.add(remainder, ring.multiply(multiplier, row))
            else:
                break
        if digit:
            return level, digit
    return ring.precision, {}


def reducible_term(ring, digit, level, pivots):
    for pivot, order, row, inverse in pivots:
        if order > level:
            continue
        for exponents, value in sorted(digit.items()):
            remaining_exponents = tuple(a - b for a, b in zip(exponents, pivot))
            if min(remaining_exponents) >= 0:
                factor = gaussian_product((-value[0], -value[1]), inverse, 3)
                return ring.term(level - order, remaining_exponents, factor), row
    return {}, {}


def check_raise_order(ring, rows, base_order, bound):
    """The first two equations force min(v_t(A),v_t(B)) at least bound."""
    threshold = base_order + 9*bound
    normalized = normalize_ninth_rows(ring, *rows)
    for row, pivot in zip(normalized, (A9, B9)):
        assert ring.initial(row)[0] == base_order
        assert set(ring.digit(row, base_order)) == {pivot}
        for level in range(base_order + 1, threshold):
            for a, b, _ in ring.digit(row, level):
                assert level + (a + b)*bound >= threshold


def check_centers(ring, rows, centers):
    for row, variable in zip(normalize_ninth_rows(ring, *rows), centers):
        _, digit = ring.initial(row)
        pivot = A9 if A9 in digit else B9
        value = gaussian_power(variable, 9, 3)
        x, y = gaussian_product(digit[pivot], value, 3)
        u, v = digit.get(CONSTANT, ZERO)
        assert ((x + u) % 3, (y + v) % 3) == ZERO


def configuration(start, axes):
    qx, qy = axes
    if start == 17:
        return (64, 189, (2, 2)) if qx == 3 < qy else (52, 153, ONE if qx == 3 else (2, 0))
    if start == 20:
        if qx == 3:
            return 61, 180, (2, 1) if qy == 3 else (1, 2)
        if qx == 9 < qy:
            return 100, 297, ONE
        return 88, 261, (1, 2) if qx == 9 else (2, 1)
    if qx >= 9:
        return 24, 23, (0, 2)
    return (28, 27, I) if qy == 3 else (24, 23, I)


def verify(start):
    count = 0
    for axes in combinations_with_replacement((3, 9, 27, 81), 2):
        precision, order, value = configuration(start, axes)
        ramification = 9 if start == 23 else 27
        integer_precision = (precision + ramification - 1) // ramification
        if integer_precision <= 3 and 81 in axes:
            continue
        moments = Moments(start, axes, precision)
        ring = moments.ring
        for high_block in range(3**(integer_precision - 1)):
            rows = moments.forms(high_block)
            expected_order = order
            if start in (17, 20) or (start == 23 and axes[0] == 3 and high_block % 3 == 0):
                expanded = Ring(3*precision, 3*ramification)
                centers = {17: ((-1, 1), (-1, 0)), 20: (I, ONE), 23: (ZERO, ZERO)}[start]
                bound = {17: 1, 20: 2, 23: 1}[start]
                leading_order = {17: 27, 20: 54, 23: 27}[start]
                check_centers(ring, rows[:2], centers)
                shifted = tuple(substitute(ring, expanded, row, centers, 0) for row in rows[:2])
                check_raise_order(expanded, shifted, leading_order, bound)
                rows = tuple(substitute(ring, expanded, row, centers, bound) for row in rows)
                ring_for_reduction = expanded
                if start == 23:
                    expected_order *= 3
            else:
                ring_for_reduction = ring
            result = reduce_last_row(ring_for_reduction, rows)
            assert result == (expected_order, {CONSTANT: value}), (start, axes, high_block, result)
            count += 1
        print(f"start {start}, axes {axes}: checked {3**(integer_precision - 1)} representatives")
    return count


if __name__ == "__main__":
    total = sum(verify(start) for start in (17, 20, 23))
    print(f"Verified {total} second-radius XY systems.")
