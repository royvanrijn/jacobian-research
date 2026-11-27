"""Exact four-row certificates for the second radius with two ternary ones."""

from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import product
from math import factorial

from verify_binary_low_digit_second_radius import form_reduction
from verify_binary_xy_second_radius import (
    C1, I, ONE, ZERO, Factorials, Moments, Ring, check_centers,
    check_raise_order, substitute,
)


SHAPES = ((1, 1), (0, 2))
CUTOFFS = {17: (27,), 20: (27, 36, 38), 23: (27, 36, 37)}


def ceiling_ninths(ninths):
    return (ninths + 8) // 9


def factorial_value(n):
    valuation = 0
    while n:
        n //= 3
        valuation += n
    return valuation


@lru_cache(None)
def free_products(modulus):
    products = [1]
    for j in range(1, modulus):
        products.append(products[-1] * (j if j % 3 else 1) % modulus)
    return products


@lru_cache(None)
def factorial_unit(n, modulus):
    products = free_products(modulus)
    unit = 1
    while n:
        blocks, residue = divmod(n, modulus)
        sign = -1 if blocks % 2 else 1
        unit = unit * sign * products[residue] % modulus
        n //= 3
    return unit


def digit_words(digit_count, gap_cap):
    """Yield the axis height and the exponent, including leading zero digits."""
    for gaps in product(range(gap_cap + 1), repeat=digit_count + 1):
        places = [sum(gaps[:j + 1]) + j for j in range(digit_count)]
        yield sum(gaps) + digit_count, sum(3**j for j in places)


@lru_cache(None)
def insertion_polynomial(height, exponent, insertions, accuracy):
    """The cross binomial as a polynomial in the axis count modulo 3^accuracy."""
    modulus = 3**accuracy
    cut = max(height - accuracy + 1, 0)
    upper = insertions * exponent // 3**cut
    axis = 3**(height - cut)
    assert upper < 3 * axis
    coefficients = [1]
    for j in range(1, upper + 1):
        power, unit = 0, j
        while unit % 3 == 0:
            power += 1
            unit //= 3
        slope = axis // 3**power * pow(unit, -1, modulus) % modulus
        if slope:
            coefficients.append(0)
            for k in range(len(coefficients) - 1, 0, -1):
                coefficients[k] = (coefficients[k] + slope * coefficients[k - 1]) % modulus
            while len(coefficients) > 1 and not coefficients[-1]:
                coefficients.pop()
    return tuple(coefficients)


@lru_cache(None)
def insertion_accuracies(cutoff):
    """Required integer precision after removing m! and u^m."""
    accuracies = []
    m = 0
    while (remaining := cutoff - 4*m - 9*factorial_value(m)) > 0:
        accuracies.append(ceiling_ninths(remaining))
        m += 1
    return tuple(accuracies)


@lru_cache(None)
def axis_signature(height, exponent, cutoff):
    return (height > 0, tuple(
        insertion_polynomial(height, exponent, m, h)
        for m, h in enumerate(insertion_accuracies(cutoff))
    ))


@lru_cache(None)
def foreign_signature(alpha, beta, cutoff):
    """The residues of F_m/m! at the effective insertion accuracies."""
    coefficients = []
    for m, h in enumerate(insertion_accuracies(cutoff)):
        modulus = 3**h
        numerator = factorial_value(m*alpha) + factorial_value(m*beta)
        denominator = m*(factorial_value(alpha) + factorial_value(beta)) + 2*factorial_value(m)
        valuation = numerator - denominator
        assert valuation >= 0
        if valuation >= h:
            coefficients.append(0)
            continue
        top = factorial_unit(m*alpha, modulus) * factorial_unit(m*beta, modulus)
        bottom = factorial_unit(m, modulus)**2 * pow(
            factorial_unit(alpha, modulus) * factorial_unit(beta, modulus), m, modulus)
        coefficients.append(3**valuation * top * pow(bottom, -1, modulus) % modulus)
    return tuple(coefficients)


def signature(x, y, cutoff):
    axes = sorted((axis_signature(*x, cutoff), axis_signature(*y, cutoff)))
    return (*axes, foreign_signature(x[1], y[1], cutoff))


def representatives(cutoff, shape):
    gap_cap = ceiling_ninths(cutoff) + 2
    profiles = {}
    for x, y in product(*(digit_words(count, gap_cap) for count in shape)):
        profiles.setdefault(signature(x, y, cutoff), (x, y))
    return profiles


def polynomial_value(polynomial, n, modulus):
    value = 0
    for coefficient in reversed(polynomial):
        value = (value*n + coefficient) % modulus
    return value


class TwoDigitMoments(Moments):
    """Build both R and DR directly from their separated axis coefficients."""

    def __init__(self, start, profile, cutoff):
        self.start = start
        self.profile = profile
        ramification = 9 if start == 23 else 27
        precision = cutoff * ramification // 9
        self.ring = ring = Ring(precision, ramification)
        self.factorials = Factorials(ring.moduli[0])
        self.root_order = 4 * ramification // 9
        self.step = {17: 1, 20: 2, 23: 1}[start]
        u_order = {17: 27, 20: 18, 23: 10}[start]
        self.b = ring.add(ring.term(coefficient=I), ring.term(self.step, (0, 1, 0)))
        k = ring.add(ring.term(self.root_order, coefficient=(-1, -1)),
                     ring.term(self.root_order + self.step, (1, 0, 0)))
        u = ring.add(ring.scale(k, (-1, 0)), ring.term(u_order, C1))
        self.previous_offsets = 27 * ((precision - 1) // ramification)
        largest_index = self.previous_offsets + start + 3
        b_powers = self.powers(self.b, largest_index)
        accuracies = insertion_accuracies(cutoff)
        u_powers = self.powers(u, len(accuracies) - 1)
        self.native, self.d_native = self.build_native(accuracies, b_powers, u_powers,
                                                       largest_index)
        self.kernel_cache = {}
        self.unit_powers = {}

    def unit_power(self, offset):
        if offset not in self.unit_powers:
            self.unit_powers[offset] = super().unit_power.__wrapped__(self, offset)
        return self.unit_powers[offset]

    def build_native(self, accuracies, b_powers, u_powers, largest_index):
        ring = self.ring
        native = [{} for _ in range(largest_index + 1)]
        d_native = [{} for _ in range(largest_index + 1)]
        pure = [1]
        for x in range(largest_index):
            next_factor = (3*x + 1) * (3*x + 2) // 2
            pure.append(pure[-1] * next_factor % ring.moduli[0])

        x_signature, y_signature, foreign = self.profile
        for m, h in enumerate(accuracies):
            if not foreign[m]:
                continue
            modulus = 3**h
            axis_rows, axis_differences = [], []
            for nonlinear, polynomials in (x_signature, y_signature):
                row = [polynomial_value(polynomials[m], x, modulus)
                       * (pure[x] if nonlinear else 1) % modulus
                       for x in range(largest_index - m + 1)]
                axis_rows.append(row)
                axis_differences.append([row[0], *((a - b) % modulus
                                         for a, b in zip(row[1:], row[:-1]))])
            for destination, (x_row, y_row) in ((native, axis_rows),
                                                (d_native, axis_differences)):
                for n in range(m, largest_index + 1):
                    terms = []
                    for y in range(n - m + 1):
                        value = x_row[n - m - y] * y_row[y] * foreign[m] % modulus
                        if value:
                            terms.append(ring.scale(b_powers[y], (value, 0)))
                    native_term = ring.scale(ring.add(*terms), (factorial(m), 0))
                    destination[n] = ring.add(destination[n], ring.multiply(native_term, u_powers[m]))
        return native, d_native

    def scalar_order(self, high_block, offset):
        base = 27 * high_block
        if base + offset < 0:
            return None, 0
        valuation, unit = self.factorials.ratio([(base, 1)], [(base + offset, 1)])
        order = self.root_order*offset + self.ring.ramification*valuation
        assert order >= 0
        return order, unit

    def forms(self, high_block):
        ring = self.ring
        rows = []
        for j in range(4):
            index = self.start + j
            native = self.native if j < 2 else self.d_native
            terms = []
            for offset in range(-self.previous_offsets, index + 1):
                order, unit = self.scalar_order(high_block, offset)
                coefficient = native[index - offset]
                if order is None or order >= ring.precision or not coefficient:
                    continue
                key = (index, offset, order)
                if key not in self.kernel_cache:
                    kernel = ring.scale(self.unit_power(offset), ONE, order)
                    self.kernel_cache[key] = ring.multiply(kernel, coefficient)
                terms.append(ring.scale(self.kernel_cache[key], (unit, 0)))
            rows.append(ring.add(*terms))
        return rows


def lift_first_two(system, rows, high_block):
    start = system.start
    source = system.ring
    if start == 23 and high_block % 3:
        h = 1 if high_block % 3 == 1 else -1
        check_centers(source, rows[:2], ((h, -h), (h, 0)))
        return source, rows
    centers = {17: ((-1, 1), (-1, 0)), 20: (I, ONE), 23: (ZERO, ZERO)}[start]
    bound = {17: 1, 20: 2, 23: 1}[start]
    leading_order = {17: 27, 20: 54, 23: 27}[start]
    check_centers(source, rows[:2], centers)
    expanded = Ring(3*source.precision, 3*source.ramification)
    shifted = [substitute(source, expanded, row, centers, 0) for row in rows[:2]]
    check_raise_order(expanded, shifted, leading_order, bound)
    return expanded, [substitute(source, expanded, row, centers, bound) for row in rows]


def filter_profile(job):
    start, cutoff, profile, high_blocks = job
    system = TwoDigitMoments(start, profile, cutoff)
    retained = []
    for high_block in high_blocks:
        ring, forms = lift_first_two(system, system.forms(high_block), high_block)
        _, value = form_reduction(ring, forms, evaluate=True)
        if value == ZERO:
            retained.append((profile, high_block))
    return retained


def verify_start(pool, start, shape):
    retained = None
    for cutoff in CUTOFFS[start]:
        period = 9 if cutoff <= 36 else 27
        jobs = []
        for profile, (x, y) in sorted(representatives(cutoff, shape).items()):
            if retained is None:
                high_blocks = list(range(period))
            else:
                old_profile = signature(x, y, previous_cutoff)
                high_blocks = [h for h in range(period) if (old_profile, h % previous_period) in retained]
            if high_blocks:
                jobs.append((start, cutoff, profile, high_blocks))
        retained = {entry for result in pool.map(filter_profile, jobs) for entry in result}
        count = sum(len(job[3]) for job in jobs)
        print(f"start {start}, shape {shape}, cutoff {cutoff}/9: "
              f"{count - len(retained)} excluded, {len(retained)} retained", flush=True)
        previous_cutoff, previous_period = cutoff, period
    assert not retained, (start, shape, len(retained))


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=4) as pool:
        for start in (17, 20, 23):
            for shape in SHAPES:
                verify_start(pool, start, shape)
