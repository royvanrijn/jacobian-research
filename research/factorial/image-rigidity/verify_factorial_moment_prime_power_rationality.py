"""Exact checks of the finite transfer, digit residue, and all-order packing."""

from fractions import Fraction
from itertools import product
from math import comb, factorial, prod


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def native_coefficient(columns, counts):
    dimension = len(columns[0])
    exponents = [sum(a[i] * m for a, m in zip(columns, counts)) for i in range(dimension)]
    denominator = prod(factorial(m) * prod(factorial(n) for n in a)**m
                       for a, m in zip(columns, counts))
    quotient, remainder = divmod(prod(factorial(n) for n in exponents), denominator)
    assert remainder == 0
    return quotient


def native_series(columns, weights, order, modulus):
    return [
        sum(native_coefficient(columns, counts) * prod(u**m for u, m in zip(weights, counts))
            for counts in compositions(n, len(columns))) % modulus
        for n in range(order + 1)
    ]


def finite_transfer(columns, weights, prime, precision, order):
    modulus = prime**precision
    maximum = max(map(max, columns))
    period = modulus
    while maximum >= prime:
        maximum //= prime
        period *= prime
    zero = (0,) * len(columns[0])
    states = {(zero, 0): 1}
    result = [1]
    for _ in range(order):
        following = {}
        for (residues, mask), value in states.items():
            for column, weight in zip(columns, weights):
                target = tuple((old + digit) % period for old, digit in zip(residues, column))
                active = mask | sum(1 << i for i, exponent in enumerate(column) if exponent)
                first = (active & -active).bit_length() - 1
                if column[first] == 0:
                    continue
                multiplicity = comb((target[first] - 1) % period, column[first] - 1)
                for i, exponent in enumerate(column):
                    if i != first:
                        multiplicity *= comb(target[i], exponent)
                key = target, active
                updated = (following.get(key, 0) + weight * multiplicity * value) % modulus
                if updated:
                    following[key] = updated
                elif key in following:
                    del following[key]
        states = following
        result.append(sum(states.values()) % modulus)
    return result


def digits(number, prime):
    result = []
    while number:
        number, digit = divmod(number, prime)
        result.append(digit)
    return result


def digit_formula(columns, weights, prime, order):
    one = [1] + [0] * order

    def multiply(left, right):
        result = [0] * (order + 1)
        for i, a in enumerate(left):
            if a:
                for j in range(order + 1 - i):
                    result[i + j] += a * right[j]
        return [n % prime for n in result]

    digit_columns = [{(i, t): digit for i, entry in enumerate(column)
                      for t, digit in enumerate(digits(entry, prime)) if digit}
                     for column in columns]
    singletons, others = [], []
    for slots, unit in zip(digit_columns, weights):
        if sum(slots.values()) == 1:
            i, t = next(iter(slots))
            singletons.append((i, t, unit))
        else:
            others.append((slots, unit))

    def vertex(i, t):
        polynomial = one.copy()
        for row, below, unit in singletons:
            if row == i and below <= t:
                degree = prime**(t - below)
                if degree <= order:
                    polynomial[degree] -= pow(unit, degree, prime)
        return [n % prime for n in polynomial]

    heights = [max((t for slots in digit_columns for row, t in slots if row == i), default=-1)
               for i in range(len(columns[0]))]
    vertices = {(i, t): vertex(i, t) for i, height in enumerate(heights)
                for t in range(height + 1)}
    numerator = [0] * (order + 1)
    for counts in product(range(prime), repeat=len(others)):
        if sum(counts) > order:
            continue
        loads = dict.fromkeys(vertices, 0)
        scalar = 1
        for (slots, unit), count in zip(others, counts):
            scalar *= pow(unit, count, prime) * pow(factorial(count), -1, prime)
            for slot, digit in slots.items():
                loads[slot] += count * digit
                scalar *= pow(factorial(digit), -count, prime)
        if max(loads.values(), default=0) >= prime:
            continue
        term = [0] * (order + 1)
        term[sum(counts)] = scalar % prime
        for slot, load in loads.items():
            term = [factorial(load) * coeff % prime for coeff in term]
            for _ in range(prime - 1 - load):
                term = multiply(term, vertices[slot])
        numerator = [(a + b) % prime for a, b in zip(numerator, term)]
    denominator = one.copy()
    for i, height in enumerate(heights):
        denominator = multiply(denominator, vertex(i, height + 1))
    result = []
    for n, coefficient in enumerate(numerator):
        contribution = sum(denominator[j] * result[n - j] for j in range(1, n + 1))
        result.append((coefficient - contribution) % prime)
    return result


def verify_cycle_poles():
    cases = ((2, 16, (0, 1, 2)), (3, 27, (0, 2)), (5, 25, (0, 1)))
    for prime, period, positions in cases:
        paths = [{0: 1}]
        for target in range(1, period + 1):
            polynomial = {}
            for position in positions:
                step = prime**position
                if target < step:
                    continue
                weight = comb(target - 1, step - 1) % prime
                for degree, coefficient in paths[target - step].items():
                    subtotal = polynomial.get(degree + 1, 0) + coefficient * weight
                    polynomial[degree + 1] = subtotal % prime
            paths.append({degree: value for degree, value in polynomial.items() if value})
        assert paths[-1] == {period // prime**position: 1 for position in positions}
    print("Three aligned singleton cycle denominators verified.")


def valuation(integer, prime):
    value = 0
    while integer % prime == 0:
        integer //= prime
        value += 1
    return value


def verify():
    samples = (
        (2, 3, ((2, 0), (1, 1), (0, 1), (3, 2)), (1, 3, 1, 2)),
        (3, 2, ((3, 0), (2, 1), (0, 1), (1, 0)), (2, 1, 3, 1)),
        (5, 2, ((1,), (5,), (2,), (7,)), (3, 2, 4, 1)),
        (3, 2, ((0, 9, 0), (2, 1, 6), (1, 0, 2)), (1, 2, 4)),
    )
    for prime, precision, columns, weights in samples:
        expected = native_series(columns, weights, 16, prime**precision)
        assert finite_transfer(columns, weights, prime, precision, 16) == expected
        assert digit_formula(columns, weights, prime, 16) == [x % prime for x in expected]

    packing_cases = (
        (2, 3, 5, ((4, 0, 0), (0, 1, 1), (4, 1, 0)), (1, 1, 1)),
        (2, 3, 6, ((4, 0), (0, 1), (4, 1)), (1, 2, 1)),
        (3, 2, 4, ((3, 0, 0), (0, 1, 1), (3, 2, 0)), (1, 1, 2)),
    )
    for prime, precision, gap, columns, weights in packing_cases:
        merged = tuple((prime**gap * a + b, *tail) for a, b, *tail in columns)
        modulus = prime**precision
        assert native_series(columns, weights, 12, modulus) == native_series(merged, weights, 12, modulus)

    carry_cases = (
        (2, ((1, 0), (3, 1)), Fraction(3, 4)),
        (2, ((2, 0), (3, 1), (2, 2), (2, 3)), Fraction(1, 4)),
        (3, ((3, 0), (3, 1), (4, 2)), Fraction(1, 4)),
    )
    checked = 0
    for prime, columns, gap in carry_cases:
        for total in range(prime, 16, prime):
            for counts in compositions(total, len(columns)):
                if counts[0] != total:
                    native_value = valuation(native_coefficient(columns, counts), prime)
                    assert native_value > gap * (total - counts[0])
                    checked += 1

    print("68 prime-power transfer coefficients and 68 explicit digit coefficients verified.")
    print(f"39 packed moment coefficients and {checked} shared-slot valuation inequalities verified.")


if __name__ == "__main__":
    verify()
    verify_cycle_poles()
