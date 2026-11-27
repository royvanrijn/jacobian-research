"""Verify the second-radius filters for total foreign digit count at least four."""

from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import comb

from verify_binary_xy_second_radius import (
    A9, B9, C1, CONSTANT, I, ONE, ZERO, Moments, Ring, check_centers,
    check_raise_order, gaussian_power, gaussian_product, normalize_ninth_rows, reducible_term,
    substitute,
)


VARIABLES = (A9, B9, C1)
AXIS_ONE = (0, 0)
AXIS_THREE = (1, 1)
CENTERS = {17: ((-1, 1), (-1, 0)), 20: (I, ONE), 23: (ZERO, ZERO)}
BOUND = {17: 1, 20: 2, 23: 1}
LEADING_ORDER = {17: 27, 20: 54, 23: 27}


def prefixes(width):
    return [(width, sum(bit * 3**j for j, bit in enumerate(bits)))
            for bits in product((0, 1), repeat=width)]


def axis_pairs(precision):
    types = [AXIS_ONE, *prefixes(precision - 1)]
    if precision == 4:
        types.insert(1, AXIS_THREE)
    for pair in combinations_with_replacement(sorted(types), 2):
        if any(width == precision - 1 for width, _ in pair):
            yield pair


def coarse_axis(axis):
    width, prefix = axis
    return (2, prefix * 3) if width == 1 else (2, prefix // 3) if width == 3 else axis


class NativeProfile(Moments):
    """The two additive components of each row in the double-insertion factor."""

    def __init__(self, start, axes, accuracy, precision=None):
        self.axis_coefficients = {}
        if precision is None:
            precision = accuracy * (9 if start == 23 else 27)
        super().__init__(start, axes, precision)
        self.components = tuple(zip(*self.native))

    def pure_axis(self, width, count):
        return 1 + 9 * comb(count + 1, 3) if width else 1

    def axis_coefficient(self, axis, count, insertions):
        key = axis, count, insertions
        if key not in self.axis_coefficients:
            width, prefix = axis
            pure = self.pure_axis(width, count)
            offset = insertions * prefix
            value = pure * comb(3**width * count + offset, offset)
            self.axis_coefficients[key] = value % self.ring.moduli[0]
        return self.axis_coefficients[key]

    def native_row(self, index, b_powers, u_powers):
        ring = self.ring
        contributions = [{}, {}, {}]
        for m in range(min(index, 2) + 1):
            terms = []
            for y in range(index - m + 1):
                x = index - m - y
                value = (self.axis_coefficient(self.axes[0], x, m)
                         * self.axis_coefficient(self.axes[1], y, m))
                terms.append(ring.scale(b_powers[y], (value, 0)))
            contributions[m] = ring.multiply(ring.add(*terms), u_powers[m])
        return ring.add(*contributions[:2]), contributions[2]

    def moment_components(self, high_block):
        base, quadratic = self.components
        self.native = base
        first = self.forms(high_block)
        self.native = quadratic
        return first, self.forms(high_block)


def add_component(ring, left, right, coefficient):
    return tuple(ring.add(row, ring.scale(other, (coefficient, 0)))
                 for row, other in zip(left, right))


def at_forced_residues(digit, pivots, ring):
    residues = []
    for _, order, row, inverse in pivots:
        constant = ring.digit(row, order).get(CONSTANT, ZERO)
        residues.append(gaussian_product(tuple(-part for part in constant), inverse, 3))
    result = ZERO
    for exponents, coefficient in digit.items():
        for exponent, residue in zip(exponents, residues):
            coefficient = gaussian_product(coefficient, gaussian_power(residue, exponent, 3), 3)
        result = tuple((a + b) % 3 for a, b in zip(result, coefficient))
    return result


def form_reduction(ring, rows, evaluate=False, full=False):
    x, y, z, remainder = rows
    x, y = normalize_ninth_rows(ring, x, y)
    pivots = []
    for exponent, row in zip(VARIABLES, (x, y, z)):
        order, digit = ring.initial(row)
        assert exponent in digit and set(digit) <= {CONSTANT, exponent}
        pivots.append((exponent, order, row, gaussian_power(digit[exponent], -1, 3)))

    for level in range(ring.precision):
        while True:
            digit = ring.digit(remainder, level)
            multiplier, row = reducible_term(ring, digit, level, pivots)
            if not multiplier:
                break
            if row is x or row is y:
                assert ring.initial(multiplier)[0] >= ring.ramification
            remainder = ring.add(remainder, ring.multiply(multiplier, row))
        if digit and not full:
            if evaluate:
                return level, at_forced_residues(digit, pivots, ring)
            assert set(digit) == {CONSTANT}
            return level, digit[CONSTANT]
    return remainder if full else (ring.precision, ZERO)


def profile_quotients(system, high_block, coefficients, evaluate=False):
    source = system.ring
    start = system.start
    base, quadratic = system.moment_components(high_block)
    if start == 23 and high_block % 3:
        for coefficient in coefficients:
            rows = add_component(source, base, quadratic, coefficient)
            level, value = form_reduction(source, rows, evaluate)
            yield coefficient, level, source.ramification, value
        return

    expanded = Ring(3 * source.precision, 3 * source.ramification)
    centers, bound = CENTERS[start], BOUND[start]
    initial_base = [substitute(source, expanded, row, centers, 0) for row in base[:2]]
    initial_quad = [substitute(source, expanded, row, centers, 0) for row in quadratic[:2]]
    final_base = [substitute(source, expanded, row, centers, bound) for row in base]
    final_quad = [substitute(source, expanded, row, centers, bound) for row in quadratic]
    for coefficient in coefficients:
        original_rows = add_component(source, base[:2], quadratic[:2], coefficient)
        check_centers(source, original_rows, centers)
        initial_rows = add_component(expanded, initial_base, initial_quad, coefficient)
        check_raise_order(expanded, initial_rows, LEADING_ORDER[start], bound)
        rows = add_component(expanded, final_base, final_quad, coefficient)
        level, value = form_reduction(expanded, rows, evaluate)
        yield coefficient, level, expanded.ramification, value


def passes_precision_three(start, axes, coefficient):
    x, y = (prefix for _, prefix in axes)
    leading = x // 3 + y // 3
    if start == 17:
        return leading == 0 and coefficient % 27 == (23 - 9 * (x + y)) % 27
    if start == 20:
        return leading == 0 and coefficient % 9 == 2
    return coefficient % 9 == (8 - 3 * leading) % 9


def passes_precision_four(start, axes, coefficient):
    top_two = [coarse_axis(axis)[1] for axis in axes]
    if start == 17:
        return False
    if start == 20:
        return not any(top_two) and coefficient % 27 == 11
    return all(prefix < 3 for prefix in top_two) and coefficient % 27 == (17 - 9 * sum(top_two)) % 27


def verify_pair(job):
    start, accuracy, axes = job
    coefficients = list(range(2, 3**accuracy, 3))
    if accuracy == 4:
        reduced_axes = tuple(map(coarse_axis, axes))
        coefficients = [f for f in coefficients
                        if passes_precision_three(start, reduced_axes, f)]
    if not coefficients:
        return Counter()
    system = NativeProfile(start, axes, accuracy)
    passes = passes_precision_three if accuracy == 3 else passes_precision_four
    counts = Counter()
    for high_block in range(9):
        for factor, order, ramification, residue in profile_quotients(
                system, high_block, coefficients):
            expected = passes(start, axes, factor)
            assert (residue == ZERO) == expected, (start, accuracy, axes, high_block,
                                                   factor, order, ramification, residue)
            counts['zero' if expected else 'excluded'] += 1
    return counts


class HigherPrecisionProfile(NativeProfile):
    def pure_axis(self, width, count):
        if not width:
            return 1
        value = 1
        for j in range(count):
            value = value * ((3*j + 1) * (3*j + 2) // 2) % 243
        return value


def verify_above_four(job):
    start, precision, axes, coefficients = job
    system = HigherPrecisionProfile(start, axes, 5, precision=precision)
    count = 0
    for high_block in range(27):
        for factor, order, ramification, residue in profile_quotients(
                system, high_block, coefficients):
            assert residue != ZERO, (start, precision, axes, high_block, factor,
                                     order, ramification, residue)
            count += 1
    return count


def stages_above_four():
    profiles = (AXIS_ONE, (3, 0), (3, 1), (3, 3), (3, 4))
    middle = []
    late_with_next_digit = []
    for x, y in combinations_with_replacement(profiles, 2):
        if y == AXIS_ONE:
            continue
        if y[1] < 3:
            middle.append((20, 114, (x, y), (11, 38, 65)))
        else:
            residue = (17 - 9 * (x[1] // 3 + y[1] // 3)) % 27
            late_with_next_digit.append((23, 37, (x, y), range(residue, 81, 27)))
    late_with_linear_axis = [(23, 41, (AXIS_ONE, (4, prefix)), (17, 44, 71))
                             for prefix in (0, 1, 3, 4)]
    return (
        ("middle start", middle),
        ("late start with a nonzero second digit", late_with_next_digit),
        ("late start with an axis of degree one", late_with_linear_axis),
    )


def verify_late_nonlinear_axes(axes):
    system = HigherPrecisionProfile(23, axes, 5, precision=42)
    next_digits = sum(prefix // 3 for _, prefix in axes)
    required = (44 - 27 * next_digits) % 81
    counts = Counter()
    for high_block in range(27):
        for factor, order, ramification, residue in profile_quotients(
                system, high_block, (17, 44, 71), evaluate=True):
            retained = high_block % 3 != 0 and factor == required
            assert (residue == ZERO) == retained, (
                axes, high_block, factor, order, ramification, residue)
            counts['retained' if retained else 'excluded'] += 1
    return counts


def translate_late_rows(ring, polynomial):
    """A' = -1-i+A_1, B' = t^2 B_1 after the first ramified substitution."""
    @lru_cache(None)
    def a_power(exponent):
        return ring.add(*(
            ring.term(exponents=(j, 0, 0),
                      coefficient=gaussian_product(
                          (comb(exponent, j), 0),
                          gaussian_power((-1, -1), exponent - j, ring.moduli[0]),
                          ring.moduli[0],
                      ))
            for j in range(exponent + 1)
        ))

    return ring.add(*(
        ring.scale(ring.multiply(a_power(a), ring.term(2*b, (0, b, c))),
                   coefficient, order)
        for (order, a, b, c), coefficient in polynomial.items()
    ))


def verify_late_residue_conflict(axes):
    system = HigherPrecisionProfile(23, axes, 5, precision=43)
    ring = system.ring
    expanded = Ring(3 * ring.precision, 3 * ring.ramification)
    factor = (44 - 27 * sum(prefix // 3 for _, prefix in axes)) % 81
    checks = 0
    for high_block in range(27):
        h = {1: 1, 2: -1}.get(high_block % 3)
        if h is None:
            continue
        rows = add_component(ring, *system.moment_components(high_block), factor)
        remainder = form_reduction(ring, rows, full=True)
        assert ring.initial(remainder) == (41, {CONSTANT: (0, -h % 3), (0, 1, 0): I})

        centers = ((h, -h), (h, 0))
        check_centers(ring, rows[:2], centers)
        shifted = tuple(substitute(ring, expanded, row, centers, 0) for row in rows[:2])
        check_raise_order(expanded, shifted, 27, 1)
        refined = tuple(substitute(ring, expanded, row, centers, 1) for row in rows[:2])
        x, y = normalize_ninth_rows(expanded, *refined)
        level, lead_x = expanded.initial(x)
        assert level == 36
        level, lead_y = expanded.initial(y)
        assert level == 36
        assert lead_y[CONSTANT] != ZERO if h == -1 else lead_y.get(CONSTANT, ZERO) == ZERO
        if h == 1:
            assert gaussian_product(lead_x[A9], gaussian_power((2, 2), 9, 3), 3) == tuple(
                -a % 3 for a in lead_x[CONSTANT])
            assert ring.digit(remainder, 42) == {(1, 1, 0): (1, 1), (1, 0, 0): (1, 1)}

            x, y = (translate_late_rows(expanded, row) for row in (x, y))
            assert expanded.initial(x) == (36, {A9: ONE})
            assert not expanded.digit(x, 37) and not expanded.digit(x, 38)
            order_39 = expanded.digit(x, 39)
            assert order_39.pop(CONSTANT) == (2, 2)
            assert all(a > 0 for a, _, _ in order_39)
            linked = expanded.add(y, expanded.scale(x, ONE, 3))
            assert expanded.initial(linked) == (
                42, {A9: (2, 1), (9, 1, 0): (2, 0)})
            assert not expanded.digit(linked, 43) and not expanded.digit(linked, 44)
            independent = {key: value for key, value in expanded.digit(linked, 45).items()
                           if key[0] == 0}
            assert independent == {(0, 1, 0): (1, 1), (0, 3, 0): (2, 1)}
        checks += 1
    return checks


def verify():
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=4) as pool:
        for accuracy in (3, 4):
            for start in (17, 20, 23):
                jobs = [(start, accuracy, axes) for axes in axis_pairs(accuracy)]
                counts = sum(pool.map(verify_pair, jobs), Counter())
                print(f"precision {accuracy}, start {start}: "
                      f"{counts['excluded']} excluded, {counts['zero']} zero remainders",
                      flush=True)
        for label, jobs in stages_above_four():
            count = sum(pool.map(verify_above_four, jobs))
            print(f"{label}: {count} excluded above precision four", flush=True)
        axes = list(combinations_with_replacement([(4, prefix) for prefix in (0, 1, 3, 4)], 2))
        counts = sum(pool.map(verify_late_nonlinear_axes, axes), Counter())
        print(f"late start with two nonlinear axes: {counts['excluded']} excluded, "
              f"{counts['retained']} retained by the residue test", flush=True)
        count = sum(pool.map(verify_late_residue_conflict, axes))
        print(f"late start with two nonlinear axes: {count} final conflicts", flush=True)


if __name__ == "__main__":
    verify()
