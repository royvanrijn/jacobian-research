#!/usr/bin/env python3
"""Supplementary native controls for binary divided-power quadratic windows.

The universal proof is written separately. This program supplies bounded checks
of the native identities, degree-bounded Hensel lifts, and even-constant reduction.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import random

if not __debug__:
    raise SystemExit('Refusing optimized Python: run without -O.')


def require(test, message):
    if not test:
        raise ValueError(message)


def pd_multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += comb(i + j, i) * x * y
    return out


def low_solution(rhs, unit, length):
    require(unit[0] % 2 == 1, 'Low divided-power coefficient is not a unit.')
    out = []
    for k in range(length):
        value = rhs[k]
        for i in range(k):
            j = k - i
            if j < len(unit):
                value -= comb(k, i) * out[i] * unit[j]
        out.append(value % 2)
    return out


def verify_product(f, g, u, precision):
    actual = pd_multiply(g, u)
    modulus = 1 << precision
    require(len(actual) == len(f), 'Product degree changed.')
    require(all((x - y) % modulus == 0 for x, y in zip(actual, f)),
            'Divided-power factor product has insufficient precision.')


def highest_bit_split(f, precision, counts):
    degree = len(f) - 1
    require(degree > 0 and degree % 2 == 0, 'Positive even degree required.')
    require(all(x % 2 == (1 if i % 2 == 0 else 0) for i, x in enumerate(f)),
            'Even-unit/odd-even divided-power parity class failed.')
    e = 1 << (degree.bit_length() - 1)
    tail = degree - e
    u = [f[e + j] % 2 for j in range(tail + 1)]
    g = low_solution([x % 2 for x in f[:e]], u, e) + [1]
    for stage in range(1, precision):
        actual = pd_multiply(g, u)
        require(all((f[i] - actual[i]) % (1 << stage) == 0
                    for i in range(degree + 1)), 'Previous lift precision failed.')
        error = [((f[i] - actual[i]) // (1 << stage)) % 2
                 for i in range(degree + 1)]
        delta_u = error[e:]
        cross = pd_multiply(g[:e], delta_u)
        rhs = [error[i] - (cross[i] if i < len(cross) else 0)
               for i in range(e)]
        delta_g = low_solution(rhs, u, e)
        for i in range(e):
            g[i] += (1 << stage) * delta_g[i]
        for j in range(tail + 1):
            u[j] += (1 << stage) * delta_u[j]
        counts['precision_lifts'] += 1
    require(g[e] == 1 and g[0] % 2 == 1, 'G normalization failed.')
    require(all(g[i] % 2 == 0 for i in range(1, e, 2)), 'G odd refinement failed.')
    require(all(x % 2 == (1 if i % 2 == 0 else 0) for i, x in enumerate(u)),
            'Recursive tail parity class failed.')
    verify_product(f, g, u, precision)
    counts['degree_bounded_splits'] += 1
    return g, u


def pd_factor(f, precision, counts):
    if len(f) == 1:
        require(f[0] % 2 == 1, 'Final constant is not a unit.')
        return f[0], []
    g, u = highest_bit_split(f, precision, counts)
    unit, rest = pd_factor(u, precision, counts)
    return unit, [g] + rest


def checked_factor(f, precision, counts):
    unit, factors = pd_factor(f, precision, counts)
    actual = [unit]
    for g in factors:
        actual = pd_multiply(actual, g)
    modulus = 1 << precision
    require(len(actual) == len(f) and
            all((x - y) % modulus == 0 for x, y in zip(actual, f)),
            'Complete binary factor product failed.')
    degree = len(f) - 1
    bits = [1 << r for r in range(degree.bit_length()) if degree >> r & 1]
    require(sorted(len(g) - 1 for g in factors) == bits, 'Wrong binary degrees.')
    return factors


def moment(m, constant):
    """Select quadratic, linear, and constant factors, then apply factorials."""
    fs = [factorial(i) for i in range(2 * m + 1)]
    return [sum(comb(m, i) * comb(m - i, b) * constant ** (m - i - b) * fs[2 * i + b]
                for b in range(m - i + 1)) for i in range(m + 1)]


def native_pd(m, constant, precision):
    coefficients = moment(m, constant)
    scale = 1 << m.bit_count()
    ordinary = []
    for i, value in enumerate(coefficients):
        value *= scale
        require(value % (1 << (2 * i)) == 0, 'Scaled moment not integral.')
        ordinary.append(value // (1 << (2 * i)))
    modulus = 1 << precision
    inverse_lead = pow(ordinary[m], -1, modulus)
    pd = []
    for k in range(m + 1):
        numerator = ordinary[m - k] * factorial(k)
        require(numerator % (1 << k) == 0, 'Native divided-power coefficient not integral.')
        pd.append((numerator // (1 << k)) * inverse_lead % modulus)
    return pd


def polynomial_add(*terms):
    length = max(map(len, terms))
    out = [0] * length
    for term in terms:
        for i, value in enumerate(term):
            out[i] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def scaled_shift(poly, scalar, shift=0):
    return [0] * shift + [scalar * x for x in poly]


def derivative_identity(m, constant, forcing_sign=1):
    f = moment(m, constant)
    prior = moment(m - 1, constant)
    derivative = [(i + 1) * f[i + 1] for i in range(m)]
    left = polynomial_add(scaled_shift(derivative, 2, 2))
    right = polynomial_add(scaled_shift(f, 2 * m, 1), scaled_shift(f, -1),
                           scaled_shift(prior, m), scaled_shift(prior, -2 * m * constant, 1),
                           [forcing_sign * constant ** m])
    require(left == right, 'Native derivative identity failed.')


def recurrence(j, constant, forcing_sign=1):
    f, prior, earlier = (moment(j, constant), moment(j - 1, constant),
                         moment(j - 2, constant))
    right = polynomial_add(scaled_shift(prior, 2 * j * (2 * j - 1), 1),
                           scaled_shift(earlier, j * (j - 1)),
                           scaled_shift(earlier, -4 * constant * j * (j - 1), 1),
                           [forcing_sign * (constant ** j + j * constant ** (j - 1))])
    require(polynomial_add(f) == right, 'Native PLUS recurrence failed.')


def f2_trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def f2_remainder(a, b):
    a = a.copy()
    while a != [0] and len(a) >= len(b):
        shift = len(a) - len(b)
        for i, value in enumerate(b):
            a[i + shift] ^= value
        f2_trim(a)
    return a


def even_constant_controls():
    count = 0
    for constant in (-4, -2, 0, 2, 4):
        rows = []
        for j in range(21):
            row = []
            for i, value in enumerate(moment(j, constant)):
                rational = Fraction(value, (1 << i) * factorial(j))
                require(rational.denominator % 2 == 1, 'Even-C coefficient not 2-integral.')
                row.append(rational.numerator % 2)
            require(row[-1] == 1, 'Reduction lost its degree.')
            if j == 0:
                expected = [1]
            elif j == 1:
                expected = [1, 1]
            else:
                expected = [0] + rows[-1]
                for i, value in enumerate(rows[-2]):
                    expected[i] ^= value
            require(row == expected, 'Even-C native reduction recurrence failed.')
            if rows:
                a, b = row, rows[-1]
                while b != [0]:
                    a, b = b, f2_remainder(a, b)
                require(a == [1], 'Consecutive reductions not coprime.')
            rows.append(row)
            count += 1
    return count


def rejected(fn):
    try:
        fn()
    except ValueError:
        return 1
    raise ValueError('Deliberately corrupted control was accepted.')


def run():
    precision = 16
    counts = {'degree_bounded_splits': 0, 'precision_lifts': 0}
    groups = []
    for m in (2, 4, 6, 8, 10, 14, 16, 18, 30, 46, 62, 90):
        for constant in (1, 3, -1, -m - 1):
            derivative_identity(m, constant)
            recurrence(m, constant)
            factors = checked_factor(native_pd(m, constant, precision), precision, counts)
            targets = []
            for g in factors:
                e = len(g) - 1
                tail = m % e
                j = e - tail - 1
                k = e - j
                require(k % 2 == 1 and g[k] % 2 == 0, 'Target odd refinement failed.')
                targets.append({'degree': e, 'target_index': j,
                                'forbidden_valuation': tail.bit_count() + 1,
                                'inherited_lower_bound': k.bit_count() + 1})
            groups.append({'m': m, 'C': constant, 'factors': targets})
    rng = random.Random(20261126)
    for m in (6, 10, 14, 30, 46, 62):
        for unused in range(20):
            f = [2 * rng.randrange(1 << (precision - 1)) + (i % 2 == 0)
                 for i in range(m + 1)]
            checked_factor(f, precision, counts)
    negative = 0
    negative += rejected(lambda: highest_bit_split([1, 1, 1], precision, counts))
    negative += rejected(lambda: highest_bit_split([1, 0, 0], precision, counts))
    negative += rejected(lambda: highest_bit_split([1, 2], precision, counts))
    negative += rejected(lambda: derivative_identity(2, 1, forcing_sign=-1))
    negative += rejected(lambda: recurrence(2, 1, forcing_sign=-1))
    negative += rejected(lambda: verify_product([1, 0, 1], [0, 0, 1], [1], 2))
    return {'status': 'PASS', 'precision_bits': precision, 'native_groups': groups,
            'native_derivative_identities': len(groups), 'native_PLUS_recurrences': len(groups),
            'random_parity_class_groups': 120, **counts,
            'even_constant_native_reductions': even_constant_controls(),
            'negative_controls': negative,
            'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Supplementary finite native and lifting checks only. The all-range proof is written separately. No formal verification, engine evaluation, or conjecture closure is inferred from a prefix.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    target = Path(args.output)
    if target.exists():
        raise SystemExit('Refusing to overwrite output: ' + str(target))
    result = run()
    with target.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({'status': result['status'], 'output': str(target),
                      'native_groups': len(result['native_groups']),
                      'degree_bounded_splits': result['degree_bounded_splits'],
                      'precision_lifts': result['precision_lifts'],
                      'negative_controls': result['negative_controls']}))


if __name__ == '__main__':
    main()