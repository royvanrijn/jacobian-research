#!/usr/bin/env python3
"""Certify the weighted tetrahedral standard gap using the T4 coefficients."""

import argparse
from functools import cache
import json
from math import comb, factorial, lcm
import platform
import sys

import flint
from flint import arb, arb_mat, ctx, fmpq

from check_tetrahedral_standard_gap import hermite, require, shifted_hermites


START = fmpq(17, 25)
GAP = fmpq(1, 1000)
# Right endpoints in units of 1/1000 and corresponding truncation degrees.
SLABS = ((770, 25), (820, 28), (850, 34), (875, 39), (890, 44),
         (900, 48), (910, 54), (915, 56), (920, 60))


def coefficient_numerators(basis, maximum):
    """Integer numerators of C_m with common denominator 2**(m+1) m!."""
    index = {item: number for number, item in enumerate(basis)}
    by_degree = []
    for m in range(1, maximum + 1):
        entries = {}
        for normal in range(0, m + 1, 2):
            for j in range(0, m - normal + 1, 2):
                i = m - normal - j
                normal_at_zero = hermite(normal, 2)[0]
                scale = (2**i * comb(m, i) * comb(m - i, j) * normal_at_zero)
                left = {(u, t): p * q * 2**t
                        for u, p in enumerate(hermite(i, 1))
                        for t, q in enumerate(hermite(j, 2))
                        if p and q and (u, t) != (0, 0)}
                right = {key: normal_at_zero * value for key, value in left.items()}
                for (u, t), value in shifted_hermites(j, normal).items():
                    for d, p in enumerate(hermite(i, 1)):
                        key = (u, t + d)
                        if p and value and not ((t + d) % 2) and key != (0, 0):
                            right[key] = right.get(key, 0) + p * value
                for x, a in left.items():
                    for y, b in right.items():
                        key = (index[x], index[y])
                        entries[key] = entries.get(key, 0) + scale * a * b
        entries = {key: value for key, value in entries.items() if value}
        require(all(entries.get((j, i)) == value for (i, j), value in entries.items()),
                f"Asymmetric Hermite coefficient {m}")
        by_degree.append(entries)
    return by_degree


def planar_numerators(basis, maximum):
    """Planar Hermite projection numerators with denominator 2**m m!."""
    index = {item: number for number, item in enumerate(basis)}
    by_degree = []
    for m in range(1, maximum + 1):
        entries = {}
        for j in range(0, m + 1, 2):
            i = m - j
            scale = 2**i * comb(m, j)
            feature = {index[u, t]: p * q * 2**t
                       for u, p in enumerate(hermite(i, 1))
                       for t, q in enumerate(hermite(j, 2))
                       if p and q and (u, t) != (0, 0)}
            for x, a in feature.items():
                for y, b in feature.items():
                    entries[x, y] = entries.get((x, y), 0) + scale * a * b
        by_degree.append(entries)
    return by_degree


def rational_upper(x):
    denominator = 10**12
    numerator = (denominator * x).upper().ceil().unique_fmpz()
    require(numerator is not None, "Unable to round the coefficient upward")
    return fmpq(numerator, denominator)


def scaled_coefficients(left, right, degree):
    sigma_left = (1 - arb(left**2)).sqrt()
    sigma_right = rational_upper((1 - arb(right**2)).sqrt())
    weights = []
    for m in range(1, degree + 1):
        critical = fmpq(m - 1, m)
        right_bound = sigma_right * right**(m - 1)
        if right**2 <= critical:
            maximum = right_bound
        elif left**2 >= critical:
            maximum = max(right_bound, rational_upper(sigma_left * arb(left**(m - 1))))
        else:
            maximum = max(right_bound, rational_upper(arb(critical).sqrt()**(m - 1) / arb(m).sqrt()))
        weights.append(maximum / (2**(m + 1) * factorial(m)))
    return weights


def planar_tail_coefficients(right, degree):
    """Negative weights for the certified deficit from the scalar tail cap."""
    sigma = (1 - arb(right**2)).sqrt()
    limit = right**degree
    partial = [fmpq(0)]
    for j in range(degree//2 + 1):
        partial.append(partial[-1] + fmpq(comb(2*j, j), 4**j) * right**(2*j))
    weights = []
    for k in range(1, degree + 1):
        j = (degree - k)//2 + 1
        value = arb(right**(k - 1)) * (1 - sigma * arb(partial[j]))
        deficit = max(fmpq(0), limit - rational_upper(value))
        weights.append(-deficit / (2**k * factorial(k)))
    return weights


def finite_matrix(coefficient_data, size, weights):
    # Use a common exact denominator before converting the result to balls.
    denominator = lcm(*(int(w.denom()) for w in weights))
    entries = [0] * (size * size)
    for table, weight in zip(coefficient_data, weights):
        multiplier = int(weight.numer()) * (denominator // int(weight.denom()))
        for (i, j), value in table.items():
            entries[i * size + j] += multiplier * value
    return arb_mat(size, size, entries) / arb(denominator)


def powers(x, order):
    values = [arb(1)]
    for _ in range(order):
        values.append(values[-1] * x)
    return values


def tilted_moments(k_minus, k_plus, order, alpha, radial):
    """Raw moments on the folded wedge, with one extra Gaussian factor."""
    p11 = 1 + 2 * (k_minus + k_plus)
    p22 = 2 + 2 * (k_minus + k_plus)
    p12 = 2 * (k_plus - k_minus)
    a = p11.sqrt()
    c = (p22 - p12 * p12 / p11).sqrt()
    slope = c / (a + p12 / a)
    beta = slope.atan()
    cosine = 1 / (1 + slope * slope).sqrt()
    sine = slope * cosine
    cos_power = powers(cosine, order + 1)
    sin_power = powers(sine, order + 1)

    @cache
    def angle(i, j):
        if i >= 2:
            return ((i - 1) * angle(i - 2, j) + cos_power[i - 1] * sin_power[j + 1]) / (i + j)
        if i == 1:
            return sin_power[j + 1] / (j + 1)
        if j >= 2:
            return ((j - 1) * angle(0, j - 2) - cosine * sin_power[j - 1]) / j
        return (beta, 1 - cosine)[j]

    x_power = powers(1 / a, order)
    y_power = powers(1 / c, order)
    mix_power = powers(-p12 / (p11 * c), order)
    prefactor = arb(2).sqrt() / (alpha * a * c)
    result = {}
    for degree in range(order + 1):
        for j in range(0, degree + 1, 2):
            i = degree - j
            integral = sum((comb(i, l) * x_power[i - l] * mix_power[l] * angle(i - l, j + l)
                            for l in range(i + 1)), arb(0))
            result[i, j] = prefactor * radial[degree] * y_power[j] * integral
    return result


def weighted_gram(basis, rho, allowance, order, alpha, radial, base):
    sigma = (1 - rho * rho).sqrt()
    rate = rho * rho / (3 * sigma * sigma)
    rates = [rate, 4 * rate / 3, 4 * rate]
    # The integrand is even only after adding each pair of reflected terms.
    moments = base.copy()
    for r in rates:
        minus = tilted_moments(r, arb(0), order, alpha, radial)
        plus = tilted_moments(arb(0), r, order, alpha, radial)
        require(minus[0, 0] > plus[0, 0], "The two folded marginal tilts must differ")
        for key in plus:
            moments[key] += (minus[key] + plus[key]) / 9
    for r in rates:
        for s in rates:
            joint = tilted_moments(r, s, order, alpha, radial)
            for key, value in joint.items():
                moments[key] += value / 81

    positive = 1 - allowance
    require(positive > 0, "Tail and claimed gap exceed the local lower bound")
    scale = (2 * alpha / arb.pi()) / positive
    mass = moments[0, 0]
    return arb_mat([[scale * (moments[i + u, j + t] - moments[i, j] * moments[u, t] / mass)
                     for u, t in basis] for i, j in basis])


def certify(slabs):
    ctx.prec = 640
    maximum = max(degree for _, degree in slabs)
    basis = [(d - j, j) for d in range(1, maximum + 1)
             for j in range(0, d + 1, 2)]
    coefficient_data = coefficient_numerators(basis, maximum)
    planar_data = planar_numerators(basis, maximum)
    u, u2, t2 = (basis.index(key) for key in [(1, 0), (2, 0), (0, 2)])
    require(coefficient_data[0] == {(u, u): 2}, "Incorrect degree-one control")
    require(coefficient_data[1] == {(u2, u2): 4, (u2, t2): 4, (t2, u2): 4, (t2, t2): 20},
            "Incorrect degree-two control")
    alpha = arb(2).sqrt().atan()
    radial = [arb(1), (arb.pi() / 2).sqrt()]
    for i in range(2, 2 * maximum + 1):
        radial.append(i * radial[-2])
    base = tilted_moments(arb(0), arb(0), 2*maximum, alpha, radial)
    require(base[0, 0].contains(1), "The untilted wedge measure does not have unit mass")
    left = START
    for endpoint, degree in slabs:
        right = fmpq(endpoint, 1000)
        require(left < right < 1, "Invalid correlation cover")
        sub_basis = [key for key in basis if sum(key) <= degree]
        weights = scaled_coefficients(left, right, degree)
        correction = planar_tail_coefficients(right, degree)
        tail = arb(right)**degree * arb.pi() / alpha
        gram = weighted_gram(sub_basis, arb(left), tail + arb(GAP), 2*degree,
                             alpha, radial, base)
        matrix = finite_matrix(coefficient_data[:degree] + planar_data[:degree],
                               len(sub_basis), weights + correction)
        power = matrix * gram
        for exponent in (2, 4, 8, 16, 32, 64):
            power = power * power
            trace = power.trace()
            if trace > 0 and trace < 1:
                break
        require(trace > 0 and trace < 1, f"No trace certificate on [{left},{right}]: {trace.str(14)}")
        yield {"left": str(left), "right": str(right), "degree": degree,
               "operator_tail_bound": tail.str(12), "trace_exponent": exponent,
               "trace": trace.str(12), "norm_bound": (trace**(arb(1) / exponent)).str(12)}
        left = right


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slabs", type=int, default=len(SLABS), choices=range(1, len(SLABS) + 1),
                        help="number of consecutive slabs to replay")
    args = parser.parse_args()
    rows = []
    for row in certify(SLABS[:args.slabs]):
        rows.append(row)
        print(f"Certified {row['left']} to {row['right']} (trace power {row['trace_exponent']})",
              flush=True, file=sys.stderr)
    print(json.dumps({"python": platform.python_version(), "python_flint": flint.__version__,
                      "precision_bits": ctx.prec,
                      "standard_gap_in_units_of_ell": str(GAP), "slabs": rows}, indent=2))


if __name__ == "__main__":
    main()
