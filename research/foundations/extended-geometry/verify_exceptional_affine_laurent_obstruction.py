#!/usr/bin/env python3
"""Replay the finite-modular affine Laurent obstructions on B3 and H3.

Requires Python 3.12 and Singular 4.3.2 (or compatible versions).
The coefficient formulas are in the field and compact extensions note.
"""

from collections import defaultdict
from functools import cache
from math import comb
from pathlib import Path
import shutil
import subprocess
import tempfile


PRIME = 131


def ratio(numerator, denominator):
    return numerator * pow(denominator, -1, PRIME) % PRIME


def multiply(left, right):
    result = defaultdict(int)
    for left_exponents, a in left.items():
        for right_exponents, b in right.items():
            exponents = tuple(i + j for i, j in zip(left_exponents, right_exponents))
            result[exponents] += a * b
    return {key: value % PRIME for key, value in result.items() if value % PRIME}


def powers(polynomial, size, highest):
    result = [{(0,) * size: 1}]
    for _ in range(highest):
        result.append(multiply(result[-1], polynomial))
    return result


@cache
def odd_product(length):
    result = 1
    for j in range(length):
        result = result * (2 * j + 1) % PRIME
    return result


def b3_moments(highest):
    p2 = powers({(1, 1, 0): 1, (1, 0, 1): 1, (0, 1, 1): 1}, 3, highest)
    p3 = powers({(1, 1, 1): 1}, 3, highest)

    @cache
    def moment(i, j):
        return sum(
            coefficient * ratio(
                odd_product(a) * odd_product(b) * odd_product(c),
                odd_product(a + b + c + 1),
            )
            for (a, b, c), coefficient in multiply(p2[i], p3[j]).items()
        ) % PRIME

    assert (moment(1, 0), moment(0, 1)) == (ratio(1, 5), ratio(1, 105))
    return moment


def h3_moments(highest):
    # The coefficients of p^-1 through p^-4 in -F(v,p)/(18662400 p^4).
    scale = ratio(1, 18662400)
    coefficients = [
        None,
        {(3,): -138240 * scale, (2,): 2073600 * scale, (1,): 4147200 * scale},
        {(3,): 348160 * scale, (2,): -46080 * scale, (1,): -82944 * scale},
        {(3,): -35840 * scale, (2,): -9216 * scale},
        {(3,): 1024 * scale},
    ]
    # Inverse square-root series from 2F Y' + F' Y = 0.
    expansion = [{(0,): 1}]
    for n in range(1, 2 * highest + 1):
        coefficient = defaultdict(int)
        for k in range(1, min(4, n) + 1):
            multiplier = ratio(-(2 * n - k), 2 * n)
            for power, value in multiply(coefficients[k], expansion[n - k]).items():
                coefficient[power] += multiplier * value
        expansion.append({
            power: value % PRIME for power, value in coefficient.items() if value % PRIME
        })

    @cache
    def moment(i, j):
        return sum(
            coefficient * ratio(pow(36, j, PRIME), 2 * (degree - j) + 1)
            for (degree,), coefficient in expansion[i + 2 * j].items()
        ) % PRIME

    # These values also follow by direct sphere integration in the x_i^2.
    assert (moment(1, 0), moment(0, 1), moment(1, 1)) == (
        ratio(-1, 21), ratio(92, 231), ratio(-172, 3927)
    )
    return moment


def centered_moments(raw):
    mean_x, mean_y = raw(1, 0), raw(0, 1)

    @cache
    def moment(i, j):
        return sum(
            comb(i, a) * comb(j, b) * pow(-mean_x, i - a, PRIME)
            * pow(-mean_y, j - b, PRIME) * raw(a, b)
            for a in range(i + 1) for b in range(j + 1)
        ) % PRIME

    assert (moment(0, 0), moment(1, 0), moment(0, 1)) == (1, 0, 0)
    return moment


def coefficient_polynomial(angular_monomials, first_parameter):
    return {
        monomial + tuple(int(i == first_parameter + j) for i in range(8)): 1
        for j, monomial in enumerate(angular_monomials)
    }


def accumulate_moment(result, polynomial, multiplier, moment):
    for (i, j, *parameters), coefficient in polynomial.items():
        result[tuple(parameters)] += multiplier * coefficient * moment(i, j)


def trinomial_equations(moment, highest):
    center = powers(coefficient_polynomial([(1, 0), (0, 1)], 0), 10, highest)
    product = powers(coefficient_polynomial(
        [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)], 2,
    ), 10, highest // 2)
    equations = {}
    for m in range(2, highest + 1):
        result = defaultdict(int)
        for k in range(m // 2 + 1):
            weight = comb(m, 2 * k) * comb(2 * k, k)
            accumulate_moment(result, multiply(center[m - 2 * k], product[k]), weight, moment)
        equations[m] = result
    return equations


def four_term_equations(moment, highest):
    center = powers(coefficient_polynomial([(1, 0), (0, 1)], 0), 10, highest)
    affine = [(0, 0), (1, 0), (0, 1)]
    positive_one = powers(coefficient_polynomial(affine, 2), 10, highest // 2)
    positive_two = powers(coefficient_polynomial(affine, 5), 10, highest // 3)
    equations = {}
    for m in range(2, highest + 1):
        result = defaultdict(int)
        for ell in range(m // 3 + 1):
            for k in range((m - 3 * ell) // 2 + 1):
                degree = m - 2 * k - 3 * ell
                weight = comb(m, degree) * comb(2 * k + 3 * ell, ell) * comb(2 * k + 2 * ell, k)
                polynomial = multiply(multiply(center[degree], positive_one[k]), positive_two[ell])
                accumulate_moment(result, polynomial, weight, moment)
        equations[m] = result
    return equations


def singular_polynomial(polynomial, names):
    terms = []
    for degrees, coefficient in polynomial.items():
        if coefficient % PRIME:
            monomial = "".join(f"*{name}^{degree}" for name, degree in zip(names, degrees) if degree)
            terms.append(str(coefficient % PRIME) + monomial)
    return "+".join(terms) or "0"


def verify(singular, label, names, equations, extra_equation, expected):
    source = [f"ring r={PRIME},({','.join(names)}),dp;"]
    generators = []
    for order, polynomial in equations.items():
        generators.append(f"f{order}")
        source.append(f"poly f{order}={singular_polynomial(polynomial, names)};")
    if extra_equation:
        generators.append("rankD")
        source.append(f"poly rankD={extra_equation};")
    source.append(f"ideal I={','.join(generators)};")
    source.append(
        """
        option(redSB);
        ideal G=std(I);
        dim(G);
        intvec powers;
        intvec e;
        int i, j, variableIndex, supportSize;
        for (i=1; i<=size(G); i++) {
            e=leadexp(G[i]);
            variableIndex=0;
            supportSize=0;
            for (j=1; j<=nvars(basering); j++) {
                if (e[j]!=0) {
                    supportSize++;
                    variableIndex=j;
                }
            }
            if (supportSize==1) {
                powers[variableIndex]=e[variableIndex];
            }
        };
        powers;
        quit;
        """
    )
    with tempfile.TemporaryDirectory() as directory:
        replay_file = Path(directory) / "verify.sing"
        replay_file.write_text("\n".join(source))
        result = subprocess.run(
            [singular, "-q", str(replay_file)], capture_output=True,
            text=True, check=True, timeout=45,
        )
    dimension, powers_line = result.stdout.strip().splitlines()
    leading_powers = tuple(int(power) for power in powers_line.split(","))
    assert dimension == "0" and leading_powers == expected, (label, result.stdout)
    print(f"{label}: dimension zero, leading pure powers {leading_powers}")


def main():
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular 4.3.2 or a compatible version is required")
    for label, moments in [("B3", b3_moments), ("H3", h3_moments)]:
        moment = centered_moments(moments(10))
        verify(
            singular, f"{label}, affine trinomial through order 10",
            ["a", "b", "d0", "d1", "d2", "d3", "d4", "d5"],
            trinomial_equations(moment, 10),
            "4*d0*d3*d5-d0*d4^2-d1^2*d5+d1*d2*d4-d2^2*d3",
            (2, 4, 6, 6, 6, 6, 7, 13),
        )
        verify(
            singular, f"{label}, constant-pole affine four-term family through order 9",
            ["a", "b", "c0", "c1", "c2", "d0", "d1", "d2"],
            four_term_equations(moment, 9),
            None, (2, 4, 4, 4, 5, 8, 8, 10),
        )


if __name__ == "__main__":
    main()
