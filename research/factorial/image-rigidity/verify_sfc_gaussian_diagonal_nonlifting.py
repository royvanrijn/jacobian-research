#!/usr/bin/env python3
"""Exact fifth jets for the four-column Gaussian nonlifting family.

Run from the repository root with .venv/bin/python.  Arithmetic is in
Q[k,a0,a1,a2,t], then Q[k,b]/(b^2+3b+15/4)[t]/(t^6).
The scale k remains an indeterminate throughout the certificate.
"""

import json
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial

import flint
import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx, fmpq_poly

LIMIT = 5
ORDERS = (1, 2, 3, 5)
RING = fmpq_mpoly_ctx.get(("a0", "a1", "a2", "k", "t"))
A0, A1, A2, K, T = RING.gens()
ZERO = (fmpq_poly([]), fmpq_poly([]))
ONE = (fmpq_poly([1]), fmpq_poly([]))
BETA = (fmpq_poly([]), fmpq_poly([1]))


def count_vectors(total, size=4):
    if size == 1:
        yield (total,)
    else:
        for head in range(total + 1):
            for tail in count_vectors(total - head, size - 1):
                yield (head,) + tail


def rational(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def scale_polynomial(value, variable):
    return sum(
        RING.constant(rational(c)) * K ** int(e[0])
        for e, c in sp.Poly(sp.expand(value), variable).terms()
    )


@lru_cache(None)
def power_sum(length, power):
    # Bernoulli(1)=+1/2: sum_{j=1}^{length*k} j^power.
    return sum(
        RING.constant(rational(sp.Rational(comb(power + 1, j), power + 1)
                               * sp.bernoulli(j)))
        * (length * K) ** (power + 1 - j)
        for j in range(power + 1)
    )


def gamma_series(order, x_length, degree, difference, limit):
    """Series in q=t^2 of the normalized finite gamma products."""
    log = [RING.constant(0)] * (limit + 1)
    for h in range(1, limit + 1):
        log[h] = sum(
            fmpq((-1) ** (h - 1) * comb(h - 1, r - 1),
                 r * (2 * order) ** r * 2 ** (h - r))
            * power_sum(x_length, r) * difference ** (h - r)
            for r in range(1, h + 1)
        )
        log[h] += (fmpq((-1) ** (h - 1), h * (2 * order) ** h)
                   * power_sum(order * degree - x_length, h))
    exp = [RING.constant(1)]
    for h in range(1, limit + 1):
        exp.append(sum(r * log[r] * exp[h - r] for r in range(1, h + 1)) / h)
    return exp


def normalized_forms(offsets, difference):
    """Exact jets in the monic cubic confluent binomial basis."""
    variable = sp.Symbol("k")
    matrix = sp.Matrix([
        [sp.prod(x * variable - j for j in range(i)) / sp.factorial(i)
         for x in offsets] for i in range(4)
    ])
    inverse = matrix.inv()
    # Clear k^3 in every coefficient; divide k^(3m) after cancellation.
    coefficients = [sum(
        scale_polynomial(variable ** 3 * inverse[j, i], variable)
        * (A0, A1, A2, RING.constant(1))[i] * T ** (3 - i)
        for i in range(4)) for j in range(4)]
    result = {}
    for order in ORDERS:
        top = 3 * order + LIMIT
        numerator = RING.constant(0)
        for alpha in count_vectors(order):
            weight = factorial(order)
            phase = RING.constant(1)
            for coefficient, multiplicity in zip(coefficients, alpha):
                weight //= factorial(multiplicity)
                phase *= coefficient ** multiplicity
            length = sum(x * n for x, n in zip(offsets, alpha))
            law = gamma_series(order, length, max(offsets), difference, top // 2)
            numerator += weight * phase * sum(c * T ** (2 * j)
                                               for j, c in enumerate(law))
        terms = []
        for exponent, coefficient in numerator.to_dict().items():
            assert exponent[4] >= 3 * order, (order, exponent, coefficient)
            if exponent[4] <= top:
                assert exponent[3] >= 3 * order
                terms.append((tuple(int(x) for x in exponent[:3])
                              + (int(exponent[3]) - 3 * order,
                                 int(exponent[4]) - 3 * order), coefficient))
        result[order] = terms
    return result


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def negative(x):
    return (-x[0], -x[1])


def multiply(x, y):
    return (x[0] * y[0] - fmpq(15, 4) * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0] - 3 * x[1] * y[1])


def inverse_constant(x):
    norm = x[0] ** 2 - 3 * x[0] * x[1] + fmpq(15, 4) * x[1] ** 2
    assert norm.degree() == 0 and norm[0]
    return ((x[0] - 3 * x[1]) / norm[0], -x[1] / norm[0])


def series_product(x, y):
    out = [ZERO] * (LIMIT + 1)
    for i, a in enumerate(x):
        for j, b in enumerate(y[:LIMIT + 1 - i]):
            out[i + j] = add(out[i + j], multiply(a, b))
    return out


def series_power(x, exponent):
    out = [ONE] + [ZERO] * LIMIT
    for _ in range(exponent):
        out = series_product(out, x)
    return out


def evaluate(form, order, coefficients):
    powers = [[series_power(a, i) for i in range(order + 1)] for a in coefficients]
    out = [ZERO] * (LIMIT + 1)
    for (i, j, k, scale, t_degree), value in form:
        term = series_product(series_product(powers[0][i], powers[1][j]), powers[2][k])
        scalar = (fmpq_poly([0] * scale + [value]), ZERO[1])
        for h in range(LIMIT + 1 - t_degree):
            out[h + t_degree] = add(out[h + t_degree], multiply(scalar, term[h]))
    return out


def solve_constant(matrix, rhs):
    augmented = [list(row) + [v] for row, v in zip(matrix, rhs)]
    for k in range(3):
        pivot = next(j for j in range(k, 3) if augmented[j][k] != ZERO)
        augmented[k], augmented[pivot] = augmented[pivot], augmented[k]
        scalar = inverse_constant(augmented[k][k])
        augmented[k] = [multiply(scalar, x) for x in augmented[k]]
        for j in range(3):
            if j != k:
                scalar = augmented[j][k]
                augmented[j] = [add(x, negative(multiply(scalar, y)))
                                for x, y in zip(augmented[j], augmented[k])]
    return [row[-1] for row in augmented]


def lift_selected(forms):
    coefficients = [[ZERO] * (LIMIT + 1), [BETA] + [ZERO] * LIMIT,
                    [ZERO] * (LIMIT + 1)]
    assert all(evaluate(forms[m], m, coefficients)[0] == ZERO for m in ORDERS)
    jacobian = []
    for m in ORDERS[:3]:
        base = evaluate(forms[m], m, coefficients)[1]
        row = []
        for k in range(3):
            perturbed = [list(x) for x in coefficients]
            perturbed[k][1] = ONE
            row.append(add(evaluate(forms[m], m, perturbed)[1], negative(base)))
        jacobian.append(row)
    for n in range(1, LIMIT + 1):
        rhs = [negative(evaluate(forms[m], m, coefficients)[n]) for m in ORDERS[:3]]
        solution = solve_constant(jacobian, rhs)
        for k in range(3):
            coefficients[k][n] = solution[k]
        assert all(evaluate(forms[m], m, coefficients)[n] == ZERO for m in ORDERS[:3])
    return evaluate(forms[5], 5, coefficients)


def direct_gamma(order, length, degree, difference, limit):
    def product(x, y):
        out = [Fraction(0)] * (limit + 1)
        for i, a in enumerate(x):
            for j, b in enumerate(y[:limit + 1 - i]):
                out[i + j] += a * b
        return out
    out = [Fraction(1)] + [Fraction(0)] * limit
    for j in range(1, length + 1):
        out = product(out, [Fraction(1), Fraction(difference, 2) + Fraction(j, 2 * order)])
    for j in range(1, order * degree - length + 1):
        out = product(out, [Fraction(1), Fraction(j, 2 * order)])
    denominator = [Fraction((-1) ** n * comb(length + n - 1, n))
                   * Fraction(difference, 2) ** n if length else Fraction(n == 0)
                   for n in range(limit + 1)]
    return product(out, denominator)


def compare_direct_products():
    comparisons = 0
    for scale in (1, 2):
        difference = scale + 2 * scale ** 3
        for order in ORDERS:
            limit = (3 * order + LIMIT) // 2
            for alpha in count_vectors(order):
                length = sum(x * n for x, n in zip((0, 3, 4, 9), alpha))
                symbolic = gamma_series(order, length, 9, K + 2 * K ** 3, limit)
                exact = direct_gamma(order, scale * length, 9 * scale, difference, limit)
                for p, q in zip(symbolic, exact):
                    value = sum(c * scale ** int(e[3]) for e, c in p.to_dict().items())
                    assert value == fmpq(q.numerator, q.denominator)
                comparisons += 1
    return comparisons


def gaussian_data():
    z, a, b, c = sp.symbols("z a b c")
    def expectation(p, m):
        return sp.expand(sum(v * sp.factorial2(j - 1) / sp.Integer(m) ** (j // 2)
                             for (j,), v in sp.Poly(sp.expand(p), z).terms() if j % 2 == 0))
    phase = z ** 3 + c * z ** 2 + b * z + a
    moments = [expectation(phase ** m, m) for m in ORDERS]
    reduced = [sp.expand(g.subs(a, -c)) for g in moments]
    assert list(sp.groebner(reduced, c, b)) == [c, b ** 2 + 3 * b + sp.Rational(15, 4)]
    # Separate chart at cubic coefficient zero: orders 1,3 force c=0,
    # and order 2 then forces b=0 (and a=0).
    lower = [expectation((a + b * z + c * z ** 2) ** m, m).subs(a, -c) for m in ORDERS]
    assert sp.expand(lower[2]) == -4 * c ** 3 / 9
    modulus = b ** 2 + 3 * b + sp.Rational(15, 4)
    def reduce(p):
        numerator, denominator = sp.fraction(sp.cancel(p))
        return sp.rem(sp.rem(numerator, modulus, b) * sp.invert(denominator, modulus, b), modulus, b)
    phase = z ** 3 + b * z
    jacobian = sp.Matrix([[reduce(m * expectation(phase ** (m - 1) * z ** j, m))
                          for j in range(3)] for m in ORDERS])
    assert reduce(jacobian[:3, :].det()) == (2 * b - 15) / 9
    u3 = reduce(-(jacobian[3, 2] - jacobian[3, 0]) / (jacobian[2, 2] - jacobian[2, 0]))
    u1 = reduce(-u3 * jacobian[2, 0] - jacobian[3, 0])
    null = [u1, 0, u3, 1]
    assert all(reduce(x) == 0 for x in sp.Matrix([null]) * jacobian)
    theta = reduce(sum(u * expectation(phase ** (m - 1) * sp.diff(phase, z), m)
                       for u, m in zip(null, ORDERS)))
    assert theta == sp.Rational(306126, 171875) - sp.Rational(1332, 6875) * b
    assert reduce(sum(u * expectation(z ** 2 * phase ** m, m)
                      for u, m in zip(null, ORDERS))) == 0
    return theta


def main():
    gaussian_data()
    residual = lift_selected(normalized_forms((0, 3, 4, 9), K + 2 * K ** 3))
    assert all(c == ZERO for c in residual[:5])
    expected = (
        fmpq_poly([0, 0, 0, fmpq(-10756503, 17187500), 0, fmpq(2774331, 8593750)]),
        fmpq_poly([0, 0, 0, fmpq(-57169, 171875), 0, fmpq(421386, 859375)]),
    )
    assert residual[5] == expected
    theta = (fmpq_poly([fmpq(306126, 171875)]), fmpq_poly([fmpq(-1332, 6875)]))
    correction = multiply(residual[5], inverse_constant(theta))
    assert 6 * correction[1] == fmpq_poly([
        0, 0, 0, fmpq(-354929845, 360419340), 0, fmpq(465440178, 360419340)])
    # A shift of one in the height difference changes the third residual
    # by -Theta/6, as the centered projection formula requires.
    shifted = lift_selected(normalized_forms((0, 3, 4, 9), K + 2 * K ** 3 + 1))
    assert shifted[3] == tuple(-x / 6 for x in theta)
    # Symmetric supports really lift: the residual must not exclude them.
    control = lift_selected(normalized_forms((0, 1, 2, 3), RING.constant(0)))
    assert all(c == ZERO for c in control)
    comparisons = compare_direct_products()
    print(json.dumps({
        "status": "PASS",
        "scope": "symbolic fifth-jet identity for every positive integer scale k; no numerical height cutoff computed",
        "orders": list(ORDERS),
        "gaussian_roots": "z^3+b*z, 4*b^2+12*b+15=0",
        "forced_difference": "k+2*k^3",
        "fifth_residual": "k^3*((5548662*k^2-10756503)+(8427720*k^2-5716900)*b)/17187500",
        "imaginary_height_correction": "k^3*(465440178*k^2-354929845)/360419340 * Im(b)",
        "direct_gamma_comparisons": comparisons,
        "symmetric_lifting_control": "passed through degree five; exact lifting follows separately by reflection",
        "python_flint": flint.__version__,
        "sympy": sp.__version__,
        "assurance": "internal exact symbolic replay; not independent review or formal verification",
    }, indent=2))


if __name__ == "__main__":
    main()
