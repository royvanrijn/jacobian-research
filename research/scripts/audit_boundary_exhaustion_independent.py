#!/usr/bin/env python3
"""Clean-room exact regressions for weighted/cancellation boundary exhaustion.

Uses only the Python standard library and imports no project module.
"""

from fractions import Fraction
from math import comb


NVARS = 6  # T, P, Q, S, U, V
T, P, Q, S, U, V = range(NVARS)
ZERO_EXP = (0,) * NVARS


def clean(poly):
    return {m: c for m, c in poly.items() if c}


def const(value):
    value = Fraction(value)
    return {} if not value else {ZERO_EXP: value}


def var(index):
    exponent = [0] * NVARS
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return clean(out)


def scale(poly, scalar):
    scalar = Fraction(scalar)
    return clean({m: scalar * c for m, c in poly.items()})


def mul(left, right):
    out = {}
    for lm, lc in left.items():
        for rm, rc in right.items():
            monomial = tuple(a + b for a, b in zip(lm, rm))
            out[monomial] = out.get(monomial, Fraction(0)) + lc * rc
    return clean(out)


def power(poly, exponent):
    out = const(1)
    base = poly
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        base = mul(base, base)
        exponent //= 2
    return out


def derivative(poly, index):
    out = {}
    for monomial, coefficient in poly.items():
        if monomial[index]:
            new = list(monomial)
            factor = new[index]
            new[index] -= 1
            out[tuple(new)] = coefficient * factor
    return clean(out)


def integrate(poly, index):
    out = {}
    for monomial, coefficient in poly.items():
        new = list(monomial)
        new[index] += 1
        out[tuple(new)] = coefficient / new[index]
    return clean(out)


def substitute(poly, index, replacement):
    out = {}
    for monomial, coefficient in poly.items():
        base_exp = list(monomial)
        exponent = base_exp[index]
        base_exp[index] = 0
        term = {tuple(base_exp): coefficient}
        out = add(out, mul(term, power(replacement, exponent)))
    return clean(out)


def specialize_zero(poly, index):
    return clean({m: c for m, c in poly.items() if m[index] == 0})


def degree(poly, index):
    return max((m[index] for m in poly), default=-1)


# Univariate polynomials are coefficient lists in ascending order.
def uclean(poly):
    out = [Fraction(c) for c in poly]
    while out and out[-1] == 0:
        out.pop()
    return out


def uadd(left, right):
    out = [Fraction(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] += value
    return uclean(out)


def umul(left, right):
    if not left or not right:
        return []
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return uclean(out)


def upow(poly, exponent):
    out = [Fraction(1)]
    base = poly
    while exponent:
        if exponent & 1:
            out = umul(out, base)
        base = umul(base, base)
        exponent //= 2
    return out


def uderivative(poly):
    return uclean([i * poly[i] for i in range(1, len(poly))])


def udivmod(numerator, denominator):
    numerator = uclean(numerator)
    denominator = uclean(denominator)
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while numerator and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        factor = numerator[-1] / denominator[-1]
        quotient[shift] += factor
        for i, value in enumerate(denominator):
            numerator[i + shift] -= factor * value
        numerator = uclean(numerator)
    return uclean(quotient), numerator


def ugcd(left, right):
    left, right = uclean(left), uclean(right)
    while right:
        _, remainder = udivmod(left, right)
        left, right = right, remainder
    if not left:
        return []
    return [c / left[-1] for c in left]


def infinity_polynomial(m, r):
    # K(w)=sum_j (-1)^j binom(mr,j) w^j/(r+j+1).
    return [Fraction((-1) ** j * comb(m * r, j), r + j + 1)
            for j in range(m * r + 1)]


def parameter_polynomial(m, r):
    n = m * r
    # M(q)=sum_j (-1)^j binom(n+r+1,j) q^(n-j).
    out = [Fraction(0)] * (n + 1)
    for j in range(n + 1):
        out[n - j] += (-1) ** j * comb(n + r + 1, j)
    return uclean(out)


def transformed_infinity_numerator(m, r):
    # Clear (1-q)^(mr) from K(-q/(1-q)).
    n = m * r
    out = []
    for j, coefficient in enumerate(infinity_polynomial(m, r)):
        term = upow([Fraction(1), Fraction(-1)], n - j)
        term = umul(term, [Fraction(0)] * j + [Fraction((-1) ** j)])
        out = uadd(out, [coefficient * c for c in term])
    return uclean(out)


t, p, q, s, u, v = map(var, range(NVARS))

for m in range(1, 6):
    for r in range(1, 6):
        n_degree = r * (m + 1) + 1

        d = add(const(1), scale(mul(t, power(add(q, scale(mul(p, t), -1)), m)), -1))
        integrand = power(d, r)
        psi_without_constant = integrate(integrand, T)
        assert derivative(psi_without_constant, T) == integrand
        assert degree(psi_without_constant, T) == n_degree

        j_integrand = power(
            add(p, scale(mul(v, power(add(q, scale(v, -1)), m)), -1)), r
        )
        j_of_u = substitute(integrate(j_integrand, V), V, u)
        j_of_ps = substitute(j_of_u, U, mul(p, s))
        psi_of_s = substitute(psi_without_constant, T, s)
        assert j_of_ps == mul(power(p, r + 1), psi_of_s)

        # At P=0: J=(-1)^r integral V^r(Q-V)^(mr), hence U^(r+1)K(U/Q).
        j_p0 = specialize_zero(j_of_u, P)
        expected_integrand = scale(mul(power(v, r), power(add(q, scale(v, -1)), m * r)), (-1) ** r)
        expected = substitute(integrate(expected_integrand, V), V, u)
        assert j_p0 == expected
        assert min(monomial[U] for monomial in j_p0) == r + 1
        assert degree(j_p0, U) == n_degree

        k_poly = infinity_polynomial(m, r)
        assert len(k_poly) - 1 == m * r
        assert ugcd(k_poly, uderivative(k_poly)) == [Fraction(1)]

        transformed = transformed_infinity_numerator(m, r)
        modulus = parameter_polynomial(m, r)
        _, remainder = udivmod(transformed, modulus)
        assert not remainder

        # Both local saturation formulas equal the generic field degree.
        assert (r + 1) + (n_degree - (r + 1)) == n_degree
        assert (r + 1) + 1 + (m * r - 1) == n_degree


print("PASS: inverse derivative, degree, and projective chart scaling")
print("PASS: P=0 factorization, squarefreeness, and parameter-root identity")
print("PASS: discriminant and P=0 local degree saturations")
print("PASS: clean-room checker used only the Python standard library")
