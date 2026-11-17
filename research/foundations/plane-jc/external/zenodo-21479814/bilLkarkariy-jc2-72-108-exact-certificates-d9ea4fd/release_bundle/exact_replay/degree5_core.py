#!/usr/bin/env python3
"""Minimal exact arithmetic for L = Q[w]/(w^5-w^4+3w^3+3w^2+26)."""

from __future__ import annotations

from dataclasses import dataclass

from flint import fmpq, fmpq_poly


MOD = fmpq_poly([26, 0, 3, 3, -1, 1])


@dataclass(frozen=True)
class L:
    p: fmpq_poly

    def __init__(self, value=0):
        if isinstance(value, L):
            poly = value.p
        elif isinstance(value, fmpq_poly):
            poly = value
        elif isinstance(value, fmpq):
            poly = fmpq_poly([value])
        elif isinstance(value, int):
            poly = fmpq_poly([value])
        elif isinstance(value, (list, tuple)):
            poly = fmpq_poly(value)
        else:
            raise TypeError(type(value))
        object.__setattr__(self, "p", poly % MOD)

    def __add__(self, other):
        return L(self.p + L(other).p)

    __radd__ = __add__

    def __neg__(self):
        return L(-self.p)

    def __sub__(self, other):
        return self + (-L(other))

    def __rsub__(self, other):
        return L(other) - self

    def __mul__(self, other):
        return L((self.p * L(other).p) % MOD)

    __rmul__ = __mul__

    def inv(self):
        if not self:
            raise ZeroDivisionError
        gcd, coefficient, _ = self.p.xgcd(MOD)
        if gcd.degree() != 0:
            raise ZeroDivisionError("nonunit in degree-five quotient")
        return L(coefficient / gcd[0])

    def __truediv__(self, other):
        return self * L(other).inv()

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inv() ** (-exponent)
        result = L(1)
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def __bool__(self):
        return not self.p.is_zero()

    def __eq__(self, other):
        return self.p == L(other).p


def decode_l(serialized):
    coefficients = [fmpq(0) for _ in range(5)]
    for degree, (numerator, denominator) in serialized.items():
        coefficients[int(degree)] = fmpq(int(numerator), int(denominator))
    return L(coefficients)


def encode_l(value):
    value = L(value)
    return {
        i: (int(value.p[i].p), int(value.p[i].q))
        for i in range(len(value.p))
        if value.p[i]
    }


class Poly:
    def __init__(self, nvars, terms=None):
        self.nvars = nvars
        self.terms = {
            tuple(m): L(c) for m, c in (terms or {}).items() if L(c)
        }

    @classmethod
    def constant(cls, nvars, value):
        value = L(value)
        return cls(nvars, {(0,) * nvars: value}) if value else cls(nvars)

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = Poly.constant(self.nvars, other)
        if other.nvars != self.nvars:
            raise ValueError("incompatible polynomial dimensions")
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            value = terms.get(monomial, L(0)) + coefficient
            if value:
                terms[monomial] = value
            else:
                terms.pop(monomial, None)
        return Poly(self.nvars, terms)

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.nvars, {m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, Poly):
            other = Poly.constant(self.nvars, other)
        if other.nvars != self.nvars:
            raise ValueError("incompatible polynomial dimensions")
        terms = {}
        for m, c in self.terms.items():
            for n, d in other.terms.items():
                monomial = tuple(a + b for a, b in zip(m, n))
                terms[monomial] = terms.get(monomial, L(0)) + c * d
        return Poly(self.nvars, terms)

    __rmul__ = __mul__

    def __bool__(self):
        return bool(self.terms)

    @property
    def degree(self):
        return max(map(sum, self.terms), default=-1)


def decode_poly(serialized, nvars=3):
    return Poly(nvars, {tuple(m): decode_l(c) for m, c in serialized.items()})


def encode_poly(poly):
    return {m: encode_l(c) for m, c in poly.terms.items()}


def specialize_first_variable_zero(poly):
    return Poly(
        poly.nvars - 1,
        {m[1:]: c for m, c in poly.terms.items() if m[0] == 0},
    )


def sign_substitution(poly, signs):
    if len(signs) != poly.nvars:
        raise ValueError("one sign is required per variable")
    return Poly(
        poly.nvars,
        {
            m: c * (-1 if sum(e for e, sign in zip(m, signs) if sign == -1) % 2 else 1)
            for m, c in poly.terms.items()
        },
    )
