"""Exact factorial/Gamma/Dirichlet bridges into the existing Gaussian producer.

This module makes no pure-null or universal-conjecture decision. All coefficients
must be exact SymPy expressions. A producer input uses the repository's existing
QQ/simple-number-field format; its certificates retain their original scope.
"""
from math import prod
import sympy as s


def factorial(poly, variables):
    return s.expand(sum(c * prod(s.factorial(a) for a in exps)
                        for exps, c in s.Poly(poly, *variables).terms()))


def gamma(poly, variables, alpha):
    if len(variables) != len(alpha):
        raise ValueError('one shape per variable required')
    return s.expand(sum(c * prod(s.rf(ai, a) for ai, a in zip(alpha, exps))
                        for exps, c in s.Poly(poly, *variables).terms()))


def dirichlet(poly, variables, alpha):
    if len(variables) != len(alpha):
        raise ValueError('one shape per variable required')
    return s.expand(sum(c * prod(s.rf(ai, a) for ai, a in zip(alpha, exps))
                        / s.rf(sum(alpha), sum(exps))
                        for exps, c in s.Poly(poly, *variables).terms()))


def homogeneous_lift(poly, variables, new_variables, degree=None):
    p = s.Poly(poly, *variables)
    if len(variables) != len(new_variables):
        raise ValueError('variable dimensions differ')
    needed = max((sum(exps) for exps, c in p.terms() if c != 0), default=0)
    degree = needed if degree is None else degree
    if not isinstance(degree, int) or degree < needed:
        raise ValueError('homogenizing degree is too small')
    radial = sum(new_variables)
    return s.expand(sum(c * prod(v**a for v, a in zip(new_variables, exps))
                        * radial**(degree-sum(exps)) for exps, c in p.terms()))


def gaussian_lift(poly, variables, block_sizes=None):
    sizes = [2] * len(variables) if block_sizes is None else list(block_sizes)
    if len(sizes) != len(variables) or any(type(k) is not int or k < 1 for k in sizes):
        raise ValueError('one positive integer block size per variable required')
    x = s.symbols('x:' + str(sum(sizes)))
    substitution, blocks, offset = {}, [], 0
    for y, k in zip(variables, sizes):
        block = x[offset:offset+k]
        blocks.append(block)
        substitution[y] = sum(v*v for v in block) / 2
        offset += k
    return s.expand(poly.subs(substitution, simultaneous=True)), x, blocks


def gaussian(poly, variables):
    return s.expand(sum(c * prod(s.factorial2(a-1) for a in exps)
                        for exps, c in s.Poly(poly, *variables).terms()
                        if not any(a % 2 for a in exps)))


def radial_monomial(exponents, block_sizes):
    """Return (Y-exponents, nonzero Reynolds factor) for an even X-monomial.

    An odd coordinate exponent gives a zero averaged marker, represented by
    (None, 0). This is an identity for every invariant phase and every power.
    """
    if any(type(a) is not int or a < 0 for a in exponents):
        raise ValueError('nonnegative integer exponents required')
    if any(type(k) is not int or k < 1 for k in block_sizes):
        raise ValueError('positive integer block sizes required')
    if len(exponents) != sum(block_sizes):
        raise ValueError('marker and block dimensions differ')
    if any(a % 2 for a in exponents):
        return None, s.Integer(0)
    offset, ys, factor = 0, [], s.Integer(1)
    for k in block_sizes:
        half = [a//2 for a in exponents[offset:offset+k]]
        b = sum(half)
        ys.append(b)
        factor *= 2**b * prod(s.rf(s.Rational(1, 2), a) for a in half) / s.rf(s.Rational(k, 2), b)
        offset += k
    return tuple(ys), s.factor(factor)


def producer_input(poly, variables, field=None):
    """Encode a Gaussian phase for gaussian_source_hom_producer.py.

    A nonrational coefficient uses the producer's symbol `a` and requires an
    explicitly supplied field spec. Float coefficients are always rejected.
    """
    if s.sympify(poly).has(s.Float):
        raise ValueError('exact coefficients required')
    terms = []
    parameter = s.Symbol('a')
    for exps, c in s.Poly(poly, *variables).terms():
        if c == 0:
            continue
        if c.free_symbols - ({parameter} if field else set()):
            raise ValueError('coefficient outside the supplied exact field')
        if field is None and not c.is_Rational:
            raise ValueError('nonrational coefficient requires an explicit field representation')
        if field is not None:
            try:
                for value in s.cancel(c).as_numer_denom():
                    s.Poly(value, parameter, domain=s.QQ)
            except (s.PolynomialError, s.CoercionFailed) as exc:
                raise ValueError('represent algebraic coefficients as rational expressions in a') from exc
        terms.append([list(exps), str(c).replace('**', '^')])
    result = {'dimension': len(variables), 'phase_terms': terms}
    if field is not None:
        result['field'] = field
    return result


def rigidity_phase(mu, variables):
    if len(mu) != len(variables):
        raise ValueError('one parameter per variable required')
    return s.expand(prod(variables) * sum(c*y for c, y in zip(mu, variables)))
