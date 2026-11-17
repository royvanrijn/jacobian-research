"""Exact local moment arithmetic at the spherical Rodrigues cubic.

Uses the charge selection rule before evaluating any Gaussian/sphere
moment. No general coefficient-support search is performed.
"""
from functools import lru_cache
from fractions import Fraction
import math

import sympy as S
from conformal_euler_cubic import z, w, t, harmonic_basis, fields

a, b = S.symbols('a b')
F0 = w+z*t-z**3/2
G1 = t*z*w-t**3/3+S.Rational(17, 21)*t+S.Rational(10, 3)*w*w
G2 = w*z*z-2*t*t*z-5*t*w+z
BASIS = tuple(f for level in harmonic_basis().values() for f in level.values())


def odd_factorial(k):
    return math.prod(range(k, 0, -2))


@lru_cache(maxsize=None)
def marked_monomial(power, i, j, k):
    """Return h(z^i w^j t^k F0^power) as a rational.

    For counts A,B,C of -z^3/2,zt,w, respectively, charge zero
    forces 2A+B=(power+j-i)/2. All surviving terms have one
    ordinary degree, so their spherical denominator is common.
    """
    if power < 0 or (power+j-i) % 2:
        return S.Integer(0)
    selected = (power+j-i)//2
    if selected < 0:
        return S.Integer(0)
    total = Fraction(0)
    degree = 2*i+3*selected+k
    for cubic in range(selected//2+1):
        quadratic = selected-2*cubic
        linear = power-selected+cubic
        if linear < 0 or (k+quadratic) % 2:
            continue
        multiplicity = math.factorial(power)//(math.factorial(cubic)*math.factorial(quadratic)*math.factorial(linear))
        numerator = multiplicity*(-1)**cubic*math.factorial(i+selected+cubic)*odd_factorial(k+quadratic-1)
        total += Fraction(numerator, 2**cubic)
    return S.Rational(total.numerator, total.denominator*odd_factorial(degree+1))


def marked(marker, power):
    if power < 0:
        return S.Integer(0)
    return S.expand(sum(coefficient*marked_monomial(power, *exponents)
                        for exponents, coefficient in S.Poly(S.expand(marker), z, w, t).terms()))


def coefficient_vector(poly):
    monomials = sorted(set().union(*(set(S.Poly(f, z, w, t).monoms()) for f in BASIS)))
    matrix = S.Matrix([[S.Poly(f, z, w, t).coeff_monomial(m) for f in BASIS] for m in monomials])
    target = S.Matrix([S.Poly(S.expand(poly), z, w, t).coeff_monomial(m) for m in monomials])
    solution = next(iter(S.linsolve((matrix, target))))
    assert S.expand(sum(c*f for c, f in zip(solution, BASIS))-poly) == 0
    return S.Matrix(solution)


@lru_cache(maxsize=1)
def linear_data():
    jacobian = S.Matrix([[n*marked(f, n-1) for f in BASIS] for n in range(1, 11)])
    pivots = tuple(jacobian.rref()[1])
    assert len(pivots) == 10
    orbit = [F0, *fields(F0)[1:]]
    directions = orbit+[G1, G2]
    kernel = S.Matrix.hstack(*(coefficient_vector(f) for f in directions))
    assert jacobian*kernel == S.zeros(10, 6)
    assert kernel.rank() == 6
    coordinates = S.Matrix.hstack(*(S.eye(16)[:, j] for j in pivots), kernel)
    assert coordinates.det() != 0
    return jacobian, pivots, jacobian[:, pivots].inv(), coordinates


def third_order():
    jacobian, pivots, inverse, coordinates = linear_data()
    G = a*G1+b*G2
    quadratic = S.Matrix([S.binomial(n, 2)*marked(G*G, n-2) for n in range(1, 11)])
    correction_vector = inverse*(-quadratic)
    correction = S.expand(sum(c*BASIS[j] for c, j in zip(correction_vector, pivots)))
    quadratic_tail = {n: S.factor(S.binomial(n, 2)*marked(G*G, n-2)) for n in range(11, 21)}
    assert not any(quadratic_tail.values())
    cubic_tail = {n: S.factor(n*(n-1)*marked(G*correction, n-2)+S.binomial(n, 3)*marked(G**3, n-3))
                  for n in range(11, 31)}
    return {'pivots': pivots, 'quadratic_correction': correction,
            'cubic_obstructions': {n: f for n, f in cubic_tail.items() if f != 0},
            'quadratic_tail': quadratic_tail}


def implicit_jets(order, progress=None, initial=None):
    """Solve the first ten moments on an affine slice transverse to the orbit.

    The slice is F0+aG1+bG2+sum u_j B_j for the ten pivot basis vectors.
    At each homogeneous parameter order, a single constant invertible
    Jacobian determines the new correction. Higher moment residuals are
    retained as obstructions, rather than silently set to zero.
    """
    jacobian, pivots, inverse, coordinates = linear_data()
    jets = {1: S.expand(a*G1+b*G2 if initial is None else initial)}
    obstructions = {}
    for degree in range(2, order+1):
        powers = {1: dict(jets)}
        for exponent in range(2, degree+1):
            powers[exponent] = {}
            for weight in range(exponent, degree+1):
                value = sum(jets[j]*powers[exponent-1].get(weight-j, 0)
                            for j in jets if j < weight)
                if value != 0:
                    powers[exponent][weight] = S.expand(value)
        markers = {k: powers[k].get(degree, S.Integer(0)) for k in range(2, degree+1)}
        residuals = {n: S.factor(sum(S.binomial(n, k)*marked(marker, n-k)
                                     for k, marker in markers.items() if k <= n))
                     for n in range(1, 10*degree+1)}
        correction = inverse*(-S.Matrix([residuals[n] for n in range(1, 11)]))
        jets[degree] = S.expand(sum(c*BASIS[j] for c, j in zip(correction, pivots)))
        tail = {n: value for n, value in residuals.items() if n > 10 and value != 0}
        obstructions[degree] = tail
        if progress is not None:
            progress(degree, jets[degree], tail)
    return jets, obstructions
