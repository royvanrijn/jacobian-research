"""Exact adapters for the conformally graded spherical cubic architecture.

Failure to find a grading is not a purity or safety decision.
"""
from functools import lru_cache
import math

import sympy as S

z, w, t = S.symbols('z w t')
W, V, Z, Y = S.symbols('W V Z Y')
q = 2*z*w+t*t


def laplace(f):
    return S.expand(2*S.diff(f, z, w)+S.diff(f, t, 2))


def harmonic_parts(f):
    polynomial = S.Poly(S.expand(f), z, w, t)
    if polynomial.total_degree() > 3:
        raise ValueError('the adapter requires degree at most three')
    pieces = {j: S.Add(*(c*z**a*w**b*t**d for (a, b, d), c in polynomial.terms()
                         if a+b+d == j)) for j in range(4)}
    trace3, trace2 = laplace(pieces[3])/10, laplace(pieces[2])/6
    return {0: S.expand(pieces[0]+trace2),
            1: S.expand(pieces[1]+trace3),
            2: S.expand(pieces[2]-q*trace2),
            3: S.expand(pieces[3]-q*trace3)}


def harmonic_basis():
    return {
        0: {0: S.Integer(1)},
        1: {-1: w, 0: t, 1: z},
        2: {-2: w*w, -1: w*t, 0: t*t-z*w, 1: z*t, 2: z*z},
        3: {-3: w**3, -2: w*w*t, -1: w*(t*t-z*w/2),
            0: t**3-3*z*w*t, 1: z*(t*t-z*w/2), 2: z*z*t, 3: z**3},
    }


def fields(f):
    return tuple(S.expand(g) for g in (
        z*S.diff(f, z)+w*S.diff(f, w)+t*S.diff(f, t),
        z*S.diff(f, z)-w*S.diff(f, w),
        z*S.diff(f, t)-t*S.diff(f, w),
        w*S.diff(f, t)-t*S.diff(f, z),
    ))


def conformal_certificate(f):
    """Solve (Eu+alpha H+beta Lplus+gamma Lminus)f=delta f.

    The canonical harmonic representative is used. The returned identity
    certifies only existence of a grading; purity is a separate premise.
    """
    parts = harmonic_parts(f)
    f = S.expand(sum(parts.values()))
    eu, h, plus, minus = fields(f)
    unknowns = S.symbols('_euler_alpha _euler_beta _euler_gamma _euler_delta', cls=S.Dummy)
    expression = S.expand(eu+unknowns[0]*h+unknowns[1]*plus+unknowns[2]*minus-unknowns[3]*f)
    equations = S.Poly(expression, z, w, t).coeffs()
    matrix, rhs = S.linear_eq_to_matrix(equations, unknowns)
    solutions = S.linsolve((matrix, rhs), unknowns)
    if solutions == S.EmptySet:
        return {'recognized': False, 'rank': matrix.rank(),
                'augmented_rank': matrix.row_join(rhs).rank(),
                'boundary': 'No conformal Euler grading; no purity or safety conclusion.'}
    solution = next(iter(solutions))
    free = {u: S.Integer(0) for u in unknowns if any(u in v.free_symbols for v in solution)}
    solution = tuple(S.simplify(v.subs(free)) for v in solution)
    assert S.expand(expression.subs(dict(zip(unknowns, solution)))) == 0
    return {'recognized': True, 'alpha': solution[0], 'beta': solution[1],
            'gamma': solution[2], 'delta': solution[3],
            'boundary': 'A verified grading identity. Purity remains a separate premise.'}


def sphere_average(f):
    value = S.Integer(0)
    for (a, b, c), coefficient in S.Poly(S.expand(f), z, w, t).terms():
        if a == b and c % 2 == 0:
            value += coefficient*math.factorial(a)*S.factorial2(c-1)/S.factorial2(2*a+c+1)
    return S.simplify(value)


def gaussian_average(f):
    value = S.Integer(0)
    for (a, b, c), coefficient in S.Poly(S.expand(f), z, w, t).terms():
        if a == b and c % 2 == 0:
            value += coefficient*math.factorial(a)*S.factorial2(c-1)
    return S.expand(value)


def shear(f, parameter):
    return S.expand(f.subs({z: z, t: t-parameter*z*z,
                           w: w+parameter*z*t-parameter**2*z**3/2}, simultaneous=True))


@lru_cache(maxsize=1)
def hopf_matrix():
    s = W*Z+V*Y
    substitution = {z: W*Y, w: 2*V*Z, t: W*Z-V*Y}
    source_basis = [(j, k, f) for j, level in harmonic_basis().items() for k, f in level.items()]
    target_basis = [W**(3-i)*V**i*Z**(3-j)*Y**j for i in range(4) for j in range(4)]
    lifts = [S.Poly(S.expand(s**(3-j)*f.subs(substitution, simultaneous=True)), W, V, Z, Y)
             for j, k, f in source_basis]
    matrix = S.Matrix([[f.coeff_monomial(v) for f in lifts] for v in target_basis])
    return matrix, source_basis, target_basis


def biform_to_phase(f):
    polynomial = S.Poly(S.expand(f), W, V, Z, Y)
    if any(a+b != 3 or c+d != 3 for (a, b, c, d), coefficient in polynomial.terms()
           if coefficient != 0):
        raise ValueError('the biform must have bidegree (3,3)')
    matrix, source_basis, target_basis = hopf_matrix()
    coefficients = matrix.inv()*S.Matrix([polynomial.coeff_monomial(v) for v in target_basis])
    return S.expand(sum(c*f for c, (j, k, f) in zip(coefficients, source_basis)))
