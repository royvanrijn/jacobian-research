"""Exact finite covers of the normalized non-null anti-halfturn locus.

The scheme-length argument is in ANTI_HALFTURN_SCHEME_LENGTH_AND_A14_GATE.md.
This module verifies coordinate identities without computing a new Groebner
basis or the actual filtered Rodrigues multiplication.
"""
from functools import lru_cache

import sympy as S

from cyclic_frame_quotient import p, r, s0, s1, s2, U
from conformal_euler_cubic import z, t, w

a, c, lam, k, d = S.symbols('a c lam k d')
D0 = S.Rational(1, 8)
FREE_TRIPLES = ((0, 0, 0), (0, 0, 2), (0, 2, 2), (2, 2, 2))
FREE_CUBIC = tuple(U[indices] for indices in FREE_TRIPLES)
ODD_TRIPLES = ((0, 0, 1), (0, 1, 2), (1, 1, 1), (1, 2, 2))


def normalized_locus():
    """Full six-coordinate Laurent presentation, prior to moment equations."""
    p_value = (c-D0**2/c**2)/a
    result = {p: p_value, r: 0, s0: a, s1: 0, s2: c}
    result.update({U[indices]: 0 for indices in ODD_TRIPLES})
    result[U[0, 1, 1]] = -a/c*U[0, 0, 0]-2*U[0, 0, 2]-p_value*U[0, 2, 2]
    result[U[1, 1, 2]] = -a/c*U[0, 0, 2]-2*U[0, 2, 2]-p_value*U[2, 2, 2]
    return result


def chart_base():
    """Scaled chart data before fixing its oriented determinant."""
    return {a: 2*lam**2, c: lam**4*k**2,
            p: lam**2*(k**2-d**2/4)/2,
            'bz': k+d/2, 'bw': k-d/2,
            'chart_D': d*k**2/2,
            'normalized_D': lam**6*d*k**2/2}


def inverse_base():
    difference = 4*D0/(a*c)
    return {'lambda_square': a/2, 'k_square': 4*c/a**2,
            'difference': difference,
            'bz': k+difference/2, 'bw': k-difference/2,
            'eta': 1/(2*k*difference)}


@lru_cache(maxsize=None)
def cubic_change():
    """Map the four quotient cubic coefficients to az,aw,ez,ew.

    The map is on the common amplitude/axis cover; its determinant is
    9*d^6*k^6*lam^20/16, a Laurent unit over Z_(1073741827).
    """
    base = chart_base()
    values = {u: 0 for u in U.values()}
    values.update({u: u for u in FREE_CUBIC})
    values[U[0, 1, 1]] = (-base[a]/base[c]*U[0, 0, 0]
                          -2*U[0, 0, 2]-base[p]*U[0, 2, 2])
    values[U[1, 1, 2]] = (-base[a]/base[c]*U[0, 0, 2]
                          -2*U[0, 2, 2]-base[p]*U[2, 2, 2])
    y = (lam*(z+w), lam**2*k*t,
         lam**3*k*((k+d/2)*z+(k-d/2)*w)/2)
    original = S.expand(sum(values[U[tuple(sorted((i, j, h)))]]*y[i]*y[j]*y[h]
                            for i in range(3) for j in range(3) for h in range(3))/lam)
    poly = S.Poly(original, z, t, w)
    coefficients = S.Matrix([poly.coeff_monomial(monomial)
                             for monomial in (z**3, w**3, z*t*t, w*t*t)])
    transform = coefficients.jacobian(FREE_CUBIC)
    return transform, original, y
