"""Exact denominator/Gauss chart, then PARI reduction without global minimization.

This defines a different finite coordinate box from the historical minimal-model
pipeline. It is an experimental mapping policy, not a coverage equivalence.
"""
from fractions import Fraction as F
from math import isqrt, lcm
from sage.all import pari
from half_lattice_pointed_sieve import linear_combination, make_chart
from pointed_quartic_search import normalize
from search_observability import transform, multiply


def mapping(model, points, centre, progress=lambda stage: None):
    progress('group_sum')
    point = linear_combination(model, points, centre['representative'])
    if point is None:
        raise ArithmeticError('finite nonzero centre required')
    x, y = point
    raw = (-3*x*x-4*model[3], -8*y, -6*x, F(0), F(1))
    progress('integral_gauss_chart')
    chart = make_chart(model, point)
    d, scale, k = chart.denominator, chart.curve_scale, chart.shift
    first = multiply((F(d)/scale, F(k)/d/scale, F(0), F(1)), chart.matrix)
    coefficients, _ = normalize(chart.coefficients)
    progress('hyperellred_only')
    polynomial = '+'.join(f'({v})*x^{i}' for i, v in enumerate(coefficients))
    result = pari('my(m,C);C=hyperellred(['+polynomial+',0],&m);[C,m]')
    second = tuple(F(str(result[1][1][i,j])) for i in range(2) for j in range(2))
    matrix = multiply(first, second)
    P = [F(str(result[0][0].polcoef(i))) for i in range(5)]
    Q = [F(str(result[0][1].polcoef(i))) for i in range(3)]
    disc = [4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0 <= i-j < 3) for i in range(5)]
    transformed = transform(raw, matrix)
    index = next(i for i, value in enumerate(disc) if value)
    ratio = transformed[index]/disc[index]
    if (ratio <= 0 or tuple(transformed) != tuple(ratio*v for v in disc)
        or isqrt(ratio.numerator)**2 != ratio.numerator
        or isqrt(ratio.denominator)**2 != ratio.denominator):
        raise ArithmeticError('exact factor-free quartic map identity failed')
    progress('exact_map_verified')
    return {'centre': centre, 'raw_coefficients': list(map(str,raw)),
            'denominator_clearing': str(lcm(*(v.denominator for v in raw))),
            'first_matrix': list(map(str,first)), 'second_matrix': list(map(str,second)),
            'matrix': list(map(str,matrix)), 'reduced_P': list(map(str,P)),
            'reduced_Q': list(map(str,Q)), 'discriminant_quartic': list(map(str,disc)),
            'square_ratio': str(ratio),
            'coordinate_policy': {'kind':'raw','matrix':list(map(str,matrix))},
            'preparation_policy': 'integral-denominator-gauss-then-hyperellred-v1',
            'hyperellminimalmodel_called': False}
