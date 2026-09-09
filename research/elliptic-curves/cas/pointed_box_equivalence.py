"""Exact projective-box key for rational height H=max(|a|,|b|).

Use only for completed pointed boxes which include infinity. The four
projective signed permutations t, -t, 1/t, -1/t preserve this height exactly.
This tests that sufficient equivalence, not every possible box coincidence.
"""
from fractions import Fraction as F


def box_key(matrix):
    a, b, c, d = map(F, matrix)
    if a*d == b*c:
        raise ArithmeticError('singular coordinate matrix')
    candidates = ((a,b,c,d), (-a,b,-c,d), (b,a,d,c), (-b,a,-d,c))
    def projective(v):
        pivot = next(x for x in v if x)
        return tuple(x/pivot for x in v)
    return min(map(projective, candidates))
