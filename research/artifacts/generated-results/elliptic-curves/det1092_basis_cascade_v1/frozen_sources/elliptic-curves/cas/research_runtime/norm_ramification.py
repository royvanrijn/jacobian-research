"""Necessary unramifiedness filter for cubic norm-projection dictionaries.

This filters coefficient words, not elliptic curves. Unremoved generators
remain UNKNOWN. No factorization, class-group data, or point oracle is used.
"""
from fractions import Fraction as Q
from math import gcd, isqrt, lcm, prod


def rational(value):
    if type(value) not in (int, str, Q):
        raise ValueError('exact integer or rational string required')
    return Q(value)


def cubic_data(coefficients, elements):
    f = list(map(rational, coefficients))
    if len(f) != 4 or f[3] != 1:
        raise ValueError('monic cubic required')
    c, b, a, _ = f
    disc = a*a*b*b - 4*b**3 - 4*a**3*c - 27*c*c + 18*a*b*c
    if not disc:
        raise ValueError('separable cubic required')
    forbidden = 2 * abs(disc.numerator) * disc.denominator
    forbidden *= prod(x.denominator for x in f)
    norms = []
    for element in elements:
        v = list(map(rational, element))
        if len(v) != 3 or not any(v):
            raise ValueError('nonzero degree-less-than-three element required')
        den = lcm(*(x.denominator for x in v))
        content = gcd(*(int(x*den) for x in v))
        forbidden *= den * abs(content)
        # Multiplication by alpha on the power basis 1, theta, theta^2.
        columns = []
        for j in range(3):
            w = [Q(0)] * j + v + [Q(0)] * (2-j)
            for k in range(4, 2, -1):
                for i in range(3):
                    w[k-3+i] -= w[k]*f[i]
            columns.append(w[:3])
        x, y, z = columns
        n = (x[0]*(y[1]*z[2]-y[2]*z[1])
             - y[0]*(x[1]*z[2]-x[2]*z[1])
             + z[0]*(x[1]*y[2]-x[2]*y[1]))
        if not n:
            raise ValueError('element must be a unit in the cubic algebra')
        norms.append(n)
        forbidden *= n.denominator
    return norms, forbidden


def isolated_remainder(norm, forbidden, other_norms):
    """Remove full prime support, including powers beyond the supplied gcd."""
    n = abs(norm.numerator)
    for support in [forbidden, *(abs(x.numerator) for x in other_norms)]:
        while (g := gcd(n, support)) > 1:
            n //= g
    return n


def preflight(coefficients, elements):
    """Prove forced-zero coordinates in EVERY unramified coefficient word.

    If i has a nonsquare norm remainder coprime to all other active norms
    and forbidden support, some good unramified prime has odd norm valuation
    only at i. Its polynomial is nonzero modulo that prime and has degree
    <3, so at least one residue component is a unit. N(alpha_i)*alpha_i
    is odd there. Thus every unramified word must have coefficient i=0.
    Remove all such coordinates and repeat. This preserves cancellations
    between generators sharing prime support, including duplicate classes.
    """
    norms, forbidden = cubic_data(coefficients, elements)
    active = list(range(len(norms)))
    rounds = []
    while active:
        rejected = []
        for i in active:
            remainder = isolated_remainder(norms[i], forbidden,
                                            [norms[j] for j in active if j != i])
            root = isqrt(remainder)
            if root*root != remainder:
                rejected.append({'index': i, 'remainder': str(remainder),
                                 'floor_square_root': str(root)})
        if not rejected:
            break
        rounds.append(rejected)
        removed = {x['index'] for x in rejected}
        active = [i for i in active if i not in removed]
    return {'status': 'PROVED_NECESSARY_COEFFICIENT_RESTRICTIONS',
            'norms': list(map(str, norms)), 'rounds': rounds,
            'input_generator_count': len(norms),
            'forced_zero_count': len(norms)-len(active),
            'unresolved_indices': active,
            'unramified_coefficient_dimension_upper_bound': len(active),
            'additional_strict_classes': 'UNKNOWN' if active else 0,
            'whole_curve_rank_decision': 'UNKNOWN'}
