"""Bounded EXACT generic shell enumeration for a new parent's V3 orbit bank.

Enumerate the closed ellipsoid w G w^t <= 10 by rational LDL, then take
parity minima. No floating qfminim completeness assertion is used. LLL is
only an exactly checked integral change of basis, supplied by the caller.
"""
from fractions import Fraction as F
from math import isqrt
from exact_parity_ellipsoid import ldl


def shell_bank(gram, transform, bound=10, node_limit=20_000_000):
    n = len(gram)
    if bound < 0 or node_limit < 1:
        raise ValueError('invalid enumeration bounds')
    if len(transform) != n or any(len(row) != n for row in transform):
        raise ValueError('invalid coordinate transport')
    L, D = ldl(gram)
    z = [0]*n
    nodes, vectors = 0, 0
    best, counts = {}, {}

    def visit(i, partial):
        nonlocal nodes, vectors
        if i < 0:
            if partial.denominator != 1:
                raise ArithmeticError('nonintegral norm')
            q = int(partial)
            vectors += 1
            counts[q] = counts.get(q, 0)+1
            w = tuple(sum(z[j]*transform[j][k] for j in range(n)) for k in range(n))
            mask = sum((v % 2) << j for j, v in enumerate(w))
            if not mask:
                return
            if next(v for v in w if v) < 0:
                w = tuple(-v for v in w)
            old = best.get(mask)
            if old is None or (q, w) < old:
                best[mask] = (q, w)
            return
        centre = -sum(L[j][i]*z[j] for j in range(i+1, n))
        r2 = (F(bound)-partial)/D[i]
        if r2 < 0:
            return
        radius = isqrt(r2.numerator//r2.denominator)+1
        a, b = centre-radius, centre+radius
        lo, hi = -((-a.numerator)//a.denominator), b.numerator//b.denominator
        for value in range(lo, hi+1):
            nodes += 1
            if nodes > node_limit:
                raise RuntimeError('EXACT_SHELL_NODE_LIMIT: no complete orbit bank published')
            total = partial+D[i]*(F(value)-centre)**2
            if total <= bound:
                z[i] = value
                visit(i-1, total)

    visit(n-1, F(0))
    return {'status': 'COMPLETE_CLOSED_ELLIPSOID', 'bound': bound, 'nodes': nodes,
            'signed_vectors_including_zero': vectors, 'norm_counts': counts,
            'rows': [{'mask': mask, 'norm': q, 'word': list(w)}
                     for mask, (q, w) in sorted(best.items()) if q in (8, 10)],
            'all_represented_parity_minima': {str(k): v[0] for k, v in sorted(best.items())},
            'scope': 'Exact minima for cosets meeting this bounded ellipsoid; not a full quotient census, covering-radius proof, or assertion that norms 8/10 are optimal point-search charts.'}
