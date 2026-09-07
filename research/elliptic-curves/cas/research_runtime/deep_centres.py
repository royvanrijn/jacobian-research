"""Finite deep-centre selection, separate from nearest-first enumeration.

Proposed representatives give upper bounds on all coset minima. Exact LDL
branch-and-bound certifies the deepest stratum: every proposed maximum must
have no shorter representative. Other upper bounds already exclude those
classes from that stratum. This avoids calling a floating CVP a proof.
"""
from fractions import Fraction
from math import isqrt

from .cvp import Hole, ldl


def norm(gram, vector):
    return sum(x*sum(g*y for g, y in zip(row, vector)) for x, row in zip(vector, gram))


def parity(vector):
    return sum((int(v) % 2) << i for i, v in enumerate(vector))


def exact_coset_minimum(gram, mask, representative, *, node_budget=2000000):
    """Exact CVP in a fixed parity, with a supplied upper bound and node cap.

    Searches strictly below the current integer norm. All pruning uses rational
    LDL costs and integer square roots; a budget stop raises, never certifies.
    """
    gram, lower, diagonal = ldl(gram)
    n = len(gram)
    if any(x.denominator != 1 for r in gram for x in r):
        raise ValueError('integral Gram required for strict integer-norm bounds')
    if not 0 <= mask < 1 << n or len(representative) != n or parity(representative) != mask:
        raise ValueError('invalid parity witness')
    if any(type(x) is not int for x in representative): raise ValueError('integral representative required')
    best = int(norm(gram, representative)); answer = tuple(representative)
    work = [0]*n; nodes = 0
    def visit(i, used):
        nonlocal best, answer, nodes
        nodes += 1
        if nodes > node_budget: raise TimeoutError('exact parity CVP node budget reached')
        if i < 0:
            if used.denominator != 1: raise ArithmeticError('nonintegral leaf norm')
            best, answer = int(used), tuple(work)
            return
        remaining = Fraction(best-1)-used
        if remaining < 0: return
        shift = sum(lower[j][i]*work[j] for j in range(i+1, n))
        # |den*z+num|^2 <= remaining*den^2/diagonal[i].
        radius2 = remaining*shift.denominator**2/diagonal[i]
        radius = isqrt(radius2.numerator//radius2.denominator)
        lo = (-radius-shift.numerator + shift.denominator-1)//shift.denominator
        hi = (radius-shift.numerator)//shift.denominator
        bit = (mask >> i) & 1
        lo += (bit-lo) % 2
        values = sorted(range(lo, hi+1, 2), key=lambda z: (abs(z+shift), z))
        for z in values:
            cost = used + diagonal[i]*(z+shift)**2
            if cost >= best: continue
            work[i] = z; visit(i-1, cost)
    visit(n-1, Fraction(0))
    if int(norm(gram, answer)) != best or parity(answer) != mask:
        raise ArithmeticError('exact CVP witness failed')
    return best, answer, nodes


def diverse_deep(rows, distances, count):
    """Farthest-first in the flat torus, with depth then mask tie breaks.

    rows are already certified deep classes. distances(a xor b) is the
    *coset minimum*, so geometry is independent of shortest-vector sign ties.
    """
    pool = {r.mask: r for r in rows}
    if len(pool) != len(rows) or not 0 <= count <= len(pool):
        raise ValueError('invalid diverse pool or exposure')
    result = []
    separation = {}
    for _ in range(count):
        if not result:
            key = max(pool, key=lambda m: (pool[m].squared_distance, -m))
        else:
            last = result[-1].mask
            for m in pool:
                separation[m] = min(separation.get(m, Fraction(10**100)), distances(m ^ last))
            key = max(pool, key=lambda m: (separation[m], pool[m].squared_distance, -m))
        result.append(pool.pop(key))
    return result
