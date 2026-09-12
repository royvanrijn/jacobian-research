"""Small exact helpers for a bounded, equation-only mod-two relation pilot."""
from math import comb, gcd


def insert(pivots, row):
    while row:
        k = row.bit_length()-1
        if k not in pivots:
            pivots[k] = row
            return True
        row ^= pivots[k]
    return False


def matrix_stats(rows, columns, canonical_count=0):
    pivots, seen, counts = {}, set(), [0]*columns
    canonical_rank = None
    for i, row in enumerate(rows):
        if row < 0 or row >> columns:
            raise ValueError('parity row exceeds factor base')
        if i == canonical_count:
            canonical_rank = len(pivots)
        insert(pivots, row)
        seen.add(row)
        v = row
        while v:
            bit = v & -v
            counts[bit.bit_length()-1] += 1
            v ^= bit
    if canonical_rank is None:
        canonical_rank = len(pivots)
    return {'columns': columns, 'rows': len(rows), 'distinct_rows': len(seen),
            'rank_mod2': len(pivots), 'deficiency': columns-len(pivots),
            'canonical_rank': canonical_rank, 'rank_gain_beyond_canonical': len(pivots)-canonical_rank,
            'zero_columns': counts.count(0), 'singleton_columns': counts.count(1),
            'zero_rows': rows.count(0), 'deficiency_fraction': (columns-len(pivots))/columns if columns else None}


def primitive_pairs(count):
    emitted, height = 0, 0
    while emitted < count:
        height += 1
        for b in range(1, height+1):
            for a in (range(-height, height+1) if b == height else (-height, height)):
                if gcd(a, b) != 1:
                    continue
                yield a, b
                emitted += 1
                if emitted == count:
                    return


def form_value(c, a, b):
    """F(a,-b) for ascending coefficients of F(X,Y)."""
    return c[3]*a**3-c[2]*a*a*b+c[1]*a*b*b-c[0]*b**3


def transform(c, matrix):
    p, q, r, s = matrix
    if abs(p*s-q*r) != 1:
        raise ValueError('nonunimodular transformation')
    out = [0]*4
    for i, value in enumerate(c):
        for j in range(i+1):
            for k in range(4-i):
                out[j+k] += value*comb(i, j)*p**j*q**(i-j)*comb(3-i, k)*r**k*s**(3-i-k)
    return out


def compose(a, b):
    p, q, r, s = a
    u, v, w, z = b
    return [p*u+q*w, p*v+q*z, r*u+s*w, r*v+s*z]


def reduce_form(c, max_steps=128):
    """Fixed finite exact height descent, deliberately NOT Julia reduction."""
    def normalized(v):
        return [-x for x in v] if v[-1] < 0 else v
    def score(v):
        return max(map(abs, v)), sum(map(abs, v)), tuple(map(abs, v)), tuple(v)
    c, total = normalized(list(c)), [1, 0, 0, 1]
    for step in range(max_steps):
        ks = {-1, 1}
        for numerator, denominator in ((-c[2], 3*c[3]), (-c[1], 3*c[0])):
            if denominator:
                k = numerator//denominator
                ks.update([k, k+1])
        moves = [[0, -1, 1, 0], [1, 0, 0, -1]]
        moves += [[1, k, 0, 1] for k in sorted(ks) if k]
        moves += [[1, 0, k, 1] for k in sorted(ks) if k]
        choices = [(normalized(transform(c, m)), m) for m in moves]
        best, move = min(choices, key=lambda item: (score(item[0]), item[1]))
        if score(best) >= score(c):
            return c, total, step, False
        c, total = best, compose(total, move)
    return c, total, max_steps, True


def strip_support(value, primorial):
    """Remove all powers supported on a fixed proven-prime product."""
    value = abs(value)
    if not value:
        raise ValueError('zero norm')
    while True:
        common = gcd(value, primorial)
        if common == 1:
            return value
        value //= common
