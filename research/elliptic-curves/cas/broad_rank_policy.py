"""Frozen broad-search policy. Scores and budgets NEVER imply a rank upper bound."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import math

NATIVE = ('074d9', '07ca9', '08234', '08f72', '103b2', '11952')
PARENTS = NATIVE + ('x1092-original', 'x1092-class1')
SALT = 'broad-rank-v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def key(*parts):
    return hashlib.sha256('/'.join(map(str, (SALT,) + parts)).encode()).hexdigest()


def addresses(offset, count):
    """Unique signed Calkin--Wilf addresses, O(count + log offset), not height-uniform."""
    require(type(offset) is int and type(count) is int and offset >= 0 and count > 0,
            'invalid address interval')
    require(offset + count <= 2**40, 'address interval exceeds declared integer bound')
    a = b = 1
    for bit in bin(offset // 2 + 1)[3:]:
        if bit == '0':
            b += a
        else:
            a += b
    for i in range(offset, offset + count):
        yield i, a if i % 2 == 0 else -a, b
        if i % 2:
            a, b = b, (2 * (a // b) + 1) * b - a


def rational_text(a, b):
    return str(a) if b == 1 else f'{a}/{b}'


def control_ids(offset, count, number):
    require(0 <= number < count, 'invalid control count')
    return sorted(range(offset, offset + count), key=lambda i: (key('control', i), i))[:number]


def select(rows, ranked, controls):
    """Controls chosen without scores; ranked arm excludes them. No outcome refills."""
    require(rows and 0 <= ranked and 0 <= controls and 0 < ranked + controls <= len(rows),
            'selection size exceeds population')
    require(len({r['index'] for r in rows}) == len(rows), 'duplicate addresses')
    blind = sorted(rows, key=lambda r: (key('control', r['index']), r['index']))[:controls]
    ids = {r['index'] for r in blind}
    top = sorted((r for r in rows if r['index'] not in ids),
                 key=lambda r: (-r['score_units'], r['model_bits'], r['index']))[:ranked]
    # Interleave proportionally so early stopping does not leave controls until last.
    order = [(r, 'ranked', (i + .5) / max(1, ranked)) for i, r in enumerate(top)]
    order += [(r, 'control', (i + .5) / max(1, controls)) for i, r in enumerate(blind)]
    order.sort(key=lambda item: (item[2], item[1], item[0]['index']))
    return [{**r, 'arm': arm, 'rescue': int(key('rescue', r['index'])[:8], 16) % 8 == 0}
            for r, arm, unused in order]


def primes(bound):
    require(type(bound) is int and 5 <= bound <= 4093, 'prime bound must be 5..4093')
    sieve = bytearray(b'\x01') * (bound + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, math.isqrt(bound) + 1):
        if sieve[p]:
            sieve[p*p::p] = b'\0' * ((bound - p*p)//p + 1)
    return [p for p in range(5, bound + 1) if sieve[p]]


def polynomial(coefficients, t):
    out = F(0)
    for c in reversed(coefficients):
        out = out*t + F(c)
    return out


def homogeneous(coefficients, a, b, degree):
    require(len(coefficients) <= degree + 1, 'coefficient degree exceeds chart')
    return sum(c * a**i * b**(degree-i) for i, c in enumerate(coefficients))


def integral_score_model(parent):
    """Exact constant weighted scaling for scoring only; never alter native packets."""
    A = list(map(F, parent['A_coefficients_low_to_high']))
    B = list(map(F, parent['B_coefficients_low_to_high']))
    require(A and B and len(A) <= 9 and len(B) <= 13, 'invalid short K3 chart')
    d = math.lcm(*(c.denominator for c in A+B))
    return [int(c*d**4) for c in A], [int(c*d**6) for c in B], d


def model_at(parent, parameter):
    t = F(parameter); d = t.denominator
    A = polynomial(parent['A_coefficients_low_to_high'], t) * d**8
    B = polynomial(parent['B_coefficients_low_to_high'], t) * d**12
    return list(map(str, (0, 0, 0, A, B)))


def j_invariant(model):
    # Search inputs are short models; exclusions can be general Weierstrass models.
    a1, a2, a3, a4, a6 = map(F, model)
    b2 = a1*a1+4*a2; b4 = a1*a3+2*a4; b6 = a3*a3+4*a6
    c4 = b2*b2-24*b4; c6 = -b2**3+36*b2*b4-216*b6
    delta = (c4**3-c6**2)/1728
    return str(c4**3/delta) if delta else None


def total_cap(rank):
    """Initial stage already gives every genuine seed amplification, even M18/M19."""
    if rank >= 27: return 8192
    if rank >= 25: return 2048
    if rank >= 23: return 1024
    if rank >= 20: return 512
    return 398


def stale_cap(rank):
    return 900 if rank >= 30 else 600 if rank >= 28 else 300 if rank >= 25 else 200


def next_stage(state, rescue, max_rounds):
    if state is None:
        return {'kind': 'fresh', 'allowance': 100, 'generic_allowance': 198,
                'bank_index': 0, 'revival': False}
    if state['status'] not in ('CERTIFIED', 'POLICY_EXHAUSTED') or state['rank'] >= 32:
        return None
    if state['round'] >= max_rounds:
        return None
    if state['rank'] == 17:
        if not rescue or state.get('rescue_used'):
            return None
        return {'kind': 'continuation', 'allowance': 100, 'generic_allowance': 198,
                'bank_index': 4 + state['round'], 'revival': True}
    budget = total_cap(state['rank']) - state['calls']
    if budget <= 0 or state['stale_calls'] >= stale_cap(state['rank']):
        return None
    n = min(100, budget)
    return {'kind': 'continuation', 'allowance': n, 'generic_allowance': n,
            'bank_index': state['round'] + 1, 'revival': False}


def score_tables(parent, bound):
    """Point counts in every P1(Fp) residue; audited independently by Legendre sums."""
    import numpy as np
    A, B, scaling = integral_score_model(parent)
    tables = []
    def valuation(n, p):
        if not n: return 10**6
        v = 0
        while n % p == 0:
            n //= p; v += 1
        return v
    for p in primes(bound):
        k = min(min(valuation(c, p) for c in A)//4,
                min(valuation(c, p) for c in B)//6)
        aa = [c//p**(4*k) % p for c in A]
        bb = [c//p**(6*k) % p for c in B]
        def ev(poly, r):
            out = 0
            for c in reversed(poly): out = (out*r+c) % p
            return out
        av = np.array([ev(aa, r) for r in range(p)] + [(aa+[0]*9)[8]], dtype=np.int64)
        bv = np.array([ev(bb, r) for r in range(p)] + [(bb+[0]*13)[12]], dtype=np.int64)
        chi = -np.ones(p, dtype=np.int64); chi[0] = 0
        for x in range(1, p): chi[x*x % p] = 1
        x = np.arange(p, dtype=np.int64)
        traces = -chi[(x[None, :]**3+av[:, None]*x[None, :]+bv[:, None]) % p].sum(axis=1)
        smooth = (4*av**3+27*bv**2) % p != 0
        checks = []
        for r in sorted({0, p, int(key('table', p)[:8], 16) % p}):
            direct = -sum(0 if (v := (z**3+int(av[r])*z+int(bv[r])) % p) == 0
                          else (1 if pow(v, (p-1)//2, p) == 1 else -1) for z in range(p))
            require(direct == int(traces[r]), 'point-count implementations disagree')
            checks.append([r, direct])
        require(all(int(t)**2 <= 4*p for t, good in zip(traces, smooth) if good), 'Hasse bound failed')
        units = [round(1e12*(2-int(t))*math.log(p)/(p+1-int(t))) if good else 0
                 for t, good in zip(traces, smooth)]
        tables.append({'prime': p, 'constant_p_scaling': k, 'traces': list(map(int, traces)),
                       'smooth': list(map(bool, smooth)), 'score_units': units, 'scalar_checks': checks})
    return {'constant_denominator_scale': str(scaling), 'tables': tables,
            'boundary': 'Coefficientwise minimization only. Singular local terms omitted, not rank exclusions.'}


def feature_rows(parent, tables, offset, count, chunk_size=4096):
    """All traces reconstructible from table + rational parameter; no duplicated vectors."""
    import numpy as np
    A, B, unused = integral_score_model(parent)
    stream = iter(addresses(offset, count))
    while True:
        import itertools
        chunk = list(itertools.islice(stream, chunk_size))
        if not chunk: break
        nums = np.array([a for i,a,b in chunk], dtype=np.int64)
        dens = np.array([b for i,a,b in chunk], dtype=np.int64)
        totals = np.zeros(len(chunk), dtype=np.int64)
        bad = [[] for unused in chunk]; bands = {}
        for tab in tables['tables']:
            p = tab['prime']; inverses = np.array([0]+[pow(i,-1,p) for i in range(1,p)], dtype=np.int64)
            residues = (nums % p) * inverses[dens % p] % p
            residues[dens % p == 0] = p
            totals += np.asarray(tab['score_units'], dtype=np.int64)[residues]
            for i in np.flatnonzero(~np.asarray(tab['smooth'], dtype=bool)[residues]): bad[int(i)].append(p)
            if p in (97, 257, 997): bands[str(p)] = totals.copy()
        for k, (i,a,b) in enumerate(chunk):
            av, bv = homogeneous(A,a,b,8), homogeneous(B,a,b,12)
            delta = -16*(4*av**3+27*bv**2)
            yield {'index': i, 'parameter': rational_text(a,b), 'score_units': int(totals[k]),
                   'scores_by_cutoff': {p:int(v[k]) for p,v in bands.items()},
                   'parameter_height': max(abs(a), b), 'model_bits': max(abs(av).bit_length(),abs(bv).bit_length()),
                   'discriminant_bits': abs(delta).bit_length(), 'nonsingular': bool(delta),
                   'local_singular_primes': bad[k], 'smooth_prime_count': len(tables['tables'])-len(bad[k])}
