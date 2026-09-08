#!/usr/bin/env sage-python
"""Bounded exact rational-span diagnosis of the unresolved split fibre.

Numerical heights propose words only. Exact rational group identities decide.
No search or proposed word is sent to the frozen population or V3 workers.
"""
import argparse
from fractions import Fraction
from pathlib import Path

from sage.all import QQ, ZZ, RealField, EllipticCurve, pari, matrix, vector
import det1092_funnel as f
from v3_warm_support import atomic, read, require, sha


def audit(seed, output):
    require(read(seed/'result.json')['status'] == 'CONIC_POINT_INDEPENDENCE_UNRESOLVED',
            'diagnostic requires an unresolved conic seed')
    paths = [seed/'m17.json', seed/'conic-points.json', seed/'result.json', Path(__file__)]
    protocol = dict(schema='det1092-split-seed-span-audit.v1', wall_seconds=120,
        rss_bytes=1024**3, precision_bits=384, denominator_bound=1000000,
        proposal_tolerance='2^-160',
        inputs={str(p.relative_to(f.ROOT)):sha(p) for p in paths},
        scope='Two retained conic witnesses only. Exact multiple identities prove rational-span membership. Failed numerical proposals remain UNKNOWN; no whole-curve upper bound.')
    atomic(output/'protocol.json', protocol, immutable=True)
    m17 = read(seed/'m17.json')
    E = EllipticCurve(QQ, list(map(QQ, m17['curve'])))
    base = [E(list(map(QQ, p))) for p in m17['points']]
    candidates = [E(list(map(QQ, p))) for p in read(seed/'conic-points.json')]
    require(len(base) == 17 and len(candidates) == 2, 'expected17 inherited points and two conic witnesses')
    RF = RealField(384)
    H = pari(E).ellheightmatrix([list(P.xy()) for P in base+candidates], precision=384)
    G = matrix(RF, 17, 17, [RF(str(H[i, j])) for i in range(17) for j in range(17)])
    results = []
    for k, P in enumerate(candidates):
        coefficients = G.solve_right(vector(RF, [RF(str(H[i, 17+k])) for i in range(17)]))
        word = [QQ(str(Fraction(str(v)).limit_denominator(1000000))) for v in coefficients]
        denominator = ZZ(1)
        for v in word:
            denominator = denominator.lcm(v.denominator())
        precise = all(abs(v-RF(w)) < RF(2)**-160 for v, w in zip(coefficients, word))
        exact = precise and sum((ZZ(denominator*w)*Q for w, Q in zip(word, base)), E(0)) == denominator*P
        results.append(dict(index=k, point=list(map(str, P.xy())),
            status='EXACT_RATIONAL_SPAN_MEMBERSHIP' if exact else 'UNKNOWN',
            word=list(map(str, word)) if exact else None,
            clearing_denominator=int(denominator) if exact else None))
    result = dict(status='COMPLETE_BOUNDED_SPAN_AUDIT', points=results,
        protocol_sha256=sha(output/'protocol.json'),
        rank_gain_of_these_conic_points=0 if all(r['status']=='EXACT_RATIONAL_SPAN_MEMBERSHIP' for r in results) else 'UNKNOWN',
        boundary='Only the rational span of the17 retained generic points and these two conic witnesses. No rank upper bound on the elliptic curve.')
    atomic(output/'result.json', result, immutable=True)
    print(result, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    audit(a.seed.resolve(), a.output.resolve())
