#!/usr/bin/env python3
"""After V3 seals, test whether its newly found subgroup lifts the fixed covers.

Character matching proposes a combination; only an exact field square root
and independent rational group-law/cover replay establish recovery.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
import verify_constructed_class_blind as exact

ROOT = exact.ROOT
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
OUT = exact.EXPERIMENT
RUN = ROOT/'artifacts/local/elliptic-curves/constructed-class-blind-v3-v2'


def read(path):
    return json.loads(path.read_text())


def write(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, sort_keys=True, indent=2)
        stream.write('\n')


def produce(output):
    from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, pari, prime_range, lcm, gcd
    started = time.monotonic()
    paths = [OUT/'covers.json', RUN/'terminal.json', RUN/'verified-v2.json',
             ROOT/'artifacts/generated-results/elliptic-curves/rank_jump_reference_strict_class_construction_inputs_v1.json']
    covers, terminal, replay, ref = [read(p) for p in paths]
    assert replay['terminal_sha256'] == exact.sha(RUN/'terminal.json')
    assert replay['status'] == 'PASS_EXACT_POINT_MAP_AND_FINITE_RANK_REPLAY'
    # No oracle or cover-worker point file is opened by this evaluator.
    R = PolynomialRing(QQ, 'z')
    f = R(covers['cubic_ascending'])
    pari.allocatemem(64000000, 268435456, silent=True)
    nf = pari.nfinit([pari(f), ref['S_finite']])
    E = EllipticCurve(list(map(QQ, terminal['curve'])))
    pts = [E(list(map(QQ, p))) for p in terminal['points']]
    n = len(pts)
    polys = [(p[0]-3)/9-R.gen() for p in pts] + [R(c['beta_ascending']) for c in covers['cases']]
    chars, width, blocks = [0]*len(polys), 0, []
    disc = f.discriminant()
    for p in prime_range(3, 5000):
        if disc % p == 0 or any(c.denominator() % p == 0 for poly in polys for c in poly):
            continue
        roots = f.change_ring(GF(p)).roots(multiplicities=False)
        if len(roots) != 3:
            continue
        vals = [[int(poly.change_ring(GF(p))(r)) for r in roots] for poly in polys]
        if any(v == 0 for row in vals for v in row):
            continue
        for i, row in enumerate(vals):
            for j, v in enumerate(row):
                chars[i] |= int(pow(v, (int(p)-1)//2, int(p)) == p-1) << (width+j)
        blocks.append({'prime': int(p), 'roots': list(map(int, roots))})
        width += 3
    A = matrix(GF(2), [[(c >> j) & 1 for j in range(width)] for c in chars[:n]])
    assert A.rank() == n
    sqrt = pari('(nf,b)->{my(y); if(nfeltissquare(nf,b,&y),y,0)}')
    rows = []
    for j, c in enumerate(covers['cases']):
        v = vector(GF(2), [(chars[n+j] >> k) & 1 for k in range(width)])
        try:
            word = list(map(int, A.transpose().solve_right(v)))
        except ValueError:
            rows.append({'column': c['column'], 'status': 'NO_FINITE_CHARACTER_MATCH'})
            continue
        point = sum((k*p for k, p in zip(word, pts)), E(0))
        assert not point.is_zero()
        X, Y = (point[0]-3)/9, point[1]/27
        beta = pari.Mod(pari(R(c['beta_ascending'])), pari(f))
        candidate = pari.Mod(pari(X-R.gen()), pari(f))/beta
        root = sqrt(nf, candidate)
        if root == 0:
            rows.append({'column': c['column'], 'status': 'FINITE_MATCH_NOT_AN_EXACT_SQUARE'})
            continue
        xi = pari.nfbasistoalg(nf, root)
        assert beta*xi**2 == pari.Mod(pari(X-R.gen()), pari(f))
        ns = QQ(c['positive_norm_square_root'])
        if ns*QQ(pari.nfeltnorm(nf, xi)) == -Y:
            xi = -xi
        assert ns*QQ(pari.nfeltnorm(nf, xi)) == Y
        coordinates = [QQ(pari.lift(xi).polcoef(i)) for i in range(3)] + [QQ(1)]
        den = lcm([x.denominator() for x in coordinates])
        integral = [ZZ(den*x) for x in coordinates]
        content = gcd(integral)
        integral = [x//content for x in integral]
        prefix = max(i+1 for i, bit in enumerate(word) if bit)
        rows.append({'column': c['column'], 'status': 'EXACT_LIFT_FROM_FRESH_V3_SUBGROUP',
            'combination_word': word, 'earliest_certified_prefix_rank': prefix,
            'primitive_cover_point': list(map(str, integral)),
            'point_on_short_model': list(map(str, point.xy())),
            'point_on_original_model': [str(X/4), str((Y-X)/8)]})
        print('V3_FIXED_COVER', c['column'], 'prefix', prefix, flush=True)
    write(output, {'schema': 'constructed-class-fresh-v3-lifts.v1',
        'status': 'PASS' if len(rows) == 2 and all(r['status'] == 'EXACT_LIFT_FROM_FRESH_V3_SUBGROUP' for r in rows) else 'UNKNOWN',
        'cases': rows, 'finite_character_rank': n, 'character_blocks': blocks,
        'wall_seconds': time.monotonic()-started,
        'bindings': {str(p.relative_to(ROOT)): exact.sha(p) for p in paths+[Path(__file__).resolve()]},
        'boundary': 'Fixed-cover evaluation after the generic-only V3 endpoint sealed. Uses only points just found by this V3 run, no previous point oracle or blind-cover-worker coordinates. Character matches alone do not prove the result.'})


def verify(source, output):
    from half_lattice_pointed_sieve import linear_combination
    started = time.monotonic()
    d, covers, run = read(source), read(OUT/'covers.json'), read(RUN/'terminal.json')
    for p, h in d['bindings'].items():
        assert exact.sha(ROOT/p) == h
    f = list(map(F, covers['cubic_ascending']))
    model = tuple(map(F, run['curve']))
    points = tuple(tuple(map(F, p)) for p in run['points'])
    by_column = {c['column']: c for c in covers['cases']}
    assert d['status'] == 'PASS' and len(d['cases']) == 2
    rows = []
    for r in d['cases']:
        c = by_column[r['column']]
        P = linear_combination(model, points, r['combination_word'])
        assert list(map(str, P)) == r['point_on_short_model']
        p = list(map(F, r['primitive_cover_point']))
        assert p[3] and all(v.denominator == 1 for v in p)
        from math import gcd
        assert gcd(*(int(v) for v in p)) == 1
        beta = list(map(F, c['beta_ascending']))
        xi = [v/p[3] for v in p[:3]]
        X, Y = (P[0]-3)/9, P[1]/27
        assert exact.mul(f, beta, exact.mul(f, xi, xi)) == [X, F(-1), F(0)]
        assert F(c['positive_norm_square_root'])*exact.norm(f, xi) == Y
        q0, q1, q2 = exact.quadrics(f, beta)
        assert exact.evaluate(q1, p) == exact.evaluate(q2, p) == 0
        assert exact.evaluate(q0, p) == X*p[3]**2
        assert Y*Y == sum(coef*X**i for i, coef in enumerate(f))
        x, y = X/4, (Y-X)/8
        assert list(map(str, [x, y])) == r['point_on_original_model']
        assert y*y+x*y == x**3+F(covers['original_curve'][3])*x+F(covers['original_curve'][4])
        rows.append({'column': r['column'], 'status': 'PASS_EXACT_FRESH_SUBGROUP_LIFT',
                     'earliest_certified_prefix_rank': r['earliest_certified_prefix_rank']})
    write(output, {'status': 'PASS', 'cases': rows, 'wall_seconds': time.monotonic()-started,
        'source_sha256': exact.sha(source), 'replayer_sha256': exact.sha(__file__),
        'boundary': 'Independent Fraction arithmetic verifies the combination of freshly discovered V3 points, source cover, field square identity and original curve equation.'})
    print('PASS both fixed covers lifted from freshly recovered V3 subgroup')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['produce', 'verify'])
    p.add_argument('--source', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    produce(a.output.resolve()) if a.mode == 'produce' else verify(a.source.resolve(), a.output.resolve())
