#!/usr/bin/env python3
"""Exact rank >=26, global minimality and pinned equation comparison, without Sage."""
import argparse
import json
from math import gcd
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = ART/'compact_r17_wide_results_v1.json'
OUTPUT = ART/'new_compact_rank26_proof_v1.json'


def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
        Path(__file__).resolve(), Path(cert.__file__).resolve(),
        ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',
        ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')}


def verify(data):
    if data['sources'] != sources():
        raise ArithmeticError('proof sources differ')
    model = tuple(map(cert.F, data['minimal_curve']))
    if any(a.denominator != 1 for a in model):
        raise ArithmeticError('model is not integral')
    inv = cert.weierstrass_invariants(model)
    if not inv['discriminant']:
        raise ArithmeticError('singular model')
    c4, c6 = inv['c4'], inv['c6']
    if any(a.denominator != 1 for a in (c4, c6)) or gcd(int(c4), int(c6)) != 27:
        raise ArithmeticError('global minimality criterion changed')
    # Any integral nonminimal model at p has p^4|c4 and p^6|c6.
    # Here gcd(c4,c6)=3^3, which is divisible by no fourth prime power.
    short = tuple(map(cert.F, data['discovery_curve']))
    b2 = inv['b2']; a1, _, a3, _, _ = model
    if short != (cert.F(0), cert.F(0), cert.F(0), -c4/48, -c6/864):
        raise ArithmeticError('short-model transport differs')
    points = [tuple(map(cert.F, p)) for p in data['points']]
    transported = [(x+b2/12, y+(a1*x+a3)/2) for x, y in points]
    if transported != [tuple(map(cert.F, p)) for p in data['discovery_points']]:
        raise ArithmeticError('point transport differs')
    proof = data['rank_certificate']
    actual = cert.checked_rank(model, points, [s['prime'] for s in proof['signatures']],
                               proof['no_rational_2_torsion_prime'])
    if len(points) != 26 or json.dumps(actual, sort_keys=True) != json.dumps(proof, sort_keys=True):
        raise ArithmeticError('independent-point proof differs')
    matches = [r['id'] for r in data['catalogue']['equations'] if cert.isomorphic(model, r['ainvs'])]
    previous = [r['address'] for r in data['previous_equations'] if cert.isomorphic(model, r['curve'])]
    if matches or previous or matches != data['icarm_matches'] or previous != data['previous_matches']:
        raise ArithmeticError('new-curve comparison failed')
    print('PASS: 26 independent points; global minimality; no match among',
          len(data['catalogue']['equations']), 'catalogue and', len(data['previous_equations']), 'earlier equations', flush=True)


def build(output):
    if output.exists():
        raise FileExistsError('preserve rank26 proof')
    source = cert.read(INPUT)
    row = next(r for r in source['curves'] if r['family'] == '07ca9' and r['parameter'] == '-2507/3068')
    short = tuple(map(cert.F, row['curve']))
    inv = cert.weierstrass_invariants(short)
    a1, a2, a3 = 1, -1, 1
    b2 = cert.F(a1*a1+4*a2)
    b4 = (b2*b2-inv['c4'])/24
    b6 = (-b2**3+36*b2*b4-inv['c6'])/216
    model = (cert.F(a1), cert.F(a2), cert.F(a3), (b4-a1*a3)/2, (b6-a3*a3)/4)
    points = []
    for x, y in (tuple(map(cert.F, p)) for p in row['points']):
        X = x-b2/12
        points.append((X, y-(a1*X+a3)/2))
    old = row['rank_certificate']
    proof = cert.checked_rank(model, points, [s['prime'] for s in old['signatures']], old['no_rational_2_torsion_prime'])
    data = {'schema': 'elliptic-curves.new-compact-rank26-proof.v1', 'sources': sources(),
        'input': {'path': str(INPUT.relative_to(ROOT)), 'sha256': cert.hashed(INPUT)},
        'family': row['family'], 'parameter': row['parameter'], 'minimal_curve': list(map(str, model)),
        'points': [list(map(str, p)) for p in points], 'rank_certificate': proof,
        'discovery_curve': row['curve'], 'discovery_points': row['points'],
        'transport': 'x_short = X - 1/4; y_short = Y + (X+1)/2',
        'minimality': 'Integral nonsingular model with gcd(c4,c6)=27; nonminimality at any prime would require that prime to the fourth power dividing this gcd.',
        'catalogue': source['catalogue'], 'previous_equations': source['previous_equations'],
        'icarm_matches': [], 'previous_matches': [],
        'claim_boundary': 'Rank at least26 and a global minimal integral model. Exact rank, conductor, and universal novelty remain unknown. No new rank28 or rank32 curve is claimed.'}
    verify(data)
    cert.write(output, data)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=OUTPUT)
    p.add_argument('--check', type=Path)
    a = p.parse_args()
    verify(cert.read(a.check)) if a.check else build(a.output)
