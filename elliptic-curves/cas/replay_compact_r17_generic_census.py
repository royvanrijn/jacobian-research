#!/usr/bin/env python3
"""Sage-free exact parity, representative norm and positive-Gram census replay."""
import argparse
from collections import Counter
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint, digest
ROOT = Path(__file__).resolve().parents[2]

def determinant(matrix):
    a = [list(row) for row in matrix]
    previous = 1
    for k in range(len(a)-1):
        pivot = a[k][k]
        if not pivot:
            raise ArithmeticError('positive principal pivot missing')
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                value = a[i][j]*pivot-a[i][k]*a[k][j]
                if value % previous:
                    raise ArithmeticError('nonexact Bareiss division')
                a[i][j] = value//previous
            a[i][k] = 0
        previous = pivot
    return a[-1][-1]

def replay(directory, family, output):
    if output.exists():
        raise FileExistsError('preserve census replay')
    protocol = cert.read(directory/'protocol.json')
    for name, h in protocol['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen census source differs')
    path = directory/family/'generic-census.json'
    data = cert.read(path)
    f = next(f for f in cert.read(spec.ATLAS)['families'] if f['family'] == family)
    if data['status'] != 'COMPLETE_DECLARED_CENSUS' or data['protocol_hash'] != digest(protocol) or data['gram'] != f['generic_height_gram']:
        raise ArithmeticError('census input binding or completeness differs')
    gram = [[F(str(q)) for q in row] for row in data['gram']]
    scale = lcm(*(q.denominator for row in gram for q in row))
    G = [[int(q*scale) for q in row] for row in gram]
    minors = [determinant([row[:i] for row in G[:i]]) for i in range(1, 18)]
    if G != [list(row) for row in zip(*G)] or min(minors) <= 0 or data['positive_leading_principal_minors'] != minors or data['integer_gram_scale'] != scale:
        raise ArithmeticError('exact positive-Gram check failed')
    if len(data['records']) != 131072:
        raise ArithmeticError('missing parity records')
    for mask, row in enumerate(data['records']):
        rep = row['representative']
        if row['mask'] != mask or len(rep) != 17 or any(type(x) is not int for x in rep):
            raise ArithmeticError('parity index or integral vector differs')
        if any((rep[j]-(mask >> j)) % 2 for j in range(17)):
            raise ArithmeticError('parity identity failed')
        nonzero = [(i, v) for i, v in enumerate(rep) if v]
        norm = sum(v*G[i][j]*w for i, v in nonzero for j, w in nonzero)
        if F(norm, scale) != F(row['norm']):
            raise ArithmeticError('exact norm identity failed')
    selected = sorted(data['records'][1:], key=lambda r: (-F(r['norm']), r['mask']))[:43]
    histogram = dict(Counter(r['norm'] for r in data['records']))
    if selected != data['selected'] or histogram != data['norm_histogram']:
        raise ArithmeticError('selected-class ordering or histogram differs')
    checkpoint(output, {'schema': 'elliptic-curves.fresh-r17-generic-census-replay.v1',
        'status': 'PASS', 'family': family, 'census_sha256': cert.hashed(path),
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (Path(__file__).resolve(), directory/'protocol.json', spec.ATLAS)},
        'parities_checked': 131072, 'norm_histogram': histogram,
        'selected_masks': [r['mask'] for r in selected],
        'claim_boundary': 'Every parity has an exact integral representative of the stated norm, on an exactly positive definite Gram. Floating CVP optimality and exact covering radius are not claimed.'})
    print('REPLAYED FRESH R17 CENSUS', family, histogram, flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--family', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    replay(a.directory.resolve(), a.family, a.output.resolve())
