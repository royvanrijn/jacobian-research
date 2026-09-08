#!/usr/bin/env python3
"""Replay all oracle translations with rational group arithmetic, without Sage/CVP."""
import argparse
from collections import Counter
from functools import lru_cache
from pathlib import Path
import certify_compact_r17_candidates as cert
from search_nagao_u42_skew_height import short_add, short_multiply
from search_observability import point_visibility
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT/'artifacts/generated-results/elliptic-curves/curve542_translated_visibility_v1.json'
INPUT = ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1/a1-fibration-04/candidate-00/result.json'

def replay(output):
    if output.exists():
        raise FileExistsError('preserve the previous replay')
    data = cert.read(ARTIFACT)
    for name, h in data['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('translation input/source changed')
    original = cert.read(INPUT)
    model = tuple(map(cert.F, data['curve']))
    basis = [tuple(map(cert.F, p)) for p in data['basis']]
    oracle = tuple(map(cert.F, data['oracle_point']))
    proof = cert.checked_rank(model, basis+[oracle])
    if proof['rank_lower_bound'] != 26:
        raise ArithmeticError('oracle direction independence failed')
    @lru_cache(None)
    def multiple(i, k):
        return short_multiply(model[3], basis[i], k)
    counts = Counter()
    for i, row in enumerate(data['proposals']):
        if len(row['word']) != 25 or row['sign'] not in (-1, 1):
            raise ArithmeticError('bad group word')
        point = oracle[0], row['sign']*oracle[1]
        for j, k in enumerate(row['word']):
            if k:
                point = short_add(model[3], point, multiple(j, k))
        if point != tuple(map(cert.F, row['point'])):
            raise ArithmeticError('independent rational group sum differs')
        if not cert.is_on_weierstrass_curve(model, point):
            raise ArithmeticError('translated point not on curve')
        visibility = point_visibility(original['charts'][row['chart_index']]['search'], point)
        if visibility != row['visibility']:
            raise ArithmeticError('exact chart visibility differs')
        counts[visibility['status']] += 1
        if (i+1) % 10 == 0:
            print('REPLAYED TRANSLATIONS', i+1, flush=True)
    if len(data['proposals']) != 86 or dict(counts) != data['status_counts']:
        raise ArithmeticError('translation census differs')
    paths = [Path(__file__).resolve(), ARTIFACT, INPUT,
             ROOT/'elliptic-curves/cas/search_nagao_u42_skew_height.py',
             ROOT/'elliptic-curves/cas/search_observability.py']
    checkpoint(output, {'schema':'elliptic-curves.curve542-translated-visibility-replay.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'status':'PASS', 'exact_group_words_checked':86, 'status_counts':dict(counts),
        'claim_boundary':'Exact oracle group translations and visibility only; no floating CVP optimality or new-curve claim.'})
    print('REPLAYED ALL 86 TRANSLATIONS', dict(counts), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    replay(p.parse_args().output)
