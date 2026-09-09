#!/usr/bin/env python3
"""Freeze curve90 M26 control/M27 discovery seeds and its exact MW16 bank.

Sage Python required. No point searches. A bounded preparation is distinct
from launching either control or discovery search.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(folder):
    from sage.all import ZZ, matrix, pari
    from research_runtime.store import checkpoint
    from memory_rank_certificate import checked_rank
    from visibility_generic_bank import enumerate_bank
    import compact_mw16_specialization as spec
    folder.mkdir(exist_ok=False)
    source_names = ('prepare_curve90_v3.py', 'visibility_generic_bank.py',
                    'visibility_lattice_fast.py', 'visibility_lattice_v2.py',
                    'compact_mw16_specialization.py', 'memory_rank_certificate.py')
    inputs = [ART/'mw16_rank27_visibility_input_v1.json', spec.ATLAS,
              LOCAL/'extended20-mw16-pari-v1/a1-fibration-01-052/result.json',
              LOCAL/'prospective-mw16-h4096-v1/a1-fibration-01/generic-census.json']
    protocol = {'schema': 'curve90-v3-preparation.v1', 'generic_rank': 16,
                'generic_gram_scale': 2, 'scaled_shells': [16, 20],
                'closed_ellipsoid_bound': 20, 'node_limit': 20000000,
                'point_searches': 0, 'family': 'a1-fibration-01', 'parameter': '-1867/270',
                'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
                'sources': {str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in source_names},
                'gate': 'Independent M26 masked visibility and exact selector replay required '
                        'before a separately frozen M27 discovery search.',
                'claim_boundary': 'Exact bounded generic shell bank and existing seed proofs only.'}
    checkpoint(folder/'protocol.json', protocol)
    data, atlas, original, census = [json.loads(p.read_text()) for p in inputs]
    proof = data['point_proof']
    family = next(f for f in atlas['families'] if f['fibration_id'] == 'a1-fibration-01')
    model = tuple(map(F, proof['discovery_curve']))
    basis = tuple(tuple(map(F, p)) for p in proof['discovery_points'])
    parent_model, generic = spec.specialize(family, '-1867/270')
    scale = F(original['family_to_curve_scale_u'])
    if (len(basis) != 27 or proof['parameter'] != '-1867/270'
            or data['initial_points'] != proof['discovery_points'][:26]
            or data['generic_points'] != proof['discovery_points'][:16]
            or tuple(v / scale**e for v, e in zip(parent_model, (1, 2, 3, 4, 6))) != model
            or tuple((x/scale**2, y/scale**3) for x, y in generic) != basis[:16]):
        raise ArithmeticError('exact native generic basis/curve identity failed')
    old = proof['rank_certificate']
    for rank in (26, 27):
        fresh = checked_rank(model, basis[:rank], [s['prime'] for s in old['signatures']],
                             old['no_rational_2_torsion_prime'])
        if rank == 27 and json.loads(json.dumps(fresh)) != old:
            raise ArithmeticError('existing rank27 proof differs')
        checkpoint(folder/f'seed-M{rank}.json', {'curve': list(map(str, model)),
            'points': [list(map(str, p)) for p in basis[:rank]], 'generic_rank': 16,
            'rank_lower_bound': rank, 'proof': fresh})
    checkpoint(folder/'withheld-point.json', {'point': list(map(str, basis[-1])),
        'boundary': 'Evaluation only after masked selection is sealed; not a selector input.'})
    gram = [[F(x) for x in row] for row in family['generic_height_gram']]
    if gram != [[F(x) for x in row] for row in census['gram']]:
        raise ArithmeticError('native generic Gram differs from historical comparison')
    if any((2*x).denominator != 1 for row in gram for x in row):
        raise ArithmeticError('unexpected generic denominator')
    G = matrix(ZZ, [[int(2*x) for x in row] for row in gram])
    U = matrix(ZZ, pari(G).qflllgram()).transpose()
    if G.nrows() != 16 or not G.is_positive_definite() or abs(U.det()) != 1:
        raise ArithmeticError('invalid generic metric or nonunimodular transport')
    checkpoint(folder/'generic-metric.json', {'gram_scale': 2,
        'gram': [list(map(int, r)) for r in G.rows()],
        'LLL': [list(map(int, r)) for r in U.rows()], 'determinant': str(G.det())})
    start = time.monotonic()
    bank = enumerate_bank([list(map(int, r)) for r in (U*G*U.transpose()).rows()],
                          [list(map(int, r)) for r in U.rows()], bound=20,
                          shells=(16, 20), node_limit=protocol['node_limit'])
    bank['wall_seconds'] = time.monotonic() - start
    for row in bank['rows']:
        w = row['word']
        if (len(w) != 16 or sum((v % 2) << i for i, v in enumerate(w)) != row['mask']
                or sum(w[i]*int(G[i,j])*w[j] for i in range(16) for j in range(16)) != row['norm']):
            raise ArithmeticError('exact exported generic word/norm differs')
    comparison = {r['mask']: int(2*F(r['norm'])) for r in census['records']
                  if F(r['norm']) in (8, 10)}
    bank['historical_census_same_shell_masks'] = comparison == {r['mask']: r['norm'] for r in bank['rows']}
    bank['generic_gram_scale'] = 2
    checkpoint(folder/'parent-bank.json', bank)
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('preparation binding changed')
    checkpoint(folder/'prepared.json', {'status': 'PASS_EXACT_MW16_BANK_AND_SEEDS',
        'protocol_sha256': sha(folder/'protocol.json'), 'point_searches': 0,
        'files': {p.name: sha(p) for p in sorted(folder.glob('*.json'))},
        'bank_rows': len(bank['rows']), 'bank_seconds': bank['wall_seconds'],
        'historical_census_same_shell_masks': bank['historical_census_same_shell_masks']})
    print('PREPARED', len(bank['rows']), bank['nodes'], bank['wall_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    prepare(parser.parse_args().folder)
