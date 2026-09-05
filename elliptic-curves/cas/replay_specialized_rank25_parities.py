#!/usr/bin/env python3
"""Exact replay of deterministic parity samples and rounded-form scheduling gates."""
import argparse
from hashlib import sha256
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint, digest
ROOT = Path(__file__).resolve().parents[2]

def replay(directory, output):
    if output.exists():
        raise FileExistsError('preserve prior sample replay')
    protocol = cert.read(directory/'protocol.json')
    for name, h in protocol['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('sample source changed')
    results = []
    for index in range(3):
        source = directory/f'input-{index}.json'
        path = directory/f'census-{index}.json'
        data, sample = cert.read(source), cert.read(path)
        if cert.hashed(source) != protocol['inputs'][index]['sha256'] or sample['input_sha256'] != cert.hashed(source):
            raise ArithmeticError('sample input changed')
        if sample['protocol_hash'] != digest(protocol) or sample['status'] != 'COMPLETE_DECLARED_SAMPLE':
            raise ArithmeticError('sample incomplete or protocol changed')
        gram = [[round(cert.F(q)*1000000) for q in row] for row in data['metric_gram']]
        if gram != sample['rounded_gram']:
            raise ArithmeticError('independent rational Gram rounding differs')
        def norm(rep):
            if len(rep) != 25 or any(type(x) is not int for x in rep):
                raise ArithmeticError('invalid integral representative')
            return sum(rep[i]*gram[i][j]*rep[j] for i in range(25) for j in range(25))
        for old in data['old_centres']:
            if norm(old['representative']) != old['metric_norm']:
                raise ArithmeticError('prior comparison metric differs')
        excluded = {r['parity'] for r in data['old_centres']}
        masks, seen, counter = [], set(excluded), 0
        while len(masks) < protocol['samples_per_curve']:
            raw = f"{protocol['seed']}|{digest(data['curve'])}|{counter}".encode()
            counter += 1
            mask = int.from_bytes(sha256(raw).digest(), 'big') & ((1 << 25)-1)
            if mask in seen or not mask >> data['generic_dimension']:
                continue
            seen.add(mask); masks.append(mask)
        if masks != sample['masks'] or [r['mask'] for r in sample['records']] != masks:
            raise ArithmeticError('deterministic sample differs')
        for row in sample['records']:
            if any((row['representative'][j]-(row['mask'] >> j)) % 2 for j in range(25)) or norm(row['representative']) != row['norm']:
                raise ArithmeticError('exact representative parity or norm differs')
        selected = sorted(sample['records'], key=lambda r:(-r['norm'], r['mask']))[:43]
        old = sorted((r['metric_norm'] for r in data['old_centres']), reverse=True)[:43]
        gate = 20*selected[21]['norm'] >= 21*old[21]
        if selected != sample['selected'] or gate != sample['point_search_gate']:
            raise ArithmeticError('frozen gate differs')
        results.append({'index':index, 'family':data['family'], 'parameter':data['parameter'],
            'census_sha256':cert.hashed(path), 'representatives_checked':len(masks),
            'old_top43_median':old[21], 'new_top43_median':selected[21]['norm'],
            'point_search_gate':gate})
    paths = [Path(__file__).resolve(), directory/'protocol.json']
    checkpoint(output, {'schema':'elliptic-curves.specialized-rank25-parity-replay.v1',
        'status':'PASS', 'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'curves':results, 'claim_boundary':'Exact sample, parity, rounded-form norms and fixed gate replay. Numerical height metric and CVP optimality remain unproved; no point search or rank gain.'})
    print('REPLAYED SPECIALIZED PARITY SAMPLES', results, flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    replay(a.directory.resolve(), a.output.resolve())
