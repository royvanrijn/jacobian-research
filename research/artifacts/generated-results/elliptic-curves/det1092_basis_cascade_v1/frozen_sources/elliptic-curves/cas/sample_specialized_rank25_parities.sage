#!/usr/bin/env sage-python
"""Fixed sample of full specialized parity classes, with no point enumeration."""
import argparse
from decimal import Decimal
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint, digest
geometry = SourceFileLoader('specialized_parity_geometry', str(CAS/'prospective_half_lattice.sage')).load_module()

def run(directory, index):
    protocol = cert.read(directory/'protocol.json')
    for name, h in protocol['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen source changed')
    input_path = directory/f'input-{index}.json'
    if cert.hashed(input_path) != protocol['inputs'][index]['sha256']:
        raise ArithmeticError('frozen input changed')
    data = cert.read(input_path)
    state = data['initial_state']['state']
    points = [tuple(map(cert.F, p)) for p in state['reductions']['points']]
    cert.checked_rank(tuple(map(cert.F, data['curve'])), points,
                      state['reductions']['primes'], state['no_two_torsion_prime'])
    if len(points) != 25:
        raise ArithmeticError('certified dimension changed')
    gram = geometry.rounded_gram([[Decimal(q) for q in row] for row in data['metric_gram']], 1000000)
    oracle = geometry.CosetOracle(gram)
    excluded = {c['parity'] for c in data['old_centres']}
    masks = []; seen = set(excluded); counter = 0
    while len(masks) < protocol['samples_per_curve']:
        raw = f"{protocol['seed']}|{digest(data['curve'])}|{counter}".encode(); counter += 1
        mask = int.from_bytes(sha256(raw).digest(), 'big') & ((1 << 25)-1)
        if mask in seen or not (mask >> data['generic_dimension']):
            continue
        seen.add(mask); masks.append(mask)
    output = directory/f'census-{index}.json'
    if output.exists():
        raise FileExistsError('preserve previous parity sample')
    result = {'schema':'elliptic-curves.specialized-rank25-parity-sample.v1',
        'protocol_hash':digest(protocol), 'input_sha256':cert.hashed(input_path),
        'curve':data['curve'], 'family':data['family'], 'parameter':data['parameter'],
        'dimension':25, 'rounded_gram':[list(row) for row in gram], 'masks':masks,
        'records':[], 'status':'RUNNING',
        'claim_boundary':'Deterministic sample of new full specialized parity classes, with floating CVP proposals and exact parity/norm checks. No CVP optimality, covering radius, point search, or rank gain.'}
    checkpoint(output, result)
    for i, mask in enumerate(masks):
        norm, rep, error = oracle.solve([(mask >> j) & 1 for j in range(25)])
        if any((rep[j]-(mask >> j)) % 2 for j in range(25)):
            raise ArithmeticError('wrong parity representative')
        result['records'].append({'mask':mask, 'representative':list(rep), 'norm':norm, 'floating_error':str(error)})
        if (i+1) % 32 == 0:
            checkpoint(output, result)
        if (i+1) % 256 == 0:
            print('SPECIALIZED PARITY SAMPLE', index, i+1, flush=True)
    selected = sorted(result['records'], key=lambda r:(-r['norm'], r['mask']))[:43]
    old = sorted((c['metric_norm'] for c in data['old_centres']), reverse=True)[:43]
    if len(old) != 43:
        raise ArithmeticError('old comparison pool missing')
    # Fixed scheduling gate: median of43 deepest sampled classes improves by5%.
    gate = 20*selected[21]['norm'] >= 21*old[21]
    result.update(status='COMPLETE_DECLARED_SAMPLE', selected=selected,
        old_top43_median=old[21], new_top43_median=selected[21]['norm'],
        old_maximum=old[0], new_maximum=selected[0]['norm'],
        point_search_gate=gate)
    checkpoint(output, result)
    print('SPECIALIZED PARITY GATE', index, gate, old[21], '->', selected[21]['norm'], flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    a = p.parse_args()
    run(a.directory.resolve(), a.index)
