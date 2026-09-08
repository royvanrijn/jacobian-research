#!/usr/bin/env sage-python
"""Fresh complete parity enumeration on the six exact compact R17 Grams.

Floating CVP proposes one representative per parity. Exact parity and norm
checks do not prove that these representatives are globally shortest.
No historical masks, curve parameters, scores or public point data are read.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from math import lcm
from pathlib import Path
import sys
from sage.all import ZZ, matrix
ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint, digest
geometry = SourceFileLoader('fresh_r17_generic_geometry', str(CAS/'prospective_half_lattice_v2.sage')).load_module()
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-r17-fresh-generic-census-v1'

def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in (
        Path(__file__).resolve(), spec.ATLAS, Path(spec.__file__).resolve(),
        CAS/'prospective_half_lattice_v2.sage', CAS/'research_runtime/store.py')}

def prepare(directory):
    if (directory/'protocol.json').exists():
        raise FileExistsError('preserve census protocol')
    atlas = cert.read(spec.ATLAS)
    checkpoint(directory/'protocol.json', {
        'schema': 'elliptic-curves.fresh-six-r17-parity-census.v1', 'sources': sources(),
        'families': [f['family'] for f in atlas['families']], 'dimension': 17,
        'classes_per_family': 131072, 'selected_nonzero_classes': 43,
        'order': 'computed representative norm descending, mask ascending',
        'wall_seconds_per_family': 300, 'rss_bytes': 1073741824, 'maximum_workers': 2,
        'checkpoint_every': 4096,
        'mathematical_gate': 'All six compact R17 models and their 102 generic sections are exact. Their earlier point pilot used inherited generic masks; a fresh enumeration removes that historical input and checks the full finite parity population before a wider balanced search.',
        'claim_boundary': 'Complete parity coverage and exact representative norms only. Floating CVP optimality, covering radius, rank jumps and rational point existence are not certified.',
        'software': {'sage': str(__import__('sage.env', fromlist=['SAGE_VERSION']).SAGE_VERSION)}})

def census(directory, family):
    protocol = cert.read(directory/'protocol.json')
    if protocol['sources'] != sources():
        raise ArithmeticError('frozen census source differs')
    f = next(f for f in cert.read(spec.ATLAS)['families'] if f['family'] == family)
    gram = [[F(str(x)) for x in row] for row in f['generic_height_gram']]
    scale = lcm(*(x.denominator for row in gram for x in row))
    integral = [[int(x*scale) for x in row] for row in gram]
    G = matrix(ZZ, integral)
    minors = [int(G[:i, :i].det()) for i in range(1, 18)]
    if G != G.transpose() or any(v <= 0 for v in minors):
        raise ArithmeticError('positive generic Gram gate failed')
    oracle = geometry.CosetOracle(integral)
    output = directory/family/'generic-census.json'
    data = cert.read(output) if output.exists() else {
        'schema': 'elliptic-curves.fresh-r17-generic-census.v1', 'protocol_hash': digest(protocol),
        'family': family, 'gram': f['generic_height_gram'], 'integer_gram_scale': scale,
        'positive_leading_principal_minors': minors, 'records': [], 'status': 'RUNNING'}
    if data['protocol_hash'] != digest(protocol) or data['gram'] != f['generic_height_gram']:
        raise ArithmeticError('census resume binding differs')
    for mask in range(len(data['records']), 131072):
        norm, rep, error = oracle.solve(tuple((mask >> j) & 1 for j in range(17)))
        if len(rep) != 17 or any((rep[j]-(mask >> j)) % 2 for j in range(17)):
            raise ArithmeticError('generic parity mismatch')
        data['records'].append({'mask': mask, 'norm': str(F(norm, scale)),
                                'representative': list(rep), 'cvp_error': error})
        if (mask+1) % 4096 == 0:
            checkpoint(output, data)
            print('FRESH R17 PARITIES', family, mask+1, flush=True)
    selected = sorted(data['records'][1:], key=lambda r: (-F(r['norm']), r['mask']))[:43]
    data.update(status='COMPLETE_DECLARED_CENSUS', selected=selected,
                norm_histogram=dict(Counter(r['norm'] for r in data['records'])))
    checkpoint(output, data)
    print('FRESH R17 GENERIC CLASSES', family, data['norm_histogram'], flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['prepare', 'census'])
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    p.add_argument('--family')
    a = p.parse_args()
    prepare(a.directory) if a.stage == 'prepare' else census(a.directory, a.family)
