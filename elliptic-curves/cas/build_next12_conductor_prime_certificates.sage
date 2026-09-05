#!/usr/bin/env sage-python
"""Generate bounded Lucas primality proof trees for one fixed discriminant."""
import argparse
from pathlib import Path
import sys
from sage.all import ZZ, prod
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/next12_rank22_exact_conductor_v1.json'

def build(output):
    if output.exists():
        raise FileExistsError('preserve prime certificates')
    original = cert.read(INPUT)
    targets = [int(p) for p, e in original['discriminant_factorization']]
    if max(p.bit_length() for p in targets) > 129:
        raise ArithmeticError('fixed 129-bit prime proof gate exceeded')
    data = {'schema': 'elliptic-curves.lucas-prime-proofs.v1', 'status': 'RUNNING',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (Path(__file__).resolve(), INPUT)},
            'targets': list(map(str, targets)), 'nodes': {}}
    checkpoint(output, data)
    def prove(p):
        if str(p) in data['nodes']:
            return
        if p == 2:
            node = {'prime': '2', 'base_case': True}
        else:
            factors = [(int(q), int(e)) for q, e in ZZ(p-1).factor(proof=True)]
            if prod(ZZ(q)**e for q, e in factors) != p-1:
                raise ArithmeticError('p-1 factorization incomplete')
            for q, e in factors:
                prove(q)
            witness = next((a for a in range(2, min(p, 100001))
                            if pow(a, p-1, p) == 1 and all(
                                ZZ(pow(a, (p-1)//q, p)-1).gcd(p) == 1 for q, e in factors)), None)
            if witness is None:
                raise ArithmeticError('fixed Lucas witness budget exhausted')
            node = {'prime': str(p), 'p_minus_one_factorization': [[str(q), e] for q, e in factors],
                    'witness': witness}
        data['nodes'][str(p)] = node
        checkpoint(output, data)
        print('LUCAS PRIME PROOF', p, flush=True)
    for p in targets:
        prove(p)
    data['status'] = 'COMPLETE_LUCAS_CERTIFICATES'
    checkpoint(output, data)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    build(p.parse_args().output.resolve())
