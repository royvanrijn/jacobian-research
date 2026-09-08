#!/usr/bin/env sage-python
"""Produce root witnesses excluding MW>=14 at old degree2 on six Mestre K3s.

All 65536 frame parity classes are accounted for. For each of the32509
classes with an isotropic lift of divisibility1, supply three independent
rational vertical roots. Completeness of PARI short-vector enumeration is
not used by the standalone verifier. One worker,120 seconds.
"""
import hashlib
import json
from pathlib import Path
import signal
import numpy as np
from sage.all import ZZ, matrix, pari

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
SOURCE = ART/'mestre_rational_ns_gram_v2.json'
OUT = ART/'mestre_degree_two_rank_gate_v1.json'
WITNESSES = ART/'mestre_degree_two_rank_gate_witnesses_v1.npz'


def choose_three(vectors):
    echelon, pivots, selected = [], [], []
    for v in vectors:
        row = v.copy()%101
        for e, p in zip(echelon, pivots):
            row = (row-row[p]*e)%101
        nz = np.flatnonzero(row)
        if len(nz):
            p = int(nz[0])
            row = row*pow(int(row[p]), -1, 101)%101
            echelon.append(row); pivots.append(p); selected.append(v)
            if len(selected) == 3:
                return np.array(selected, dtype=np.int16), pivots
    raise AssertionError('No three independent witnesses: leave this class unresolved')


def produce():
    assert not OUT.exists() and not WITNESSES.exists(), 'Preserve certificates'
    source = json.loads(SOURCE.read_text())
    r = source['rows'][0]
    G = matrix(ZZ, r['rational_NS_Gram'])
    B = matrix(ZZ, r['U_split_change'])
    H = -(B*G*B.transpose())[2:,2:]
    U = matrix(ZZ, pari(H).qflllgram())
    assert abs(U.det()) == 1
    V = U.transpose()*B[2:,:]
    L = -V*G*V.transpose()
    gram = np.array(L.rows(), dtype=np.int64)
    short = pari(L).qfminim(8, 1000000, 2)
    vectors = np.array(matrix(ZZ, short[2]).transpose().rows(), dtype=np.int64)
    norms = np.sum((vectors@gram)*vectors, axis=1)
    old_roots = vectors[norms == 2]
    norm8 = vectors[norms == 8]
    powers = 2**np.arange(16, dtype=np.int64)
    labels = (norm8%2)@powers
    bins = {}
    for label, v in zip(labels, norm8):
        bins.setdefault(int(label), []).append(v)
    words = ((np.arange(65536, dtype=np.int64)[:,None]>>np.arange(16))&1)
    pair = words@gram
    norm = np.sum(pair*words, axis=1)
    divisibility = np.gcd.reduce(np.column_stack((np.full(65536,2),norm//4,pair)), axis=1)
    valid = (norm%4 == 0)&(divisibility == 1)
    masks = np.flatnonzero(valid)
    assert len(masks) == 32509
    selected, minors = [], []
    for n, mask in enumerate(masks):
        old = 2*old_roots[(old_roots@pair[mask])%2 == 0]
        available = list(old)+bins.get(int(mask), [])
        witness, minor = choose_three(available)
        selected.append(witness); minors.append(minor)
        if n%4096 == 0:
            print('WITNESSES', n, len(masks), flush=True)
    np.savez_compressed(WITNESSES, masks=masks.astype(np.uint16),
                        root_numerators=np.array(selected,dtype=np.int16),
                        nonzero_minor_columns=np.array(minors,dtype=np.uint8))
    def rows(m): return [[int(c) for c in row] for row in m.rows()]
    result = {'schema': 'mestre.degree-two-rank-gate.v1',
              'status': 'ALL_OLD_DEGREE_TWO_Q_JACOBIAN_FIBRATIONS_HAVE_MW_RANK_AT_MOST13',
              'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), SOURCE]},
              'old_fibre': list(map(int,B[0])), 'hyperbolic_partner': list(map(int,B[1])),
              'frame_basis': rows(V), 'frame_gram': rows(L),
              'parity_classes': 65536, 'isotropic_divisibility_one_classes': len(masks),
              'witness_roots_per_class': 3, 'root_rank_lower_bound': 3,
              'rational_Picard_rank': 18, 'arithmetic_MW_rank_upper_bound': 13,
              'witness_file': str(WITNESSES.relative_to(ROOT)),
              'witness_sha256': hashlib.sha256(WITNESSES.read_bytes()).hexdigest(),
              'protocol': {'seconds': 120, 'workers': 1, 'short_vector_norm_bound': 8,
                           'short_vector_storage_cap': 1000000, 'rank_witness_prime': 101},
              'boundary': 'All Q-defined Jacobian elliptic fibrations of intersection2 with the original fibre on the six certified determinant468 parents. Higher old degree, other surfaces and all302 parents remain open. No nef representative, equation, generic section basis or302 parameter is constructed. Three exhibited roots per class suffice; short-vector enumeration completeness is not required.'}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('PASS_PRODUCED',len(masks),flush=True)


if __name__ == '__main__':
    signal.alarm(120)
    produce()
