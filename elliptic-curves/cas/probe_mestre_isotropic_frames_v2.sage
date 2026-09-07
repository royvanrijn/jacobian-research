#!/usr/bin/env sage-python
"""Fixed 192-class lattice gate at old degrees3,4,5 on the determinant468 rational NS lattice.

No equation or nefness claim. All geometric roots and their Galois-fixed
span are checked, so rational roots alone cannot overestimate MW rank.
120 seconds, one worker, at most1024 deterministic proposals per degree.
"""
import argparse
import hashlib
import json
from pathlib import Path
import signal
from sage.all import QQ, ZZ, matrix, vector, pari, gcd, xgcd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'artifacts/generated-results/elliptic-curves/mestre_rational_ns_gram_v2.json'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/mestre_isotropic_frames_v2.json'
CHECKPOINT = ROOT/'artifacts/local/elliptic-curves/mestre-isotropic-frames-v2.jsonl'


def rows(m):
    return [[int(c) for c in r] for r in m.rows()]


def frame(G, f, z):
    C = matrix(ZZ, [f*G, z*G]).right_kernel_matrix()
    H = -C*G*C.transpose()
    U = matrix(ZZ, pari(H).qflllgram())
    assert abs(U.det()) == 1
    C = U.transpose()*C
    H = -C*G*C.transpose()
    assert H.is_positive_definite()
    result = pari(H).qfminim(2, 50000, 2)
    V = matrix(ZZ, result[2])
    assert int(result[0]) == 2*V.ncols()
    assert all(v*H*v == 2 for v in V.columns())
    return C, H, V.transpose()


def section_class(G, f):
    pair = f*G
    value = ZZ(0)
    z = vector(ZZ, G.nrows())
    for i, a in enumerate(pair):
        d, s, t = xgcd(value, a)
        z *= s
        z[i] += t
        value = d
    assert value == 1 and z*G*f == 1
    z -= ((z*G*z+2)//2)*f
    assert z*G*z == -2 and z*G*f == 1
    return z


def compute(check=False):
    source = json.loads(SOURCE.read_text())
    row = source['rows'][0]
    assert all(r['rational_NS_Gram'] == row['rational_NS_Gram'] for r in source['rows'])
    G = matrix(ZZ, row['rational_NS_Gram'])
    B = matrix(ZZ, row['U_split_change'])
    L = -(B*G*B.transpose())[2:,2:]
    U = matrix(ZZ, pari(L).qflllgram())
    assert abs(U.det()) == 1
    V = U.transpose()*B[2:,:]
    L = -V*G*V.transpose()
    glue = row['eligible_index_two_geometric_glues'][0]
    GG = matrix(ZZ, glue['geometric_Gram'])
    action = matrix(ZZ, glue['Galois_action'])
    change = matrix.identity(QQ, 19)
    change[glue['replace_basis_index']] = vector(QQ, glue['fixed_numerator']+[1])/2
    embedding = change.inverse()[:18,:]
    assert embedding*GG*embedding.transpose() == G
    assert embedding*action == embedding
    records = []
    seen = set()
    if not check:
        CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
        assert not CHECKPOINT.exists(), 'Preserve prior checkpoint'
    for degree in (3,4,5):
        seen = set()
        for proposal in range(1024):
            digest = hashlib.sha256(('mestre468-isotropic-v2:'+str(degree)+':'+str(proposal)).encode()).digest()
            word = tuple(int(c)%degree for c in digest[:16])
            if not any(word) or word in seen:
                continue
            seen.add(word)
            w = vector(ZZ, word)
            norm = w*L*w
            if norm % (2*degree):
                continue
            f = (norm//(2*degree))*B[0]+degree*B[1]+w*V
            assert f*G*f == 0 and f*G*B[0] == degree
            if gcd(list(f*G)) != 1:
                continue
            z = section_class(G, f)
            cf, hf, rf = frame(G, f, z)
            fg, zg = vector(ZZ, f*embedding), vector(ZZ, z*embedding)
            cg, hg, rg = frame(GG, fg, zg)
            assert hf.det() == 468 and hg.det() == GG.det() == 468
            geometric_roots = rg*cg
            # Averaging spans the invariant part of the full geometric root space.
            sums = geometric_roots*(matrix.identity(ZZ,19)+action)
            fixed_rank = sums.rank()
            rank = 16-fixed_rank
            assert fixed_rank >= rf.rank()
            record = {'proposal': proposal, 'degree': degree, 'word': list(word), 'lift_norm': int(norm),
                      'fibre_class': list(map(int,f)), 'section_class': list(map(int,z)),
                      'rational_frame_basis': rows(cf), 'rational_frame_gram': rows(hf),
                      'rational_root_representatives': rows(rf),
                      'geometric_frame_basis': rows(cg), 'geometric_frame_gram': rows(hg),
                      'geometric_root_representatives': rows(rg),
                      'geometric_root_rank': int(rg.rank()),
                      'rational_root_rank': int(rf.rank()),
                      'fixed_geometric_root_rank': int(fixed_rank),
                      'arithmetic_MW_rank_if_realized': int(rank)}
            records.append(record)
            if not check:
                with CHECKPOINT.open('a') as out:
                    out.write(json.dumps(record, sort_keys=True)+'\n')
            print('FRAME', len(records)-1, degree, proposal, 'ranks', rf.rank(), rg.rank(), fixed_rank, 'MW', rank, flush=True)
            if len(records) == 64*(degree-2):
                break
        assert len(records) == 64*(degree-2), 'Incomplete declared roster'
    return {'schema': 'mestre.isotropic-frames.v2', 'status': 'PASS_FINITE_LATTICE_GATE',
            'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), SOURCE]},
            'protocol': {'classes': 192, 'proposal_cap': 1024, 'seconds': 120, 'workers': 1,
                         'old_fibre_degrees': [3,4,5], 'signed_root_storage_cap': 100000},
            'reduced_frame_basis_in_old_NS': rows(V), 'reduced_frame_gram': rows(L),
            'rational_to_geometric_embedding': rows(embedding), 'records': records,
            'rank_histogram': {str(r): sum(v['arithmetic_MW_rank_if_realized'] == r for v in records) for r in sorted({v['arithmetic_MW_rank_if_realized'] for v in records})},
            'boundary': 'A finite isotropic-class gate, not a classification of fibrations or a302 parent. Nef representatives, equations, section coordinates and302 specialization are not constructed. Geometric roots and their Galois-fixed span are included; a rational root census alone does not determine arithmetic MW rank.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--check', action='store_true'); args = p.parse_args()
    signal.alarm(120)
    result = compute(args.check)
    if args.check:
        assert result == json.loads(OUT.read_text())
    else:
        assert not OUT.exists(), 'Preserve certificate'
        OUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('PASS', result['rank_histogram'], flush=True)
