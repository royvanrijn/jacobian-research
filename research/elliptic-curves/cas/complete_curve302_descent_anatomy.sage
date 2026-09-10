#!/usr/bin/env sage-python
"""Complete the five retained nonunit residues, only at primes 47,67,89."""
import json
import hashlib
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, pari, prod

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / 'artifacts/local/elliptic-curves/curve302-descent-anatomy-v1'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
DEST = ART / 'curve302_descent_anatomy_v1.json'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def compute():
    old = read(WORK / 'result.json')
    h = read(WORK / 'half-ideals.json')
    a = read(ART / 'rank_jump_curve302_strict_constructor_arithmetic_v1.json')
    R = PolynomialRing(QQ, 'z')
    f = R(h['cubic_ascending'])
    nf = pari.nfinit([pari(f), a['S_finite']])
    gs = [pari.Mod(pari(R([QQ(P['a']), -QQ(P['d'])**2])), pari(f)) for P in h['points']]
    betas = [prod(g for bit, g in zip(w, gs) if bit) for w in h['words']]
    repairs = []
    M = []
    for c in old['columns']:
        bits = []
        for i, e in enumerate(c['entries']):
            if 'bit' in e:
                bits.append(e['bit'])
                continue
            p = ZZ(e['gcd'])
            assert p in [47, 67, 89] and p.is_prime(proof=True) and f.discriminant() % p
            N = ZZ(c['norm'])
            exponent = N.valuation(p)
            cofactor = N // p**exponent
            H = pari(matrix(QQ, c['coprime_ideal']))
            Ps = [P for P in pari.idealprimedec(nf, p) if pari.idealval(nf, H, P) > 0]
            assert len(Ps) == 1
            P = Ps[0]
            assert P[2] == P[3] == 1 and pari.idealval(nf, H, P) == exponent
            beta = betas[i]
            coords = pari.nfalgtobasis(nf, beta)
            residue = ZZ(coords[0] - H[0, 1]*coords[1] - H[0, 2]*coords[2]) % cofactor
            assert residue.gcd(cofactor) == 1
            symbol = int(pari.kronecker(residue, cofactor))
            v = int(pari.idealval(nf, beta, P))
            assert v % 2 == 0
            local_bit = int(pari.nfislocalpower(nf, P, beta, 2) == 0)
            bit = int(symbol == -1) ^ ((int(exponent) % 2)*local_bit)
            tc = pari.nfalgtobasis(nf, pari.Mod(pari(R.gen()), pari(f)))
            root = int((tc[0] - H[0, 1]*tc[1] - H[0, 2]*tc[2]) % p)
            assert pari.idealval(nf, pari.Mod(pari(R.gen()-root), pari(f)), P) > 0
            repairs.append({'column': c['index'], 'character': i, 'prime': int(p),
                            'root': root, 'exponent': int(exponent), 'cofactor': str(cofactor),
                            'cofactor_residue': str(residue), 'cofactor_symbol': symbol,
                            'valuation': v, 'local_bit': local_bit, 'bit': bit})
            bits.append(bit)
        M.append(bits)
    A = matrix(GF(2), M).transpose()
    assert len(repairs) == 5
    rank = int(A.rank())
    result = {'schema': 'curve302.descent-anatomy.complete.v1', 'status': 'PASS_EXACT_ARTIN_MATRIX',
              'bindings': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), WORK/'result.json', WORK/'half-ideals.json']},
              'half_ideal_packet': h, 'initial_artin_packet': old,
              'repairs': repairs, 'artin_matrix': [[int(v) for v in row] for row in A.rows()],
              'artin_rank': rank, 'right_kernel_words': [list(map(int, v)) for v in A.right_kernel().basis()],
              'selected_ideal_indices': list(map(int, A.pivots())),
              'ordinary_half_ideal_rank_interval': [max(8, rank), 10],
              'unit_kernel_dimension_interval': [0, min(2, 10-rank)],
              'full_S_class_rank': 'UNKNOWN; at least 10',
              'full_selmer_dimension': '21 + c_S = 31 + (c_S - 10)',
              'full_rank': 'UNKNOWN; at least 31'}
    if DEST.exists():
        assert read(DEST) == result
    else:
        with DEST.open('x') as out:
            json.dump(result, out, indent=2, sort_keys=True)
            out.write('\n')
    print('PASS complete matrix: rank', rank, 'kernel', result['right_kernel_words'])


if __name__ == '__main__':
    compute()
