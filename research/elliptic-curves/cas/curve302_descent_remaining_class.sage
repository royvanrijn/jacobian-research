#!/usr/bin/env sage-python
"""Test the single Artin-kernel half ideal against six inherited characters."""
import json
import hashlib
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, matrix, pari, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
INPUT = ART / 'curve302_descent_anatomy_v1.json'
ARITH = ART / 'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
ANCHORS = ART / 'rank_jump_generic_only_class_anchors_v1.json'
OUTPUT = ART / 'curve302_descent_remaining_class_v1.json'


def read(p):
    return json.loads(p.read_text())


def compute():
    d, a, anchors = read(INPUT), read(ARITH), read(ANCHORS)['cases'][1]
    R = PolynomialRing(QQ, 'z')
    f = R(a['cubic_ascending'])
    nf = pari.nfinit([pari(f), a['S_finite']])
    assert d['artin_rank'] == 9
    word, = d['right_kernel_words']
    J = pari.idealhnf(nf, 1)
    for bit, H in zip(word, d['half_ideal_packet']['half_ideals']):
        if bit:
            J = pari.idealmul(nf, J, pari(matrix(QQ, H)))
    H, alpha = pari.idealred(nf, [J, 1])
    assert pari.idealmul(nf, H, alpha) == J
    good, removed = H, []
    for p in a['S_finite']:
        for j, P in enumerate(pari.idealprimedec(nf, p)):
            e = int(pari.idealval(nf, H, P))
            assert e >= 0
            if e:
                removed.append((p, j, P, e))
                good = pari.idealmul(nf, good, pari.idealpow(nf, P, -e))
    N = ZZ(pari.idealnorm(nf, good))
    assert N > 0 and all(N % p for p in a['S_finite'])
    assert good[0, 0] == N and good[1, 1] == good[2, 2] == 1
    gs = [pari.Mod(pari(R(c['beta_ascending'])), pari(f)) for c in a['generic_classes']]
    entries = []
    for mask0 in anchors['generic_coefficient_masks']:
        mask = int(mask0)
        beta = prod(g for i, g in enumerate(gs) if mask >> i & 1)
        coords = pari.nfalgtobasis(nf, beta)
        residue = ZZ(coords[0]-good[0, 1]*coords[1]-good[0, 2]*coords[2]) % N
        if residue.gcd(N) != 1:
            entries.append({'mask': mask0, 'status': 'UNKNOWN_NONUNIT', 'gcd': str(residue.gcd(N))})
            continue
        symbol = int(pari.kronecker(residue, N))
        bits = []
        for p, j, P, e in removed:
            bit = int(pari.nfislocalpower(nf, P, beta, 2) == 0)
            bits.append({'p': p, 'j': j, 'exponent': e, 'bit': bit})
        bit = (int(symbol == -1) + sum(t['bit']*t['exponent'] for t in bits)) % 2
        entries.append({'mask': mask0, 'residue': str(residue), 'symbol': symbol, 'S_terms': bits, 'bit': bit})
    def record(M):
        return [[str(M[i, j]) for j in range(3)] for i in range(3)]
    detected = any(e.get('bit') == 1 for e in entries)
    result = {'schema': 'curve302.descent-remaining-class.v1',
              'bindings': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__), INPUT, ARITH, ANCHORS]},
              'strict_word': word, 'reduced_ideal': record(H), 'coprime_ideal': record(good), 'norm': str(N),
              'multiplier': [str(pari.lift(pari.nfbasistoalg(nf, alpha)).polcoef(i)) for i in range(3)],
              'entries': entries, 'remaining_class_detected': detected,
              'ordinary_half_ideal_rank_interval': [10 if detected else 9, 10],
              'unit_kernel_dimension_interval': [0, 0 if detected else 1],
              'boundary': 'The extra characters are ordinary unramified but not S-split, so every removed S-prime contribution is restored. A zero vector does not prove principality.'}
    if OUTPUT.exists():
        assert read(OUTPUT) == result
    else:
        with OUTPUT.open('x') as out:
            json.dump(result, out, indent=2, sort_keys=True)
            out.write('\n')
    print('PASS remaining-class test', [e.get('bit', e.get('status')) for e in entries],
          'unit kernel', result['unit_kernel_dimension_interval'])


if __name__ == '__main__':
    compute()
