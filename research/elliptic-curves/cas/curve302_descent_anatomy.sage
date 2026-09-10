#!/usr/bin/env sage-python
"""Exact retrospective half-ideal/Artin anatomy; no BNF or point search.

Run under timeout 120 sage -python ... --output NEW_DIRECTORY.
Every completed column is checkpointed. An absent/partial column is UNKNOWN.
The ten character and ideal words are frozen by the prior local filtration.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, pari, prod
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/curve302-descent-anatomy-v1'
ARITH = ART / 'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
FILTRATION = ART / 'curve302_recovered_quotient_local_filtration_v1.json'
sys.path.insert(0, str(Path(__file__).parent))
import icarm_curve302 as source


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(M):
    return [[str(M[i, j]) for j in range(3)] for i in range(3)]


def write(path, data):
    with path.open('x') as out:
        json.dump(data, out, indent=2, sort_keys=True)
        out.write('\n')


def compute(out):
    out.mkdir(parents=True, exist_ok=True)
    a, filt = read(ARITH), read(FILTRATION)
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in
                [Path(__file__), ARITH, FILTRATION, Path(source.__file__)]}
    protocol = {'schema': 'curve302.descent-anatomy.protocol.v1',
                'bindings': bindings, 'scope': 'Retrospective ten strict words; exact ideals and Artin entries only.',
                'bounds': {'wall_seconds': 120, 'workers': 1, 'pari_stack_max_bytes': 536870912,
                           'strict_words': 10, 'point_searches': 0, 'bnf_calls': 0},
                'failure_semantics': 'Missing Artin entries or unreturned computations remain UNKNOWN.'}
    if (out / 'protocol.json').exists():
        assert read(out / 'protocol.json') == protocol
    else:
        write(out / 'protocol.json', protocol)
    pari.allocatemem(64000000, 536870912, silent=True)
    R = PolynomialRing(QQ, 'z')
    f = R(a['cubic_ascending'])
    nf = pari.nfinit([pari(f), a['S_finite']])
    assert list(map(str, nf.nf_get_zk())) == a['maximal_order_basis']
    assert str(nf.disc()) == a['field_discriminant']
    S = [(p, j, P) for p in a['S_finite'] for j, P in enumerate(pari.idealprimedec(nf, p))]
    gammas, gcds, ys, points = [], [], [], []
    for point in source.POINTS:
        assert source.on_curve(point)
        x, y = map(QQ, point)
        X, Y = 4*x, 8*y + 4*x + 4
        assert Y**2 == f(X)
        d = ZZ(X.denominator()).sqrt()
        assert d in ZZ and d*d == X.denominator()
        aa, bb = ZZ(X*d*d), ZZ(Y*d**3)
        g = pari.Mod(pari(R([aa, -d*d])), pari(f))
        assert pari.nfeltnorm(nf, g) == bb*bb
        gammas.append(g)
        gcds.append(pari.idealadd(nf, bb, g))
        ys.append(abs(bb))
        points.append({'a': str(aa), 'b': str(bb), 'd': str(d)})
    words = filt['local_filtration_mod_2']['strict_kernel_public_words']
    assert len(words) == 10 and matrix(GF(2), words).rank() == 10
    betas = [prod(g for bit, g in zip(w, gammas) if bit) for w in words]
    ideals = []
    for k, (word, beta) in enumerate(zip(words, betas)):
        J = pari.idealhnf(nf, 1)
        for bit, I in zip(word, gcds):
            if bit:
                J = pari.idealmul(nf, J, I)
        for p, j, P in S:
            v = int(pari.idealval(nf, beta, P))
            assert v % 2 == 0
            e = v//2 - int(pari.idealval(nf, J, P))
            if e:
                J = pari.idealmul(nf, J, pari.idealpow(nf, P, e))
        assert pari.idealpow(nf, J, 2) == pari.idealhnf(nf, beta)
        assert pari.idealnorm(nf, J) == prod(y for bit, y in zip(word, ys) if bit)
        ideals.append(J)
    half_packet = {'schema': 'curve302.strict-half-ideals.v1', 'bindings': bindings,
                   'cubic_ascending': a['cubic_ascending'], 'maximal_order_basis': a['maximal_order_basis'],
                   'field_discriminant': a['field_discriminant'], 'signature': list(map(int, nf.nf_get_sign())),
                   'points': points, 'words': words,
                   'half_ideals': [record(J) for J in ideals],
                   'software': {'sage': version, 'pari': str(pari.version())}}
    hp = out / 'half-ideals.json'
    if hp.exists():
        assert read(hp) == half_packet
    else:
        write(hp, half_packet)
    print('PASS ten half-ideal square identities', flush=True)
    columns = []
    for k, J in enumerate(ideals):
        path = out / ('column-%02d.json' % k)
        if path.exists():
            columns.append(read(path))
            continue
        H, alpha = pari.idealred(nf, [J, 1])
        assert pari.idealmul(nf, H, alpha) == J
        good, removed = H, []
        for p, j, P in S:
            e = int(pari.idealval(nf, H, P))
            assert e >= 0
            if e:
                good = pari.idealmul(nf, good, pari.idealpow(nf, P, -e))
                removed.append({'p': p, 'j': j, 'exponent': e})
        N = ZZ(pari.idealnorm(nf, good))
        assert N > 0 and N % 2 and all(N % p for p in a['S_finite'])
        cyclic = good[0, 0] == N and good[1, 1] == good[2, 2] == 1
        entries = []
        for beta in betas:
            if not cyclic:
                entries.append({'status': 'UNKNOWN_NONCYCLIC_RING'})
                continue
            c = pari.nfalgtobasis(nf, beta)
            residue = ZZ(c[0] - good[0, 1]*c[1] - good[0, 2]*c[2]) % N
            if residue.gcd(N) != 1:
                entries.append({'status': 'UNKNOWN_NONUNIT', 'gcd': str(residue.gcd(N))})
                continue
            symbol = int(pari.kronecker(residue, N))
            assert symbol in [-1, 1]
            entries.append({'residue': str(residue), 'symbol': symbol, 'bit': int(symbol == -1)})
        col = {'index': k, 'reduced_ideal': record(H),
               'multiplier': [str(pari.lift(pari.nfbasistoalg(nf, alpha)).polcoef(j)) for j in range(3)],
               'removed_S_factors': removed, 'coprime_ideal': record(good),
               'norm': str(N), 'cyclic': bool(cyclic), 'entries': entries}
        write(path, col)
        columns.append(col)
        print('COLUMN', k, 'known entries', sum('bit' in e for e in entries), flush=True)
    complete = [c for c in columns if all('bit' in e for e in c['entries'])]
    rank = int(matrix(GF(2), [[e['bit'] for e in c['entries']] for c in complete]).rank())
    result = {'schema': 'curve302.descent-anatomy.v1', 'bindings': bindings,
              'half_ideals_sha256': sha(hp), 'strict_dimension': 10,
              'complete_columns': len(complete), 'artin_rank_lower_bound': rank,
              'ordinary_half_ideal_rank_interval': [max(8, rank), 10],
              'unit_kernel_dimension_interval': [0, min(2, 10-rank)],
              'S_class_rank_lower_bound': 10, 'full_S_class_rank': 'UNKNOWN',
              'full_selmer_dimension': '21 + c_S = 31 + (c_S - 10)',
              'full_rank': 'UNKNOWN; at least 31',
              'boundary': 'Artin rank is not CT rank. Lower bounds are not class-group completeness. The full Selmer identity uses derivative reciprocity and strict class-field identification; see proof note.',
              'columns': columns}
    dest = out / 'result.json'
    if dest.exists():
        assert read(dest) == result
    else:
        write(dest, result)
    print('RESULT Artin rank', rank, 'unit kernel', result['unit_kernel_dimension_interval'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, default=LOCAL)
    compute(p.parse_args().output)
