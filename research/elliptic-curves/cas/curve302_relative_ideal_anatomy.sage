#!/usr/bin/env sage-python
"""Full inherited even-valuation ideal image versus the frozen strict302 block.

Use all unramified characters in the displayed rational31 space, no new
characters from an unbounded search. Fixed18 ideal columns. Checkpoint every
column; only factor nonunit gcds supported on primes <=1009. Run timeout120s.
"""
import hashlib
import json
from pathlib import Path
from sage.all import AA, QQ, ZZ, GF, PolynomialRing, matrix, pari, prod, prime_range

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
WORK = ROOT/'artifacts/local/elliptic-curves/curve302-relative-ideal-anatomy-v1'
ARITH = ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
STRICT = ART/'curve302_descent_anatomy_v1.json'
GEN = [ART/('det1092_generic_virtual_units_v1/basis-%02d.json' % i) for i in range(8)]
DEST = ART/'curve302_relative_ideal_anatomy_v1.json'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in matrix(QQ, M).rows()]


def save(path, data):
    if path.exists():
        assert read(path) == data
    else:
        with path.open('x') as out:
            json.dump(data, out, indent=2, sort_keys=True)
            out.write('\n')


def compute():
    WORK.mkdir(parents=True, exist_ok=True)
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ARITH, STRICT, *GEN]}
    protocol = {'schema': 'curve302.relative-ideal-anatomy.protocol.v1', 'bindings': bindings,
                'scope': 'Retrospective all8 inherited even-valuation classes and the ten frozen strict classes; all ordinary unramified characters in the displayed31 space.',
                'bounds': {'wall_seconds': 120, 'workers': 1, 'columns': 18,
                           'character_dimension_max': 31, 'nonunit_gcd_trial_prime_bound': 1009,
                           'new_points': 0, 'new_parameters': 0, 'bnf_calls': 0},
                'failure_semantics': 'Incomplete factor support, noncyclic quotients and undetected ideal relations remain UNKNOWN.'}
    save(WORK/'protocol.json', protocol)
    a, s = read(ARITH), read(STRICT)
    R = PolynomialRing(QQ, 'z')
    f = R(a['cubic_ascending'])
    nf = pari.nfinit([pari(f), a['S_finite']])
    assert list(map(str, nf.nf_get_zk())) == a['maximal_order_basis']
    h = s['half_ideal_packet']
    gs = [pari.Mod(pari(R([QQ(p['a']), -QQ(p['d'])**2])), pari(f)) for p in h['points']]
    places = [(p, j, P) for p in a['S_finite'] for j, P in enumerate(pari.idealprimedec(nf, p))]
    constraints, vals, unit_checks = [], [], []
    for p, j, P in places:
        v = [int(pari.idealval(nf, g, P)) for g in gs]
        vals.append(v)
        constraints.append([x % 2 for x in v])
        if p != 2:
            continue
        e = int(P[2])
        assert int(P[3]) == 1
        pi = pari(2) if e == 1 else pari.nfbasistoalg(nf, pari.idealappr(nf, P))
        assert pari.idealval(nf, pi, P) == 1
        for k in range(1, 2*e+1):
            unit = 1+pi**k
            bits = [int(pari.nfhilbert(nf, g, unit, P) == -1) for g in gs]
            constraints.append(bits)
            unit_checks.append({'prime_index': j, 'power': k, 'unit': str(unit), 'bits': bits})
    roots = f.roots(AA, multiplicities=False)
    polynomials = [R([QQ(pari.lift(g).polcoef(i)) for i in range(3)]) for g in gs]
    signs = [[int(g(root)<0) for g in polynomials] for root in roots]
    constraints.extend(signs)
    C = matrix(GF(2), constraints)
    words = C.right_kernel().basis_matrix()
    betas = [prod(g for bit, g in zip(word, gs) if bit) for word in words]
    character_packet = {'constraints': [[int(x) for x in row] for row in C],
                        'bad_prime_valuations': vals, 'dyadic_units': unit_checks,
                        'real_signs': signs, 'words': [[int(x) for x in row] for row in words],
                        'dimension': int(words.nrows()), 'valuation_rank': int(matrix(GF(2), vals).rank()),
                        'even_dimension': 31-int(matrix(GF(2), vals).rank())}
    for beta in betas:
        for p, j, P in places:
            assert pari.idealval(nf, beta, P) % 2 == 0
            if p == 2:
                assert pari.nfislocalpower(nf, P, beta, 2) or pari.nfislocalpower(nf, P, 5*beta, 2)
    save(WORK/'characters.json', character_packet)
    print('CHARACTERS',len(betas),'EVEN',character_packet['even_dimension'],flush=True)
    columns = []
    roster = [('generic', i, pari(matrix(QQ, read(path)['ideal']))) for i, path in enumerate(GEN)]
    roster += [('strict', i, pari(matrix(QQ, c['reduced_ideal']))) for i, c in enumerate(s['initial_artin_packet']['columns'])]
    for index, (kind, source_index, H) in enumerate(roster):
        path = WORK/('column-%02d.json' % index)
        if path.exists():
            columns.append(read(path))
            continue
        good, bad = H, []
        for p, j, P in places:
            e = int(pari.idealval(nf, H, P))
            assert e >= 0
            if e:
                good = pari.idealmul(nf, good, pari.idealpow(nf, P, -e))
                bad.append((p, j, P, e))
        N = ZZ(pari.idealnorm(nf, good))
        assert N > 0 and N % 2 and all(N % p for p in a['S_finite'])
        cyclic = good[0, 0] == N and good[1, 1] == good[2, 2] == 1
        entries = []
        for beta in betas:
            if not cyclic:
                entries.append({'status': 'UNKNOWN_NONCYCLIC'})
                continue
            cs = pari.nfalgtobasis(nf, beta)
            residue = ZZ(cs[0]-good[0, 1]*cs[1]-good[0, 2]*cs[2]) % N
            gcd = residue.gcd(N)
            missing = gcd
            factors = []
            for p in prime_range(3, 1010):
                if missing % p:
                    continue
                while missing % p == 0:
                    missing //= p
                factors.append(int(p))
            if missing != 1:
                entries.append({'status': 'UNKNOWN_NONUNIT_SUPPORT', 'gcd': str(gcd), 'unfactored': str(missing)})
                continue
            cofactor = N
            repairs = []
            for p in factors:
                exponent = int(N.valuation(p))
                cofactor //= ZZ(p)**exponent
                Ps = [P for P in pari.idealprimedec(nf, p) if pari.idealval(nf, good, P) > 0]
                assert len(Ps) == 1
                P = Ps[0]
                assert P[2] == P[3] == 1 and pari.idealval(nf, good, P) == exponent
                v = int(pari.idealval(nf, beta, P))
                assert v % 2 == 0
                bit = int(not pari.nfislocalpower(nf, P, beta, 2))
                repairs.append({'p': p, 'exponent': exponent, 'valuation': v, 'bit': bit})
            reduced_residue = residue % cofactor
            assert reduced_residue.gcd(cofactor) == 1
            symbol = int(pari.kronecker(reduced_residue, cofactor))
            assert symbol in [-1, 1]
            bad_terms = [{'p': p, 'j': j, 'exponent': e,
                          'bit': int(not pari.nfislocalpower(nf, P, beta, 2))} for p, j, P, e in bad]
            bit = (int(symbol == -1) + sum(t['exponent']*t['bit'] for t in repairs+bad_terms)) % 2
            entries.append({'residue': str(residue), 'gcd': str(gcd), 'cofactor': str(cofactor),
                            'symbol': symbol, 'repairs': repairs, 'bad_terms': bad_terms, 'bit': bit})
        col = {'index': index, 'kind': kind, 'source_index': source_index,
               'ideal': rows(H), 'good_ideal': rows(good), 'norm': str(N), 'cyclic': bool(cyclic),
               'bad_parts': [{'p': p, 'j': j, 'exponent': e} for p, j, P, e in bad], 'entries': entries}
        save(path, col)
        columns.append(col)
        print('COLUMN',index,'COMPLETE',sum('bit' in e for e in entries),'/',len(entries),flush=True)
    complete = all(all('bit' in e for e in c['entries']) for c in columns)
    result = {'schema': 'curve302.relative-ideal-anatomy.v1', 'bindings': bindings,
              'protocol': protocol, 'characters': character_packet, 'columns': columns, 'complete': complete}
    if complete:
        M = matrix(GF(2), [[e['bit'] for e in c['entries']] for c in columns]).transpose()
        rg, rv, ru = map(int, [M[:, :8].rank(), M[:, 8:].rank(), M.rank()])
        result.update(artin_matrix=[[int(x) for x in row] for row in M],
                      generic_detected_rank=rg, strict_detected_rank=rv, union_detected_rank=ru,
                      detected_relative_dimension=ru-rg,
                      joint_kernel_words=[[int(x) for x in row] for row in M.right_kernel().basis()])
        print('RANKS',rg,rv,ru,'RELATIVE',ru-rg,'KERNEL',result['joint_kernel_words'],flush=True)
    save(DEST, result)


if __name__ == '__main__':
    compute()
