#!/usr/bin/env sage-python
"""Replay existing results via integral-monic PARI fields, never fill UNKNOWN.

Use a bounded outer timeout. Known factors are replay witnesses, not hints given
to the original profiling workers. This shares Sage/PARI with the worker and is
not a separate full descent implementation.
"""
import argparse
from pathlib import Path
from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, pari
import historical_external_arithmetic as lane


def verify(out):
    plan = lane.check(out, complete=True)
    records, bindings = [], {}
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    for row in plan['rows']:
        key = row['curve_key']
        paths = {mode: out / mode / f'{key}.json' for mode in ('inputs', 'base', 'local')}
        bindings.update({str(p.relative_to(out)): lane.sha(p) for p in paths.values()})
        base, local = lane.read(paths['base']), lane.read(paths['local'])
        record = {'curve_key': key, 't': row['t'], **lane.labels(lane.HISTORICAL),
                  'base': base['status'], 'local': local['status']}
        if base['status'] != 'PASS':
            assert base['status'].startswith('UNKNOWN')
            assert local['status'].startswith('UNKNOWN'), 'LOCAL PASS without BASE requires explicit replay'
            records.append(record)
            continue
        original = EllipticCurve(QQ, row['ainvs'])
        E = EllipticCurve(QQ, list(map(QQ, base['minimal_ainvs'])))
        iso = original.isomorphism_to(E)
        assert original.discriminant() == iso.tuple()[0]**12 * E.discriminant()
        assert E.minimal_model().ainvs() == E.ainvs()
        assert str(E.discriminant()) == base['minimal_discriminant']
        assert E.c4()**3 - E.c6()**2 == 1728*E.discriminant()
        f = R(list(map(QQ, base['two_division_cubic'])))
        assert f.monic() == (4*x**3 + E.b2()*x*x + 2*E.b4()*x + E.b6()).monic()
        assert str(f.discriminant()) == base['two_division_cubic_discriminant']
        assert f.is_irreducible() and base['rational_2torsion_rank'] == 0
        if local['status'] == 'PASS':
            factors = [(ZZ(p), int(e)) for p, e in local['minimal_discriminant_factorization']]
            assert len({p for p, e in factors}) == len(factors)
            product = ZZ(1)
            for p, e in factors:
                assert p.is_prime(proof=True) and e > 0
                product *= p**e
            assert product == abs(E.discriminant())
            primes = [p for p, e in factors]
            pari.addprimes(primes)
            # X=4x changes the worker's rational monic equation to an integral one.
            F = x**3 + E.b2()*x*x + 8*E.b4()*x + 16*E.b6()
            nf = pari.nfinit([pari(F), sorted(set([ZZ(2)] + primes))])
            assert pari.nfcertify(nf) == []
            assert str(nf.disc()) == local['field_discriminant']
            assert list(map(int, nf.nf_get_sign())) == local['field_signature']
            d = abs(ZZ(nf.disc()))
            ramified = list(map(ZZ, local['field_ramified_primes']))
            assert len(set(ramified)) == len(ramified) == local['field_ramified_prime_count']
            for p in ramified:
                assert p.is_prime(proof=True) and d % p == 0
                while d % p == 0:
                    d //= p
            assert d == 1
            assert {ZZ(r['p']) for r in local['local_data']} == set(primes)
            assert len(local['local_data']) == len(primes)
            phi_m, phi_a, additive_splits = [], [], []
            conductor, root = ZZ(1), -1
            for r in local['local_data']:
                p = ZZ(r['p'])
                ld = pari.elllocalred(pari(E), p)
                exponent = int(ld[0])
                assert exponent == r['conductor_exponent']
                assert E.discriminant().valuation(p) == r['v_delta']
                assert str(E.local_data(p).kodaira_symbol()) == r['kodaira']
                local_root = int(pari.ellrootno(pari(E), p))
                assert local_root == r['root_number']
                root *= local_root
                conductor *= p**exponent
                if exponent == 1:
                    assert r['kind'] == 'multiplicative'
                    if r['v_delta'] % 2 == 0:
                        phi_m.append(str(p))
                else:
                    assert r['kind'] == 'additive'
                    phi_a.append(str(p))
                    n = len(pari.idealprimedec(nf, p))
                    assert n == r['cubic_prime_count']
                    additive_splits.append(n)
            assert phi_m == local['phi_m'] and phi_a == local['phi_a']
            assert len(phi_m) == local['phi_m_count'] and len(phi_a) == local['phi_a_count']
            assert additive_splits == local['additive_split_counts']
            u = 1 if E.discriminant() < 0 else 2
            n = len(phi_m) + sum(v-1 for v in additive_splits)
            assert (u, n, u+n) == (local['bk_u'], local['bk_n'], local['bk_local_term'])
            assert str(conductor) == local['conductor']
            assert root == local['root_number'] == int(pari.ellrootno(pari(E)))
            record.update(bk_local_term=u+n, certified_maximal_order='PASS', local_replay='PASS')
        else:
            assert local['status'].startswith('UNKNOWN')
        records.append(record)
        print(f"HISTORICAL_REPLAY|t={row['t']}|BASE=PASS|LOCAL={local['status']}", flush=True)
    lane.save(out / 'ARITHMETIC_REPLAY.json', {
        'schema': 'elliptic-curves.historical-external-arithmetic-replay.v1',
        'status': 'PASS_EXISTING_CHECKPOINT_REPLAY', 'rows': records,
        'checkpoint_sha256': bindings, 'checker_sha256': lane.sha(Path(__file__)),
        'boundary': 'Exact BASE identities, certified integral-monic maximal orders and reductions for completed LOCAL rows. Shared Sage/PARI, not independent full descent. UNKNOWN unchanged. No point searches or class groups.'})
    lane.check(out, complete=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=lane.DEFAULT_OUTPUT)
    args = ap.parse_args()
    verify(args.output.resolve())
