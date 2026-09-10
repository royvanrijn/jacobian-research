#!/usr/bin/env sage-python
"""Read-only lattice, cyclic-ring and Hensel replay of the 302 descent block.

The half ideals are checked without idealred/idealadd/idealpow construction.
Integer Jacobi and simple-root Hensel calculations replace the Artin producer.
Existing exact local-filtration and ordinary-character checkers are rerun.
"""
import argparse
import hashlib
import json
import runpy
import sys
from pathlib import Path
from sage.all import AA, QQ, ZZ, GF, NumberField, PolynomialRing, matrix, vector, pari, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
sys.path.insert(0, str(ROOT/'elliptic-curves/rank-jump'))
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
from verify_strict_half_ideals import check_lattice
from verify_half_ideal_artin import jacobi
import icarm_curve302 as curve


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def require_bindings(d, embedded=None):
    for name, expected in d['bindings'].items():
        if embedded and name in embedded:
            raw = (json.dumps(embedded[name], indent=2, sort_keys=True)+'\n').encode()
            assert hashlib.sha256(raw).hexdigest() == expected
        else:
            assert sha(ROOT/name) == expected, name


def hensel_bit(f, beta, p, root, valuation):
    assert f(root) % p == 0 and f.derivative()(root) % p
    modulus = ZZ(p)
    for _ in range(valuation):
        step = (-ZZ(f(root)//modulus) * ZZ(f.derivative()(root)).inverse_mod(p)) % p
        root += modulus*step
        modulus *= p
        assert f(root) % modulus == 0
    value = ZZ(beta(root)) % modulus
    assert value.valuation(p) == valuation and valuation % 2 == 0
    return int(pow(int(value//ZZ(p)**valuation), (p-1)//2, p) == p-1)


def verify(negative_controls=False):
    manifest = ART/'curve302_descent_anatomy_manifest_v1.json'
    if manifest.exists():
        require_bindings(read(manifest))
    d = read(ART/'curve302_descent_anatomy_v1.json')
    last = read(ART/'curve302_descent_remaining_class_v1.json')
    h, old = d['half_ideal_packet'], d['initial_artin_packet']
    require_bindings(d, {'artifacts/local/elliptic-curves/curve302-descent-anatomy-v1/result.json': old,
                         'artifacts/local/elliptic-curves/curve302-descent-anatomy-v1/half-ideals.json': h})
    require_bindings(last)
    require_bindings(h)
    assert h['bindings'] == old['bindings']
    filtration_path = ART/'curve302_recovered_quotient_local_filtration_v1.json'
    previous = runpy.run_path(str(ROOT/'elliptic-curves/cas/verify_curve302_recovered_quotient_local_filtration.sage'))
    assert previous['build']() == read(filtration_path)
    unramified = runpy.run_path(str(ROOT/'elliptic-curves/rank-jump/verify_generic_only_class_anchors.py'))
    assert unramified['compute']() == read(ART/'rank_jump_generic_only_class_anchors_verification_v1.json')
    a = read(ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json')
    R = PolynomialRing(QQ, 'z')
    f = R(h['cubic_ascending'])
    assert f.is_irreducible() and f == R(a['cubic_ascending'])
    remaining = abs(ZZ(f.discriminant()))
    for p in a['S_finite']:
        assert ZZ(p).is_prime(proof=True)
        while remaining % p == 0:
            remaining //= p
    assert remaining == 1
    nf = pari.nfinit([pari(f), a['S_finite']])
    assert list(map(str, nf.nf_get_zk())) == h['maximal_order_basis']
    K = NumberField(f, 'theta')
    theta = K.gen()
    basis = [K(R(b)) for b in h['maximal_order_basis']]
    B = matrix(QQ, [list(b) for b in basis]).transpose()
    Bi = B.inverse()
    assert matrix(QQ, 3, 3, lambda i, j: (basis[i]*basis[j]).trace()).det() == ZZ(h['field_discriminant'])
    roots = f.roots(AA, multiplicities=False)
    assert len(roots) == 3
    gs, ys = [], []
    assert len(h['points']) == len(curve.POINTS) == 31
    for P, public in zip(h['points'], curve.POINTS):
        aa, bb, dd = (ZZ(P[k]) for k in ['a', 'b', 'd'])
        x, y = map(QQ, public)
        assert (aa/dd**2, bb/dd**3) == (4*x, 8*y+4*x+4)
        assert curve.on_curve(public) and aa.gcd(dd) == bb.gcd(dd) == 1
        g = aa-dd*dd*theta
        assert g.norm() == bb*bb
        gs.append(g)
        ys.append(abs(bb))
    assert h['words'] == read(filtration_path)['local_filtration_mod_2']['strict_kernel_public_words']
    betas = [prod(g for bit, g in zip(w, gs) if bit) for w in h['words']]
    Hs = [matrix(QQ, H) for H in h['half_ideals']]
    for word, beta, H in zip(h['words'], betas, Hs):
        norm = prod(y for bit, y in zip(word, ys) if bit)
        assert all(R(list(beta))(r) > 0 for r in roots)
        check_lattice(K, basis, H, beta, norm)
    if negative_controls:
        for alter in [2*Hs[0], matrix(QQ, Hs[0])]:
            if alter == Hs[0]:
                alter[0, 2] += 1
            try:
                check_lattice(K, basis, alter, betas[0], Hs[0].det())
            except AssertionError:
                pass
            else:
                raise AssertionError('Altered half ideal accepted')

    def multiplication(x):
        return matrix(QQ, [list(Bi*vector(QQ, list(x*b))) for b in basis]).transpose()

    def same_lattice(X, Y):
        C = X.inverse()*Y
        assert all(v in ZZ for v in C.list()) and abs(C.det()) == 1

    def cyclic_ring(H, N):
        assert H[0, 0] == H.det() == N and H[1, 1] == H[2, 2] == 1
        assert all(H[i, j] == 0 for i in range(3) for j in range(i))
        assert N > 0 and all(N % p for p in a['S_finite'])
        residues = vector(QQ, [1, -H[0, 1], -H[0, 2]])
        for i in range(3):
            for j in range(3):
                cs = Bi*vector(QQ, list(basis[i]*basis[j]))
                assert all(c in ZZ for c in cs)
                assert (residues.dot_product(cs)-residues[i]*residues[j]) % N == 0
        return residues

    def value(beta, residues, N):
        cs = Bi*vector(QQ, list(beta))
        assert all(c in ZZ for c in cs)
        return ZZ(residues.dot_product(cs)) % N

    def check_removal(reduced, good, terms):
        target = pari(good)
        for term in terms:
            p, j, e = term['p'], term['j'], term['exponent']
            P = pari.idealprimedec(nf, p)[j]
            assert int(pari.idealval(nf, pari(reduced), P)) == e > 0
            target = pari.idealmul(nf, target, pari.idealpow(nf, P, e))
        same_lattice(reduced, matrix(QQ, target))

    repairs = {(r['column'], r['character']): r for r in d['repairs']}
    cols = []
    for j, col in enumerate(old['columns']):
        assert col['index'] == j
        reduced = matrix(QQ, col['reduced_ideal'])
        alpha = K(list(map(QQ, col['multiplier'])))
        same_lattice(Hs[j], multiplication(alpha)*reduced)
        good, N = matrix(QQ, col['coprime_ideal']), ZZ(col['norm'])
        check_removal(reduced, good, col['removed_S_factors'])
        residues = cyclic_ring(good, N)
        bits = []
        for i, (beta, entry) in enumerate(zip(betas, col['entries'])):
            val = value(beta, residues, N)
            if 'bit' in entry:
                symbol = jacobi(int(val), int(N))
                assert symbol in [-1, 1] and str(val) == entry['residue'] and symbol == entry['symbol']
                bit = int(symbol == -1)
                assert bit == entry['bit']
            else:
                r = repairs[j, i]
                p, e = r['prime'], r['exponent']
                assert val.gcd(N) == p and N.valuation(p) == e
                cofactor = N//ZZ(p)**e
                assert str(cofactor) == r['cofactor']
                residue = value(beta, residues, cofactor)
                assert str(residue) == r['cofactor_residue']
                symbol = jacobi(int(residue), int(cofactor))
                assert symbol == r['cofactor_symbol'] and symbol in [-1, 1]
                assert int(value(theta, residues, p)) == r['root']
                local = hensel_bit(f, R(list(beta)), p, ZZ(r['root']), r['valuation'])
                assert local == r['local_bit']
                bit = int(symbol == -1) ^ ((e % 2)*local)
                assert bit == r['bit']
            bits.append(bit)
        cols.append(bits)
    A = matrix(GF(2), cols).transpose()
    assert [list(map(int, row)) for row in A] == d['artin_matrix']
    assert A.rank() == d['artin_rank'] == 9
    kernel = A.right_kernel().basis_matrix()
    assert [list(map(int, row)) for row in kernel] == d['right_kernel_words'] == [last['strict_word']]

    # Verify the kernel product ideal by its square identity, avoiding the
    # producer's product-of-ideals construction.
    beta = prod(b for bit, b in zip(last['strict_word'], betas) if bit)
    alpha = K(list(map(QQ, last['multiplier'])))
    H = matrix(QQ, last['reduced_ideal'])
    norm = prod(H.det() for bit, H in zip(last['strict_word'], Hs) if bit)
    recovered = multiplication(alpha)*H
    check_lattice(K, basis, recovered, beta, norm)
    good, N = matrix(QQ, last['coprime_ideal']), ZZ(last['norm'])
    terms = last['entries'][0]['S_terms']
    check_removal(H, good, terms)
    residues = cyclic_ring(good, N)
    generic = [K(R(c['beta_ascending'])) for c in a['generic_classes']]
    anchor = read(ART/'rank_jump_generic_only_class_anchors_v1.json')['cases'][1]
    assert [e['mask'] for e in last['entries']] == anchor['generic_coefficient_masks']
    tested = []
    for entry in last['entries']:
        mask = int(entry['mask'])
        beta = prod(g for i, g in enumerate(generic) if mask >> i & 1)
        val = value(beta, residues, N)
        symbol = jacobi(int(val), int(N))
        assert str(val) == entry['residue'] and symbol == entry['symbol'] and symbol in [-1, 1]
        assert [(t['p'], t['j'], t['exponent']) for t in entry['S_terms']] == [(7, 0, 1)]
        P = pari.idealprimedec(nf, 7)[0]
        # The removed ideal is unramified of degree one; its simple root is
        # recovered directly from its lattice and checked by Hensel lifting.
        PH = matrix(QQ, pari.idealhnf(nf, P))
        pres = vector(QQ, [1, -PH[0, 1], -PH[0, 2]])
        root = value(theta, pres, ZZ(7))
        poly = R(list(beta))
        valuation = int(pari.idealval(nf, pari.Mod(pari(poly), pari(f)), P))
        local = hensel_bit(f, poly, 7, root, valuation)
        assert local == entry['S_terms'][0]['bit']
        bit = int(symbol == -1) ^ local
        assert bit == entry['bit']
        tested.append(bit)
    assert tested == [1, 1, 0, 0, 0, 0] and last['remaining_class_detected']
    assert last['ordinary_half_ideal_rank_interval'] == [10, 10]
    assert last['unit_kernel_dimension_interval'] == [0, 0]
    # Derivative reciprocity removes the one boundary bit left by the local
    # point product. The sign test is exact in the algebraic real field.
    derivative = -f.discriminant()*f.derivative()
    assert [int(derivative(r) < 0) for r in sorted(roots)] == [1, 0, 1]
    point_signs = set(tuple(int((4*QQ(x)-r) < 0) for r in sorted(roots)) for x, y in curve.POINTS)
    assert point_signs == {(0, 0, 0), (0, 1, 1)}
    loc = read(filtration_path)['local_filtration_mod_2']
    assert loc['local_product_dimension'] == 22 and loc['joint_local_dimensions']['D'] == 21
    print('PASS 10 half-ideal lattice identities; 100 Artin entries; rank 9 plus detected kernel;')
    print('PASS ordinary/narrow half-ideal dimension 10, unit kernel 0; Selmer dimension 21+c_S.')
    if negative_controls:
        print('PASS two altered half ideals rejected.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--negative-controls', action='store_true')
    verify(p.parse_args().negative_controls)
