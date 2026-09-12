#!/usr/bin/env python3
"""Independent SymPy replay of the fixed complete-word obstruction certificate.

Run with `sage -python` to use the repository's bundled SymPy. No producer
module is imported. No ideals, carriers, points or alternative primes are
searched. The certificate supplies a finite-field Bezout identity.
"""
import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = GEN/'frozen_dependency_continuation_v1'


def read(p):
    return json.loads(p.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Write a new independent-replay receipt exclusively.')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    started = time.process_time()
    definition, result = read(OUT/'definition.json'), read(OUT/'result.json')
    assert hashlib.sha256((OUT/'definition.json').read_bytes()).hexdigest() == result['definition_sha256']
    sources = {}
    for path, digest in definition['bindings'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
        sources[path] = read(ROOT/path)
    old = sources['artifacts/generated-results/elliptic-curves/rank_jump_reference_additional_strict_verification_v1.json']
    original_gamma = sources['artifacts/generated-results/elliptic-curves/rank_jump_reference_strict_class_construction_inputs_v1.json']
    family = sources['artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/geometry.json']['family']
    block = sources['artifacts/generated-results/elliptic-curves/marked_two_class_propagation_v1/frozen-block.json']
    assert definition['factor_labels'] == old['class_representatives'][0]['factor_labels']
    ids = [a['index'] for a in definition['factor_labels'] if a['kind'] == 'projected_atom']
    corrections = [a['index'] for a in definition['factor_labels'] if a['kind'] == 'generic']
    assert len(ids) == len(set(ids)) == 1676 and corrections == [3, 5, 6, 7, 9, 10]
    assert corrections == definition['generic_correction_indices']
    assert definition['column'] == 6 and definition['distinguished_atom'] == min(ids) == 2
    old_atoms = {a['index']: a for a in old['atoms']}
    assert definition['atoms'] == [{'index': i, 'alpha_ascending': old_atoms[i]['alpha_ascending'],
                                     'norm_at_control': old_atoms[i]['norm']} for i in ids]
    assert definition['control_compaction_circuit'] == next(c for c in block['circuits'] if c['column'] == 6)
    t, theta = sp.symbols('t theta')
    def qp(cs):
        return sp.Poly.from_list([sp.Rational(c) for c in reversed(cs)], t, domain=sp.QQ)
    u, t0 = sp.Rational(definition['short_scale_u']), sp.Rational(definition['control_parameter'])
    assert u == sp.Rational(289, 2) and t0 == sp.Rational(3, 17)
    A, B = [qp(family[k+'_coefficients_low_to_high']) for k in ('A', 'B')]
    c0, c1, one, other_one = map(qp, definition['field_cubic_ascending_in_theta'])
    assert one == other_one == qp(['1'])
    assert c1 == (u**4*A+27).mul_ground(sp.Rational(1, 81))
    assert c0 == (u**6*B+3*u**4*A+27).mul_ground(sp.Rational(1, 729))
    f0 = [c0.eval(t0), c1.eval(t0), sp.Rational(1), sp.Rational(1)]
    assert f0 == list(map(sp.Rational, original_gamma['cubic_ascending']))
    assert f0 == list(map(sp.Rational, definition['field_specialization_cubic_ascending']))
    irred_p = result['field_specialization_irreducibility_prime']
    assert sp.isprime(irred_p)
    assert sp.Poly.from_list(list(reversed(f0)), theta, modulus=irred_p).is_irreducible
    # A monic factorization over Q(t) would specialize to a monic factorization
    # of f_0; irreducibility at the control proves K_t is a field.
    N = qp(definition['distinguished_norm_primitive_ascending'])
    assert N.degree() == 12 and N.LC() > 0
    assert all(c.q == 1 for c in N.all_coeffs()) and N.content() == 1
    a, b, zero = map(sp.Rational, old_atoms[2]['alpha_ascending'])
    f = theta**3+theta**2+c1.as_expr()*theta+c0.as_expr()
    exact_norm = sp.Poly(sp.resultant(f, a+b*theta, theta), t, domain=sp.QQ)
    assert exact_norm == sp.Rational(definition['distinguished_norm_constant'])*N
    p = result['prime']
    assert p == 1000003 and sp.isprime(p)
    def fp(q):
        if not isinstance(q, sp.Poly):
            q = qp(q)
        return sp.Poly.from_list([int(c.p)*pow(int(c.q), -1, p) % p for c in q.all_coeffs()], t, modulus=p)
    np, cp0, cp1 = fp(N), fp(c0), fp(c1)
    assert np.degree() == 12 and np == fp(result['primitive_norm_mod_p'])
    disc = c1*c1-4*c1**3-4*c0-27*c0*c0+18*c1*c0
    product = (np.diff()*fp(disc)).rem(np)
    for i in ids:
        aa, bb, zz = map(sp.Rational, old_atoms[i]['alpha_ascending'])
        assert not zz and bb
        norm_at_control = aa**3-aa*aa*bb+c1.eval(t0)*aa*bb*bb-c0.eval(t0)*bb**3
        assert norm_at_control == sp.Rational(old_atoms[i]['norm']) != 0
        if i != 2:
            aa, bb = [int(q.p)*pow(int(q.q), -1, p) % p for q in (aa, bb)]
            mod_norm = cp1.mul_ground(aa*bb*bb)-cp0.mul_ground(bb**3)+(aa**3-aa*aa*bb)
            product = (product*mod_norm).rem(np)
    assert len(definition['generic_sections']) == 16
    for j, packet in enumerate(definition['generic_sections']):
        assert packet['index'] == j
        section, old_gamma = family['sections'][j], original_gamma['generic_classes'][j]
        coords = []
        for name in ('X', 'Y'):
            den = qp(section[name]['denominator_coefficients_low_to_high'])
            assert den.degree() == 0 and den.LC()
            coords.append(qp(section[name]['numerator_coefficients_low_to_high']).mul_ground(1/den.LC()))
        x, y = coords
        assert y*y == x*x*x+A*x+B
        z, w, s = qp(packet['z_ascending']), qp(packet['w_ascending']), sp.Rational(packet['scale'])
        assert z == (u*u*x-3).mul_ground(sp.Rational(1, 9)) and w == y.mul_ground(u**3/27)
        assert w*w == z**3+z*z+c1*z+c0
        assert s > 0 and sp.sqrt(s).is_Rational
        assert [s*z.eval(t0), -s, 0] == list(map(sp.Rational, old_gamma['beta_ascending']))
        assert s**3*w.eval(t0)**2 == sp.Rational(old_gamma['norm'])
        product = (product*fp(w)).rem(np)
    assert product == fp(result['joint_product_remainder_mod_p'])
    left, right = [fp(result[k]) for k in ('bezout_norm_multiplier_mod_p', 'bezout_product_multiplier_mod_p')]
    assert left*np+right*product == fp(['1'])
    assert result['status'] == 'ODD_HORIZONTAL_VALUATION_OBSTRUCTION'
    receipt = {
        'schema': 'marked-dependency.independent-replay.v1', 'status': 'VERIFIED',
        'verifier': str(Path(__file__).relative_to(ROOT)),
        'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'definition_sha256': result['definition_sha256'],
        'result_sha256': hashlib.sha256((OUT/'result.json').read_bytes()).hexdigest(),
        'sympy_version': sp.__version__, 'atoms_checked': len(ids), 'generic_sections_checked': 16,
        'norm_identity': 'Independent symbolic resultant for atom 2',
        'coprimality_certificate': 'Complete product recomputed; supplied Bezout identity equals 1 over F_1000003[t]',
        'cpu_seconds': time.process_time()-started,
        'peak_rss_bytes': 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'boundary': 'Checks the complete specified family-element continuation and its good-divisor valuation obstruction. It inherits the earlier number-field strictness and compaction circuit; it does not replay ideal factorization, search for points or conclude anything about other continuations.'}
    if args.write:
        with (OUT/'independent-replay.json').open('x') as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
