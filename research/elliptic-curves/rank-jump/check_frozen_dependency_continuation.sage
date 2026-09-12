#!/usr/bin/env sage
"""Bounded compatibility check for one complete frozen principal dependency.

The definition is written before the modular test. Failure retains that
definition with UNKNOWN. This never discovers new atoms or varies a budget.
"""
import hashlib
import json
import resource
import time
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, gcd, lcm, prime_range

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT/'artifacts/generated-results/elliptic-curves'
BANK = GEN/'marked_two_class_propagation_v1'
OUT = GEN/'frozen_dependency_continuation_v1'
PROTOCOL = Path(__file__).with_name('FROZEN_DEPENDENCY_CONTINUATION_PROTOCOL.json')


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def coefficients(poly):
    return list(map(str, poly.list())) or ['0']


def main():
    protocol = read(PROTOCOL)
    resource.setrlimit(resource.RLIMIT_CPU, (protocol['cpu_seconds'], protocol['cpu_seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (protocol['address_space_bytes'], protocol['address_space_bytes']))
    started = time.process_time()
    OUT.mkdir(exist_ok=False)
    save(OUT/'protocol.json', protocol)
    source_paths = [BANK/'geometry.json', BANK/'frozen-block.json',
                    GEN/'rank_jump_reference_additional_strict_verification_v1.json',
                    GEN/'rank_jump_reference_strict_class_construction_inputs_v1.json', PROTOCOL]
    bindings = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    family = read(source_paths[0])['family']
    block, old, generic = map(read, source_paths[1:4])
    for path in source_paths[2:4]:
        assert bindings[str(path.relative_to(ROOT))] == block['source_bindings'][str(path.relative_to(ROOT))]
    labels = old['class_representatives'][0]['factor_labels']
    atom_ids = [a['index'] for a in labels if a['kind'] == 'projected_atom']
    corrections = [a['index'] for a in labels if a['kind'] == 'generic']
    assert len(set(atom_ids)) == len(atom_ids) == protocol['required_projected_atom_count'] == 1676
    assert corrections == [3, 5, 6, 7, 9, 10]
    atoms = {a['index']: a for a in old['atoms']}
    chosen = min(atom_ids)
    assert chosen == 2 and all(QQ(atoms[i]['alpha_ascending'][2]) == 0 for i in atom_ids)
    R = PolynomialRing(QQ, 't')
    t = R.gen()
    t0, u = QQ(3)/17, QQ(block['generic_control_transport']['short_scale_u'])
    assert u == QQ(289)/2
    A = R(family['A_coefficients_low_to_high'])
    B = R(family['B_coefficients_low_to_high'])
    c1, c0 = (u**4*A+27)/81, (u**6*B+3*u**4*A+27)/729
    f0 = [c0(t0), c1(t0), QQ(1), QQ(1)]
    assert f0 == list(map(QQ, old['equation']['cubic_ascending'])) == list(map(QQ, generic['cubic_ascending']))
    gamma_packet = []
    for section, gamma in zip(family['sections'], generic['generic_classes']):
        def section_poly(key):
            v = section[key]
            den = R(v['denominator_coefficients_low_to_high'])
            assert den.degree() == 0 and den
            return R(v['numerator_coefficients_low_to_high'])/den[0]
        x, y = section_poly('X'), section_poly('Y')
        assert y*y == x*x*x+A*x+B
        z, w = (u*u*x-3)/9, u**3*y/27
        assert w*w == z**3+z*z+c1*z+c0
        g = list(map(QQ, gamma['beta_ascending']))
        scale = -g[1]
        assert scale > 0 and scale.is_square() and g == [scale*z(t0), -scale, QQ(0)]
        assert scale**3*w(t0)**2 == QQ(gamma['norm'])
        gamma_packet.append({'index': int(section['basis_index']), 'scale': str(scale),
                             'z_ascending': coefficients(z), 'w_ascending': coefficients(w)})
    norms = {}
    for i in atom_ids:
        a, b, zero = map(QQ, atoms[i]['alpha_ascending'])
        assert b and not zero
        # Norm(a+b*theta)=-b^3*f_t(-a/b).
        norms[i] = a**3-a*a*b+c1*a*b*b-c0*b**3
        assert norms[i](t0) == QQ(atoms[i]['norm']) != 0
    n = norms[chosen]
    ZR = PolynomialRing(ZZ, 't')
    integer_n = ZR(n*lcm([x.denominator() for x in n.list()]))
    primitive = integer_n // gcd(integer_n.list())
    if primitive.leading_coefficient() < 0:
        primitive = -primitive
    constant = n.leading_coefficient()/primitive.leading_coefficient()
    assert n == constant*R(primitive) and primitive.degree() == 12
    definition = {
        'schema': 'marked-dependency.coefficientwise-definition.v1', 'bindings': bindings,
        'column': 6, 'control_parameter': str(t0), 'short_scale_u': str(u),
        'field_cubic_ascending_in_theta': [coefficients(c0), coefficients(c1), ['1'], ['1']],
        'field_specialization_cubic_ascending': list(map(str, f0)),
        'coordinate_transport': 'z=(u^2*Xcompact-3)/9=4*x_original; w=u^3*Ycompact/27=8*y_original+4*x_original',
        'element_formula': 'alpha_i(t)=a_i+b_i*theta_t; pi_i(t)=Norm(alpha_i(t))*alpha_i(t); gamma_j(t)=scale_j*(z_j(t)-theta_t)',
        'dependency_formula': 'beta_6(t)=product of all listed pi_i(t), times gamma_3*gamma_5*gamma_6*gamma_7*gamma_9*gamma_10',
        'norm_square_root_formula': 'product_i Norm(alpha_i(t))^2 * product_j scale_j^(3/2)*w_j(t), with j the fixed generic corrections',
        'factor_labels': labels, 'generic_correction_indices': corrections,
        'atoms': [{'index': i, 'alpha_ascending': atoms[i]['alpha_ascending'],
                   'norm_at_control': atoms[i]['norm']} for i in atom_ids],
        'generic_sections': gamma_packet, 'distinguished_atom': chosen,
        'distinguished_norm_primitive_ascending': coefficients(primitive),
        'distinguished_norm_constant': str(constant),
        'control_compaction_circuit': next(c for c in block['circuits'] if c['column'] == 6),
        'boundary': 'The complete pre-compaction word specializes factor-by-factor to the original certificate. Its previous square-equivalence circuit to compact column 6 is inherited, not recomputed. No identification of specialized prime ideals across t is made.'
    }
    save(OUT/'definition.json', definition)
    status = {'schema': 'marked-dependency.continuation-check.v1', 'status': 'UNKNOWN',
              'definition_sha256': hashlib.sha256((OUT/'definition.json').read_bytes()).hexdigest()}
    try:
        p = ZZ(protocol['coprimality_prime'])
        assert p.is_prime()
        Rp = PolynomialRing(GF(p), 't')
        N = Rp(primitive)
        assert N.degree() == primitive.degree()
        discriminant = c1*c1-4*c1**3-4*c0-27*c0*c0+18*c1*c0
        product = (N.derivative()*Rp(discriminant)) % N
        for i in atom_ids:
            if i != chosen:
                product = product * Rp(norms[i]) % N
        # All sixteen generic corrections are units at these divisors.
        for g in gamma_packet:
            product = product * Rp(g['w_ascending']) % N
        g, left, right = N.xgcd(product)
        assert g == 1
        irreducible_prime = next(q for q in prime_range(5, 102)
                                if all(x.denominator() % q for x in f0)
                                and PolynomialRing(GF(q), 'z')(f0).is_irreducible())
        status.update({
            'status': 'ODD_HORIZONTAL_VALUATION_OBSTRUCTION',
            'prime': int(p), 'norm_degree': int(N.degree()), 'distinguished_atom': chosen,
            'field_specialization_irreducibility_prime': int(irreducible_prime),
            'primitive_norm_mod_p': coefficients(N), 'joint_product_remainder_mod_p': coefficients(product),
            'bezout_norm_multiplier_mod_p': coefficients(left),
            'bezout_product_multiplier_mod_p': coefficients(right),
            'joint_product_definition': "N'(t)*disc_theta(f_t)*product_(i != 2) Norm(alpha_i(t))*product_(j=0..15) w_j(t), reduced modulo N(t)",
            'checked_atoms': len(atom_ids), 'checked_other_atom_coprimalities': len(atom_ids)-1,
            'checked_generic_sections': len(gamma_packet),
            'geometric_valuation_pattern': [2, 1, 1],
            'conclusion': 'The specified beta_6(t) has odd valuation on the two residual sheets over every irreducible factor of the distinguished norm. These are good-reduction divisors. Thus this word is not a rational-point Kummer class over Q(t), even after a square or generic-point correction.',
            'boundary': 'Only this coefficientwise continuation is obstructed. The soluble strict class at 3/17, other continuations, arithmetic incidence conditions and ramified base changes remain unaffected.'})
    except (AssertionError, ValueError, ZeroDivisionError, StopIteration) as error:
        status['failure'] = type(error).__name__+': '+str(error)
    status['cpu_seconds'] = time.process_time()-started
    status['peak_rss_bytes'] = 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    save(OUT/'result.json', status)
    print(json.dumps({k: v for k, v in status.items() if 'mod_p' not in k}, indent=2))


if __name__ == '__main__':
    main()
