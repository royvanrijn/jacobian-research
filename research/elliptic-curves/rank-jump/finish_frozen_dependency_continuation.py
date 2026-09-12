#!/usr/bin/env python3
"""Finish the fixed compatibility calculation using Fraction polynomial arithmetic.

Sage's setup failed before the modular test. This uses the original atom,
prime, complete dependency and memory bound. It writes new evidence files
exclusively and performs no search for a replacement witness.
"""
import hashlib
import json
import math
import resource
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT/'artifacts/generated-results/elliptic-curves'
BANK = GEN/'marked_two_class_propagation_v1'
OUT = GEN/'frozen_dependency_continuation_v1'
PROTOCOL = Path(__file__).with_name('FROZEN_DEPENDENCY_CONTINUATION_PROTOCOL.json')


def read(p):
    return json.loads(p.read_text())


def save(p, obj):
    with p.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a or [0]


def add(a, b):
    c = [0]*max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return trim(c)


def scale(a, b):
    return trim([b*x for x in a])


def mul(a, b):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def ev(a, t):
    s = 0
    for x in reversed(a):
        s = s*t+x
    return s


def residue(a, p):
    return trim([Q(x).numerator*pow(Q(x).denominator, -1, p) % p for x in a])


def divmodp(a, b, p):
    a, b = residue(a, p), residue(b, p)
    assert b != [0]
    q = [0]*max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        k, c = len(a)-len(b), a[-1]*pow(b[-1], -1, p) % p
        q[k] = c
        for i, x in enumerate(b):
            a[i+k] = (a[i+k]-c*x) % p
        a = trim(a)
    return trim(q), a


def xgcdp(a, b, p):
    x, xx, y, yy = [1], [0], [0], [1]
    while b != [0]:
        q, r = divmodp(a, b, p)
        a, b = b, r
        x, xx = xx, residue(add(x, scale(mul(q, xx), -1)), p)
        y, yy = yy, residue(add(y, scale(mul(q, yy), -1)), p)
    c = pow(a[-1], -1, p)
    return residue(scale(a, c), p), residue(scale(x, c), p), residue(scale(y, c), p)


def norm(a, b, c0, c1):
    return add([a**3-a*a*b], add(scale(c1, a*b*b), scale(c0, -b**3)))


def prime(p):
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p)+1))


def main():
    protocol = read(PROTOCOL)
    assert read(OUT/'protocol.json') == protocol
    resource.setrlimit(resource.RLIMIT_CPU, (protocol['cpu_seconds'],)*2)
    resource.setrlimit(resource.RLIMIT_AS, (protocol['address_space_bytes'],)*2)
    started = time.process_time()
    assert not (OUT/'definition.json').exists() and not (OUT/'result.json').exists()
    failed_source = Path(__file__).with_name('check_frozen_dependency_continuation.sage')
    save(OUT/'setup-failure.json', {
        'stage': 'Sage polynomial-ring coercion while preparing the distinguished norm',
        'failure': 'ImportError: failed to map segment from shared object under the 1 GiB address-space bound',
        'mathematical_modular_test_reached': False,
        'failed_source': str(failed_source.relative_to(ROOT)),
        'failed_source_sha256': hashlib.sha256(failed_source.read_bytes()).hexdigest(),
        'repair': 'Use plain Python for the same fixed arithmetic, with the same 1 GiB address-space ceiling. No alternate atom, prime, word or parameter is tested.',
        'prior_setup_cpu_seconds': 'UNMEASURED; excluded from the component timing below'
    })
    paths = [BANK/'geometry.json', BANK/'frozen-block.json',
             GEN/'rank_jump_reference_additional_strict_verification_v1.json',
             GEN/'rank_jump_reference_strict_class_construction_inputs_v1.json', PROTOCOL]
    bindings = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    family, block, old, generic = read(paths[0])['family'], read(paths[1]), read(paths[2]), read(paths[3])
    for path in paths[2:4]:
        assert bindings[str(path.relative_to(ROOT))] == block['source_bindings'][str(path.relative_to(ROOT))]
    labels = old['class_representatives'][0]['factor_labels']
    ids = [a['index'] for a in labels if a['kind'] == 'projected_atom']
    corrections = [a['index'] for a in labels if a['kind'] == 'generic']
    assert len(ids) == len(set(ids)) == protocol['required_projected_atom_count'] == 1676
    assert corrections == [3, 5, 6, 7, 9, 10]
    atoms = {a['index']: a for a in old['atoms']}
    chosen = min(ids)
    assert chosen == 2
    t0, u = Q(3, 17), Q(block['generic_control_transport']['short_scale_u'])
    assert u == Q(289, 2)
    A, B = [list(map(Q, family[k+'_coefficients_low_to_high'])) for k in ('A', 'B')]
    c1 = scale(add(scale(A, u**4), [27]), Q(1, 81))
    c0 = scale(add(add(scale(B, u**6), scale(A, 3*u**4)), [27]), Q(1, 729))
    f0 = [ev(c0, t0), ev(c1, t0), Q(1), Q(1)]
    assert f0 == list(map(Q, old['equation']['cubic_ascending'])) == list(map(Q, generic['cubic_ascending']))
    gammas = []
    assert len(family['sections']) == len(generic['generic_classes']) == 16
    for section, gamma in zip(family['sections'], generic['generic_classes']):
        def section_poly(key):
            v = section[key]
            den = trim(list(map(Q, v['denominator_coefficients_low_to_high'])))
            assert len(den) == 1 and den[0]
            return scale(list(map(Q, v['numerator_coefficients_low_to_high'])), 1/den[0])
        x, y = section_poly('X'), section_poly('Y')
        assert mul(y, y) == add(add(mul(mul(x, x), x), mul(A, x)), B)
        z, w = scale(add(scale(x, u*u), [-3]), Q(1, 9)), scale(y, u**3/27)
        assert mul(w, w) == add(add(add(mul(mul(z, z), z), mul(z, z)), mul(c1, z)), c0)
        g = list(map(Q, gamma['beta_ascending']))
        s = -g[1]
        assert s > 0 and math.isqrt(s.numerator)**2 == s.numerator and math.isqrt(s.denominator)**2 == s.denominator
        assert g == [s*ev(z, t0), -s, 0] and s**3*ev(w, t0)**2 == Q(gamma['norm'])
        gammas.append({'index': section['basis_index'], 'scale': str(s),
                       'z_ascending': list(map(str, z)), 'w_ascending': list(map(str, w))})
    norms = {}
    for i in ids:
        a, b, zero = map(Q, atoms[i]['alpha_ascending'])
        assert b and not zero
        norms[i] = norm(a, b, c0, c1)
        assert ev(norms[i], t0) == Q(atoms[i]['norm']) != 0
    n = norms[chosen]
    denominator = math.lcm(*(q.denominator for q in n))
    integers = [int(q*denominator) for q in n]
    content = math.gcd(*integers)*(1 if integers[-1] > 0 else -1)
    primitive = [q//content for q in integers]
    constant = Q(content, denominator)
    assert n == scale(primitive, constant) and len(primitive) == 13
    definition = {
        'schema': 'marked-dependency.coefficientwise-definition.v1', 'bindings': bindings,
        'column': 6, 'control_parameter': str(t0), 'short_scale_u': str(u),
        'field_cubic_ascending_in_theta': [list(map(str, c0)), list(map(str, c1)), ['1'], ['1']],
        'field_specialization_cubic_ascending': list(map(str, f0)),
        'coordinate_transport': 'z=(u^2*Xcompact-3)/9=4*x_original; w=u^3*Ycompact/27=8*y_original+4*x_original',
        'element_formula': 'alpha_i(t)=a_i+b_i*theta_t; pi_i(t)=Norm(alpha_i(t))*alpha_i(t); gamma_j(t)=scale_j*(z_j(t)-theta_t)',
        'dependency_formula': 'beta_6(t)=product of all listed pi_i(t), times gamma_3*gamma_5*gamma_6*gamma_7*gamma_9*gamma_10',
        'norm_square_root_formula': 'product_i Norm(alpha_i(t))^2 * product_j scale_j^(3/2)*w_j(t), with j the fixed generic corrections',
        'factor_labels': labels, 'generic_correction_indices': corrections,
        'atoms': [{'index': i, 'alpha_ascending': atoms[i]['alpha_ascending'],
                   'norm_at_control': atoms[i]['norm']} for i in ids],
        'generic_sections': gammas, 'distinguished_atom': chosen,
        'distinguished_norm_primitive_ascending': list(map(str, primitive)),
        'distinguished_norm_constant': str(constant),
        'control_compaction_circuit': next(c for c in block['circuits'] if c['column'] == 6),
        'boundary': 'The complete pre-compaction word specializes factor-by-factor to the original certificate. Its previous square-equivalence circuit to compact column 6 is inherited, not recomputed. No identification of specialized prime ideals across t is made.'
    }
    save(OUT/'definition.json', definition)
    result = {'schema': 'marked-dependency.continuation-check.v1', 'status': 'UNKNOWN',
              'definition_sha256': hashlib.sha256((OUT/'definition.json').read_bytes()).hexdigest()}
    try:
        p = protocol['coprimality_prime']
        assert prime(p)
        N = residue(primitive, p)
        assert len(N) == len(primitive)
        derivative = [i*x % p for i, x in enumerate(N)][1:]
        disc = add(add(mul(c1, c1), scale(mul(mul(c1, c1), c1), -4)),
                   add(scale(c0, -4), add(scale(mul(c0, c0), -27), scale(mul(c1, c0), 18))))
        product = divmodp(mul(derivative, residue(disc, p)), N, p)[1]
        for i in ids:
            if i != chosen:
                product = divmodp(mul(product, residue(norms[i], p)), N, p)[1]
        for g in gammas:
            product = divmodp(mul(product, residue(g['w_ascending'], p)), N, p)[1]
        gcd, left, right = xgcdp(N, product, p)
        assert gcd == [1]
        assert residue(add(mul(left, N), mul(right, product)), p) == [1]
        irred_prime = next(q for q in range(5, 102) if prime(q)
                           and all(x.denominator % q for x in f0)
                           and all(ev(residue(f0, q), a) % q for a in range(q)))
        result.update({
            'status': 'ODD_HORIZONTAL_VALUATION_OBSTRUCTION', 'prime': p, 'norm_degree': len(N)-1,
            'distinguished_atom': chosen, 'field_specialization_irreducibility_prime': irred_prime,
            'primitive_norm_mod_p': N, 'joint_product_remainder_mod_p': product,
            'bezout_norm_multiplier_mod_p': left, 'bezout_product_multiplier_mod_p': right,
            'joint_product_definition': "N'(t)*disc_theta(f_t)*product_(i != 2) Norm(alpha_i(t))*product_(j=0..15) w_j(t), reduced modulo N(t)",
            'checked_atoms': len(ids), 'checked_other_atom_coprimalities': len(ids)-1,
            'checked_generic_sections': len(gammas), 'geometric_valuation_pattern': [2, 1, 1],
            'conclusion': 'The specified beta_6(t) has odd valuation on the two residual sheets over every irreducible factor of the distinguished norm. These are good-reduction divisors. Thus this word is not a rational-point Kummer class over Q(t), even after a square or generic-point correction.',
            'boundary': 'Only this coefficientwise continuation is obstructed. The soluble strict class at 3/17, other continuations, arithmetic incidence conditions and ramified base changes remain unaffected.'})
    except (AssertionError, ValueError, ZeroDivisionError, StopIteration) as error:
        result['failure'] = type(error).__name__+': '+str(error)
    result['cpu_seconds'] = time.process_time()-started
    result['peak_rss_bytes'] = 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    save(OUT/'result.json', result)
    print(json.dumps({k: v for k, v in result.items() if 'mod_p' not in k}, indent=2))


if __name__ == '__main__':
    main()
