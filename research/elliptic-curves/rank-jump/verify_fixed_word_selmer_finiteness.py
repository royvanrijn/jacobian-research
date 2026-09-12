#!/usr/bin/env python3
"""Independent SymPy replay of the finite-prime binding; no producer imports.

The exact integer resultants remain Sylvester-determinant circuits. This
replay independently checks every nonzero residue and every scalar guard.
Faltings's theorem is a cited mathematical input, not machine-formalized here.
"""
import argparse
import hashlib
import json
import math
import resource
import time
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/generated-results/elliptic-curves'
PRIOR = BASE/'frozen_dependency_continuation_v1'
OUT = BASE/'fixed_word_selmer_finiteness_v1'


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    started = time.process_time()
    binding, definition, prior = read(OUT/'sigma-binding.json'), read(PRIOR/'definition.json'), read(PRIOR/'result.json')
    for path, digest in binding['sources'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    for path, digest in definition['bindings'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    assert binding['column'] == definition['column'] == 6
    assert prior['definition_sha256'] == hashlib.sha256((PRIOR/'definition.json').read_bytes()).hexdigest()
    old_atoms = read(BASE/'rank_jump_reference_additional_strict_verification_v1.json')
    assert definition['factor_labels'] == old_atoms['class_representatives'][0]['factor_labels']
    lookup = {a['index']: a for a in old_atoms['atoms']}
    ids = [a['index'] for a in definition['factor_labels'] if a['kind'] == 'projected_atom']
    assert len(ids) == len(set(ids)) == len(definition['atoms']) == 1676
    assert definition['atoms'] == [{'index': i, 'alpha_ascending': lookup[i]['alpha_ascending'],
                                    'norm_at_control': lookup[i]['norm']} for i in ids]
    assert definition['generic_correction_indices'] == [3, 5, 6, 7, 9, 10]
    t, theta, m, n = sp.symbols('t theta m n')
    def qp(cs):
        return sp.Poly.from_list([sp.Rational(q) for q in reversed(cs)], t, domain=sp.QQ)
    N = qp(binding['distinguished_norm_primitive_ascending'])
    assert N == qp(definition['distinguished_norm_primitive_ascending'])
    assert N.degree() == 12 and N.LC() > 0 and N.content() == 1
    assert all(c.q == 1 for c in N.all_coeffs())
    c0, c1, one, other_one = map(qp, definition['field_cubic_ascending_in_theta'])
    assert one == other_one == qp(['1'])
    f = theta**3+theta**2+c1.as_expr()*theta+c0.as_expr()
    disc = sp.Poly(sp.discriminant(f, theta), t, domain=sp.QQ)
    a2, b2, zero = map(sp.Rational, lookup[2]['alpha_ascending'])
    assert b2 and not zero and definition['distinguished_atom'] == 2
    lam = sp.Rational(definition['distinguished_norm_constant'])
    assert sp.Poly(sp.resultant(f, a2+b2*theta, theta), t, domain=sp.QQ) == lam*N
    t0 = sp.Rational(definition['control_parameter'])
    f0 = [c0.eval(t0), c1.eval(t0), sp.Rational(1), sp.Rational(1)]
    assert f0 == list(map(sp.Rational, definition['field_specialization_cubic_ascending']))
    raw = [('N_derivative', 'derivative', None, N.diff()),
           ('cubic_discriminant', 'discriminant', None, disc)]
    for atom in definition['atoms']:
        a, b, zero = map(sp.Rational, atom['alpha_ascending'])
        assert not zero and b
        norm = c1.mul_ground(a*b*b)-c0.mul_ground(b**3)+(a**3-a*a*b)
        assert norm.eval(t0) == sp.Rational(atom['norm_at_control']) != 0
        if atom['index'] != 2:
            raw.append((f"atom_norm_{atom['index']:05d}", 'other_atom_norm', atom['index'], norm))
    for j, section in enumerate(definition['generic_sections']):
        assert section['index'] == j
        z, w, s = qp(section['z_ascending']), qp(section['w_ascending']), sp.Rational(section['scale'])
        assert s > 0 and sp.sqrt(s).is_Rational
        assert w*w == z**3+z*z+c1*z+c0
        raw.append((f'generic_w_{j:02d}', 'generic_ordinate', j, w))
    assert len(raw) == len(binding['polynomials']) == binding['resultant_count'] == 1693
    p = binding['nonvanishing_prime']
    assert p == prior['prime'] == 1000003 and sp.isprime(p)
    def fp(q):
        if not isinstance(q, sp.Poly):
            q = qp(q)
        return sp.Poly.from_list([int(c.p)*pow(int(c.q), -1, p) % p for c in q.all_coeffs()], t, modulus=p)
    np = fp(N)
    assert np.degree() == 12
    product = fp(['1'])
    resultants_product, denominator_product = 1, 1
    denominators = []
    for (label, kind, index, q), row in zip(raw, binding['polynomials']):
        assert (row['id'], row['kind'], row['source_index']) == (label, kind, index)
        D, integer_q = q.clear_denoms(convert=True)
        assert D == int(row['denominator_multiplier'])
        assert integer_q == qp(row['integer_coefficients_ascending']).set_domain(sp.ZZ)
        assert integer_q.degree() == row['degree']
        mod_q = fp(integer_q)
        fixed_degree_factor = pow(int(np.LC()) % p, row['degree']-mod_q.degree(), p)
        residue = int(np.resultant(mod_q))*fixed_degree_factor % p
        assert residue == row['resultant']['nonzero_residue_mod_p'] != 0
        assert row['resultant']['exact_definition'] == 'det Sylvester(N,Q)'
        assert row['resultant']['sylvester_matrix_size'] == 12+row['degree']
        assert row['resultant']['actual_right_degree_mod_p'] == mod_q.degree()
        assert row['resultant']['integer_value_expanded'] is False
        resultants_product = resultants_product * residue % p
        denominator_product = denominator_product * int(D) % p
        denominators.append(int(D))
        product = (product*mod_q).rem(np)
    scaled = binding['scaled_prior_bezout']
    assert denominator_product == scaled['denominator_product_mod_p']
    assert product == fp(scaled['product_remainder_mod_p'])
    assert product == fp(prior['joint_product_remainder_mod_p']).mul_ground(denominator_product)
    left, right = fp(scaled['norm_multiplier_mod_p']), fp(scaled['product_multiplier_mod_p'])
    assert left == fp(prior['bezout_norm_multiplier_mod_p'])
    assert right == fp(prior['bezout_product_multiplier_mod_p']).mul_ground(pow(denominator_product, -1, p))
    assert left*np+right*product == fp(['1'])
    # Independently verify that C makes every necessary scalar a p-adic unit
    # at all primes outside its support. Extra guard factors are harmless.
    guard_values = [int(g['value']) for g in binding['scalar_unit_guards']]
    assert all(v > 1 for v in guard_values) and len(set(guard_values)) == len(guard_values)
    C = math.prod(guard_values)
    assert C == int(binding['scalar_product_C'])
    required = [2, 3, abs(int(N.LC()))]+denominators
    for q in (lam, b2, sp.Rational(definition['short_scale_u'])):
        assert q
        required += [abs(int(q.p)), int(q.q)]
    for coefficients in definition['field_cubic_ascending_in_theta']:
        required += [int(sp.Rational(q).q) for q in coefficients]
    for atom in definition['atoms']:
        required += [int(sp.Rational(q).q) for q in atom['alpha_ascending']]
    for section in definition['generic_sections']:
        s = sp.Rational(section['scale'])
        required += [abs(int(s.p)), int(s.q)]
        for name in ('z_ascending', 'w_ascending'):
            required += [int(sp.Rational(q).q) for q in section[name]]
    assert all(C % v == 0 for v in required)
    assert resultants_product == binding['product_of_resultants_mod_p']
    assert C % p * resultants_product % p == binding['exceptional_integer_M_mod_p'] != 0
    assert binding['exceptional_integer_M'] == 'C * product_(all listed Q) det Sylvester(N,Q)'
    homogeneous = sum(N.nth(j)*m**j*n**(12-j) for j in range(13))
    assert sp.cancel(n**12*N.as_expr().subs(t, m/n)-homogeneous) == 0
    assert sp.Poly(homogeneous, n).nth(0) == N.LC()*m**12
    assert binding['twist_genus'] == (N.degree()-2)//2 == 5
    receipt = {
        'schema': 'fixed-word.finiteness-independent-replay.v1', 'status': 'VERIFIED',
        'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sigma_binding_sha256': hashlib.sha256((OUT/'sigma-binding.json').read_bytes()).hexdigest(),
        'prior_definition_sha256': hashlib.sha256((PRIOR/'definition.json').read_bytes()).hexdigest(),
        'independent_engine': 'SymPy '+sp.__version__, 'full_atom_specializations_checked': 1676,
        'generic_section_identities_checked': 16, 'integer_polynomial_normalizations_checked': len(raw),
        'independent_modular_resultants_checked': len(raw), 'M_mod_1000003': binding['exceptional_integer_M_mod_p'],
        'scalar_guard_count': len(guard_values), 'homogeneous_denominator_identity': 'VERIFIED SYMBOLICALLY',
        'cpu_seconds': time.process_time()-started,
        'peak_rss_bytes': 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'boundary': 'Verifies the finite-prime binding and exact polynomial hypotheses. The written finiteness deduction invokes Faltings. It does not expand/factor resultants, enumerate Sigma or twists, test new specializations, replay ideal factorizations, or use column7.'}
    if args.write:
        with (OUT/'independent-replay.json').open('x') as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
