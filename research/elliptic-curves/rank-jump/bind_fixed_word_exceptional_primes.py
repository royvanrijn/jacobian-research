#!/usr/bin/env python3
"""Bind an exact finite exceptional-prime set for the fixed column-6 word.

Integer resultants are specified as exact Sylvester determinants of expanded
integer polynomials. Their nonzero residues are checked, but the determinants
and their product are not expanded or factored. No twist or parameter search.
"""
import argparse
import hashlib
import json
import math
import resource
import time
from fractions import Fraction as Q
from pathlib import Path
import finish_frozen_dependency_continuation as arithmetic

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/generated-results/elliptic-curves'
PRIOR = BASE/'frozen_dependency_continuation_v1'
DEFAULT = BASE/'fixed_word_selmer_finiteness_v1'


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def resultant_mod_p(f, g, original_g_degree, p):
    """Reduction of the fixed-degree integer Sylvester determinant.

The degree of f is preserved. If the degree of g drops, the leading
coefficient of f supplies the missing factors in the fixed-degree formula.
"""
    f, g = arithmetic.residue(f, p), arithmetic.residue(g, p)
    assert g != [0]
    result = pow(f[-1], original_g_degree-(len(g)-1), p)
    while len(g) > 1:
        m, n = len(f)-1, len(g)-1
        remainder = arithmetic.divmodp(f, g, p)[1]
        if remainder == [0]:
            return 0
        result = result * (-1 if m*n % 2 else 1) * pow(g[-1], m-(len(remainder)-1), p) % p
        f, g = g, remainder
    return result * pow(g[0], len(f)-1, p) % p


def raw_polynomials(definition):
    N = list(map(Q, definition['distinguished_norm_primitive_ascending']))
    c0, c1 = [list(map(Q, v)) for v in definition['field_cubic_ascending_in_theta'][:2]]
    add, mul, scale = arithmetic.add, arithmetic.mul, arithmetic.scale
    yield 'N_derivative', 'derivative', None, [i*x for i, x in enumerate(N)][1:]
    disc = add(add(mul(c1, c1), scale(mul(mul(c1, c1), c1), -4)),
               add(scale(c0, -4), add(scale(mul(c0, c0), -27), scale(mul(c1, c0), 18))))
    yield 'cubic_discriminant', 'discriminant', None, disc
    for atom in definition['atoms']:
        if atom['index'] != 2:
            a, b, zero = map(Q, atom['alpha_ascending'])
            assert b and not zero
            yield f"atom_norm_{atom['index']:05d}", 'other_atom_norm', atom['index'], arithmetic.norm(a, b, c0, c1)
    for section in definition['generic_sections']:
        yield f"generic_w_{section['index']:02d}", 'generic_ordinate', section['index'], list(map(Q, section['w_ascending']))


def scalar_guards(definition, polynomial_rows):
    """Every retained guard is an explicit nonzero integer, without factoring."""
    factors = {}
    def include(n, reason):
        n = abs(int(n))
        assert n
        if n != 1:
            factors.setdefault(n, set()).add(reason)
    def denominator(q, reason):
        include(Q(q).denominator, reason)
    def unit(q, reason):
        q = Q(q)
        assert q
        include(q.numerator, reason+':numerator')
        include(q.denominator, reason+':denominator')
    include(2, 'odd-residue-characteristic and discriminant(E)=16*discriminant(f)')
    include(definition['distinguished_norm_primitive_ascending'][-1], 'leading_coefficient_of_N')
    unit(definition['distinguished_norm_constant'], 'n_2=lambda*N')
    unit(definition['short_scale_u'], 'coordinate_scale_u')
    include(3, 'coordinate_transport_divisors_9_and_27')
    for j, coefficients in enumerate(definition['field_cubic_ascending_in_theta']):
        for k, q in enumerate(coefficients):
            denominator(q, f'field_coefficient:{j}:{k}')
    for atom in definition['atoms']:
        for k, q in enumerate(atom['alpha_ascending']):
            denominator(q, f"atom_coefficient:{atom['index']}:{k}")
        if atom['index'] == 2:
            unit(atom['alpha_ascending'][1], 'distinguished_atom_slope_b_2')
    for section in definition['generic_sections']:
        j = section['index']
        unit(section['scale'], f'generic_scale:{j}')
        for key in ('z_ascending', 'w_ascending'):
            for k, q in enumerate(section[key]):
                denominator(q, f'generic_coefficient:{j}:{key}:{k}')
    for row in polynomial_rows:
        include(row['denominator_multiplier'], 'denominator_clearing:'+row['id'])
    return [{'value': str(n), 'reasons': sorted(reasons)} for n, reasons in sorted(factors.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    started = time.process_time()
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = {'schema': 'fixed-word.finiteness-binding-protocol.v1', 'column': 6,
                'cpu_seconds': 60, 'address_space_bytes': 1024**3, 'worker_count': 1,
                'purpose': 'Bind the finite-prime lemma using only the retained polynomial identities and the complete fixed word.',
                'resultants': 'Exact integer Sylvester-determinant circuits; no integer expansion or factorization.',
                'new_atoms': 0, 'parameter_searches': 0, 'twist_searches': 0,
                'column_7': 'untouched', 'on_failure': 'Preserve completed evidence; no alternate prime, word or increased bound.'}
    save(args.output/'protocol.json', protocol)
    definition, prior = read(PRIOR/'definition.json'), read(PRIOR/'result.json')
    assert hashlib.sha256((PRIOR/'definition.json').read_bytes()).hexdigest() == prior['definition_sha256']
    assert definition['column'] == 6 and definition['distinguished_atom'] == 2
    assert len(definition['atoms']) == 1676
    assert prior['status'] == 'ODD_HORIZONTAL_VALUATION_OBSTRUCTION'
    assert len(definition['generic_sections']) == 16
    N = list(map(int, definition['distinguished_norm_primitive_ascending']))
    assert len(N) == 13 and N[-1] > 0 and math.gcd(*N) == 1
    p = prior['prime']
    assert p == 1000003 and arithmetic.prime(p) and N[-1] % p
    np = arithmetic.residue(N, p)
    rows, product, denominator_product = [], [1], 1
    resultant_product = 1
    for label, kind, index, raw in raw_polynomials(definition):
        raw = list(map(Q, arithmetic.trim(raw)))
        D = math.lcm(*(q.denominator for q in raw))
        coefficients = [int(D*q) for q in raw]
        assert coefficients[-1] and D % p
        degree = len(coefficients)-1
        residue = resultant_mod_p(np, coefficients, degree, p)
        assert residue, label
        rows.append({'id': label, 'kind': kind, 'source_index': index,
                     'denominator_multiplier': str(D), 'integer_coefficients_ascending': list(map(str, coefficients)),
                     'degree': degree, 'normalization': 'Q=D*q; no primitive-content or squareclass normalization',
                     'resultant': {'exact_definition': 'det Sylvester(N,Q)', 'sylvester_matrix_size': 12+degree,
                                   'integer_value_expanded': False, 'nonzero_residue_mod_p': residue,
                                   'actual_right_degree_mod_p': len(arithmetic.residue(coefficients, p))-1}})
        product = arithmetic.divmodp(arithmetic.mul(product, arithmetic.residue(coefficients, p)), np, p)[1]
        denominator_product = denominator_product * D % p
        resultant_product = resultant_product * residue % p
    assert len(rows) == 1693
    old_product = prior['joint_product_remainder_mod_p']
    assert product == arithmetic.residue(arithmetic.scale(old_product, denominator_product), p)
    left = prior['bezout_norm_multiplier_mod_p']
    right = arithmetic.residue(arithmetic.scale(prior['bezout_product_multiplier_mod_p'], pow(denominator_product, -1, p)), p)
    assert arithmetic.residue(arithmetic.add(arithmetic.mul(left, np), arithmetic.mul(right, product)), p) == [1]
    guards = scalar_guards(definition, rows)
    C = math.prod(int(g['value']) for g in guards)
    assert C % p
    source_paths = [PRIOR/'definition.json', PRIOR/'result.json', PRIOR/'independent-replay.json',
                    Path(__file__), Path(arithmetic.__file__)]
    result = {
        'schema': 'fixed-word.exceptional-prime-binding.v1', 'status': 'FINITE_PRIME_SET_BOUND',
        'column': 6, 'sources': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
        'distinguished_norm_primitive_ascending': list(map(str, N)), 'norm_degree': 12,
        'scalar_unit_guards': guards, 'scalar_product_C': str(C),
        'polynomials': rows, 'resultant_count': len(rows), 'nonvanishing_prime': p,
        'product_of_resultants_mod_p': resultant_product,
        'exceptional_integer_M_mod_p': C % p * resultant_product % p,
        'exceptional_integer_M': 'C * product_(all listed Q) det Sylvester(N,Q)',
        'Sigma': 'The rational primes dividing abs(M). M is specified by the exact integer determinant circuit above, not by a factored or expanded integer.',
        'expanded_resultants': False, 'factored_resultants': False, 'enumerated_Sigma': False,
        'scaled_prior_bezout': {'denominator_product_mod_p': denominator_product,
                                'product_remainder_mod_p': product, 'norm_multiplier_mod_p': left,
                                'product_multiplier_mod_p': right},
        'homogeneous_norm': 'N_hom(m,n)=sum_(j=0..12) N_j*m^j*n^(12-j)=n^12*N(m/n)',
        'local_necessary_condition': 'For every p outside Sigma, a locally soluble eligible specialization has even v_p(N_hom(m,n)).',
        'twist_containment': 'All locally soluble eligible parameters lie in the images of delta*y^2=N(t), for signed squarefree delta supported on Sigma.',
        'twist_genus': 5,
        'boundary': 'This is an exact finite-prime binding, not a computation of Selmer membership, twist rational points, the finite parameter set, a useful numerical bound, or an obstruction to parameter-dependent constructors. Generic-ordinate and N-derivative resultants enlarge Sigma conservatively. Only column6 is used.',
        'cpu_seconds': time.process_time()-started,
        'peak_rss_bytes': 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    save(args.output/'sigma-binding.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'resultant_count', 'nonvanishing_prime', 'exceptional_integer_M_mod_p', 'twist_genus', 'cpu_seconds', 'peak_rss_bytes')}, indent=2))


if __name__ == '__main__':
    main()
