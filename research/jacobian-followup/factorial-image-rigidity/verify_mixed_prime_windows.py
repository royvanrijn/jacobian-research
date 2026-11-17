"""Exact controls for mixed-prime lifting, kernels and square-divisor profiles."""
from pathlib import Path
from fractions import Fraction
from math import comb, factorial
import hashlib
import json
import sympy as s
import moment_bridges as b
import prime_power_windows as p
import mixed_prime_windows as m

HERE = Path(__file__).resolve().parent
counts = {}
examples = []


def check(group, value):
    if not value:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


for cofactor, prime in ((1, 2), (2, 3), (3, 2), (4, 3), (6, 5), (3, 3)):
    for insertion_power in (1, 2):
        for n in (1, 2, 5, 7):
            for shape in sorted({1, prime**insertion_power}):
                certificate = m.lifting_certificate(n, cofactor, prime, insertion_power, shape)
                check('literal_insertion_preserves_hull_and_scaled_residuals', certificate['status'] == 'newton_data_lifting_identity_verified')
                base = m.square_divisor_degrees(certificate['base_edges'], prime)
                lifted = m.square_divisor_degrees(certificate['lifted_edges'], prime)
                check('all_residual_degree_budgets_preserved_by_insertion', base['possible_degrees'] == lifted['possible_degrees'])
                shift = certificate['parameter_scale_denominator_exponent']
                d = certificate['exponent']
                for j in range(n+1):
                    multiplier = Fraction(int(s.rf(shape, d*j)), factorial(cofactor*j)*prime**(shift*j))
                    check('coefficient_insertion_multiplier_is_geometric_unit', p.valuation(multiplier, prime) == 0
                          and p.residue(multiplier, prime) == (-1)**(shift*j) % prime)

T, x = s.symbols('T x')
families = [(2, 1, 3), (2, 2, 7), (3, 1, 5), (3, 2, 7), (4, 1, 5), (1, 2, 3)]
for cofactor, k, prime in families:
    family = m.indexed_family_certificate(cofactor, k, prime)
    check('finite_kernel_accepts_unbounded_a_r_family', family['status'] == 'unbounded_indexed_family_certified')
    examples.append(family)
    for r in (1, 2):
        n = k*prime**r
        e = prime**r
        c = (e-1)//(prime-1)
        coefficients = [comb(n, j)*factorial(cofactor*j) for j in range(n+1)]
        rows = p.newton_residuals(coefficients, prime)
        check('single_edge_literal_vertices', len(rows) == 1 and rows[0]['left'] == [0, 0]
              and rows[0]['right'] == [n, cofactor*k*c]
              and rows[0]['e'] == e and rows[0]['h'] == cofactor*c)
        expected = [comb(k, j)*factorial(cofactor*j)*(-1)**(cofactor*c*j) % prime for j in range(k+1)]
        check('single_edge_literal_small_kernel_residual', rows[0]['residual_coefficients'] == expected)
        for j, coefficient in enumerate(coefficients):
            excess = p.valuation(coefficient, prime)-Fraction(cofactor*j*c, e)
            check('all_kernel_edge_lattice_points_and_strict_off_edge_points', excess >= 0 and
                  ((excess == 0) == (j % e == 0)))
        # Literal insertion checks are deliberately small controls of the formula.
        if n <= 27:
            for a in (1, 2):
                shape = prime**a
                lifted = m.lifting_certificate(n, cofactor, prime, a, shape)
                check('kernel_family_literal_lifted_residuals_squarefree', all(row['residual_squarefree'] for row in lifted['lifted_edges']))
                if n <= 14 and cofactor*prime**a <= 25:
                    f = p.moment_polynomial(n, cofactor*prime**a, shape, x)
                    previous = p.moment_polynomial(n-1, cofactor*prime**a, shape, x)
                    check('separate_kernel_family_rational_gcd', s.gcd(f, previous).degree() == 0)

bad_kernel = m.indexed_family_certificate(2, 2, 5)
check('bad_kernel_is_unknown_not_counterexample', bad_kernel['status'] == 'residual_kernel_not_squarefree_no_decision'
      and bad_kernel['counterexample_claimed'] is False and bad_kernel['theorem'] is None)
examples.append(bad_kernel)

for prime, expression, expected in ((2, T**2+T+1, True), (2, T**2+1, False),
                                    (3, T**2-1, False), (3, T**2+1, True),
                                    (5, T+2, True), (5, (T+1)**3, False)):
    check('separate_frobenius_factor_irreducibility_controls',
          m._irreducible_verified(s.Poly(expression, T, modulus=prime), prime) == expected)

square = s.Poly((1+x)**4*(x*x+x+1)**2, x)
coefficients = [int(square.nth(j)) for j in range(square.degree()+1)]
for prime in (2, 3, 5, 7):
    rows = p.newton_residuals(coefficients, prime)
    profile = m.square_divisor_degrees(rows, prime)
    check('known_square_divisor_degree_survives_every_local_profile', 4 in profile['possible_degrees'])
    check('every_exported_residual_factor_irreducibility_verified', all(row['factor_irreducibility_verified'] for row in profile['branches']))

# Independent Sturm controls for the all-shape real-place theorem.
for shape in (s.Rational(1,2), s.Integer(1), s.Rational(5,3)):
    for degree in (2, 6):
        for n in range(1, 8):
            actual = s.Poly(sum(s.binomial(n,j)*s.rf(shape,degree*j)*x**j for j in range(n+1)), x, domain=s.QQ)
            previous = s.Poly(sum(s.binomial(n-1,j)*s.rf(shape,degree*j)*x**j for j in range(n)), x, domain=s.QQ)
            check('exact_adjacent_derivative_identity_for_actual_Gamma_moments',
                  actual-actual.diff()*s.Poly(x/n,x)-previous == 0)
            check('Sturm_real_root_count_matches_real_place_theorem',
                  actual.count_roots(-s.oo,s.oo) == n % 2)
            check('no_repeated_real_root_in_separate_exact_gcd',
                  s.gcd(actual,actual.diff()).count_roots(-s.oo,s.oo) == 0)
            if n % 2:
                check('odd_index_real_root_is_negative_and_derivative_has_no_real_root',
                      actual.count_roots(-s.oo,0)==1 and actual.diff().count_roots(-s.oo,s.oo)==0)
            if n <= 3:
                check('first_three_indices_squarefree_for_positive_shapes', s.gcd(actual,actual.diff()).degree()==0)

raw_seven = m.mixed_window_certificate(7,6,[3],use_real_parity=False)
real_seven = m.mixed_window_certificate(7,6,[3])
check('one_prime_repeated_residual_plus_real_place_certifies_squarefree',
      raw_seven['status']=='unresolved_by_supplied_local_data'
      and raw_seven['possible_square_divisor_degrees']==[0,1]
      and real_seven['local_possible_square_divisor_degrees']==[0,1]
      and real_seven['possible_square_divisor_degrees']==[0]
      and real_seven['status']=='adjacent_window_certified')
check('real_place_filter_is_recorded_separately_from_local_data',
      real_seven['real_place_parity_used'] and not raw_seven['real_place_parity_used']
      and real_seven['profiles']==raw_seven['profiles'])
check('generic_rational_square_divisor_test_does_not_drop_odd_degrees',
      1 in m.square_divisor_degrees(p.newton_residuals([1,2,1],3),3)['possible_degrees'])
examples.append(real_seven)

point_rows = []
for n in range(2, 33):
    certificate = m.mixed_window_certificate(n, 6, [2, 3])
    point_rows.append(certificate)
    f = p.moment_polynomial(n, 6, 1, x)
    previous = p.moment_polynomial(n-1, 6, 1, x)
    check('d6_finite_control_rational_adjacent_gcd', s.gcd(f, previous).degree() == 0)
    if certificate['status'] == 'adjacent_window_certified':
        check('accepted_degree_intersection_is_zero_only', certificate['possible_square_divisor_degrees'] == [0])
    else:
        check('real_parity_leaves_genuine_even_degree_survivor', certificate['possible_square_divisor_degrees']==[0,2]
              and certificate['real_place_parity_used'])
        check('unresolved_intersection_is_not_counterexample', certificate['counterexample_claimed'] is False
              and certificate['ordinary_pair_profile'] is None)
check('first_order_two_prime_method_has_explicit_survivors',
      [row['terminal_index'] for row in point_rows if row['status'] != 'adjacent_window_certified'] == [21, 23])
for n, expected in ((7, [[0, 2], [0, 1]]), (22, [[0, 8], [0, 1, 2, 3, 4]])):
    certificate = next(row for row in point_rows if row['terminal_index'] == n)
    check('neither_prime_alone_suffices_but_intersection_does', certificate['status'] == 'adjacent_window_certified'
          and [row['possible_degrees'] for row in certificate['profiles']] == expected)
    examples.append(certificate)
for n in (21, 23):
    certificate = m.mixed_window_certificate(n, 6, [2, 3, 5])
    check('third_prime_closes_supplied_finite_survivor_only', certificate['status'] == 'adjacent_window_certified'
          and certificate['all_starting_indices_claimed'] is False)
    examples.append(certificate)

ys = [s.Symbol('Y')]
for n, degree, shape in ((5, 6, s.Rational(1, 2)), (4, 10, s.Rational(1, 3)),
                         (6, 12, 2), (5, 15, 1)):
    certificate = m.mixed_window_certificate(n, degree, [2, 3, 5, 7], shape=shape)
    f = p.moment_polynomial(n, degree, 1, x) if shape == 1 else s.Poly(
        sum(s.binomial(n, j)*s.rf(shape, degree*j)*x**j for j in range(n+1)), x)
    direct = b.gamma(s.expand((1+x*ys[0]**degree)**n), ys, [shape])
    check('actual_rational_gamma_coefficients_in_local_certificate', s.expand(direct-f.as_expr()) == 0)
    if certificate['status'] == 'adjacent_window_certified':
        check('separate_rational_gamma_squarefree_gcd', s.gcd(f, f.diff()).degree() == 0)
    check('no_universal_gamma_shape_claim_from_point_input', certificate['all_starting_indices_claimed'] is False)
    examples.append(certificate)

for n, degree in ((7, 6), (5, 10), (4, 12)):
    certificate = m.mixed_window_certificate(n, degree, [2, 3, 5, 7])
    check('ordinary_pair_control_is_certified_first', certificate['status'] == 'adjacent_window_certified')
    A, B = s.symbols('A B')
    moments = [sum(s.binomial(k, j)*s.factorial(degree*j)*A**(k-j)*B**j
                   for j in range(k+1)) for k in (n-1, n)]
    gb = s.groebner(moments, A, B, order='grevlex', domain=s.QQ)
    leading = [poly.monoms(order=gb.order)[0] for poly in gb.polys]
    standard = [(i, total-i) for total in range(2*n-1) for i in range(total+1)
                if not any(i >= u and total-i >= v for u, v in leading)]
    profile = certificate['ordinary_pair_profile']
    check('ordinary_pair_length_and_top_degree', len(standard) == profile['length']
          and max(i+j for i, j in standard) == profile['highest_nonzero_degree'])
    check('ordinary_sharp_nilpotence_containment', all(gb.reduce(A**i*B**(2*n-2-i))[1] == 0 for i in range(2*n-1))
          and any(gb.reduce(A**i*B**(2*n-3-i))[1] != 0 for i in range(2*n-2)))

# Actual odd-block Gaussian transfer for a finite half-shape certificate.
certificate = m.mixed_window_certificate(5, 6, [2, 3, 5, 7], shape=s.Rational(1, 2))
check('half_shape_pair_has_a_local_certificate', certificate['status'] == 'adjacent_window_certified')
phase = 1+s.I*ys[0]**6
lift, xs, _ = b.gaussian_lift(phase, ys, [1])
for k in (4, 5):
    check('half_shape_original_gaussian_moments_match', b.gaussian(s.expand(lift**k), xs)
          == b.gamma(s.expand(phase**k), ys, [s.Rational(1, 2)]))

for args in ((2, 2, 4), (2, 2, 3), (0, 1, 3)):
    try:
        m.indexed_family_certificate(*args)
    except ValueError:
        check('unsupported_small_kernel_inputs_rejected', True)
    else:
        raise AssertionError('unsupported kernel accepted')
for args in ((3, 2, 3, 1, 4), (3, 2, 4, 1, 1)):
    try:
        m.lifting_certificate(*args)
    except ValueError:
        check('unsupported_lifting_shape_or_prime_rejected', True)
    else:
        raise AssertionError('unsupported lifting input accepted')
for invalid_flag in (0, 1, 'yes', None):
    try:
        m.mixed_window_certificate(7,6,[3],use_real_parity=invalid_flag)
    except ValueError:
        check('real_place_control_requires_explicit_boolean', True)
    else:
        raise AssertionError('nonboolean real-place flag accepted')

for shape in (0, 0.5, s.sqrt(2)):
    try:
        m.mixed_window_certificate(3, 6, shape=shape)
    except ValueError:
        check('unsupported_point_shape_fails_closed', True)
    else:
        raise AssertionError('unsupported point shape accepted')

probe = {'date': '2026-11-16', 'scope': 'finite terminal indices2..32 at exponent6, primes2/3 plus the proved real-place parity restriction',
         'all_starting_indices_claimed': False,
         'unresolved_indices': [row['terminal_index'] for row in point_rows if row['status'] != 'adjacent_window_certified'],
         'certificates': point_rows}
(HERE/'mixed_prime_d6_probe.json').write_text(json.dumps(probe, indent=2)+'\n')
receipt = {'date': '2026-11-16', 'status': 'pass', 'checks': counts, 'total': sum(counts.values()),
           'scope': 'Finite exact insertion, single-edge kernel, residual factor/Frobenius, degree-intersection plus real-place parity/Sturm, original Gamma/Gaussian and ordinary-pair controls. Universal unbounded insertion/indexed-family statements use the written proof. Finite d6 survivors remain explicitly UNKNOWN under their supplied two-prime data.',
           'files': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                     ('MIXED_PRIME_RESIDUAL_WINDOWS.tex', 'mixed_prime_windows.py', 'verify_mixed_prime_windows.py',
                      'prime_power_windows.py', 'moment_bridges.py')},
           'independent_replay': False, 'formal_verification': False, 'external_review': False,
           'novelty_claim': False}
(HERE/'mixed_prime_window_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'mixed_prime_window_examples.json').write_text(json.dumps(examples, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
