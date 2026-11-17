"""Exact controls for the all-index prime-power two-term SFC proof."""
from pathlib import Path
from fractions import Fraction
from math import factorial, comb
import hashlib
import json
import sympy as s
import moment_bridges as bridge
import prime_power_windows as p

HERE = Path(__file__).resolve().parent
counts = {}
examples = []


def check(group, value):
    if not value:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


def digits(n, prime):
    result = []
    while n:
        result.append(n % prime)
        n //= prime
    return result or [0]


x, T, Z = s.symbols('x T Z')
for prime in (2, 3, 5, 7, 11):
    for n in (1, 2, 3, 4, 6, 8, 9, 12, 17, 20, 24):
        gamma = prime-1
        for unit_kind in (0, 1):
            units = [Fraction(gamma**k*(1+prime*(k*k+n)), 1 if unit_kind == 0 else 1+prime*(k+1))
                     for k in range(n+1)]
            coeffs = [Fraction(factorial(n), factorial(n-k))*units[k] for k in range(n+1)]
            actual = p.newton_residuals(coeffs, prime)
            expected = p.expected_edges(n, prime, gamma)
            check('generic_geometric_unit_deformation_edges_and_residuals', actual == expected)
            check('generic_deformation_residuals_squarefree', all(row['residual_squarefree'] for row in actual))
            # Verify every coefficient's relation to every actual edge line.
            for row in actual:
                (left, vl), (right, vr) = row['left'], row['right']
                slope = Fraction(vr-vl, right-left)
                for k in range(left, right+1):
                    excess = p.valuation(coeffs[k], prime)-vl-(k-left)*slope
                    check('edge_points_above_line_with_exact_lattice_equalities', excess >= 0 and
                          ((excess == 0) == ((k-left) % row['e'] == 0)))
            if unit_kind == 0 and n <= 12:
                polynomial = s.Poly(sum(s.Rational(c.numerator, c.denominator)*x**k for k, c in enumerate(coeffs)), x)
                check('independent_generic_deformation_rational_gcd', s.gcd(polynomial, polynomial.diff()).degree() == 0)

for prime in (2, 3, 5, 7, 11):
    for degree in range(prime):
        E = s.Poly(sum(s.Rational(1, factorial(j))*T**j for j in range(degree+1)), T, domain=s.QQ)
        # Construct finite-field coefficients without reducing rational denominators incorrectly.
        reduced = s.Poly.from_dict({m: p.residue(Fraction(c), prime) for m, c in E.terms()}, T, modulus=prime)
        check('truncated_exponential_separability_below_prime', s.gcd(reduced, reduced.diff()).degree() == 0)
        if degree:
            previous = s.Poly(sum(pow(factorial(j), -1, prime)*T**j for j in range(degree)), T, modulus=prime)
            check('truncated_exponential_derivative_identity', reduced.diff() == previous)

for prime in (2, 3, 5, 7):
    for n in (0, 1, prime-1, prime, prime+1, prime**2+2):
        value = factorial(n)
        valuation = p.factorial_valuation(n, prime)
        unit = value//prime**valuation
        expected = (-1)**valuation
        for digit in digits(n, prime):
            expected *= factorial(digit)
        check('factorial_unit_digit_formula', unit % prime == expected % prime)

degrees = [2, 4, 8, 3, 9, 5, 25, 7, 49]
indices = [1, 2, 3, 4, 5, 6, 8, 9, 10, 14, 17, 19, 20, 24]
for degree in degrees:
    prime, exponent = p.prime_power(degree)
    b = (degree-1)//(prime-1)
    for n in indices:
        for shape in sorted({1, min(2, degree), degree}):
            certificate = p.polynomial_certificate(n, degree, shape)
            check('literal_prime_power_gamma_newton_certificate', certificate['window'] == 2)
            for k in (0, 1, n//2, n):
                value = s.rf(shape, degree*k)/s.factorial(degree*k)
                check('gamma_shape_multiplier_is_unit_one', p.valuation(Fraction(value), prime) == 0 and
                      p.residue(Fraction(value), prime) == 1)
            if n <= 10 and shape in (1, degree):
                polynomial = p.moment_polynomial(n, degree, shape, x)
                check('independent_prime_power_squarefree_gcd', s.gcd(polynomial, polynomial.diff()).degree() == 0)
                if n > 1:
                    previous = p.moment_polynomial(n-1, degree, shape, x)
                    check('independent_adjacent_moment_gcd', s.gcd(polynomial, previous).degree() == 0)
                    check('adjacent_derivative_identity', n*previous == n*polynomial-s.Poly(x, x)*polynomial.diff())
    examples.append(p.polynomial_certificate(14, degree, degree))

for n, degree, shape in ((5, 3, 1), (7, 4, 4), (8, 9, 3), (6, 5, 5)):
    f = p.moment_polynomial(n, degree, shape, x)
    previous = p.moment_polynomial(n-1, degree, shape, x)
    u, v, gcd = s.gcdex(f, previous)
    check('explicit_adjacent_bezout_identity', u*f+v*previous == gcd and gcd.as_expr() == 1)

# A repeated square divisor is detected; failure of a residual test alone
# does not imply that a characteristic-zero polynomial has a repeated root.
repeated = p.newton_residuals([1, 2, 1], 3)
check('repeated_factor_residual_countercontrol', repeated[0]['residual_squarefree'] is False and
      s.gcd(s.Poly((1+x)**2, x), s.Poly((1+x)**2, x).diff()).degree() == 1)
wild_repeated = p.newton_residuals([9, 0, 0, -6, 0, 0, 1], 3)
check('wild_slope_repeated_residual_detected', wild_repeated[0]['e'] == 3 and not wild_repeated[0]['residual_squarefree'])
sufficient_only = p.newton_residuals([1, 0, 1], 2)
check('failed_residual_test_is_not_non_squarefree_decision', not sufficient_only[0]['residual_squarefree'] and
      s.gcd(s.Poly(1+x*x, x), s.Poly(1+x*x, x).diff()).degree() == 0)
wild_simple = p.polynomial_certificate(3, 3)['newton_edges'][0]
residual = s.Poly(sum(c*T**j for j, c in enumerate(wild_simple['residual_coefficients'])), T, modulus=3)
check('test_T_residual_not_inseparable_Z_power', residual.diff() != s.Poly(0, T, modulus=3)
      and s.Poly(residual.as_expr().subs(T, Z**3), Z, modulus=3).diff().is_zero)

y = s.Symbol('y')
for degree, values in ((2, [1, 2]), (4, [2, -1, 0, 2]), (9, [1, 3, 0, 3])):
    marker = sum(s.Rational(value, factorial(e))*y**e for e, value in enumerate(values))
    certificate = p.marked_class_certificate(degree, values)
    check('accepted_low_degree_marker_class', certificate['status'] == 'two_consecutive_marked_moments_cannot_both_vanish')
    prime, exponent = p.prime_power(degree)
    b = (degree-1)//(prime-1)
    mass = sum(values)
    for n in (1, 2, 5, 8):
        moments = [bridge.factorial(s.expand(marker*y**(degree*k)), [y])/mass for k in range(n+1)]
        polynomial = s.Poly(sum(s.binomial(n, k)*moments[k]*x**k for k in range(n+1)), x)
        scaled = [Fraction(polynomial.nth(k))/prime**(b*k) for k in range(n+1)]
        check('literal_marked_newton_certificate', p.newton_residuals(scaled, prime) == p.expected_edges(n, prime, (-1)**b % prime))
        check('independent_marked_parameter_squarefreeness', s.gcd(polynomial, polynomial.diff()).degree() == 0)
        direct = bridge.factorial(s.expand(marker*(1+x*y**degree)**n), [y])/mass
        check('marked_coefficients_are_actual_original_moments', s.expand(direct-polynomial.as_expr()) == 0)
    examples.append(certificate)

for degree in (2, 3, 4, 5):
    for marker_degree in range(degree):
        values = [0]*marker_degree+[1]
        check('every_monomial_marker_below_exponent_is_accepted', p.marked_class_certificate(degree, values)['window'] == 2)
for degree, values in ((2, [1, 1]), (3, [1, 2]), (4, [1, -1])):
    result = p.marked_class_certificate(degree, values)
    check('cancelling_or_zero_mass_marker_is_outside_only', result['status'] == 'outside_unit_mass_marker_criterion' and result['period_is_zero'] is None)

u = s.Symbol('u')
base = 1+s.I*y**3
coefficient = 2-s.I
f = [bridge.factorial(s.expand(base**n), [y])/s.factorial(n) for n in range(7)]
g = [bridge.factorial(s.expand((base+coefficient*u)**n), [y, u])/s.factorial(n) for n in range(7)]
for n in range(1, 7):
    check('independent_linear_adjunction_original_moments', s.expand(g[n]-coefficient*g[n-1]-f[n]) == 0)

for shape in (1, 2, 3):
    phase = 1+s.I*y**3
    lift, xs, _ = bridge.gaussian_lift(phase, [y], [2*shape])
    for n in range(3):
        check('gamma_gaussian_original_phase_identity', bridge.gaussian(s.expand(lift**n), xs)
              == bridge.gamma(s.expand(phase**n), [y], [shape]))
        check('odd_gaussian_marker_has_separate_zero_branch', bridge.gaussian(s.expand(xs[0]*lift**n), xs) == 0)

for args in ((6, 1), (10, 1), (3, 4), (4, 0), (1, 1), (4.0, 1)):
    try:
        p.class_certificate(*args)
    except ValueError:
        check('outside_exponent_or_shape_fails_closed', True)
    else:
        raise AssertionError('unsupported class accepted')
for args in ((3, [1, 0, 0, 1]), (4, [0.5, 1])):
    try:
        p.marked_class_certificate(*args)
    except ValueError:
        check('outside_marker_degree_or_inexact_input_rejected', True)
    else:
        raise AssertionError('unsupported marker accepted')

# Preserve ordinary pair ideals, including nilpotents and the sharp exponent.
A, B = s.symbols('A B')
for degree, shape in ((2, 1), (3, 1), (4, 4), (9, 3)):
    for n in range(1, 6):
        moments = [sum(s.binomial(k, j)*s.rf(shape, degree*j)*A**(k-j)*B**j
                       for j in range(k+1)) for k in (n, n+1)]
        groebner = s.groebner(moments, A, B, order='grevlex', domain=s.QQ)
        leading = [poly.monoms(order=groebner.order)[0] for poly in groebner.polys]
        standard = [(i, j) for total in range(2*n+1) for i in range(total+1)
                    for j in [total-i] if not any(i >= u and j >= v for u, v in leading)]
        check('ordinary_pair_length_is_complete_intersection_product', len(standard) == n*(n+1))
        check('ordinary_pair_top_degree_and_dimension', max(i+j for i, j in standard) == 2*n-1
              and sum(i+j == 2*n-1 for i, j in standard) == 1)
        check('all_degree_2n_coefficients_vanish', all(groebner.reduce(A**i*B**(2*n-i))[1] == 0
                                                     for i in range(2*n+1)))
        check('universal_phase_nilpotence_exponent_is_sharp', any(groebner.reduce(A**i*B**(2*n-1-i))[1] != 0
                                                                 for i in range(2*n)))
        check('pair_profile_export_matches_actual_quotient', p.polynomial_certificate(n, degree, shape)['ordinary_pair_profile']['length'] == len(standard))

receipt = {'date': '2026-11-16', 'status': 'pass', 'checks': counts, 'total': sum(counts.values()),
           'scope': 'Finite exact controls for geometric-unit Newton residuals, edge lattice points, prime-power/Gamma/marker coefficients, separate gcd/Bezout identities, wild/repeated residual controls, original moments and Gaussian/disjoint-variable transfers. Ordinary Groebner controls check pair lengths, degree2n containments and degree2n-1 survivors. Universal all-n,p,a claims use the written proof.',
           'files': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                     ('PRIME_POWER_TWO_TERM_SFC.tex', 'prime_power_windows.py', 'verify_prime_power_windows.py', 'moment_bridges.py')},
           'independent_replay': False, 'formal_verification': False, 'external_review': False,
           'novelty_claim': False}
(HERE/'prime_power_window_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'prime_power_window_examples.json').write_text(json.dumps(examples, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
