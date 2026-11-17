"""Finite exact controls for all-index window identities and source boundaries."""
from pathlib import Path
import hashlib
import itertools
import json
import sympy as s
import moment_bridges as b
import moment_windows as w
import two_angular as a

HERE = Path(__file__).resolve().parent
counts = {}
examples = []


def check(group, result):
    if not result:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


cases = [([1, -1], [1, 1], 1), ([1, -1], [1, 1], 2),
         ([1, 1, -1], [1, 2, 1], 3), ([0, 2], [1, 2], 2),
         ([1, s.I], [2, 1], 2), ([s.sqrt(2), -s.sqrt(2)], [1, 1], 2),
         ([1, 1, -1], [s.Rational(1, 2), s.Rational(1, 2), 1], 2),
         ([0, 1], [s.Rational(1, 2), 2], 1)]
for nodes, shapes, power in cases:
    ys = s.symbols('y:'+str(len(nodes)))
    linear = sum(c*y for c, y in zip(nodes, ys))
    for marker in (s.S.One, ys[0]-ys[-1], ys[0]**2+2*ys[-1]+s.I):
        cert = w.linear_power_window(nodes, shapes, power, marker, ys)
        actual = w.coefficients(cert, 5)
        for m in range(6):
            expected = b.gamma(s.expand(marker*linear**(power*m)), ys, shapes)/s.factorial(power*m)
            check('rational_section_vs_direct_original_gamma_moments', s.simplify(actual[m]-expected) == 0)
        check('cancelled_rational_function_is_proper_or_zero', cert['numerator'] == 0 or
              s.degree(cert['numerator'], cert['variable']) < s.degree(cert['denominator'], cert['variable']))
        if marker == 1:
            check('unit_sequence_is_nonzero', cert['status'] == 'nonzero_in_every_certified_window')
            check('unit_window_within_active_total_shape', cert['window'] <= sum(shapes[i] for i, c in enumerate(nodes) if c != 0))
        examples.append(w.export(cert))

# A root-of-unity example attains the complete consecutive-window bound.
root = (-1+s.I*s.sqrt(3))/2
for power, period in ((1, 3), (2, 3), (3, 1)):
    cert = w.linear_power_window([1, root, s.conjugate(root)], [1, 1, 1], power)
    t = cert['variable']
    check('cyclotomic_section_exact_sharp_example', s.expand(cert['denominator']-(1-t**period)) == 0 and cert['numerator'] == 1)
    check('cyclotomic_window', cert['window'] == period)

ys = s.symbols('y:2')
zero = w.linear_power_window([1, -1], [1, 1], 2, ys[0]-ys[1], ys)
check('nonzero_marker_can_have_identically_zero_selected_sequence', zero['status'] == 'identically_zero_marked_sequence')
phase_zero = w.linear_power_window([0, 0], [s.Rational(1, 2), 2], 3, ys[0], ys)
check('zero_phase_has_separate_finite_prefix', phase_zero['status'] == 'zero_phase'
      and phase_zero['constant_marked_moment'] == s.Rational(1, 2) and phase_zero['window'] is None)
check('zero_marker_is_certified_identically', w.linear_power_window([1], [1], marker=0)['status'] == 'identically_zero_marked_sequence')

# Every original simplex marker is retained by homogenizing on sum(T)=1.
ts = s.symbols('t:3')
S = sum(ts)
linear = ts[0]-ts[1]+s.I*ts[2]
original_marker = ts[0]-ts[1]**2
lifted_marker = ts[0]*S-ts[1]**2
cert = w.linear_power_window([1, -1, s.I], [1, 2, 1], 2, lifted_marker, ts)
coefficients = w.coefficients(cert, 4)
for m in range(5):
    gamma_marked = coefficients[m]*s.factorial(2*m)
    simplex = b.dirichlet(s.expand(original_marker*linear**(2*m)), ts, [1, 2, 1])
    check('all_original_simplex_markers_by_homogenization', s.simplify(gamma_marked-s.rf(4, 2*m+2)*simplex) == 0)

linear = ys[0]-ys[1]
marker = ys[0]**2+3*ys[1]
cert = w.linear_power_window([1, -1], [1, 1], 2, marker, ys)
lift, xs, blocks = b.gaussian_lift(linear**2, ys)
gaussian_marker, _, _ = b.gaussian_lift(marker, ys)
for m in range(4):
    expected = b.gaussian(s.expand(gaussian_marker*lift**m), xs)
    check('gaussian_original_marker_translation', s.simplify(expected-s.factorial(2*m)*w.coefficients(cert, m)[m]) == 0)

# The affine draft typo is false; the proof's longer window identity is exact.
y = s.Symbol('y')
check('rocks_affine_display_counterexample', b.factorial(y-1, [y]) == 0 and s.Poly(y-1, y).is_zero is False)
check('rocks_correct_affine_window_countercontrol', b.factorial((y-1)**2, [y]) == 1)
for nodes, shapes in [([1, -1], [1, 1]), ([2, 2], [1, 2])]:
    ys = s.symbols('v:2')
    L = sum(c*v for c, v in zip(nodes, ys))
    z = s.Symbol('z')
    P = s.Poly(s.prod((1-c*z)**alpha for c, alpha in zip(nodes, shapes)), z)
    constant = s.I
    values = [b.gamma(s.expand((constant+L)**k), ys, shapes)/s.factorial(k) for k in range(10)]
    for k in range(P.degree(), 10):
        result = sum(P.nth(j)*values[k-j] for j in range(P.degree()+1))
        check('affine_window_inhomogeneous_recurrence', s.expand(result-constant**k/s.factorial(k)) == 0)

# Independent Gamma addition contributes a degree-a difference polynomial.
v, u, z = s.symbols('v u z')
base = v**2+s.I*v+1
for alpha in (1, 2, 3):
    coefficient = 2+s.I
    f = [b.gamma(s.expand(base**k), [v], [1])/s.factorial(k) for k in range(9)]
    g = [b.gamma(s.expand((base+coefficient*u)**k), [v, u], [1, alpha])/s.factorial(k) for k in range(9)]
    D = s.Poly((1-coefficient*z)**alpha, z)
    for k in range(alpha, 9):
        check('new_independent_gamma_difference_identity', s.expand(sum(D.nth(j)*g[k-j] for j in range(alpha+1))-f[k]) == 0)

for prime in (2, 3, 5):
    for beta in (1, 2, 3):
        slope = (prime**beta-1)//(prime-1)
        for k in range(6):
            value = s.factorial(prime**beta*k)/(prime**(slope*k)*s.factorial(k))
            check('prime_power_factorial_ratio_integer', value.is_Integer is True)
            check('wilson_scaled_factorial_unit', int(value) % prime == (-1)**(slope*k) % prime)

# Full homogeneous coefficient reduction, not an F_p-rational point scan.
x0, x1 = s.symbols('x0 x1')
for prime, beta in ((2, (1, 2)), (3, (1, 1)), (5, (1, 2))):
    slopes = [(prime**e-1)//(prime-1) for e in beta]
    for n in range(1, 6):
        scaled = sum(s.factorial(prime**beta[0]*j)*s.factorial(prime**beta[1]*(n-j)) /
                     (s.factorial(j)*s.factorial(n-j)*(-prime)**(slopes[0]*j+slopes[1]*(n-j)))
                     *x0**j*x1**(n-j) for j in range(n+1))
        check('common_prime_normalized_coefficients_integral', all(c.is_Integer for c in s.Poly(scaled, x0, x1).coeffs()))
        check('common_prime_reduction_is_complete_homogeneous', all(int(c) % prime == 1 for c in s.Poly(scaled, x0, x1).coeffs()))

for degree in range(1, 6):
    polynomials = [sum(s.binomial(n, j)*s.factorial(degree*j)*z**j for j in range(n+1)) for n in range(9)]
    for n in range(1, 9):
        check('two_term_adjacent_derivative_identity', s.expand(n*polynomials[n-1]-n*polynomials[n]+z*s.diff(polynomials[n], z)) == 0)
        repeated = s.gcd(polynomials[n], s.diff(polynomials[n], z))
        adjacent = s.gcd(polynomials[n], polynomials[n-1])
        check('adjacent_common_roots_equal_repeated_roots', s.Poly(repeated, z).monic() == s.Poly(adjacent, z).monic())

for bad_nodes, bad_shapes, bad_power in [([1], [s.Rational(1, 2)], 1), ([1], [0], 1),
                                        ([1.0], [1], 1), ([s.pi], [1], 1), ([1], [1], 0)]:
    try:
        w.linear_power_window(bad_nodes, bad_shapes, bad_power)
    except ValueError:
        check('unsupported_window_inputs_fail_closed', True)
    else:
        raise AssertionError('unsupported input accepted')

# Retain the failed blanket mod-five upgrade of the cubic phase.
def mod5(poly):
    return s.Poly.from_dict({m: int(c.p)*pow(int(c.q), -1, 5) % 5
                            for m, c in s.Poly(poly, a.h, domain=s.QQ).terms()}, a.h, domain=s.GF(5))


modular_exceptions = []
for exponents in itertools.product(range(5), repeat=3):
    cert = a.normal_form(a.barycentric_amplitude(exponents))
    ordinary_gcd = s.Poly(s.gcd(cert['a'], cert['b']), a.h)
    check('bounded_cubic_rational_gcd_not_universal_evidence', ordinary_gcd.degree() == 0)
    gcd = s.gcd(mod5(cert['a']), mod5(cert['b']))
    if gcd.degree() > 0:
        modular_exceptions.append({'exponents': list(exponents), 'common_factor_mod5': str(gcd.as_expr())})
check('blanket_all_residue_mod5_upgrade_is_false', len(modular_exceptions) == 12)
matrices = []
for amplitude, determinant in ((a.u, 3), (a.v, 3), (a.w, 1)):
    matrix = s.Matrix([[mod5(a.normal_form(a.X**j*amplitude**5)[key]).as_expr()
                        for j in (0, 1)] for key in ('a', 'b')])
    check('fifth_power_matrices_invertible_despite_residue_exceptions', mod5(matrix.det()).as_expr() == determinant if determinant <= 2 else mod5(matrix.det()).as_expr() == -2)
    matrices.append([[str(value) for value in row] for row in matrix.tolist()])
probe = {'date': '2026-11-16', 'scope': 'only the125 residue triples0..4 in each exponent',
         'rational_gcds_in_tested_box': 'all constant; no all-shape extrapolation',
         'mod5_exception_count': len(modular_exceptions), 'mod5_exceptions': modular_exceptions,
         'fifth_power_matrices_mod5': matrices,
         'conclusion': 'The all-residue nonzero-vector premise fails, so invertible fifth-power matrices alone do not upgrade the current strict-negative-valuation theorem.'}
(HERE/'two_angular_residue_probe.json').write_text(json.dumps(probe, indent=2)+'\n')

receipt = {'date': '2026-11-16', 'status': 'pass', 'checks': counts, 'total': sum(counts.values()),
           'scope': 'Finite exact controls for rational section identities, original Gamma/simplex/Gaussian markers, affine windows, independent-variable differences, common-prime congruences, adjacent two-term roots and failed cubic mod-five upgrade. All-index assertions use the written proofs.',
           'files': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                     ('SFC_WINDOWS_AND_ROCKS_AUDIT.tex', 'moment_windows.py', 'verify_moment_windows.py',
                      'moment_bridges.py', 'two_angular.py')},
           'independent_replay': False, 'formal_verification': False, 'external_review': False,
           'novelty_claim': False}
(HERE/'moment_window_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'moment_window_examples.json').write_text(json.dumps(examples, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
