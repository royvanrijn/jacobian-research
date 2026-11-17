"""Exact controls for coupled E-period reduction and the full block limit."""
from pathlib import Path
from math import prod
import hashlib
import json
import mpmath as mp
import sympy as s
import moment_bridges as b
import coupled_block_periods as c

HERE = Path(__file__).resolve().parent
counts, examples = {}, []


def check(group, value):
    if not value:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


def cube_integral(polynomial, variables):
    value = s.expand(polynomial)
    for variable in variables:
        primitive = s.Poly(value, variable).integrate()
        value = s.expand(primitive.eval(1)-primitive.eval(0))
    return s.simplify(value)


u = s.symbols('u:4')
z = s.Symbol('z')
fixtures = [(u[:1], 1), (u[:1], u[0]), (u[:1], 1+u[0]**4),
            (u[:2], 1), (u[:2], u[0]**2*u[1]),
            (u[:2], (1-u[0])**2*(1-u[1])),
            (u[:2], s.sqrt(2)*u[0]+s.I*u[1]**3),
            (u[:2], u[0]-u[1]),
            (u[:3], 1+u[0]**2*u[1]**2*u[2]**2),
            (u[:4], u[0]**3*(1-u[1])*(1+u[2]*u[3]))]
for variables, marker in fixtures:
    reduction = c.cube_reduction(marker, variables, z)
    t, W = reduction['moment_parameter'], reduction['rational_moment']
    maximum_shift = max((p['shift'] for p in reduction['partial_fractions']), default=1)
    ceiling = 9+maximum_shift
    functions = [s.Integer(1), sum(z**n/s.factorial(n) for n in range(ceiling+1))]
    functions += [sum(z**n/(s.factorial(n)*n**j) for n in range(1, ceiling+1))
                  for j in range(1, len(variables))]
    represented = s.expand(sum(a*f for a, f in zip(reduction['coefficients'], functions)))
    for k in range(9):
        direct = cube_integral(marker*prod(variables)**k, variables)
        check('literal_rational_moment_against_direct_cube_integration', s.simplify(W.subs(t, k)-direct) == 0)
        check('E_basis_Laurent_reduction_Taylor_coefficients',
              s.simplify(represented.coeff(z, k)-direct/s.factorial(k)) == 0)
    for k in range(-maximum_shift, 0):
        check('removable_origin_poles_cancel_in_actual_period', s.simplify(represented.coeff(z, k)) == 0)
    for scale in (0, 1, s.I):
        cert = c.period_certificate(marker, variables, scale)
        check('exact_period_decision_never_implies_full_Gamma_zero_sequence',
              not cert['full_Gamma_marked_sequence_zero_claimed']
              and not cert['zero_period_alone_decides_whole_sequence'])
    examples.append(c.period_certificate(marker, variables, 1))

for r in range(1, 8):
    variables = s.symbols('v:'+str(r))
    reduction = c.cube_reduction(1, variables, z)
    expected = [-1/z, 1/z] if r == 1 else [s.Integer(0)]*r+[1/z]
    check('unit_coupled_period_is_J_rminus1_over_z', all(s.simplify(a-bb) == 0
          for a, bb in zip(reduction['coefficients'], expected)))
    for scale in (0, 1, s.I, s.sqrt(2), 1+2*s.I):
        cert = c.period_certificate(1, variables, scale)
        check('all_rank_unit_period_certificates_nonzero', cert['status'] == 'nonzero_period_certified')
        check('unit_period_transcendence_only_at_nonzero_algebraic_scale', cert['transcendence_certified'] == (scale != 0))
    # Verify the actual triangular system against independent Taylor polynomials.
    functions = [s.Integer(1), sum(z**n/s.factorial(n) for n in range(13))]
    functions += [sum(z**n/(s.factorial(n)*n**j) for n in range(1, 13)) for j in range(1, r)]
    matrix = s.zeros(r+1)
    matrix[1, 1] = 1
    if r > 1:
        matrix[2, 0], matrix[2, 1] = -1/z, 1/z
    for j in range(2, r):
        matrix[j+1, j] = 1/z
    for row, function in enumerate(functions):
        difference = s.expand(s.diff(function, z)-sum(matrix[row, col]*f for col, f in enumerate(functions)))
        check('actual_E_function_system_Taylor_identity', all(difference.coeff(z, k) == 0 for k in range(10)))
    check('E_system_has_only_finite_pole_zero', all(s.denom(entry) in (1, z) for entry in matrix))

for r in (2, 3):
    variables = u[:r]
    q = prod(variables)
    marker = 1-2*variables[0]+variables[0]*(1-variables[0])*prod(variables[1:])
    check('actual_Stokes_zero_period_identity', s.expand(1-2*variables[0]
          +s.diff(q, variables[0])*variables[0]*(1-variables[0])-marker) == 0)
    cert = c.period_certificate(marker, variables, 1)
    check('specialized_zero_period_does_not_mean_whole_zero_sequence', cert['status'] == 'zero_period_certified'
          and not cert['cube_marked_power_sequence_identically_zero'] and cube_integral(marker, variables) != 0)
    check('specialized_zero_marker_period_is_not_identically_zero_function',
          c.period_certificate(marker, variables, 2)['status'] == 'nonzero_period_certified')
    examples.append(cert)
algebraic = c.period_certificate(u[0], u[:1], 1)
check('positive_marker_can_have_nonzero_algebraic_period', algebraic['status'] == 'nonzero_period_certified'
      and algebraic['algebraic_value'] == '1' and not algebraic['transcendence_certified'])
antisymmetric = c.period_certificate(u[0]-u[1], u[:2], s.I)
check('exact_cube_zero_sequence_identity_is_distinguished', antisymmetric['cube_marked_power_sequence_identically_zero']
      and antisymmetric['status'] == 'zero_period_certified')
examples.append(algebraic)

y = s.symbols('y:6')
blocks, exponents = ((0, 1), (2, 3)), (2, 1)
F0 = (y[0]+y[1])**2*(y[2]+y[3])
phase = (s.sqrt(2)+s.I)*F0+2*y[0]*y[2]+s.I*y[1]+1
marker = y[0]**2+s.I*y[3]-1
for shapes in ((1, 1, 1, 1), (s.Rational(1, 4),)*2+(s.Rational(1, 6),)*2,
               (s.Rational(1, 2), s.Rational(3, 2), s.Rational(2, 3), s.Rational(1, 3))):
    data = c.block_limit_data(phase, y[:4], blocks, exponents, shapes, marker)
    check('literal_block_top_and_scalar_recognized', data['status'] == 'block_limit_identity_certified'
          and s.expand(data['scalar']-s.sqrt(2)-s.I) == 0)
    T = data['angular_variables']
    check('actual_angular_next_layer_and_marker_top', s.simplify(data['angular_phase']
          -T[0]*T[2]/(s.sqrt(2)+s.I)) == 0 and data['angular_marker'] == 4*T[0]**2)
    radii = s.symbols('s:2')
    substitution = {y[i]: radii[j]*T[i] for j, block in enumerate(blocks) for i in block}
    for n in range(4):
        polynomial = s.expand(marker*phase**n)
        split = s.expand(polynomial.subs(substitution, simultaneous=True))
        for block in blocks:
            split = b.dirichlet(split, tuple(T[i] for i in block), tuple(shapes[i] for i in block))
        split = b.gamma(split, radii, data['block_shapes'])
        direct = b.gamma(polynomial, y[:4], shapes)
        check('full_Gamma_Dirichlet_block_split_actual_moments', s.simplify(split-direct) == 0)
    examples.append(c.block_limit_certificate(phase, y[:4], blocks, exponents, shapes, marker))
check('positive_small_total_shapes_explicitly_in_scope',
      c.block_limit_data(phase, y[:4], blocks, exponents, (s.Rational(1, 4),)*4)['block_shapes'] == (s.Rational(1, 2),)*2)
check('wrong_actual_top_returns_no_FC_decision',
      c.block_limit_certificate(phase+y[0]**3, y[:4], blocks, exponents, (1,)*4)['status']
      == 'outside_block_radial_top_no_decision')
# Coordinate-volume and probability normalizations differ for block size three.
for n in range(1, 4):
    pure = (y[0]+y[1]+y[2])*(y[3]+y[4])
    actual = b.factorial(pure**n, y[:5])
    check('source_coordinate_volume_normalization_retained', actual/(s.factorial(n+2)*s.factorial(n+1)) == s.Rational(1, 2))
    check('same_actual_moment_has_probability_limit_one', actual/(s.rf(3, n)*s.rf(2, n)) == 1)
S1, S2 = y[0]+y[1], y[2]+y[3]
separated_G = 3*s.I*S2**2*y[0]**2+2*(1+s.I)*S1**3*y[2]
data = c.block_limit_data(S1**3*S2**2+separated_G, y[:4], blocks, (3, 2), (1,)*4)
T = data['angular_variables']
check('source_block_separable_angular_product_step', s.expand(data['angular_phase'].subs(
      {T[1]: 1-T[0], T[3]: 1-T[2]}, simultaneous=True)
      -s.I*T[0]**2-(1+s.I)*T[2]) == 0)

for degrees in ((2,), (2, 1), (2, 2, 1)):
    variables = y[:2*len(degrees)]
    D = sum(degrees)
    lower = 1+s.I*variables[0] if D >= 3 else s.Integer(1)
    cert = c.coupled_family_certificate(degrees, 1+s.I, lower, variables=variables)
    F = s.sympify(cert['phase'])
    check('coupled_family_actual_angular_identity_and_nonvanishing', cert['unit_moments_eventually_all_nonzero']
          and cert['unit_period']['status'] == 'nonzero_period_certified'
          and not cert['SFC_claimed'] and not cert['all_complex_scale_nonzero_limit_claimed'])
    top = prod((variables[2*j]+variables[2*j+1])**d for j, d in enumerate(degrees))
    G = c._homogeneous(s.Poly(F, *variables), D-1)
    for power in (2, 3):
        powered = s.Poly(s.expand(F**power), *variables)
        check('literal_block_family_power_top', s.expand(c._homogeneous(powered, power*D)-top**power) == 0)
        check('literal_block_family_power_next_layer', s.expand(c._homogeneous(powered, power*D-1)
              -power*top**(power-1)*G) == 0)
        remainder = s.expand(powered.as_expr()-top**power-power*top**(power-1)*G)
        lifted = c.coupled_family_certificate(tuple(power*d for d in degrees), 1+s.I,
                                               remainder, variables=variables)
        check('powers_preserve_same_coupled_scale_and_family', s.expand(s.sympify(lifted['phase'])-F**power) == 0)
    examples.append(cert)

# Inactive variables: exact original moments, active-degree markers and true couplings.
x, a, z, w = s.symbols('x a z w')
variables = (x, a, z, w)
Q2 = s.I*x*x+(2+s.I)*x*a+3*x*z+z*z+z*w+s.sqrt(2)*a*w+w*w
cubic = x**3+Q2+x+z+w+1
for shapes in ((1,)*4, (s.Rational(1, 4), s.Rational(2, 3), s.Rational(1, 2), s.Rational(5, 4))):
    for H in (1, z**4, x*z*z+z**7, z-shapes[2], x*z+x*x*w):
        cert = c.monomial_top_certificate(cubic, variables, (3, 0, 0, 0), shapes, H)
        check('inactive_cubic_proposed_exponent_Q2_at_3_over_27', s.simplify(s.sympify(cert['unit_limit_exponent'])-s.I/3) == 0)
        check('inactive_indices_and_zero_exponents_are_retained', cert['active_indices'] == [0]
              and cert['inactive_indices'] == [1, 2, 3] and not cert['every_variable_occurs_in_top'])
        hp = s.Poly(H, *variables)
        active_degree = max(powers[0] for powers, coeff in hp.terms())
        independent_top = sum(coeff*s.prod(v**p for v, p in zip(variables, powers))
                              for powers, coeff in hp.terms() if powers[0] == active_degree)
        independent_amplitude = b.gamma(independent_top.subs(x, 3), variables[1:], shapes[1:])
        check('inactive_marker_uses_active_degree_and_actual_Gamma_expectation',
              cert['marker_degree'] == active_degree and s.simplify(s.sympify(cert['leading_marker_amplitude'])-independent_amplitude) == 0)
        data = c.block_limit_data(cubic, variables, ((0,),), (3,), shapes, H)
        for n in range(3):
            literal = s.expand(H*cubic**n)
            first_inactive = b.gamma(literal, variables[1:], shapes[1:])
            split = b.gamma(first_inactive, (x,), shapes[:1])
            check('inactive_Gamma_Fubini_matches_full_original_moments', s.simplify(split-b.gamma(literal, variables, shapes)) == 0)
        examples.append(cert)
# Inactive indices need not be last, and highest-degree extra terms must be rejected.
permuted = c.monomial_top_certificate(cubic, (z, x, a, w), (0, 3, 0, 0), (1,)*4)
check('inactive_order_is_arbitrary', permuted['active_indices'] == [1] and s.sympify(permuted['unit_limit_exponent']) == s.I/3)
check('inactive_top_degree_term_is_outside_the_theorem',
      c.monomial_top_certificate(cubic+z**3, variables, (3, 0, 0, 0), (1,)*4)['status']
      == 'outside_block_radial_top_no_decision')
binary = c.monomial_top_certificate(x*x*a+Q2, variables, (2, 1, 0, 0), (1,)*4)
check('inactive_binary_top_has_the_right_normalization_exponent',
      s.simplify(s.sympify(binary['unit_limit_exponent'])-Q2.subs({x:2, a:1, z:0, w:0})/4) == 0)
linear = c.monomial_top_certificate(2*x+3*s.I, (x, z), (1, 0), (s.Rational(1, 2), s.Rational(1, 4)), z*z)
check('degree_one_and_unchanged_inactive_marker_mass', s.sympify(linear['unit_limit_exponent']) == 3*s.I/2
      and s.sympify(linear['leading_marker_amplitude']) == s.Rational(5, 16))
zero_leading = c.monomial_top_certificate(x**3+x*z, (x, z), (3, 0), (1, 1), z-1)
check('inactive_zero_mean_marker_can_have_nonzero_actual_moments',
      not zero_leading['marked_eventual_nonzero_certified'] and b.factorial((z-1)*(x**3+x*z), (x,z)) == 1
      and not zero_leading['zero_leading_amplitude_decides_whole_sequence'])
# Two angular blocks, subunit active shape totals, and two inactive shapes.
all_vars = y[:4]+(z,w)
shapes = (s.Rational(1,4),)*2+(s.Rational(1,6),)*2+(s.Rational(1,3),s.Rational(1,2))
F = (y[0]+y[1])**2*(y[2]+y[3])+2*y[0]*y[2]+z*y[1]+w*w+z*w+y[0]+w
H = y[0]*z*z+y[1]*w**3+z**7
info = c.block_limit_data(F, all_vars, blocks, exponents, shapes, H)
T = info['angular_variables']
check('coupled_inactive_next_layer_disappears_exactly', info['angular_phase'] == T[0]*T[2])
check('coupled_inactive_marker_retains_Gamma_moments', info['marker_degree'] == 1
      and s.expand(info['angular_marker']-s.Rational(8,9)*T[0]-s.Rational(15,4)*T[1]) == 0)
radii = s.symbols('radius:2')
substitution = {y[i]: radii[j]*T[i] for j, block in enumerate(blocks) for i in block}
for n in range(4):
    direct = b.gamma(H*F**n, all_vars, shapes)
    split = s.expand((H*F**n).subs(substitution, simultaneous=True))
    for block in blocks:
        split = b.dirichlet(split, tuple(T[i] for i in block), tuple(shapes[i] for i in block))
    split = b.gamma(split, radii+(z,w), info['block_shapes']+shapes[4:])
    check('coupled_inactive_full_Gamma_Dirichlet_split', s.simplify(split-direct) == 0)
K = z*y[1]+w*w+z*w
extended = c.coupled_family_certificate((2,1), s.I, y[0]+w, H, y[:4], (z,w), shapes[4:], K)
base = c.coupled_family_certificate((2,1), s.I, variables=y[:4])
check('coupled_exporter_integrates_inactive_marker_before_period_decision',
      extended['marker_degree'] == 1 and extended['inactive_expectation_taken']
      and extended['unit_period']['specialized_coefficients'] == base['unit_period']['specialized_coefficients'])
for power in (2,3):
    actual = s.sympify(extended['phase'])
    top = F0
    G = c._homogeneous(s.Poly(actual, *all_vars), 2)
    next_K = power*top**(power-1)*K
    remainder = s.expand(actual**power-top**power-power*top**(power-1)*G)
    lift = c.coupled_family_certificate((2*power,power), s.I, remainder, variables=y[:4],
                                       inactive_variables=(z,w), inactive_shapes=shapes[4:], coupling=next_K)
    check('coupled_inactive_family_has_exact_power_closure', s.expand(s.sympify(lift['phase'])-actual**power) == 0)
# The inactive Gaussian lift includes its original odd and even marker moments.
lifted, xs_inactive, _ = b.gaussian_lift(x**3+x*z+s.I*x*x, (x,z))
for n in range(3):
    check('inactive_original_Gaussian_marker_is_preserved',
          s.simplify(b.gaussian(xs_inactive[2]**2*lifted**n, xs_inactive)
                     -b.factorial(z*(x**3+x*z+s.I*x*x)**n, (x,z))) == 0)
    check('inactive_original_odd_Gaussian_marker_average_is_zero',
          b.gaussian(xs_inactive[2]*lifted**n, xs_inactive) == 0)

# Original Gaussian markers, including exact zero averages.
for shapes in ((1, 1, 1), (s.Rational(1, 4), s.Rational(1, 3), s.Rational(1, 2))):
    phase3 = y[0]*y[1]*y[2]+s.I*y[0]**2+2*y[1]*y[2]+y[0]-1
    singleton = c.monomial_top_certificate(phase3, y[:3], (1, 1, 1), shapes, 1+y[0])
    check('single_top_monomial_all_complex_lower_terms', singleton['unit_moments_eventually_all_nonzero']
          and s.sympify(singleton['unit_limit_exponent']) == 2+s.I
          and singleton['leading_marker_amplitude'] == '1')
    check('single_top_monomial_marker_leading_zero_is_not_sequence_zero',
          not c.monomial_top_certificate(phase3, y[:3], (1, 1, 1), shapes,
                                         y[0]-y[1])['marked_eventual_nonzero_certified'])
    examples.append(singleton)
check('wrong_positive_exponents_do_not_match_a_top_with_inactive_variables',
      c.monomial_top_certificate(y[0]*y[1]+y[2], y[:3], (1, 1, 1), (1, 1, 1))['status']
      == 'outside_block_radial_top_no_decision')

cert = c.coupled_family_certificate((2, 1), s.I, 1+y[0], variables=y[:4])
F = s.sympify(cert['phase'])
gaussian_phase, xs, gaussian_blocks = b.gaussian_lift(F, y[:4])
for original_marker in (xs[0]**2+s.I*xs[4]**2-1, xs[0]*xs[2], xs[0]**2-xs[1]**2):
    average = s.Integer(0)
    for powers, amplitude in s.Poly(original_marker, *xs).terms():
        radial, factor = b.radial_monomial(powers, (2,)*4)
        if factor:
            average += amplitude*factor*prod(v**a for v, a in zip(y[:4], radial))
    for n in range(3):
        check('original_Gaussian_markers_match_actual_Gamma_moments',
              s.simplify(b.gaussian(original_marker*gaussian_phase**n, xs)
                         -b.gamma(average*F**n, y[:4], (1,)*4)) == 0)
    leading = c.coupled_family_certificate((2, 1), s.I, 1+y[0], average, y[:4])
    check('leading_marker_zero_is_not_full_sequence_decision',
          not leading['zero_leading_marker_period_decides_full_sequence'])
    if average == 0:
        check('original_odd_and_angular_zero_averages_retained', all(
            b.gaussian(original_marker*gaussian_phase**n, xs) == 0 for n in range(3)))
check('actual_Gaussian_dimension_is_four_per_angular_block', len(xs) == cert['Gaussian_coordinate_count'] == 8)

# A nonzero homogeneous marker can have zero LEADING period and nonzero moments.
rho = (y[0]+y[1])/4+(y[2]+y[3])/2
H = s.expand(rho**3-y[0]*rho**2+y[0]*y[2]*rho/2-y[0]**2*y[2]/4)
leading_zero = c.coupled_family_certificate((2, 1), 1, marker=H, variables=y[:4])
check('nonzero_marker_zero_leading_period_is_not_full_sequence_identity',
      leading_zero['leading_marker_period']['status'] == 'zero_period_certified'
      and not leading_zero['whole_marked_sequence_zero_certified']
      and not leading_zero['zero_leading_marker_period_decides_full_sequence']
      and b.factorial(H*s.sympify(leading_zero['phase']), y[:4]) != 0)
examples.append(leading_zero)

# An exact finite formula checks normalization before a numerical limit diagnostic.
plain = c.coupled_family_certificate((2, 1), 1, variables=y[:4])
F = s.sympify(plain['phase'])
for n in range(1, 5):
    finite = sum(s.binomial(n, k)*2**k/s.rf(2*n+2-k, k)/(k+1)**2 for k in range(n+1))
    direct = b.factorial(F**n, y[:4])/(s.rf(2, 2*n)*s.rf(2, n))
    check('literal_normalized_moment_formula_before_limit', s.simplify(finite-direct) == 0)
xi = s.Symbol('xi')
bare = (y[0]+y[1])**2*(y[2]+y[3])+2*xi*y[0]*y[2]
mu1, mu2 = b.factorial(bare, y[:4]), b.factorial(bare**2, y[:4])
check('ordinary_first_two_moment_unit_witness_in_bare_subfamily', s.expand(mu2-(8*xi+24)*mu1) == 432)
mp.mp.dps = 50
diagnostics = []
for r in (1, 2, 3):
    target = mp.fsum(mp.mpf(1)/(mp.factorial(k)*(k+1)**r) for k in range(160))
    samples = []
    for n in (8, 32, 96):
        finite = sum(s.binomial(n, k)*2**k/s.rf(2*n+2-k, k)/(k+1)**r for k in range(n+1))
        value = mp.mpf(str(s.N(finite, 52)))
        samples.append(dict(index=n, normalized_moment=mp.nstr(value, 25),
                            absolute_error=mp.nstr(abs(value-target), 12)))
    diagnostics.append(dict(angular_dimension=r, scale=1, unit_limit=mp.nstr(target, 25), samples=samples,
                            role='numerical diagnostic only; not universal nonvanishing evidence'))

invalid = [lambda: c.period_certificate(1, u[:2], 0.1),
           lambda: c.period_certificate(1, u[:2], 2*s.pi*s.I),
           lambda: c.period_certificate(s.exp(1), u[:2], 1),
           lambda: c.period_certificate(s.Symbol('unbound'), u[:2], 1),
           lambda: c.cube_reduction(1, ()),
           lambda: c.block_limit_data(phase, y[:4], ((0, 1), (1, 2, 3)), (2, 1), (1,)*4),
           lambda: c.block_limit_data(phase, y[:4], blocks, exponents, (0, 1, 1, 1)),
           lambda: c.block_limit_data(phase, y[:4], blocks, exponents, (0.5, 1, 1, 1)),
           lambda: c.coupled_family_certificate((1, 2), 1),
           lambda: c.coupled_family_certificate((2, 1), 1, y[0]**2, variables=y[:4]),
           lambda: c.coupled_family_certificate((True, 1), 1),
           lambda: c.block_limit_data(x+z, (x,z), (), (), (1,1)),
           lambda: c.block_limit_data(x+z, (x,z), ((2,),), (1,), (1,1)),
           lambda: c.monomial_top_certificate(x+z, (x,z), (0,0), (1,1)),
           lambda: c.monomial_top_certificate(x, (x,z), (1,), (1,1)),
           lambda: c.monomial_top_certificate(x, (x,z), (1,-1), (1,1)),
           lambda: c.coupled_family_certificate((2,1), 1, variables=y[:4], inactive_variables=(y[0],)),
           lambda: c.coupled_family_certificate((2,1), 1, variables=y[:4], inactive_variables=(z,), coupling=z**3),
           lambda: c.coupled_family_certificate((2,1), 1, variables=y[:4], inactive_variables=(z,), coupling=y[0]**2),
           lambda: c.coupled_family_certificate((2,1), 1, variables=y[:4], inactive_variables=(z,), inactive_shapes=(0,))]
for operation in invalid:
    try:
        operation()
    except ValueError:
        check('outside_exact_domain_fails_closed', True)
    else:
        raise AssertionError('invalid input accepted')

draft_receipt = json.loads((HERE/'../../jacobian/cross-programme-connections/check_inactive_variable_factorial_extension_20261116.json').read_text())
check('retained_inactive_draft_50_transcription_checks_pass',
      draft_receipt['exact_checks'] == 50 and draft_receipt['all_exact_checks_pass'])

files = ['COUPLED_BLOCK_PERIODS.tex', 'coupled_block_periods.py', 'verify_coupled_block_periods.py',
         'moment_bridges.py', '../../jacobian/cross-programme-connections/check_inactive_variable_factorial_extension_20261116.py',
         '../../jacobian/cross-programme-connections/check_inactive_variable_factorial_extension_20261116.json', 'sources/long-eit-vi-current.tex', 'sources/beukers-math0405549.html']
receipt = dict(status='pass', checks=sum(counts.values()), groups=counts,
               files={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
               universal_proofs='COUPLED_BLOCK_PERIODS.tex',
               mathematical_block_shape_domain='all positive real', executable_block_shape_domain='positive rational',
               coupled_period_scale_domain='all algebraic scales; zero handled separately',
               coupled_family_complex_scale_FC=True, coupled_family_complex_scale_eventually_all_nonzero_claimed=False,
               SFC_claimed=False, independent_replay=False, external_review=False, formal_verification=False,
               novelty_claim=False, inactive_variable_extension_checked=True, general_FC_and_SFC_closed=False)
(HERE/'coupled_block_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'coupled_block_examples.json').write_text(json.dumps(dict(certificates=examples,
    limit_diagnostics=diagnostics), indent=2)+'\n')
print(json.dumps(receipt, indent=2))
