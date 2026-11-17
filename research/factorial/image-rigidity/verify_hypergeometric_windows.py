"""Exact controls for the written all-index Gamma/Dirichlet window proof."""
from pathlib import Path
from math import prod
import hashlib
import json
import sympy as s
import moment_bridges as bridges
import hypergeometric_windows as h

HERE = Path(__file__).resolve().parent
counts, examples = {}, []


def check(group, value):
    if not value:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


def euler_operator(operator, polynomial, x):
    value, result = polynomial, s.Integer(0)
    for j in range(operator.degree()+1):
        result += operator.nth(j)*value
        value = s.expand(x*s.diff(value, x))
    return s.expand(result)


def ordinary_coefficients(P, B, x):
    order = max(P.degree(), B.degree())
    return [s.expand(x**j*sum(s.functions.combinatorial.numbers.stirling(k, j, kind=2)
                             *(P.nth(k)-x*B.nth(k)) for k in range(j, order+1)))
            for j in range(order+1)]


y = s.symbols('y:3')
x, u = s.symbols('x u')
fixtures = [
    h.make_model((6,), (1,), variables=y[:1]),
    h.make_model((6,), (s.Rational(13, 2),), y[0]**2-s.Rational(195, 4), y[:1]),
    h.make_model((2, 1, 0), (s.Rational(1, 2), s.Rational(3, 2), s.Rational(2, 3)),
                 y[0]**2+s.I*y[1]-3*y[2]/2+1, y),
    h.make_model((1, 1), (1, 1), y[0]-y[1], y[:2]),
    h.make_model((1,), (1,), y[0]**3-9*y[0]**2+18*y[0]-6, y[:1]),
    h.make_model((2,), (s.Rational(1, 2),), s.sqrt(2)*y[0]+s.I*s.sqrt(3), y[:1]),
    h.make_model((6, 0), (s.Rational(1, 2), s.Rational(1, 3)), variables=y[:2], law='dirichlet'),
    h.make_model((2, 1, 0), (1, s.Rational(3, 2), s.Rational(2, 3)),
                 1+y[0]-s.I*y[1]*y[2], y, 'dirichlet'),
    h.make_model((1, 1), (1, 1), y[0]-y[1], y[:2], 'dirichlet'),
    h.make_model((2, 0), (1, 1), y[0]+y[1]-1, y[:2], 'dirichlet'),
    h.make_model((1, 1), (1, 1), (y[0]-y[1])**2, y[:2], 'dirichlet'),
    h.make_model((2, 0), (1, 1), y[0]-s.Rational(1, 2), y[:2], 'dirichlet'),
    h.make_model((1, 1, 1), (s.Rational(1, 2), s.Rational(3, 4), 2),
                 1-2*y[0]+s.sqrt(2)*y[1]**2, y, 'dirichlet'),
    h.make_model((1, 0), (s.Rational(1, 2), s.Rational(3, 2)), y[0]+y[1], y[:2], 'dirichlet'),
]

for model in fixtures:
    direct = bridges.gamma if model.law == 'gamma' else bridges.dirichlet
    Z = prod(v**d for v, d in zip(model.variables, model.exponents))
    for k in range(7):
        value = direct(model.marker*Z**k, model.variables, model.shapes)
        check('marker_polynomial_against_direct_monomial_integration', s.simplify(value-h.coefficient(model, k)) == 0)
    for n in (0, 1, 3, 6):
        value = direct(model.marker*(2+s.I*Z)**n, model.variables, model.shapes)
        check('affine_phase_against_direct_expansion', s.simplify(value-h.moment(model, n, 2, s.I)) == 0)
    certificate = h.window_certificate(model, 4)
    examples.append(certificate)
    if model.R.is_zero:
        check('zero_marker_sequence_identity_retains_index_zero', certificate['status'] == 'marked_zero_sequence_identity'
              and certificate['includes_index_zero'] and all(h.moment(model, n, 2, s.I) == 0 for n in range(8)))
        continue
    L = model.degree+model.marker_degree+1
    check('all_starting_index_window_is_not_two_term_sfc', certificate['window_length'] == L
          and certificate['terminal'] == 4+L-1 and certificate['all_starting_indices']
          and certificate['two_term_sfc_claimed'] is False)
    for n in sorted({0, 1, 3, L+1}):
        f = h.moment_polynomial(model, n, x).as_expr()
        P, B = h.operator_polynomials(model, n)
        check('literal_euler_differential_equation', s.simplify(euler_operator(P, f, x)-x*euler_operator(B, f, x)) == 0)
        value = f
        for j in range(min(n, 4)+1):
            check('adjacent_moments_equal_falling_euler_derivatives',
                  s.simplify(value-s.ff(n, j)*h.moment_polynomial(model, n-j, x).as_expr()) == 0)
            value = s.expand((n-j)*value-x*s.diff(value, x))
        if n >= model.marker_degree:
            check('moment_polynomial_nonzero_after_marker_degree', f != 0)
    n = L+3
    P, B = h.operator_polynomials(model, n)
    a = ordinary_coefficients(P, B, x)
    check('exact_differential_order', len(a) == L+1)
    qD, re = model.Q.LC(), model.R.LC()
    if model.law == 'gamma':
        check('Gamma_highest_ordinary_derivative_coefficient', s.simplify(a[L]-qD*re*x**(L+1)) == 0)
    else:
        star = -s.Rational(model.degree**model.degree, qD)
        check('Dirichlet_highest_ordinary_derivative_coefficient',
              s.simplify(a[L]-re*x**L*(model.degree**model.degree+qD*x)) == 0)
        check('Dirichlet_extra_singularity_simple', s.diff(a[L], x).subs(x, star) != 0)
        indicial = L-1-a[L-1].subs(x, star)/s.diff(a[L], x).subs(x, star)
        check('Dirichlet_indicial_exponent_from_ordinary_coefficients', s.simplify(indicial-n-model.beta) == 0)
        check('Dirichlet_indicial_exponent_exceeds_degree_for_all_N', model.beta > 0)
        for order in (L-1, L, L+2):
            test = (x-star)**order
            image = s.expand(sum(aj*s.diff(test, x, j) for j, aj in enumerate(a)))
            leading = s.Poly(image.subs(x, star+u), u).nth(order-L+1)
            expected = s.ff(order, L-1)*s.diff(a[L], x).subs(x, star)*(order-n-model.beta)
            check('singular_root_leading_term_exact', s.simplify(leading-expected) == 0)

laguerre = fixtures[4]
check('short_zero_prefix_does_not_mean_whole_zero_sequence',
      s.expand(laguerre.R.as_expr()-laguerre.t*(laguerre.t-1)*(laguerre.t-2)) == 0
      and all(h.moment_polynomial(laguerre, n, x).is_zero for n in range(3))
      and not h.moment_polynomial(laguerre, 3, x).is_zero)
check('Dirichlet_marker_cancellation_reduces_window', fixtures[10].homogenization_degree == 2
      and fixtures[10].marker_degree == 1 and h.window_certificate(fixtures[10])['window_length'] == 4)
check('Dirichlet_homogenized_constant_marker_same_moments', all(
    h.moment(fixtures[-1], n, 1, 2) == h.moment(h.make_model((1, 0), fixtures[-1].shapes,
    variables=y[:2], law='dirichlet'), n, 1, 2) for n in range(7)))

mass_zero = h.make_model((6,), (s.Rational(1, 2),), y[0]-s.Rational(1, 2), y[:1])
check('constant_phase_zero_mass_exception_is_explicit', not mass_zero.R.is_zero
      and h.coefficient(mass_zero, 0) == 0 and all(h.moment(mass_zero, n, 3, 0) == 0 for n in range(7)))
check('zero_constant_term_phase_uses_integer_R_roots', all(
    h.moment(mass_zero, n, 0, s.I) == s.I**n*h.coefficient(mass_zero, n) for n in range(7)))

# Original Gaussian markers are averaged before using the Gamma model.
sizes = (1, 2)
gaussian_phase, xs, blocks = bridges.gaussian_lift(y[0]**2*y[1], y[:2], sizes)
markers = [xs[0]**2*xs[1]**2+xs[0]*xs[1]+s.I*xs[2]**4-3,
           xs[0]*xs[1], xs[1]**2-xs[2]**2, xs[0]**2-1]
for marker in markers:
    averaged = s.Integer(0)
    for powers, c in s.Poly(marker, *xs).terms():
        radial, factor = bridges.radial_monomial(powers, sizes)
        if factor:
            averaged += c*factor*prod(v**j for v, j in zip(y[:2], radial))
    model = h.make_model((2, 1), (s.Rational(1, 2), 1), s.expand(averaged), y[:2])
    for n in range(6):
        direct = bridges.gaussian(marker*(1+s.I*gaussian_phase)**n, xs)
        check('original_Gaussian_marker_matches_exported_sequence', s.simplify(direct-h.moment(model, n, 1, s.I)) == 0)
    if marker in markers[1:3]:
        check('odd_and_angular_Gaussian_markers_remain_exact_zero_sequences', model.R.is_zero)

# Ordinary ideals: a regular subideal has the advertised Hilbert series;
# the full window ideal can be strictly smaller in quotient dimension.
U, V = s.symbols('U V')
ideal_fixtures = [
    (h.make_model((1,), (s.Rational(1, 2),)), 2),
    (h.make_model((6,), (s.Rational(1, 2),)), 1),
    (h.make_model((2, 1), (s.Rational(1, 2), 2)), 2),
    (h.make_model((2, 1), (1, 1), law='dirichlet'), 2),
]
for model, start in ideal_fixtures:
    cert = h.ordinary_ideal_certificate(model, start)
    examples.append(cert)
    N, D = start+model.degree, model.degree
    moments = [h.moment(model, j, U, V) for j in range(start, N+1)]
    value = cert['selected_integer']
    G = s.expand(sum(value**(j-1)*U**(D-j)*moments[j] for j in range(1, D+1)))
    check('regular_subideal_is_actual_window_combination', s.expand(G.subs({U: 1, V: x})-s.sympify(cert['dehomogenized_second'])) == 0)
    check('regular_subideal_projective_gcd', s.gcd(moments[0], G) == 1)
    gb = s.groebner(moments, U, V, domain=s.QQ)
    ci = s.groebner([moments[0], G], U, V, domain=s.QQ)
    exponent = cert['window_nilpotence_exponent_upper_bound']
    for j in range(exponent+1):
        monomial = U**j*V**(exponent-j)
        check('ordinary_window_ideal_contains_full_power', gb.reduce(monomial)[1] == 0)
        check('ordinary_regular_subideal_contains_full_power', ci.reduce(monomial)[1] == 0)
    leading = [p.LM(order='lex').exponents for p in ci.polys]
    dims = []
    for total in range(exponent):
        dims.append(sum(not any(i >= a and total-i >= b for a, b in leading) for i in range(total+1)))
    check('regular_subideal_exact_length_and_top_degree', sum(dims) == start*N and dims[-1] == 1)
    check('window_ideal_bound_does_not_claim_sharpness_or_CI', cert['sharpness_claimed'] is False
          and cert['window_ideal_is_complete_intersection_claimed'] is False)
    if D == 6 and start == 1:
        check('full_window_quotient_strictly_smaller_than_selected_CI', gb.reduce(V**2)[1] == 0
              and ci.reduce(V**2)[1] != 0)

for sizes, exponents, accepted in [((2,), (6,), True), ((1,), (6,), False),
                                  ((2, 2), (2, 3), True), ((2, 1), (1, 6), False),
                                  ((3, 2), (4, 2), True), ((1,), (1,), True)]:
    for extra in (0, 2):
        cert = h.block_sum_sfc_certificate(sizes, exponents, extra)
        check('SFC_support_budget_is_explicit', (cert['status'] == 'expanded_block_sum_sfc_certified') == accepted)
        check('no_SFC_for_entire_phase_powers_or_Gaussian_linear_adjunction', not cert['powers_of_entire_phase_claimed']
              and not cert['Gaussian_linear_adjunction_claimed'])
        if not extra:
            variables = s.symbols('v:'+str(sum(sizes)))
            offset, phase = 0, s.Integer(1)
            for size, d in zip(sizes, exponents):
                phase *= sum(variables[offset:offset+size])**d
                offset += size
            check('actual_expanded_support_count', len(s.Poly(s.expand(1+phase), *variables).terms())
                  == cert['actual_expanded_support_when_a_b_nonzero'])
        examples.append(cert)

invalid = [lambda: h.make_model((0,), (1,)),
           lambda: h.make_model((2,), (0,)),
           lambda: h.make_model((2,), (0.5,)),
           lambda: h.make_model((2,), (s.sqrt(2),)),
           lambda: h.make_model((2,), (1,), law='dirichlet'),
           lambda: h.make_model((2,), (1,), s.exp(1)),
           lambda: h.make_model((2,), (1,), 0.1),
           lambda: h.make_model((2,), (1,), s.Symbol('unbound')),
           lambda: h.make_model((True,), (1,)),
           lambda: h.window_certificate(fixtures[0], 0),
           lambda: h.moment(fixtures[0], 1, 0.1, 1),
           lambda: h.ordinary_ideal_certificate(mass_zero),
           lambda: h.block_sum_sfc_certificate((0,), (2,))]
for operation in invalid:
    try:
        operation()
    except ValueError:
        check('invalid_or_outside_executable_domain_fails_closed', True)
    else:
        raise AssertionError('invalid input accepted')

files = ['HYPERGEOMETRIC_MONOMIAL_WINDOWS.tex', 'hypergeometric_windows.py',
         'verify_hypergeometric_windows.py', 'moment_bridges.py']
receipt = dict(status='pass', checks=sum(counts.values()), groups=counts,
               files={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
               universal_proofs='HYPERGEOMETRIC_MONOMIAL_WINDOWS.tex',
               mathematical_shape_domain='all positive real shapes', executable_shape_domain='positive rational',
               independent_replay=False, external_review=False, formal_verification=False, novelty_claim=False,
               full_two_term_exponent_six_SFC_proved=False, whole_goal_complete=False)
(HERE/'hypergeometric_window_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'hypergeometric_window_examples.json').write_text(json.dumps(examples, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
