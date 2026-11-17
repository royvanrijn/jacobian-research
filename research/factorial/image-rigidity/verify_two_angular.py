"""Exact implementation checks; universal quantifiers use the written proof."""
from pathlib import Path
import hashlib
import itertools
import json
import sympy as s
import moment_bridges as bridge
import two_angular as t

HERE = Path(__file__).resolve().parent
counts = {}


def check(group, result):
    if not result:
        raise AssertionError(group)
    counts[group] = counts.get(group, 0)+1


X, Y, h = t.X, t.Y, t.h
for i in range(9):
    for j in range(6):
        poly = X**i*Y**j
        cert = t.normal_form(poly)
        check('monomial_stokes_reconstruction', t.stokes_residual(poly, cert) == 0)
        # Independent classical even-Y formula and X recurrence without primitives.
        ys = (0 if j % 2 else (-h/2)**(j//2)*s.factorial2(j-1))
        xs = [s.S.One, X]
        for k in range(2, i+1):
            xs.append(s.expand(xs[k-2]-(h*s.Rational(k-2, 3)*xs[k-3] if k > 2 else 0)))
        check('independent_tensor_normal_form', s.expand(cert['a']+cert['b']*X-xs[i]*ys) == 0)

for i in range(4):
    for j in range(4):
        p = X**i*Y**j
        for exact in (h*s.diff(p, X)+s.diff(t.q, X)*p,
                      h*s.diff(p, Y)+s.diff(t.q, Y)*p):
            rem = sum(h**k*(lambda c: c['a']+c['b']*X)(t.normal_form(c))
                      for (k,), c in s.Poly(exact, h).terms())
            check('twisted_images_have_zero_remainder', s.expand(rem) == 0)

for exponents in itertools.product(range(4), repeat=3):
    poly = t.barycentric_amplitude(exponents)
    cert = t.normal_form(poly)
    check('weighted_stokes_reconstruction', t.stokes_residual(poly, cert) == 0)
    for critical in (1, -1):
        value = s.expand(cert['a']+critical*cert['b']).subs(h, 0)
        expected = poly.subs({X: critical, Y: 0})
        check('jacobian_specialization_at_both_criticals', value == expected)
    unit = s.Rational(cert['a'].subs(h, 0)+cert['b'].subs(h, 0))
    check('all_shape_h0_5adic_unit_control', t._valuation5(unit) == 0)
    check('normal_form_is_5_integral', all(t._valuation5(c) in (None, 0)
          or t._valuation5(c) > 0
          for key in ('a', 'b') for c in s.Poly(cert[key], h).all_coeffs()))

q0, q1 = t.normal_form(t.q), t.normal_form(t.q*X)
check('connection_first_row', (q0['a'], q0['b']) == (-5*h/6, -2))
check('connection_second_row', (q1['a'], q1['b']) == (-2, -7*h/6))
z = s.Symbol('z')
A = s.Matrix([[-s.Rational(5, 6)/z, -2], [-2, -s.Rational(7, 6)/z]])
e = s.Matrix([[1, 0]])
de = e*A
d2e = de.diff(z)+de*A
check('cyclic_scalar_equation', s.simplify(d2e+2*de/z+(s.Rational(5, 36)/z**2-4)*e) == s.zeros(1, 2))
check('cyclic_coordinate_is_rank_two', s.det(s.Matrix.vstack(e, de)) == -2)
rho = s.Symbol('rho')
check('riccati_origin_residues', set(s.solve(rho*(rho-1)+s.Rational(5, 36), rho)) == {s.Rational(1, 6), s.Rational(5, 6)})
for epsilon in (2, -2):
    b0 = s.Symbol('b0')
    ansatz = epsilon+b0/z
    residual = s.diff(ansatz, z)+ansatz**2-4+s.Rational(5, 36)/z**2
    check('riccati_infinity_1_over_z_forces_zero', s.expand(residual).coeff(z, -1) == 2*epsilon*b0)

check('strict_interior_minimum_identity', s.expand(t.q+2-(X-1)**2*(X+2)-Y**2) == 0)
check('two_angular_coordinates_are_essential', s.det(s.hessian(t.q, (X, Y))) == 12*X)
uu, vv = s.symbols('uu vv')
map_x, map_y = 2*(1-uu-vv), vv-uu
check('triangle_jacobian_is_four', s.det(s.Matrix([map_x, map_y]).jacobian([uu, vv])) == -4)
check('barycentric_inverse', all(s.expand(a.subs({X: map_x, Y: map_y})-b) == 0
                               for a, b in zip((t.u, t.v, t.w), (uu, vv, 1-uu-vv))))
check('boundary_phase', s.expand(t.q.subs(Y, 1-X/2)) == X**3+X**2/4-4*X+1)

for value in (101*s.I/5, s.sqrt(2)/5, s.sqrt(s.Rational(1, 5)), (s.sqrt(2)+s.I)/5, 1/s.Integer(125)):
    cert = t.integer_shapes_certificate(value)
    check('negative_place_certified_by_minimal_polynomial', cert['status'] == 'nonzero_for_all_integer_block_shapes')
for value in (101*s.I, s.sqrt(5), 1, s.Rational(1, 3)):
    cert = t.integer_shapes_certificate(value)
    check('integral_scale_is_outside_uniform_test_only', cert['status'] == 'outside_uniform_arithmetic_criterion' and cert['period_is_zero'] is None)
for bad in (0.2, s.pi, s.Symbol('xi')):
    try:
        t.integer_shapes_certificate(bad)
    except ValueError:
        check('inexact_or_non_algebraic_arithmetic_input_rejected', True)
    else:
        raise AssertionError('unsupported arithmetic input accepted')
check('zero_scale_separate', t.integer_shapes_certificate(0)['status'] == 'nonzero_constant_phase_limit')

examples = []
for shapes in ([1, 1, 1], [2, 3, 4], [4, 2, 1]):
    cert = t.weighted_certificate(shapes, 101*s.I/5)
    check('actual_specialized_remainder_nonzero', cert['status'] == 'nonzero_transcendental_period')
    examples.append(cert)
boundary = t.amplitude_certificate(1+2*Y**2, 1)
check('exact_amplitude_retains_boundary_no_zero_decision', boundary['status'] == 'zero_remainder_no_period_decision'
      and boundary['period_is_zero'] is None and boundary['boundary_term_retained'])
check('positive_exact_amplitude', s.Poly(1+2*Y**2, Y).all_coeffs() == [2, 0, 1])
examples.append(boundary)

ys = s.symbols('y:3')
S = sum(ys)
G = 8*ys[2]**3-6*ys[2]*S**2+(ys[1]-ys[0])**2*S
for shapes in ([1, 1, 1], [2, 1, 3], [3, 2, 1]):
    a, b, c = (value-1 for value in shapes)
    weight = t.barycentric_amplitude((a, b, c))
    normalization = s.gamma(sum(shapes))/(4*s.prod(s.gamma(value) for value in shapes))
    for m in range(4):
        integrand = s.expand(weight*t.q**m)
        int_y = s.integrate(integrand, (Y, -(1-X/2), 1-X/2))
        triangle = s.integrate(int_y, (X, 0, 2))
        angular = bridge.dirichlet(s.expand(G**m), ys, shapes)
        check('triangle_dirichlet_normalization', s.expand(normalization*triangle-angular) == 0)
        check('homogeneous_gamma_radial_identity', bridge.gamma(s.expand(G**m), ys, shapes)
              == s.rf(sum(shapes), 3*m)*angular)

for shapes, blocks in [([1, 1, 2, 1], [[0], [1], [2, 3]]),
                       ([2, 1, 1, 3], [[0, 1], [2], [3]])]:
    source = s.symbols('b:4')
    sums = [sum(source[j] for j in block) for block in blocks]
    radial = sum(source)
    cubic = 8*sums[2]**3-6*sums[2]*radial**2+(sums[1]-sums[0])**2*radial
    grouped = [sum(shapes[j] for j in block) for block in blocks]
    for m in range(3):
        check('gamma_grouping_preserves_original_variables', bridge.gamma(s.expand(cubic**m), source, shapes)
              == bridge.gamma(s.expand(G**m), ys, grouped))
    cert = t.family_certificate(source, blocks, shapes, 4, 101*s.I/5,
                                source[0]*source[3]+s.I*source[1]+2)
    check('family_constructs_arbitrary_lower_couplings', cert['status'] == 'eventually_every_pure_moment_nonzero')
    examples.append({k: str(value) if isinstance(value, s.Basic) else value for k, value in cert.items()})

lift, xs, blocks = bridge.gaussian_lift(G, ys, [2, 2, 2])
for m in range(3):
    check('three_even_block_gaussian_bridge', bridge.gaussian(s.expand(lift**m), xs)
          == bridge.gamma(s.expand(G**m), ys, [1, 1, 1]))

# Individual half shapes may aggregate to three integer block shapes.
source = s.symbols('c:6')
radial = sum(source)
group_sums = [source[0]+source[1], source[2]+source[3], source[4]+source[5]]
cubic = 8*group_sums[2]**3-6*group_sums[2]*radial**2+(group_sums[1]-group_sums[0])**2*radial
for m in range(3):
    check('half_shapes_aggregate_to_integer_blocks', bridge.gamma(s.expand(cubic**m), source, [s.Rational(1, 2)]*6)
          == bridge.gamma(s.expand(G**m), ys, [1, 1, 1]))

scale = 101*s.I/5
phase = S**4+4*scale*G+ys[0]*ys[1]+s.I
for power in (2, 3):
    remainder = s.expand(phase**power-S**(4*power)-4*power*scale*S**(4*power-4)*G)
    check('power_closure_preserves_scale_and_remainder_bound', s.Poly(remainder, *ys).total_degree() <= 4*power-2)

for bad_args in [([1, 1, s.Rational(1, 2)], 1), ([1, 0, 2], 1), ([1, 2], 1)]:
    try:
        t.weighted_certificate(*bad_args)
    except ValueError:
        check('unsupported_block_shapes_rejected', True)
    else:
        raise AssertionError('unsupported shape vector accepted')
try:
    t.family_certificate(ys, [[0], [1], [2]], [1, 1, 1], 4, scale, ys[0]**3)
except ValueError:
    check('next_layer_cannot_be_hidden_in_lower_terms', True)
else:
    raise AssertionError('degree-three lower perturbation accepted')

receipt = {'date': '2026-11-16', 'status': 'pass', 'checks': counts,
           'total': sum(counts.values()),
           'scope': 'Finite exact Stokes identities, free normal form, rank-two connection/Riccati controls, all-shape h0 units, arithmetic scales, actual boundaries, Gamma/Dirichlet/Gaussian constants and power closure. All-exponent conclusions use the written proof.',
           'files': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                     ('TWO_ANGULAR_COORDINATES.tex', 'two_angular.py', 'verify_two_angular.py', 'moment_bridges.py')},
           'independent_replay': False, 'formal_verification': False, 'external_review': False,
           'novelty_claim': False}
(HERE/'two_angular_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
(HERE/'two_angular_examples.json').write_text(json.dumps(examples, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
