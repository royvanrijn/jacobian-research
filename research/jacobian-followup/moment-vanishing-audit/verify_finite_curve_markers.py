"""Selected exact controls, not a verifier of the universal written proofs.

Run from the repository root with .venv/bin/python. No moment-prefix,
degree/support search, numerical Riemann map, or exceptional purity inference.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
t, u, z, v = s.symbols('t u z v', real=True)
p = t**2-s.Rational(1, 2)+s.I*(t**3-t)
x, y = s.re(p).expand(), s.im(p).expand()
assert s.degree(p, t) == 3 and p.subs(t, -1) == p.subs(t, 1) == s.Rational(1, 2)
rx = s.cancel((x-x.subs(t, v))/(t-v))
ry = s.cancel((y-y.subs(t, v))/(t-v))
assert rx == t+v and ry == t**2+t*v+v**2-1
assert s.expand(ry.subs(v, -t)) == t**2-1
assert s.gcd(s.diff(x, t), s.diff(y, t)) == 1
assert s.expand(y**2-t**2*(1-t**2)**2) == 0
assert x.subs(t, 0) == -s.Rational(1, 2)
area = s.integrate(2*s.sqrt(u)*(1-u), (u, 0, 1))
signed_area = s.integrate((x*s.diff(y, t)-y*s.diff(x, t))/2, (t, -1, 1))
assert area == s.Rational(8, 15) and abs(signed_area) == area
assert 0 < s.Rational(1, 2) < 1  # the origin satisfies the strict domain inequalities

# Two selected nonmonic module-reduction controls. Written proof supplies all d.
controls = []
for phase, q in [
    (p, (2+s.I)*t**14-3*t**11+t**7+s.I*t-5),
    (3*t**2+2*t-7, t**11+(2-s.I)*t**6-4*t+1),
]:
    degree = s.degree(phase, t)
    rem = s.rem(q, phase-z, t, domain=s.QQ_I.poly_ring(z))
    assert s.degree(rem, t) < degree
    assert s.expand(rem.subs(z, phase)-q) == 0
    coeffs = [s.expand(rem).coeff(t, j) for j in range(degree)]
    assert all(s.denom(s.cancel(a)).free_symbols == set() for a in coeffs)
    controls.append({'phase_degree': int(degree), 'input_degree': int(s.degree(q, t)),
                     'basis_size': int(degree), 'reconstruction_exact': True,
                     'coefficients_polynomial_in_phase': True})

# The nilpotent conclusion, including a nonreduced algebra, is intentional.
epsilon = s.Symbol('epsilon')
basis = [1, epsilon, epsilon**2]
def frobenius(a):
    return s.rem(a, epsilon**3, epsilon).expand().coeff(epsilon, 2)
gram = s.Matrix([[frobenius(a*b) for b in basis] for a in basis])
assert gram.det() == -1
assert all(frobenius(a*epsilon**3) == 0 for a in basis)

# Finite reduced algebra: coefficient-of-t^2 pairing, t^3=2.
def finite_coefficient(a):
    return s.rem(a, t**3-2, t).expand().coeff(t, 2)
finite_gram = s.Matrix([[finite_coefficient(t**(a+b)) for b in range(3)] for a in range(3)])
assert finite_gram.det() == -1
assert finite_coefficient(t**2) == 1
assert s.rem(t**3-2, t**3-2, t) == 0  # exact tail recurrence a_(m+3)=2*a_m

# Boundary control: Laurent CT is nondegenerate, but the algebra is not finite
# over C[z]. A selected Laurent marker has an eventually zero z-tail.
marker = { -5: s.Rational(2), -1: s.Rational(-3), 0: s.Rational(7), 4: s.Rational(1) }
tail = [str(marker.get(-m, 0)) for m in range(8)]
assert tail == ['7', '-3', '0', '0', '0', '2', '0', '0']

receipt = {
    'date': '2026-11-15', 'completed': True,
    'checker_sha256': sha(Path(__file__)),
    'python_version': platform.python_version(), 'sympy_version': s.__version__,
    'loop_polynomial': str(p), 'loop_degree': 3,
    'only_distinct_collision_parameters_on_interval': ['-1,1'],
    'real_derivative_components_coprime': True,
    'loop_area': str(area), 'loop_signed_area': str(signed_area),
    'origin_strictly_inside_loop': True,
    'selected_free_module_reductions': controls,
    'nilpotent_frobenius_pairing_determinant': str(gram.det()),
    'finite_reduced_pairing_determinant': str(finite_gram.det()),
    'finite_reduced_nonzero_tail': 'a_(3k+2)=2^k by t^3=2',
    'laurent_boundary_marker_tail_control': tail,
    'laurent_finiteness_hypothesis': False,
    'all_degree_marker_theorem': 'INTERNAL WRITTEN:finite module+nondegenerate pairing implies safety iff nilpotence',
    'full_exceptional_thin_pure_phase_consequence': 'CONDITIONAL MZ FAILURE:some fixed h^j,1<=j<deg p,has infinitely many nonzero marked moments',
    'abstract_loop_all_order_purity': 'INTERNAL WRITTEN:Riemann map+Jordan boundary extension+mean-value property',
    'abstract_loop_fixed_marker_list': ['t', 't^2'],
    'abstract_loop_which_marker_identified': False,
    'abstract_loop_is_exceptional_Haar_counterexample': False,
    'abstract_loop_Haar_realization_proved': False,
    'universal_algebra_analytic_and_Haar_proofs_machine_verified': False,
    'new_support_degree_or_moment_prefix_campaign': False,
    'exceptional_MZ_decisions': 'OPEN on both models',
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'finite_curve_marker_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
