#!/usr/bin/env python3
"""Declared exact controls for QFC1; universal arguments are written proofs."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


from math import comb, factorial
from itertools import product
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h = sp.symbols('h', nonzero=True)
u, v, w, z, eps = sp.symbols('u v w z eps')
variables = (u, v, w)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()




# Taylor coefficient extraction over the Q-algebra Q[h].
def at_zero_by_derivatives(poly, coords):
    poly = sp.Poly(poly, *coords)
    out = 0
    for alpha in product(*(range(poly.degree(t)+1) for t in coords)):
        term = poly.as_expr()
        factor = 1
        for t, exponent in zip(coords, alpha):
            term = sp.diff(term, t, exponent)
            factor *= (-t)**exponent/sp.factorial(exponent)
        out += factor*term
    return sp.expand(out)


test_poly = (h+1)*u**4*v**2+3*h**2*u*v**3+v**5-7*u+2*h+5
assert at_zero_by_derivatives(test_poly, (u, v)) == test_poly.subs({u: 0, v: 0})
coefficient_checks = []
for alpha, coefficient in sp.Poly(test_poly, u, v).terms():
    derivative = sp.diff(test_poly, u, alpha[0], v, alpha[1])
    extracted = at_zero_by_derivatives(derivative, (u, v))
    extracted /= factorial(alpha[0])*factorial(alpha[1])
    assert sp.expand(extracted-coefficient) == 0
    coefficient_checks.append({'exponent': list(alpha), 'coefficient': str(coefficient)})


# A finite presentation with two compatible translation connections.
U = sp.Matrix([[1, u], [0, 1]])
V = sp.Matrix([[1, 0], [v, 1]])
A = U*sp.diag(h, h**2)*V
assert sp.expand(A.det()) == h**3
assert sp.expand(A[0, 0]-u*A[1, 0]) == h
for t in (u, v):
    gamma_target = -U.diff(t)*U.inv()
    gamma_source = V.inv()*V.diff(t)
    assert sp.simplify(A.diff(t)+gamma_target*A-A*gamma_source) == sp.zeros(2)
    identity = (sp.eye(2)+eps*gamma_target)*(A+eps*A.diff(t))-A*(sp.eye(2)+eps*gamma_source)
    for entry in identity:
        assert sp.expand(entry).coeff(eps, 0) == 0
        assert sp.expand(entry).coeff(eps, 1) == 0
    assert sp.diff(A.det(), t) == 0
    for entry in A.diff(t):
        assert sp.rem(sp.Poly(entry, h), sp.Poly(h, h)).as_expr() == 0


# The finite escaping module is h-divisible but cannot have a u-connection.
q = 1-h*u
assert sp.expand(q-u*sp.diff(q, u)) == 1
assert q.subs(u, 1/h) == 0
assert (h*u).subs(u, 1/h) == 1
geometric_prefix = sum(h**j*u**j for j in range(7))
assert sp.expand(q*geometric_prefix-(1-h**7*u**7)) == 0


# Closed polynomial two-form in three variables, with the radial primitive.
potential = [h*(u*v**2+w**2), h*(u**2*w+v), h*(v**2*u+w**3)]
curvature = [
    [sp.diff(potential[j], variables[i])-sp.diff(potential[i], variables[j])
     for j in range(3)]
    for i in range(3)
]
assert sp.expand(
    sp.diff(curvature[1][2], u)+sp.diff(curvature[2][0], v)+sp.diff(curvature[0][1], w)
) == 0
primitive = [
    sp.expand(sum(
        variables[i]*sp.integrate(
            z*curvature[i][j].subs({u: z*u, v: z*v, w: z*w}, simultaneous=True),
            (z, 0, 1)
        )
        for i in range(3)
    ))
    for j in range(3)
]
for i in range(3):
    assert sp.expand(primitive[i].subs(h, 0)) == 0
    for j in range(3):
        assert sp.expand(
            sp.diff(primitive[j], variables[i])-sp.diff(primitive[i], variables[j])
            -curvature[i][j]
        ) == 0


# Independent small PBW implementation: normal order x1,x2,p1,p2, [p_i,x_j]=h delta_ij.
ZERO = {}
ONE = {(0, 0, 0, 0): sp.Integer(1)}


def clean(op):
    return {key: value for key, old in op.items() if (value := sp.expand(old)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, 0)+value
    return clean(out)


def scale(op, value):
    return clean({key: value*coefficient for key, coefficient in op.items()})


def mul(left, right):
    out = {}
    for (ix, iy, ip, iq), c in left.items():
        for (jx, jy, jp, jq), e in right.items():
            for k in range(min(ip, jx)+1):
                for ell in range(min(iq, jy)+1):
                    key = (ix+jx-k, iy+jy-ell, ip+jp-k, iq+jq-ell)
                    coefficient = (
                        c*e*comb(ip, k)*(factorial(jx)//factorial(jx-k))
                        *comb(iq, ell)*(factorial(jy)//factorial(jy-ell))*h**(k+ell)
                    )
                    out[key] = out.get(key, 0)+coefficient
    return clean(out)


def comm(left, right):
    return add(mul(left, right), scale(mul(right, left), -1))


def power(op, exponent):
    out = ONE
    for _ in range(exponent):
        out = mul(out, op)
    return out


def degree(op):
    return max((sum(key) for key in op), default=-1)


def classical(op):
    return clean({key: coefficient.subs(h, 0) for key, coefficient in op.items()})


x1 = {(1, 0, 0, 0): sp.Integer(1)}
x2 = {(0, 1, 0, 0): sp.Integer(1)}
p1 = {(0, 0, 1, 0): sp.Integer(1)}
p2 = {(0, 0, 0, 1): sp.Integer(1)}
positions = [add(x1, scale(power(p1, 2), h**2)), add(x2, scale(power(p2, 2), h**2))]
extra = scale(mul(power(positions[0], 2), positions[1]), h**2)
mates = [scale(p1, -1), add(scale(p2, -1), extra)]
assert classical(positions[0]) == x1
assert classical(positions[1]) == x2
assert classical(mates[0]) == scale(p1, -1)
assert classical(mates[1]) == scale(p2, -1)
assert comm(*positions) == ZERO
for i in range(2):
    for j in range(2):
        assert comm(positions[j], mates[i]) == (scale(ONE, h) if i == j else ZERO)
C = comm(*mates)
assert C == scale(mul(*positions), -2*h**3)
b2 = scale(extra, -1)
flat = [mates[0], add(mates[1], b2)]
assert flat[1] == scale(p2, -1)
assert comm(*flat) == ZERO
for i in range(2):
    for j in range(2):
        assert comm(positions[j], flat[i]) == (scale(ONE, h) if i == j else ZERO)
d1, d2 = map(degree, positions)
e1, e2 = map(degree, mates)
bound = e1+e2-2+d1
assert degree(C) <= e1+e2-2
assert degree(b2) <= degree(C)+d1 <= bound
assert degree(flat[1]) <= max(e2, bound)
assert degree(C) == 4 and degree(b2) == 6 and degree(flat[1]) == 1


# The complete Poisson-centralizer and module-finiteness arguments are written.
# These controls replay the parity basis and derivation identities on declared polynomials.
a, b, c, d = sp.symbols('a b c d')
sample = 3*a**7*b**4-a**2*b**5+11*a*b+4
parity_parts = {(i, j): 0 for i in range(2) for j in range(2)}
for (i, j), coefficient in sp.Poly(sample, a, b).terms():
    parity_parts[(i % 2, j % 2)] += coefficient*c**(i//2)*d**(j//2)
rebuilt = sum(
    a**i*b**j*part.subs({c: a*a, d: b*b})
    for (i, j), part in parity_parts.items()
)
assert sp.expand(rebuilt-sample) == 0
symbol_sample = u**3*a**2+u*v*b+a*b**2+v**4
for t, p in [(u, a), (v, b)]:
    bracket_with_square = sp.diff(symbol_sample, p)*sp.diff(p*p, t)-sp.diff(symbol_sample, t)*sp.diff(p*p, p)
    assert sp.expand(bracket_with_square+2*p*sp.diff(symbol_sample, t)) == 0


# A finite polynomial Keller shear: actual inverse and all Weyl relations.
k_positions = [add(x1, power(x2, 3)), x2]
k_mates = [scale(p1, -1), add(scale(p2, -1), scale(mul(power(x2, 2), p1), 3))]
assert add(k_positions[0], scale(power(k_positions[1], 3), -1)) == x1
assert comm(*k_positions) == ZERO and comm(*k_mates) == ZERO
for i in range(2):
    for j in range(2):
        assert comm(k_positions[j], k_mates[i]) == (scale(ONE, h) if i == j else ZERO)


# Birationality alone does not give module finiteness: (x,xy) has an infinite fiber.
assert sp.Matrix([[sp.diff(u, u), sp.diff(u, v)],
                  [sp.diff(u*v, u), sp.diff(u*v, v)]]).det() == u
assert sp.cancel((u*v)/u-v) == 0
assert (u*v).subs(u, 0) == 0


hash_names = ['FINITE_CENTRALIZERS_AND_COMPLETION.md', 'check_finite_centralizer_completion.py', 'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md', 'BOUNDED_REALIZATION.md']
receipt = {
    'claim': 'QFC1', 'status': 'PASS',
    'Taylor_projector_and_coefficients': coefficient_checks,
    'dual_number_Fitting_control': {
        'presentation': str(A), 'Fitt0': 'h^3', 'Fitt1': '(h)',
        'two_compatible_connections': 'PASS', 'semilinear_equivalence_mod_epsilon2': 'PASS',
    },
    'escaping_support_control': {
        'module': 'R0[u]/(1-h*u)', 'nonzero_and_h_divisible': True,
        'annihilator_derivative_unit_identity': 'q-u*qprime=1',
        'coordinate_connection_possible': False,
        'formal_inverse_prefix_length': 7,
    },
    'three_variable_closed_curvature_control': {
        'closed': True, 'finite_radial_primitive': [str(x) for x in primitive],
        'preserves_classical_tuple': True,
    },
    'nonpositional_finite_Weyl_control': {
        'classical_graph': 'identity', 'Long_graph': False,
        'positions': ['x1+h^2*p1^2', 'x2+h^2*p2^2'],
        'initial_curvature': '-2*h^3*F1*F2', 'correction': '-h^2*F1^2*F2',
        'all_final_Weyl_relations': 'PASS', 'leading_symbol_module_basis_rank': 4,
        'position_degrees': [d1, d2], 'initial_mate_degrees': [e1, e2],
        'curvature_degree': degree(C), 'correction_degree': degree(b2),
        'proved_correction_bound_in_control': bound, 'corrected_second_mate_degree': degree(flat[1]),
    },
    'parity_basis_control': {'basis': ['1', 'p1', 'p2', 'p1*p2'], 'sample_reconstruction': 'PASS'},
    'Keller_shear_control': {'full_relations_and_polynomial_inverse': 'PASS', 'noninvertible': False},
    'birational_nonfinite_control': {
        'map': '(u,u*v)', 'fraction_field_degree': 1,
        'fiber_over_zero': 'K[v]', 'Jacobian': 'u', 'Keller': False,
    },
    'external_dependencies': [
        {'source': 'https://stacks.math.columbia.edu/tag/07ZA',
         'used': 'Fitting base change, presentation invariance and annihilator facts; differential ideal proof is given directly.'},
        {'source': 'https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27',
         'used': 'Remark3.3(c), Theorem3.4 and the fundamental group exact sequence: finite etale affine-space covers in characteristic zero are constant.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Declared finite Taylor, presentation, escaping support, radial homotopy, '
        'PBW, parity and birational controls only. Universal centralizer exhaustion, '
        'Fitting stability, integral and generic finiteness equivalence, all-rank '
        'flat completion, full filtered certificate and Keller/JC boundary use '
        'written proofs with the stated external covering theorem. No actual '
        'Long candidate satisfying the finiteness and individual-mate hypotheses '
        'is constructed or certified. No finite sweep proves module finiteness.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(_check_output(HERE/'finite_centralizer_controls.json')).write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QFC1', 'status': 'PASS', 'all_final_Weyl_relations': True,
    'correction_degree': degree(b2), 'corrected_second_mate_degree': degree(flat[1]),
    'radial_primitive_in_three_variables': True, 'Long_candidate_constructed': False,
}, indent=2))
