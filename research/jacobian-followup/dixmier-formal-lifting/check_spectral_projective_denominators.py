#!/usr/bin/env python3
"""Declared QPD1 controls; universal Morita clearing uses written dimension proofs."""
from pathlib import Path
import hashlib
import json
import sympy as sp
import check_hamiltonian_idealizer as hi

HERE = Path(__file__).resolve().parent
u, v = hi.u, hi.v
xi = sp.symbols('xi')
mul, comm, zero = hi.mul, hi.comm, hi.zero


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['hamiltonian_idealizer_controls.json', 'projective_frame_controls.json',
             'taylor_centralizer_controls.json', 'rational_spectral_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, digest in prior['hashes'].items():
        assert sha(HERE/rel) == digest, rel


m1, m2 = u, v+1/u
H0, L0 = u*v, v*v-2/u**2
multipliers = [
    (u**2, (u**2, 0), (mul(v, u)+1, 0)),
    (u**3, (u**3, 0), (mul(v, u**2)+u, 0)),
    (H0, (-1, u**2), (v**2, 0)),
    (L0, (-v**2, 2*u*v), (0, v**2)),
]
memberships = 0
for r, coefficients1, coefficients2 in multipliers:
    for m, coefficients in [(m1, coefficients1), (m2, coefficients2)]:
        zero(mul(m, r)-mul(coefficients[0], m1)-mul(coefficients[1], m2))
        memberships += 1


F = u*u
theta = u
T = -v/(2*u)
actual_mate_numerator = -H0/2
zero(comm(F, theta))
zero(comm(F, T)-1)
zero(mul(F, theta)-u**3)
zero(mul(F, T)-actual_mate_numerator)
zero(comm(F, actual_mate_numerator)-F)
zero(mul(u**3, u**3)-mul(F, mul(F, F)))
assert hi.degree_symbol(theta)[0] == 1
assert hi.degree_symbol(T)[0] == 0
assert sp.diff(hi.degree_symbol(T)[1], v) != 0

T_square = mul(T, T)
square_numerator = (mul(H0, H0)-2*H0)/4
zero(comm(F, T_square)-2*T)
zero(comm(F, comm(F, T_square))-2)
zero(comm(F, comm(F, comm(F, T_square))))
zero(mul(mul(F, F), T_square)-square_numerator)
zero(mul(mul(F, F), 1/u)-u**3)

# Scalar spectral shift preserves the exact linear clearing witness.
a = hi.a
shifted_F = F+a
zero(comm(shifted_F, actual_mate_numerator)-(shifted_F-a))
cubic_F, cubic_T = u**3, -v/(3*u**2)
zero(comm(cubic_F, cubic_T)-1)
zero(mul(cubic_F, cubic_T)+H0/3)
zero(comm(cubic_F, -H0/3)-cubic_F)


def conjugate(q):
    return mul(mul(1/u, q), u)


# The fractional frame M->M*u gives the original nonfree projective ideal.
new_m1, new_m2 = mul(m1, u), mul(m2, u)
zero(new_m1-u*u)
zero(new_m2-(u*v+2))
for r, coefficients1, coefficients2 in multipliers:
    transformed = conjugate(r)
    for m, coefficients in [(new_m1, coefficients1), (new_m2, coefficients2)]:
        zero(mul(m, transformed)-mul(coefficients[0], new_m1)-mul(coefficients[1], new_m2))
        memberships += 1
zero(conjugate(H0)-(H0+1))
zero(conjugate(F)-F)
zero(comm(conjugate(F), conjugate(T))-1)
zero(mul(conjugate(F), conjugate(T))-conjugate(actual_mate_numerator))
zero(conjugate(T)-(-v/(2*u)-1/(2*u*u)))

# Declared semigroup controls support the written all-degree cusp proof.
cusp_controls = []
for degree in [2, 3, 4, 5, 7, 8]:
    odd = degree % 2
    even_power = (degree-3*odd)//2
    assert even_power >= 0
    zero((u*u)**even_power*(u**3)**odd-u**degree)
    cusp_controls.append(degree)

# Division-ring degree controls are not a finite proof about all algebraic elements.
degree_controls = []
for q in [1/u, v/u, v*v-2/u**2, sp.Integer(3)]:
    degree, symbol = hi.degree_symbol(q)
    degree_controls.append({'operator': str(q), 'degree': degree, 'symbol': str(symbol)})
assert hi.degree_symbol(1/u)[0] == -1
assert hi.degree_symbol(v/u)[0] == 0
assert sp.diff(hi.degree_symbol(v/u)[1], v) != 0

# Rank two has nonconstant units, outside the rank-one unit theorem.
matrix_unit = sp.Matrix([[1, u], [0, 1]])
matrix_inverse = sp.Matrix([[1, -u], [0, 1]])
assert matrix_unit*matrix_inverse == sp.eye(2)
assert matrix_inverse*matrix_unit == sp.eye(2)

hash_names = [
    'SPECTRAL_DENOMINATORS_AND_PROJECTIVE_CENTRALIZERS.md',
    'check_spectral_projective_denominators.py',
    'HAMILTONIAN_IDEALIZER_CRITERION.md', 'check_hamiltonian_idealizer.py',
    'hamiltonian_idealizer_controls.json', 'PROJECTIVE_FRAME_BOUNDARY.md',
    'projective_frame_controls.json', 'TAYLOR_BRANCHES_AND_CENTRALIZER_RIGIDITY.md',
    'taylor_centralizer_controls.json', 'ETALE_CENTRALIZERS_AND_GALOIS_COMPLETION.md',
    'FINITE_CENTRALIZERS_AND_COMPLETION.md',
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
    'RATIONAL_SPECTRAL_OBSTRUCTIONS.md', 'rational_spectral_controls.json',
    'BOUNDED_REALIZATION.md',
]
receipt = {
    'claim': 'QPD1', 'status': 'PASS',
    'charge_one_projective_control': {
        'module_generators': ['u', 'v+u^-1'], 'nonfree': 'retained QPF1 written theorem',
        'right_multiplier_memberships_in_two_fractional_frames': memberships,
        'position': 'u^2', 'rational_commutant': 'u', 'commutant_spectral_denominator': 'F',
        'cleared_commutant': 'u^3', 'rational_mate': '-v/(2*u)',
        'mate_spectral_denominator': 'F', 'cleared_mate': '-u*v/2',
        'spectral_class_nonzero': 'retained QHI1 all-partner exclusion and QPD1 field argument',
        'exact_spectral_annihilator': '(xi), written ideal proof',
        'centralizer': 'L[u^2,u^3], written all-degree proof',
        'rational_centralizer': 'L(u)', 'fraction_field_degree_over_positions': 2,
        'actual_radical_certificate': '(u^3)^2=F*(u^2)^2',
    },
    'finite_adjoint_orbit_control': {
        'operator': 'T^2', 'successive_commutators': ['2*T', '2', '0'],
        'spectral_denominator': 'F^2', 'cleared_numerator': '(H0^2-2*H0)/4',
        'all_identities': 'PASS',
    },
    'shifted_and_cubic_scalar_controls': {
        'shifted_position': 'u^2+a', 'spectral_denominator': 'xi-a',
        'cubic_position': 'u^3', 'cubic_denominator': 'xi', 'identities': 'PASS',
    },
    'fractional_frame_transport': {
        'g': 'u', 'new_module_generators': ['u^2', 'u*v+2'],
        'endomorphisms': 'q -> u^-1*q*u', 'spectral_clearing_preserved': True,
        'scalar_free_frame_constructed': False,
    },
    'cusp_declared_semigroup_degrees': cusp_controls,
    'division_symbol_controls': degree_controls,
    'higher_rank_unit_boundary': {'matrix': str(matrix_unit), 'nonconstant_unit': True},
    'source_boundary_correction': {
        'claim': 'QST1', 'change': 'Later spectral clearing supplies the rational-mate/torsion converse; prior residue proofs and scope remain unchanged.',
    },
    'external_dependencies': [
        {'source': 'https://arxiv.org/html/math/0401263',
         'used': 'Theorems1.2,1.4,3.3 give module dimension bounds; Lemma3.7 is the scalar spectral-clearing model, extended here through rank-one Morita equivalence.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Declared charge-one rational Ore identities, actual multiplier memberships, '
        'a finite adjoint orbit, exact spectral numerators, cusp semigroup samples, '
        'fractional transport, scalar shift/cubic examples, division symbols and '
        'a matrix-unit boundary only. Scalar idealizer units reuse QHI1. Universal '
        'Morita denominator clearing, rational centralizer saturation, torsion iff '
        'rational individual admission, coefficient-curve application, algebraicity '
        'and horizontal constants, and projective birational/Galois/radical transfer '
        'use written proofs with explicit external dimension results. No denominator '
        'degree bound, actual full Long carrier or all-representative nonlifting '
        'is certified. Torsion is not actual vanishing or full flat completion.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'spectral_projective_denominator_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QPD1', 'status': 'PASS', 'projective_multiplier_memberships': memberships,
    'rational_mate_cleared_by_F': True, 'T_square_cleared_by_F_square': True,
    'actual_class_zero': False, 'Long_carrier_constructed': False,
}, indent=2))
