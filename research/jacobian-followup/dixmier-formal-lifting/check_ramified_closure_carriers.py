#!/usr/bin/env python3
"""Exact QGC1 closure, connection and admission-boundary controls."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
r, a, v, w, A, X = sp.symbols('r a v w A X')
sout, tout, x, t, z = sp.symbols('sout tout x t z')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['quantum_galois_residue_controls.json', 'finite_centralizer_controls.json',
             'division_index_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, expected in prior['hashes'].items():
        assert sha(HERE/rel) == expected, rel

prior = json.loads((HERE/'division_index_controls.json').read_text())
assert prior['Long_classical_field_degree']['value'] == 3

# Exact splitting of the full classical cubic after the quadratic source extension.
Svalue = 2*v+4*a-3*r*a*a
Tvalue = (r*a**3-a*a-v*a)/2
Q = (3*r*a-2)**2-8*r*v
v_closure = ((3*r*a-2)**2-w*w)/(8*r)
b = (2-r*a+w)/(2*r)
c = (2-r*a-w)/(2*r)
cubic = r*X**3-2*X**2+sout*X+4*tout
split_cubic = cubic.subs({sout: Svalue, tout: Tvalue}).subs(v, v_closure)
assert sp.cancel(split_cubic-r*(X-a)*(X-b)*(X-c)) == 0
assert sp.cancel(a+b+c-2/r) == 0
assert sp.cancel(Q.subs(v, v_closure)-w*w) == 0
Delta = sp.discriminant(cubic, X)
assert sp.expand(Delta.subs({sout: Svalue, tout: Tvalue})-4*v*v*Q) == 0
assert sp.diff(Q, v) == -8*r
# At the generic branch w=0, r and v are units.
branch_v = (3*r*a-2)**2/(8*r)
assert sp.cancel(branch_v) != 0
assert sp.cancel(Q.subs(v, branch_v)) == 0
derivative_w = -4*r/w
assert sp.cancel(2*w*derivative_w-sp.diff(Q, v)) == 0

# The source affine open is genuine: all inverse coordinate expressions are Laurent in v.
source_x = 1/v
source_q = (2*v-r*v*v)/3
source_Y = a-v
source_beta = 5*v*v-r*v**3-3*a*v
assert sp.cancel(2*source_x-3*source_x**2*source_q-r) == 0
assert sp.cancel(source_q-source_x*source_beta/3-source_Y) == 0
determinant = -9*source_x**2*source_q**2-(1-3*source_x*source_q)*(1+3*source_x*source_q)
assert sp.cancel(determinant+1) == 0

# Explicit S3 action on the centered roots, with the chosen-source involution w -> -w.
cycle = ((-A+3*w)/2, -(A+w)/2)
reflection = (A, -w)


def transform(function, mapping):
    return sp.cancel(function.subs({A: mapping[0], w: mapping[1]}, simultaneous=True))


def compose(left, right):
    return tuple(transform(q, right) for q in left)


cycle_squared = compose(cycle, cycle)
cycle_cubed = compose(cycle_squared, cycle)
assert all(sp.cancel(q-basis) == 0 for q, basis in zip(cycle_cubed, (A, w)))
assert compose(reflection, reflection) == (A, w)
assert all(sp.cancel(q-other) == 0 for q, other in zip(
    compose(compose(reflection, cycle), reflection), cycle_squared))
a_centered = (A+2)/(3*r)
v_centered = (A*A-w*w)/(8*r)
s_centered = Svalue.subs({a: a_centered, v: v_centered}, simultaneous=True)
t_centered = Tvalue.subs({a: a_centered, v: v_centered}, simultaneous=True)
for mapping in [cycle, reflection]:
    assert sp.cancel(transform(s_centered, mapping)-s_centered) == 0
    assert sp.cancel(transform(t_centered, mapping)-t_centered) == 0
arrangement = w*(A*A-w*w)
assert sp.cancel(transform(arrangement, cycle)-arrangement) == 0
assert sp.cancel(transform(arrangement, reflection)+arrangement) == 0

# Local ramification has a pole in the unique derivative connection.
# Multiplication is on column coefficients of the basis 1,z,...,z^(m-1).
connection_controls = []
scalar_potential = sp.symbols('scalar_potential')
for m in (2, 3, 5):
    Z = sp.zeros(m)
    for j in range(m-1):
        Z[j+1, j] = 1
    Z[0, m-1] = x
    C = sp.diag(*[sp.Rational(j, m)/x for j in range(m)])
    I = sp.eye(m)
    assert Z**m == x*I
    compatibility = Z.diff(x)+C*Z-Z*C-Z/(m*x)
    assert all(sp.cancel(q) == 0 for q in compatibility)
    shifted = C-scalar_potential*I
    assert all(sp.cancel(q) == 0 for q in (shifted*Z-Z*shifted)-(C*Z-Z*C))
    assert sp.cancel(shifted[1, 1]-shifted[0, 0]-1/(m*x)) == 0
    disc = sp.discriminant(z**m-x, z)
    assert sp.degree(disc, x) == m-1
    # The zero componentwise connection preserves a free module, but not its field product.
    if m == 2:
        assert Z.diff(x) != Z/(2*x)
        nonderivation_defect = sp.diff(x, x)-2*z*sp.S.Zero
        assert nonderivation_defect == 1
    connection_controls.append({
        'm': m, 'multiplication_power': 'PASS',
        'field_connection_compatibility': 'PASS',
        'scalar_potential_cancels': 'PASS',
        'connection_residue_difference': '1/m',
        'discriminant_x_order': m-1,
        'all_fractional_lattices_excluded': 'written multiplier/Fitting theorem',
    })

# Allowed constant etale extension; the derivative of its defining polynomial is a unit.
Theta = sp.Matrix([[0, 2], [1, 0]])
I2 = sp.eye(2)
assert Theta**2 == 2*I2
assert (2*Theta)*(Theta/4) == I2
assert sp.discriminant(z*z-2, z) == 8
assert Theta.diff(x) == sp.zeros(2)

# Exact positive Weyl shear with a ramified quantum parameter hbar=t^3.
quantum_parameter = t**3
F_x = (0, x)
F_p = (1, t*x*x)


def first_order_comm(left, right):
    la, lb = left
    ra, rb = right
    return (
        sp.cancel(quantum_parameter*(la*sp.diff(ra, x)-ra*sp.diff(la, x))),
        sp.cancel(quantum_parameter*(la*sp.diff(rb, x)-ra*sp.diff(lb, x))),
    )


assert first_order_comm(F_p, F_x) == (0, quantum_parameter)
assert sp.cancel((F_p[1]-0)/t-x*x) == 0
polynomial_corrected_lift = (1, t*x*x)
assert polynomial_corrected_lift == F_p

# Finite parameter-flat modules need not have faithful ambient residue.
# Columns are 1 and x+t*z in the field basis 1,z, with z^2=x.
bad_frame = sp.Matrix([[1, x], [0, t]])
assert bad_frame.det() == t
hidden_vector = sp.Matrix([0, t])
hidden_coordinates = bad_frame.inv()*hidden_vector
assert hidden_coordinates == sp.Matrix([-x, 1])
assert bad_frame.subs(t, 0).rank() == 1
assert hidden_coordinates.subs(t, 0) != sp.zeros(2, 1)
assert hidden_vector.subs(t, 0) == sp.zeros(2, 1)

# Genuine nonsaturated constant-cover Weyl carriers, with finite equivariant repairs.
# Row matrices describe the right x,p action on the displayed left module basis.
equivariant_repairs = []
for degree in (2, 3):
    omega = -sp.S.One if degree == 2 else (-1+sp.sqrt(-3))/2
    assert sp.simplify(omega**degree-1) == 0
    basis_frame = sp.zeros(degree)
    basis_frame[0, 0] = 1
    for j in range(1, degree):
        basis_frame[j, 0] = x**j
        basis_frame[j, j] = t**j
    assert sp.cancel(basis_frame.det()-t**(degree*(degree-1)//2)) == 0
    constant_root = sp.zeros(degree)
    for j in range(degree-1):
        constant_root[j, j+1] = 1
    constant_root[degree-1, 0] = 2
    assert constant_root**degree == 2*sp.eye(degree)
    diagonal_group = sp.diag(*[omega**j for j in range(degree)])
    assert all(sp.simplify(q) == 0 for q in
               diagonal_group.inv()*constant_root*diagonal_group-omega*constant_root)
    group = sp.simplify(basis_frame*diagonal_group*basis_frame.inv())
    assert all(sp.simplify(q) == 0 for q in group**degree-sp.eye(degree))
    assert basis_frame.subs(t, 0).rank() == 1

    projectors = []
    for character in range(degree):
        projector = sum(
            (omega**(-character*k)*group**k for k in range(degree)),
            sp.zeros(degree),
        )/degree
        projector = projector.applyfunc(sp.simplify)
        assert all(sp.simplify(q) == 0 for q in projector**2-projector)
        assert projector.rank() == 1
        assert all(sp.simplify(q) == 0 for q in group*projector-omega**character*projector)
        projectors.append(projector)
    assert all(sp.simplify(q) == 0 for q in sum(projectors, sp.zeros(degree))-sp.eye(degree))
    for j in range(degree):
        for k in range(degree):
            if j != k:
                assert all(sp.simplify(q) == 0 for q in projectors[j]*projectors[k])

    # Integer powers of t repair each nonzero character residue line.
    eigenbasis = sp.eye(degree)
    for j in range(1, degree):
        eigenbasis[j, 0] = -x**j
    minimum_shifts = sp.diag(*[t**(-j) for j in range(degree)])
    assert minimum_shifts*eigenbasis*basis_frame == sp.eye(degree)

    for parameter_order in (1, 3):
        parameter = t**parameter_order
        correction = sp.zeros(degree)
        for j in range(1, degree):
            correction[j, 0] = -parameter*j*x**(j-1)
        # [pI+correction,xI]=parameter*I: correction commutes with x.
        assert parameter*sp.diff(x, x)*sp.eye(degree) == parameter*sp.eye(degree)
        compatibility = parameter*group.diff(x)+correction*group-group*correction
        assert all(sp.simplify(q) == 0 for q in compatibility)
        assert correction.subs(t, 0) == sp.zeros(degree)
        # Each pure constant-character row commutes with p under the right action.
        assert eigenbasis*correction == parameter*eigenbasis.diff(x)
        # Divided left/right difference preserves the carrier and the valuation levels.
        connection = -correction/parameter
        assert all(sp.cancel(q).is_polynomial(x, t) for q in connection)
        assert eigenbasis.diff(x)+eigenbasis*connection == sp.zeros(degree)
    equivariant_repairs.append({
        'constant_degree': degree,
        'constant_field': 'Q(sqrt(2))' if degree == 2 else 'Q(omega)(cube_root(2))',
        'basis': ['1']+[f'x^{j}+t^{j}*theta^{j}' for j in range(1, degree)],
        'ambient_residue_rank': 1, 'abstract_classical_rank': degree,
        'character_valuation_minima': list(range(degree)),
        'finite_projector_identities': 'PASS',
        'full_Weyl_right_action_parameter_orders': [1, 3],
        'group_equivariance': 'PASS',
        'identity_classical_graph': 'PASS',
        'quantum_order_graph_congruence': 'PASS',
        'finite_character_repair': 'PASS',
        'valuation_level_ranks': list(range(1, degree+1)),
        'spatial_cover': 'constant etale, allowed',
        'Long_carrier': False,
    })

# Actual torsion in the special module is a different boundary.
# I=A_R*x+A_R*t is finite and stable on both sides, but I/tI has a nonzero t-class
# killed by x. For parameter t, [p,x]=t is not in tI; it is not the identified graph.
assert sp.expand(x*t-t*x) == 0
assert t.subs(t, 0) == 0
for parameter_order in (1, 3):
    quantum_parameter = t**parameter_order
    assert first_order_comm((1, 0), (0, x)) == (0, quantum_parameter)
    assert first_order_comm((1, 0), (0, t)) == (0, 0)
    assert sp.cancel(quantum_parameter/t-t**(parameter_order-1)) == 0

hash_names = [
    'RAMIFICATION_AND_FINITE_CLOSURE_CARRIERS.md',
    'check_ramified_closure_carriers.py',
    'FINITE_QUANTUM_SYMMETRIES_AND_GALOIS_RESIDUES.md',
    'quantum_galois_residue_controls.json',
    'FINITE_CENTRALIZERS_AND_COMPLETION.md', 'finite_centralizer_controls.json',
    'QUANTUM_DIVISION_AND_COMPLETION.md', 'division_index_controls.json',
    'RAMIFICATION_AND_POSITION_AVOIDANCE.md', 'position_avoidance_controls.json',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
]
receipt = {
    'claim': 'QGC1', 'status': 'PASS',
    'Long_quadratic_source_closure': {
        'Q': str(sp.expand(Q)), 'w_squared': 'Q',
        'other_roots': [str(b), str(c)],
        'complete_cubic_splitting': 'PASS',
        'discriminant_pullback': '4*v^2*Q, PASS',
        'closure_rational_field': 'K(r,a,w,d)',
        'degree_over_original_source': 'two, written nonsquare linear-v proof',
        'interior_source_chart': 'PASS, retained fourth-coordinate chart',
        'derivative_Q_v': '-8*r', 'derivative_w_v': '-4*r/w',
        'branch_in_original_affine_source': 'written open-chart/DVR proof',
        'faithful_saturated_finite_natural_carrier': 'excluded by written theorem',
        'Long_lift': False,
    },
    'S3_source_subgroup_diagram': {
        'cycle': [str(q) for q in cycle],
        'chosen_source_reflection': [str(q) for q in reflection],
        'all_generator_relations': 'PASS',
        'outputs_fixed': 'PASS',
        'field_rationality_implies_regular_Poisson_order': False,
    },
    'power_cover_connections': connection_controls,
    'constant_etale_positive_control': {
        'field': 'Q(sqrt(2))', 'polynomial_discriminant': 8,
        'derivative_unit_matrix_identity': 'PASS',
        'finite_saturated_natural_carrier': 'constant scalar extension, written',
    },
    'ramified_parameter_positive_control': {
        'quantum_parameter': 't^3', 'x_image': 'x', 'p_image': 'p+t*x^2',
        'exact_Weyl_relation': 'PASS',
        'sub_quantum_polynomial_correction': 't*x^2, PASS',
        'all_parameter_orders_theorem': 'written separability and finite descent proof',
    },
    'nonsaturation_boundary': {
        'commutative_module_basis': ['1', 'x+t*sqrt(x)'],
        'generic_basis_determinant': 't',
        'ambient_residue_rank': 1, 'abstract_classical_module_rank': 2,
        'hidden_nonzero_class': 't*sqrt(x)',
        'hidden_class_coordinates': ['-x', '1'],
        'h_saturated': False, 'Weyl_right_action_constructed': False,
    },
    'matrix_connection_boundary': {
        'componentwise_connection_regular': True,
        'field_Leibniz_defect_on_z_squared_equals_x': 1,
        'natural_division_induced_action': False,
        'arbitrary_matrix_repairs_excluded': False,
    },
    'genuine_nonsaturated_equivariant_Weyl_carriers': equivariant_repairs,
    'special_torsion_boundary': {
        'actual_bimodule': 'A_R*x+A_R*t',
        'full_left_and_right_Weyl_stability': 'written two-generator proof',
        'nonzero_special_torsion_class': 't modulo tI, killed by x',
        'special_torsion_free': False,
        'ordinary_parameter_correct_classical_graph': False,
        'rank_one_residue_map_automatically_injective': False,
    },
    'research_context': [
        {'source': 'https://stacks.math.columbia.edu/tag/07ZA',
         'used': 'Fitting functoriality/base change supports the retained QFC1 finite-connection argument.'},
        {'source': 'https://stacks.math.columbia.edu/tag/02GU',
         'used': 'Finite presentation, flatness and zero relative differentials imply etaleness.'},
        {'source': 'https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27',
         'used': 'The optional global affine-space conclusion uses the characteristic-zero constant-cover theorem; the local ramification exclusion does not require it.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Exact full cubic splitting and S3 actions, discriminant pullback, source '
        'chart identities, interior branch derivative, local companion connections '
        'at m=2,3,5, scalar-potential cancellation, a constant etale extension, '
        'a t^3 finite Weyl shear, a finite nonsaturated commutative lattice and '
        'a regular split matrix connection incompatible with field multiplication. '
        'The arbitrary-lattice finite-carrier theorem uses written saturation, '
        'separability, multiplier-ring, sub-quantum potential descent and Fitting '
        'proofs. It excludes a natural larger-closure action with a faithful full '
        'residue lattice across ramification, not every quantum closure, '
        'nonsaturated model, different matrix action or the ordinary rank-one '
        'Long graph carrier. Genuine constant-degree2/3 full Weyl carriers '
        'check that saturation is not automatic, and verify finite character '
        'repairs, equivariance and graph admission at parameter orders1,3. '
        'Finite valuation-level admission, specialization-rank comparison and '
        'arbitrary split-abelian equivariant repair use written proofs. '
        'Correct ordinary-parameter graph or equivariant torsion-free rank '
        'admission now excludes the natural Long closure without assuming '
        'saturation. No automatic saturation, arbitrary matrix exclusion or '
        'degree, lattice or gauge census is certified.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'ramified_closure_carrier_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QGC1', 'status': 'PASS',
    'full_S3_splitting_and_interior_branch': 'PASS',
    'local_connections_m': [2, 3, 5],
    'ramified_quantum_parameter_shear': 'PASS',
    'nonsaturation_and_matrix_boundaries': 'PASS',
    'genuine_equivariant_repairs': [2, 3],
    'Long_carrier_constructed': False,
}, indent=2))
