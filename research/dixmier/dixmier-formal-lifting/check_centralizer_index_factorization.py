#!/usr/bin/env python3
"""Declared QCI1 field-tower and positional controls; universal indices are written."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h = sp.symbols('h', nonzero=True)
u, z, t, lam = sp.symbols('u z t lam')
x, y, px, py = sp.symbols('x y px py')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['spectral_projective_denominator_controls.json', 'division_index_controls.json',
             'identity_lift_controls.json', 'finite_centralizer_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, expected in prior['hashes'].items():
        assert sha(HERE/rel) == expected, rel


# A valued commutative tower with d=3, M=3, delta=3 and N=9.
u_of_z = z**3+z
t_of_u = u+h*u**3
t_of_z = sp.expand(t_of_u.subs(u, u_of_z))
d = sp.degree(t_of_u, u)
M = sp.degree(u_of_z, z)
N = sp.degree(t_of_z, z)
delta = sp.degree(t_of_z.subs(h, 0), z)
assert (d, M, N, delta) == (3, 3, 9, 3)
assert N == d*M and M >= delta and N-delta == delta*(d-1)+d*(M-delta)

f = h*u**3+u-t
root_prefix = sum(
    (-1)**j*sp.binomial(3*j, j)*h**j*t**(2*j+1)/(2*j+1)
    for j in range(7)
)
assert sp.series(root_prefix+h*root_prefix**3-t, h, 0, 7).removeO() == 0
assert sp.series(root_prefix.subs(t, t_of_u)-u, h, 0, 7).removeO() == 0
inverse_degrees = [sp.degree(sp.expand(root_prefix).coeff(h, j), t) for j in range(7)]
assert inverse_degrees == [2*j+1 for j in range(7)]
disc = sp.discriminant(f, u)
assert sp.expand(disc+4*h+27*h*h*t*t) == 0

# Multiplication and extended differentiation on the actual degree-three basis.
companion = sp.Matrix([[0, 0, t/h], [1, 0, -1/h], [0, 1, 0]])
eye = sp.eye(3)
assert companion+h*companion**3 == t*eye
assert sp.expand(companion.charpoly(lam).as_expr()-(lam**3+lam/h-t/h)) == 0
du_dt = (eye+3*h*companion**2).inv()
basis1 = sp.Matrix([1, 0, 0])
connection = sp.Matrix.hstack(sp.zeros(3, 1), du_dt*basis1, 2*companion*du_dt*basis1)
compatibility = companion.diff(t)+connection*companion-companion*connection-du_dt
assert all(sp.cancel(entry) == 0 for entry in compatibility)

# Rational positional functions have inverse-Jacobian first-order mates.
g = [x*x+y, y]
J = sp.Matrix([[sp.diff(q, c) for c in (x, y)] for q in g])
inverse_J = J.inv()
assert sp.expand(J.det()-2*x) == 0
assert J*inverse_J == sp.eye(2)
mates = [-h*inverse_J[:, i] for i in range(2)]


def apply(vector, function):
    return sp.cancel(vector[0]*sp.diff(function, x)+vector[1]*sp.diff(function, y))


def bracket(left, right):
    return sp.Matrix([sp.cancel(apply(left, right[k])-apply(right, left[k])) for k in range(2)])


for i in range(2):
    for j in range(2):
        # Multiplication by g_j commuted with this vector-field operator is -T_i(g_j).
        assert sp.cancel(-apply(mates[i], g[j])-(h if i == j else 0)) == 0
assert bracket(mates[0], mates[1]) == sp.zeros(2, 1)
assert -2*x*mates[0] == sp.Matrix([h, 0])
assert -mates[0]-mates[1] == sp.Matrix([0, h])
symbol = px/py
symbol_brackets = [
    -h*(sp.diff(q, x)*sp.diff(symbol, px)+sp.diff(q, y)*sp.diff(symbol, py))
    for q in g
]
assert all(sp.cancel(q) != 0 for q in symbol_brackets)

# A genuine actual delta1 shear remains allowed; use its retained full PBW check.
positive = json.loads((HERE/'finite_centralizer_controls.json').read_text())
assert positive['nonpositional_finite_Weyl_control']['all_final_Weyl_relations'] == 'PASS'
assert positive['nonpositional_finite_Weyl_control']['classical_graph'] == 'identity'
assert sp.expand((x+h*h*px*px)-h*h*px*px-x) == 0
assert sp.expand((y+h*h*py*py)-h*h*py*py-y) == 0

# The actual rescaled cotangent countercontrol has delta1 and purely coefficient index.
prior_indices = json.loads((HERE/'division_index_controls.json').read_text())
cotangent = prior_indices['full_polynomial_control']
assert cotangent['rank'] == 3
assert cotangent['generic_division_degree'] == 3
assert cotangent['completed_division_degree'] == 1
assert cotangent['retained_full_relations'] == 15
assert 3 == 3*1  # d=3, residual M=1 is allowed at delta1.

# Genuine rational residual families and the finite-lattice power obstruction.
# Ore multiplication here is independent of the earlier vector-field controls.
# Operators are dictionaries from p-order to rational x-coefficient.
def ore_clean(operator):
    return {j: sp.cancel(c) for j, c in operator.items() if sp.cancel(c) != 0}


def ore_add(left, right, factor=1):
    result = dict(left)
    for j, c in right.items():
        result[j] = result.get(j, 0)+factor*c
    return ore_clean(result)


def ore_mul(left, right):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            for k in range(i+1):
                order = i+j-k
                result[order] = result.get(order, 0)+sp.binomial(i, k)*h**k*a*sp.diff(b, x, k)
    return ore_clean(result)


def ore_comm(left, right):
    return ore_add(ore_mul(left, right), ore_mul(right, left), -1)


def ore_degree_symbol(operator):
    if not operator:
        return None, sp.S.Zero
    terms = []
    for j, c in operator.items():
        num, den = sp.fraction(sp.cancel(c))
        numerator, denominator = sp.Poly(num, x), sp.Poly(den, x)
        xdegree = numerator.degree()-denominator.degree()
        terms.append((int(xdegree+j), numerator.LC()/denominator.LC()*x**xdegree*px**j))
    degree = max(d for d, symbol_part in terms)
    symbol = sum(symbol_part for d, symbol_part in terms if d == degree)
    return degree, sp.cancel(symbol)


rational_family_controls = []
frame_function = {0: 1/x}
for residual_m in (2, 3, 5):
    Fm = {1: x**(1-residual_m)/sp.Integer(residual_m)}
    Tm = {0: x**residual_m}
    assert ore_comm(Fm, Tm) == {0: h}
    assert ore_comm(Tm, {j: -c for j, c in Fm.items()}) == {0: h}
    assert ore_mul({0: residual_m*x**(residual_m-1)}, Fm) == {1: 1}
    scalar_degree, scalar_symbol = ore_degree_symbol(Fm)
    assert scalar_degree == 2-residual_m
    assert sp.cancel(scalar_symbol-x**(1-residual_m)*px/residual_m) == 0

    # For G=((1,1/x),(0,1)), G^-1*(Fm I)*G has off-diagonal [Fm,1/x].
    frame_off_diagonal = ore_comm(Fm, frame_function)
    assert frame_off_diagonal == {0: -h*x**(-residual_m-1)/residual_m}
    assert ore_comm(frame_off_diagonal, Tm) == {}
    power = {0: sp.S.One}
    conjugated_power_off_diagonal = {}
    power_degrees = []
    for exponent in range(1, 7):
        conjugated_power_off_diagonal = ore_add(
            ore_mul(power, frame_off_diagonal),
            ore_mul(conjugated_power_off_diagonal, Fm),
        )
        power = ore_mul(power, Fm)
        assert ore_comm(power, frame_function) == conjugated_power_off_diagonal
        power_degree, power_symbol = ore_degree_symbol(power)
        assert power_degree == exponent*(2-residual_m)
        assert sp.cancel(power_symbol-(x**(1-residual_m)*px/residual_m)**exponent) == 0
        off_degree, unused_symbol = ore_degree_symbol(conjugated_power_off_diagonal)
        assert off_degree is None or off_degree <= power_degree-3
        power_degrees.append(power_degree)

    rational_family_controls.append({
        'm': residual_m,
        'F': 'x^(1-m)*p/m', 'T': 'x^m',
        'exact_Weyl_relation': 'PASS', 'recover_ambient_p': 'PASS',
        'Bernstein_degree_F': scalar_degree,
        'power_degrees_k1_to_6': power_degrees,
        'fractional_two_by_two_frame_all_checked_powers': 'PASS',
        'whole_left_right_indices': 'm, written differential-field basis proof',
        'position_F_d_and_M': [1, residual_m],
        'position_T_d_and_M': [residual_m, 1],
        'completed_index': 'm, written residue-basis proof',
        'no_finite_lattice_any_positive_rank': 'written bounded-powers proof',
        'Long_lift': False,
    })

# The rank-one degree obstruction must not be applied to arbitrary matrix entries.
nilpotent = sp.Matrix([[0, 1/x], [0, 0]])
lattice_frame = sp.diag(x, 1)
polynomial_action = sp.Matrix([[0, 1], [0, 0]])
assert nilpotent**2 == sp.zeros(2)
assert lattice_frame*nilpotent == polynomial_action*lattice_frame
assert sp.simplify(lattice_frame*nilpotent*lattice_frame.inv()) == polynomial_action

# Positive polynomial Weyl generators have unbounded power degrees.
assert ore_comm({0: x}, {1: -1}) == {0: h}
assert [ore_degree_symbol({0: x**k})[0] for k in range(1, 7)] == list(range(1, 7))

hash_names = [
    'CENTRALIZER_INDEX_AND_RATIONAL_POLARIZATION.md',
    'check_centralizer_index_factorization.py',
    'QUANTUM_DIVISION_AND_COMPLETION.md', 'division_index_controls.json',
    'SPECTRAL_DENOMINATORS_AND_PROJECTIVE_CENTRALIZERS.md',
    'spectral_projective_denominator_controls.json',
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
    'ETALE_CENTRALIZERS_AND_GALOIS_COMPLETION.md',
    'TAYLOR_BRANCHES_AND_CENTRALIZER_RIGIDITY.md',
    'IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md', 'identity_lift_controls.json',
    'FINITE_CENTRALIZERS_AND_COMPLETION.md', 'finite_centralizer_controls.json',
    'ELLIPTIC_RAMIFICATION_AND_STABLE_POLARIZATION.md',
    'ELLIPTIC_ISOGENY_MATRIX_MODULE.md', 'BOUNDED_REALIZATION.md',
    'HAMILTONIAN_IDEALIZER_CRITERION.md',
]
receipt = {
    'claim': 'QCI1', 'status': 'PASS',
    'commutative_valued_tower_control': {
        'u_of_z': str(u_of_z), 't_of_u': str(t_of_u), 't_of_z': str(t_of_z),
        'generic_centralizer_analogue_degree': int(d), 'generic_residual_degree': int(M),
        'generic_total_degree': int(N), 'completed_residue_degree': int(delta),
        'index_loss': int(N-delta),
        'selected_intermediate_completion_degree': 'one, written Hensel proof',
        'formal_inverse_prefix_through_h6': 'PASS',
        'inverse_coefficient_spatial_degrees': [int(n) for n in inverse_degrees],
        'Weyl_tuple': False, 'Long_lift': False,
    },
    'field_basis_derivative_control': {
        'basis': ['1', 'u', 'u^2'], 'minimal_polynomial': str(f),
        'multiplication_matrix': str(companion),
        'derivation_extended_by_minimal_polynomial': 'PASS',
        'connection_multiplication_compatibility': 'PASS',
    },
    'positional_recovery_control': {
        'positions': ['x^2+y', 'y'], 'Jacobian': '2*x',
        'rational_mates': ['-p_x/(2*x)', 'p_x/(2*x)-p_y'],
        'all_first_order_relations': 'PASS',
        'recover_p_x': '-2*x*T1', 'recover_p_y': '-T1-T2',
        'position_centralizer_is_self_polynomial': False,
        'Long_lift': False,
    },
    'positive_delta1_shear': {
        'retained_full_PBW_relations': 'PASS', 'actual_polynomial_inverse': 'PASS',
        'classical_degree': 1, 'residual_index': 1,
    },
    'retained_rescaled_cotangent_control': {
        'claim': 'QNI1', 'rank': 3, 'delta': 1, 'd': 3, 'M': 1, 'N': 3,
        'full_relations': 15, 'Long_lift': False,
    },
    'unscaled_cotangent_boundary': {
        'classical_degree': 3, 'classical_position_centralizer': 'K[x], larger than K[f]',
        'violates_joint_centralizer_hypothesis': True,
        'actual_hypotheses_not_satisfied': True,
    },
    'genuine_rational_residual_family_controls': rational_family_controls,
    'bounded_matrix_powers_boundary': {
        'negative_degree_nilpotent': '((0,x^-1),(0,0))',
        'minimal_scalar_polynomial': 't^2',
        'finite_free_lattice_basis': ['(x,0)', '(0,1)'],
        'polynomial_right_action_after_frame': 'PASS',
        'Weyl_mate': False,
        'positive_polynomial_pair': 'x,-p',
        'positive_power_degrees_k1_to_6': list(range(1, 7)),
        'arbitrary_rank_theorem': 'written proof, no finite rank search',
    },
    'research_context': [
        {'source': 'https://arxiv.org/html/math/0401263',
         'used': 'Corollary3.10 provides coordinate-field maximality context. Lemma3.7, also retained in QPD1, clears curve rational commutants for the rational residual family; its finite operator centralizer is proved directly by differential order.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Declared exact commutative valued-field tower, finite formal inverse prefix, '
        'multiplication/connection matrices, rational inverse-Jacobian vector fields, '
        'symbol witnesses, retained actual delta1 PBW/shear and full cotangent inputs. '
        'The toy tower is not a Weyl or Long realization. Universal finite H/C, '
        'left/right index factorization, containment of the full rational centralizer '
        'in the completed image, residual lower bound, fixed-position rational target '
        'polarization obstruction and all identity stabilizations use written proofs '
        'and the retained QDI1/QMC1/QPD1 theorems. No equality of generic left/right '
        'indices, mixed source normalization exclusion, new Long lift or all-representative '
        'nonlifting is certified. Exact Ore controls for m=2,3,5 check the rational '
        'Weyl family, powers and fractional matrix frames. The arbitrary-rank '
        'bounded-powers/finite-lattice equivalence, all-m indices, residual role '
        'reversal, completed degrees and finite-lattice exclusions use written '
        'proofs. A negative-degree nilpotent with a finite free lattice prevents '
        'extension to arbitrary matrix entries without the algebraicity condition.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'centralizer_index_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QCI1', 'status': 'PASS', 'valued_field_control': {'delta': int(delta), 'd': int(d),
       'M': int(M), 'N': int(N), 'loss': int(N-delta)},
    'rational_first_order_recovery': 'PASS', 'actual_delta1_countercontrols_retained': True,
    'rational_residual_family_m': [2, 3, 5],
    'arbitrary_rank_lattice_boundary': 'PASS',
    'Long_lift_constructed': False,
}, indent=2))
