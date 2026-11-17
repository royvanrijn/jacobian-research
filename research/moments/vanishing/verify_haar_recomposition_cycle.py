"""Selected exact controls for the written primitive-cycle transfer theorem.

These do not compute exceptional Betti monodromy or verify universal topology.
"""
from pathlib import Path
import hashlib
import json
import sympy as s


def sum_powers(A, k):
    return sum((A**v for v in range(k)), s.zeros(A.rows))


def cyclic(T, e):
    n = T.rows
    A = s.zeros(e*n)
    for j in range(e-1):
        A[(j+1)*n:(j+2)*n, j*n:(j+1)*n] = s.eye(n)
    A[:n, (e-1)*n:e*n] = T
    assert A**e == s.diag(*([T]*e))
    return A


def eig1(A):
    q = s.symbols('q')
    cp = A.charpoly(q).as_expr()
    multiplicity = 0
    while s.rem(cp, q-1, q) == 0:
        cp = s.div(cp, q-1, q)[0]
        multiplicity += 1
    if multiplicity == 0:
        return s.zeros(A.rows)
    p = s.rem(cp*s.invert(cp, (q-1)**multiplicity, q),
              cp*(q-1)**multiplicity, q)
    E = s.zeros(A.rows)
    for (v,), coeff in s.Poly(p, q).terms():
        E += coeff*A**v
    assert E**2 == E and E*A == A*E
    return E


# Universal nilpotent-jet expression: S_k(1+N) has constant term k.
k = s.symbols('k', positive=True, integer=True)
N = s.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
Sk = k*s.eye(3) + k*(k-1)/2*N + k*(k-1)*(k-2)/6*N**2
assert s.factor(Sk.det()) == k**3
for value in range(1, 6):
    assert Sk.subs(k, value) == sum_powers(s.eye(3)+N, value)

# Cyclic ramification blocks retain both eigenvalue-1 Jordan motion and
# other eigenvalues. Primitive full-return motion follows e turns, not one.
T = s.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, -1]])
E = s.diag(1, 1, 0)
blocks = []
for e in (2, 3):
    A = cyclic(T, e)
    P = eig1(A)
    S = s.diag(*([E]*e))
    assert A*S == S*A and S*P == P
    d = s.Matrix(list(range(1, A.rows+1)))
    full_return = sum_powers(A, e)*d
    L = 6
    D = sum_powers(A, L)*d
    for j in range(e):
        selection = slice(3*j, 3*j+3)
        assert D[selection, :] == sum_powers(T, L//e)*full_return[selection, :]
    assert P*(s.eye(A.rows)-S) == s.zeros(A.rows)
    assert P*sum_powers(A, L)*P == sum_powers(A, L)*P
    blocks.append({'branch_order': e, 'common_loop_power': L,
                   'primitive_loop_count': L//e,
                   'full_return_and_spectral_inclusion_exact': True})

# Without eventual purity the primitive integer motion can be nonzero even
# if the composite eigenvalue-1 variation is zero. This is an ABSTRACT block.
J = s.Matrix([[1, 1], [0, 1]])
A = cyclic(J, 2)
B = s.eye(2)+(J-s.eye(2))/2
assert B**2 == J
v = s.Matrix([0, 1])
d = (-B*v).col_join(v)
assert eig1(A)*d == s.zeros(4, 1)
power_motion = (s.eye(4)+A)*d
assert power_motion == s.Matrix([s.Rational(1, 2), 0, -s.Rational(1, 2), 0])
assert power_motion != s.zeros(4, 1)

# The finite unit-jet formula is independent of the degree of H. Its first
# four jets suffice here although H has a freely varying degree-five term.
z, b1, b2, b5 = s.symbols('z b1 b2 b5')
H = b1*z+b2*z**2+b5*z**5
unit_moments = [1, 2, -3, 5]
jet_values = []
for m in range(1, 6):
    value = s.expand(sum(a*s.expand(H**m).coeff(z, j)
                         for j, a in enumerate(unit_moments)))
    jet_values.append(str(value))
    if m > 3:
        assert value == 0
assert jet_values[:3] == ['2*b1 - 3*b2', '-3*b1**2 + 10*b1*b2', '5*b1**3']

# Actual trace-square control: the primitive root 0 is regular but the
# composite is non-pure; zero variation alone is not a purity test.
catalan = [s.binomial(2*m, m)/(m+1) for m in range(7)]
assert catalan == [1, 1, 2, 5, 14, 42, 132]

out = {
    'completed': True, 'sympy_version': s.__version__,
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'symbolic_unipotent_sum_determinant': str(s.factor(Sk.det())),
    'cyclic_ramification_controls': blocks,
    'eventual_purity_cannot_be_dropped_countercontrol_exact': True,
    'abstract_countercontrol_power_motion': [str(x) for x in power_motion],
    'finite_unit_jet_control': jet_values,
    'trace_square_Catalan_moments': [str(x) for x in catalan],
    'universal_topology_and_all_degree_classification_machine_verified': False,
    'actual_exceptional_primitive_monodromy_computed': False,
    'primitive_purity_reflection_claimed': False,
    'uniform_Haar_orbit_condition_proved_for_every_exceptional_factor': False,
    'abstract_matrix_controls_are_exceptional_phase_realizations': False,
    'all_coefficient_ring_topological_promotion_claimed': False,
    'full_exceptional_MZ_cases': 'OPEN on BOTH quotients',
    'independent_replay': False, 'formal_verification': False,
    'external_review': False, 'new_degree_support_or_prefix_campaign': False}
Path(__file__).with_name('haar_recomposition_cycle_verification.json').write_text(
    json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
