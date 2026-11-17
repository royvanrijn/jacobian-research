#!/usr/bin/env python3
"""Exact declared controls for QSC1's finite spectral/curvature framework."""
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h = sp.symbols('h', nonzero=True)
x, y, px, py = sp.symbols('x y px py')
r, v, a, d, s = sp.symbols('r v a d s')
spectral = sp.symbols('spectral', nonzero=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for receipt_name in ['adler_admission_controls.json', 'matrix_carrier_controls.json']:
    prior = json.loads((HERE/receipt_name).read_text())
    assert prior['status'] == 'PASS', receipt_name
    for name, digest in prior['hashes'].items():
        assert sha(HERE/name) == digest, name


# The actual inverse localized Long chart: all expressions are Laurent in v.
Fgrad = 3*v-sp.Rational(5, 3)*r*v*v+sp.Rational(1, 3)*r*r*v**3
Ggrad = -sp.Rational(5, 12)*v**4+sp.Rational(2, 15)*r*v**5
X = 1/v
Q = (2*v-r*v*v)/3
Z = -3*(d-Ggrad)/v**2
P = -v**3*a+2*(r*v-1)*(d-Ggrad)+v**3*Fgrad
assert sp.cancel(X*(2-3*X*Q)-r) == 0
assert sp.cancel(-X**3*P-sp.Rational(2, 3)*X*(1-3*X*Q)*Z+Fgrad-a) == 0
assert sp.cancel(-Z/(3*X*X)+Ggrad-d) == 0


def poisson(f, g):
    return sp.cancel(
        sp.diff(f, d)*sp.diff(g, r)-sp.diff(f, r)*sp.diff(g, d)
        +(sp.diff(f, a)*sp.diff(g, v)-sp.diff(f, v)*sp.diff(g, a))/v
    )


assert poisson(P, X) == 1
assert poisson(Z, Q) == 1
for f, g in [(P, Z), (P, Q), (Z, X), (X, Q)]:
    assert poisson(f, g) == 0
g = 4*a-3*r*a*a
Sclass = g+2*v
assert poisson(r, Sclass) == 0
assert poisson(Sclass, a) == -2/v
assert sp.Poly(s-g, a).nth(1) == -4
for denominator in [r, s, r+s]:
    assert sp.Poly(denominator, r, s, a).gcd(sp.Poly(s-g, r, s, a)).total_degree() == 0
test_local = a**3/(s-g)**2+r*s
as_chart = test_local.subs(s, Sclass)
assert sp.cancel(
    poisson(Sclass, as_chart).subs(v, (s-g)/2)
    +4*sp.diff(test_local, a)/(s-g)
) == 0


# Exact finite rank-two Weyl PBW algebra; normal order x,y,px,py.
ZERO = {}
ONE = {(0, 0, 0, 0): sp.Integer(1)}


def clean(op):
    return {key: z for key, value in op.items() if (z := sp.cancel(value)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, 0)+value
    return clean(out)


def scale(op, value):
    return clean({key: value*x for key, x in op.items()})


def mul(left, right):
    out = {}
    for (ix, iy, ip, iq), c in left.items():
        for (jx, jy, jp, jq), e in right.items():
            for k in range(min(ip, jx)+1):
                for ell in range(min(iq, jy)+1):
                    key = (ix+jx-k, iy+jy-ell, ip+jp-k, iq+jq-ell)
                    value = (
                        c*e*comb(ip, k)*(factorial(jx)//factorial(jx-k))
                        *comb(iq, ell)*(factorial(jy)//factorial(jy-ell))*h**(k+ell)
                    )
                    out[key] = out.get(key, 0)+value
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


def mono(key, value=1):
    return {key: sp.sympify(value)}


Xop = mono((1, 0, 0, 0))
Yop = mono((0, 1, 0, 0))
Pop = mono((0, 0, 1, 0))
Qop = mono((0, 0, 0, 1))
assert comm(Pop, Xop) == scale(ONE, h)
assert comm(Qop, Yop) == scale(ONE, h)

# A globally admitted finite nonpositional identity deformation.
F1 = add(Xop, scale(Qop, h**2))
F2 = add(Yop, scale(Pop, h**2))
assert comm(F1, F2) == ZERO
F12 = mul(F1, F2)
T1 = add(scale(Pop, -1), scale(F12, 2))
T2 = add(scale(Qop, -1), power(F1, 2), scale(power(F2, 2), 3))
for i, F in enumerate([F1, F2]):
    for j, T in enumerate([T1, T2]):
        assert comm(F, T) == (scale(ONE, h) if i == j else ZERO)
assert comm(T1, T2) == ZERO

# The two actual Koszul maps on a finite noncentral test operator.
test = add(mul(power(Xop, 2), Qop), mul(Yop, power(Pop, 2)))
d0 = [comm(F1, test), comm(F2, test)]
d1d0 = add(comm(F1, d0[1]), scale(comm(F2, d0[0]), -1))
assert d1d0 == ZERO

# Right evaluation of a polynomial in two central spectral parameters.
# Spectral polynomial key (u1 power,u2 power), value a finite PBW operator.
def spectral_add(*polys):
    out = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = add(out.get(key, ZERO), value)
    return {key: value for key, value in out.items() if value}


def spectral_mul(left, right):
    out = {}
    for (i, j), value in left.items():
        for (k, ell), other in right.items():
            key = (i+k, j+ell)
            out[key] = add(out.get(key, ZERO), mul(value, other))
    return {key: value for key, value in out.items() if value}


def evaluate(poly):
    out = ZERO
    for (i, j), value in poly.items():
        out = add(out, mul(mul(value, power(F1, i)), power(F2, j)))
    return out


L1 = {(0, 0): F1, (1, 0): scale(ONE, -1)}
L2 = {(0, 0): F2, (0, 1): scale(ONE, -1)}
spectral_test = {(1, 0): Xop, (0, 2): Pop, (1, 1): mul(Yop, Qop)}
assert evaluate(spectral_mul(spectral_test, L1)) == ZERO
assert evaluate(spectral_mul(spectral_test, L2)) == ZERO
assert spectral_add(spectral_mul(L1, L2),
                    {key: scale(value, -1) for key, value in spectral_mul(L2, L1).items()}) == {}

# Universal spectral torsion: every fixed fiber is exact, the global class is not.
Sconst = power(Pop, 2)
torsion_witness = scale(mul(Xop, Pop), sp.Rational(1, 2))
assert comm(Sconst, torsion_witness) == scale(Sconst, h)
generic_witness = scale(mul(Xop, Pop), 1/(2*spectral))
assert comm(Sconst, generic_witness) == scale(Sconst, h/spectral)
zero_witness = add(
    scale(power(Xop, 2), 1/(2*h)),
    scale(mul(power(Xop, 3), Pop), -1/(3*h**2)),
)
assert comm(Sconst, zero_witness) == add(
    scale(ONE, h), scale(mul(power(Xop, 2), Sconst), -2/h)
)

# Finite curvature correction and its choice-independent class.
linear1 = scale(Pop, -1)
linear2 = add(scale(Qop, -1), mul(power(F1, 2), F2))
C = comm(linear1, linear2)
assert C == scale(F12, -2*h)
b1 = ZERO
b2 = scale(mul(power(F1, 2), F2), -1)
assert add(C, comm(linear1, b2), scale(comm(linear2, b1), -1)) == ZERO
assert comm(add(linear1, b1), add(linear2, b2)) == ZERO


def mmul(left, right):
    rows, middle, columns = len(left), len(right), len(right[0])
    return [[add(*[mul(left[i][k], right[k][j]) for k in range(middle)])
             for j in range(columns)] for i in range(rows)]


def madd(left, right):
    return [[add(x, y) for x, y in zip(row, other)] for row, other in zip(left, right)]


def mscale(matrix, value):
    return [[scale(op, value) for op in row] for row in matrix]


def mcomm(left, right):
    return madd(mmul(left, right), mscale(mmul(right, left), -1))


def diag(op):
    return [[op, ZERO], [ZERO, op]]


def mtrace(matrix):
    return scale(add(*[matrix[i][i] for i in range(len(matrix))]), sp.Rational(1, len(matrix)))


I = diag(ONE)
Mzero = [[ZERO, ZERO], [ZERO, ZERO]]
# A genuine polynomial GL2(B) gauge with nontrivial off-diagonal connection terms.
G = [[add(ONE, F12), F1], [F2, ONE]]
Ginv = [[ONE, scale(F1, -1)], [scale(F2, -1), add(ONE, F12)]]
assert mmul(G, Ginv) == I and mmul(Ginv, G) == I
U1 = mmul(mmul(Ginv, diag(T1)), G)
U2 = mmul(mmul(Ginv, diag(T2)), G)
assert U1 != diag(T1) and U2 != diag(T2)
assert mcomm(U1, U2) == Mzero
for i, F in enumerate([F1, F2]):
    for j, U in enumerate([U1, U2]):
        assert mcomm(diag(F), U) == (mscale(I, h) if i == j else Mzero)
average1, average2 = mtrace(U1), mtrace(U2)
assert average1 == T1 and average2 == T2
assert comm(average1, average2) == ZERO
entry_bound = max(degree(op) for matrix in [U1, U2] for row in matrix for op in row)
trace_bound = max(degree(average1), degree(average2))
assert trace_bound <= entry_bound
B1 = madd(U1, mscale(diag(average1), -1))
B2 = madd(U2, mscale(diag(average2), -1))
assert mtrace(B1) == ZERO and mtrace(B2) == ZERO
for matrix in [B1, B2]:
    for row in matrix:
        for op in row:
            assert comm(F1, op) == ZERO and comm(F2, op) == ZERO
assert mtrace(mcomm(B1, B2)) == ZERO

# A nonconstant finite projective corner with identity classical graph.
W = [[ONE, scale(Xop, h)], [ZERO, ONE]]
Winv = [[ONE, scale(Xop, -h)], [ZERO, ONE]]
e11 = [[ONE, ZERO], [ZERO, ZERO]]
e = mmul(mmul(W, e11), Winv)
assert mmul(e, e) == e
assert e != e11


def corner(op):
    return mmul(mmul(W, [[op, ZERO], [ZERO, ZERO]]), Winv)


for i, F in enumerate([F1, F2]):
    for j, T in enumerate([T1, T2]):
        assert mcomm(corner(F), corner(T)) == (mscale(e, h) if i == j else Mzero)
assert mcomm(corner(T1), corner(T2)) == Mzero

# Noncommutative-centralizer boundary: flat outer matrices need not trace to a flat pair.
innerA = sp.Matrix([[0, 1], [0, 0]])
innerB = sp.Matrix([[0, 0], [1, 0]])
zero2 = sp.zeros(2)
outer1 = sp.diag(innerA, zero2)
outer2 = sp.diag(zero2, innerB)
assert outer1*outer2-outer2*outer1 == sp.zeros(4)
averaged_curvature = (innerA*innerB-innerB*innerA)/4
assert averaged_curvature != sp.zeros(2)

# Escaping-boundary curvature: a PBW-compatible twisted Ore algebra.
qx, qy = 1-h*x, 1-h*y
curvature = h**2/(qx*qy)
double_residue = sp.residue(sp.residue(curvature, x, 1/h), y, 1/h)
assert double_residue == 1
test_function = x*x*y+1/(qx*qy**2)
# The p2,p1,f critical overlap reduces to commutation of base derivations.
assert sp.cancel(
    h**2*sp.diff(test_function, x, y)-h**2*sp.diff(test_function, y, x)
) == 0
formal_correction = -sp.log(qx)/qy
assert sp.cancel(sp.diff(formal_correction, x)-h/(qx*qy)) == 0
series_controls = []
for N in range(1, 6):
    coefficient = sum(x**k*y**(N-k)/sp.Integer(k) for k in range(1, N+1))
    expected_derivative = sum(x**j*y**(N-1-j) for j in range(N))
    assert sp.expand(sp.diff(coefficient, x)-expected_derivative) == 0
    assert sp.Poly(coefficient, x, y).total_degree() == N
    assert sp.Poly(coefficient, x, y).coeff_monomial(x**N) == sp.Rational(1, N)
    series_controls.append({'h_order': N, 'spatial_degree': N})
truncated = sum(
    h**N*sum(x**k*y**(N-k)/sp.Integer(k) for k in range(1, N+1))
    for N in range(1, 6)
)
assert sp.series(curvature-h*sp.diff(truncated, x), h, 0, 7).removeO() == 0

hash_names = [
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md', 'check_spectral_curvature.py',
    'QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md', 'phase_carrier_controls.json',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md', 'matrix_carrier_controls.json',
    'CENTRALIZER_NORMAL_FORM.md', 'ALGEBRAIC_INFINITY_LIMITS.md',
    'ADLER_PERIODS_AND_FINITE_ADMISSION.md', 'adler_admission_controls.json',
    'BOUNDED_REALIZATION.md',
]
receipt = {
    'claim': 'QSC1', 'status': 'PASS',
    'Long_classical_chart': {
        'inverse_ring_formulas': 'PASS',
        'canonical_brackets': 'PASS',
        'r_kernel_removes_d': 'written all-degree kernel proof',
        's_kernel_on_r_s_a': '-4/(s-g)*partial_a',
        'unit_a_linear_coefficient': -4,
        'joint_polynomial_centralizer': 'written proof K[R,S]',
    },
    'nonpositional_quantum_control': {
        'classical_graph': 'identity',
        'positions': ['x+h^2*py', 'y+h^2*px'],
        'all_four_scalar_relations': 'PASS',
        'global_Long_deformation_constructed': False,
    },
    'Koszul_d1_d0_and_right_evaluation': 'PASS',
    'universal_spectral_torsion_control': {
        'S': 'p^2', 'global_class_nonzero': 'QVL1 written constant-coefficient theorem',
        'annihilated_by_spectral_u': True,
        'nonzero_fiber_witness': 'x*p/(2*lambda)',
        'zero_fiber_witness': 'x^2/(2*h)-x^3*p/(3*h^2)',
        'all_fiber_translation_classes_zero': True,
        'global_translation_class_zero': False,
    },
    'finite_curvature_correction': {
        'curvature': '-2*h*F1*F2', 'b1': '0', 'b2': '-F1^2*F2',
        'flat_after_finite_correction': True,
    },
    'scalar_position_matrix_completion': {
        'size': 2, 'gauge': '[[1+F1*F2,F1],[F2,1]]',
        'nontrivial_offdiagonal_entries': True,
        'full_matrix_relations': 'PASS',
        'normalized_trace_is_flat_scalar_pair': True,
        'entry_spatial_degree_bound': entry_bound,
        'trace_spatial_degree_bound': trace_bound,
    },
    'nonconstant_projective_corner_identity_control': 'PASS',
    'noncommutative_centralizer_boundary': {
        'outer_constant_matrices_commute': True,
        'average_curvature': str(averaged_curvature),
        'average_pair_flat': False,
    },
    'escaping_boundary_curvature': {
        'algebra': 'twisted Ore algebra, not the finite Weyl algebra',
        'PBW_overlap_compatibility': 'PASS',
        'curvature': 'h^2/((1-h*x)*(1-h*y))',
        'double_residue': str(double_residue),
        'formal_correction': '-log(1-h*x)/(1-h*y)',
        'series_controls': series_controls,
        'formal_flatness_through_h6': 'PASS',
        'no_rational_flat_correction': 'written iterated-residue proof',
    },
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'research_context': [
        {
            'source': 'https://arxiv.org/html/2604.16636',
            'used': 'Canonical first-order cohomology classes as related context; no identification with the spectral classes is asserted.',
        },
        {
            'source': 'https://arxiv.org/abs/1904.03128',
            'used': 'Graph-bimodule formulation as related context; direct spectral and curvature calculations are supplied.',
        },
        {
            'source': 'https://stacks.math.columbia.edu/tag/09CN',
            'used': 'Principal-parts sequence and connection viewpoint as context; the relative noncommutative Koszul resolution is proved directly.',
        },
    ],
    'coverage_boundary': (
        'Declared Long chart identities and kernel witnesses, a globally finite '
        'nonpositional identity deformation, Koszul/right-evaluation controls, '
        'universal versus fiber spectral torsion, finite curvature correction, '
        'nontrivial scalar-position matrix gauges and their traces, an identity '
        'projective corner, noncommutative-centralizer boundary and escaping '
        'curvature model. Universal Long joint-centralizer exhaustion, all '
        'commuting deformations, carrier faithfulness, exact spectral/curvature '
        'criteria, all matrix sizes, identity stabilization and nonzero '
        'boundary curvature use written proofs. The twisted Ore boundary is '
        'not a finite Weyl or Long construction. No full Long carrier or '
        'nonzero obstruction for every Long deformation is certified.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'spectral_curvature_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QSC1', 'status': 'PASS',
    'matrix_entry_bound': entry_bound, 'scalar_trace_bound': trace_bound,
    'all_fibers_zero_global_class_nonzero': True,
    'escaping_curvature_double_residue': str(double_residue),
}, indent=2))
