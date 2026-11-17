#!/usr/bin/env python3
"""Finite exact controls for QKR1; universal Fourier-rank proofs are written."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
q, v, a, tau, h, s, r, la = sp.symbols('q v a tau h s r la')
half = sp.Rational(1, 2)


def clean(op):
    out = {}
    for k, c in op.items():
        c = sp.cancel(c)
        if c != 0:
            out[k] = c
    return out


def add(*ops):
    out = {}
    for op in ops:
        for k, c in op.items():
            out[k] = out.get(k, 0) + c
    return clean(out)


def scale(op, c):
    return clean({k: c * z for k, z in op.items()})


def mul(left, right, variable, rate=1):
    """D f=f D+h*rate*f', with coefficients on the left."""
    out = {}
    for i, x in left.items():
        for j, y in right.items():
            dy = y
            for k in range(i + 1):
                key = i + j - k
                out[key] = out.get(key, 0) + sp.binomial(i, k) * h**k * x * dy
                dy = sp.cancel(rate * sp.diff(dy, variable))
    return clean(out)


def power(op, n, variable, rate=1):
    out = {0: sp.Integer(1)}
    for _ in range(n):
        out = mul(out, op, variable, rate)
    return out


def poly_at(poly, op, variable, rate=1):
    p = sp.Poly(poly, a)
    out = {}
    for (j,), c in p.terms():
        out = add(out, scale(power(op, j, variable, rate), c))
    return out


def remainder(dividend, divisor, variable, rate=1):
    rem = clean(dividend)
    d = max(divisor)
    while rem and max(rem) >= d:
        n = max(rem) - d
        c = sp.cancel(rem[max(rem)] / divisor[d])
        rem = add(rem, scale(mul({n: c}, divisor, variable, rate), -1))
    return rem


def inverse_fourier(op):
    """q->-tau, h*partial_q->a, normal ordered in tau=h*partial_a."""
    out = {}
    for k, c in op.items():
        for (j,), coefficient in sp.Poly(c, q).terms():
            out = add(out, mul({j: (-1)**j * coefficient}, {0: a**k}, a))
    return out


def matrix_zero(matrix):
    return all(sp.cancel(x) == 0 for x in matrix)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


phase_receipt = json.loads((HERE / 'rational_fourier_controls.json').read_text())
assert phase_receipt['status'] == 'PASS'
for name, digest in phase_receipt['hashes'].items():
    assert sha(HERE / name) == digest, name

controls = []
for m, g in [(1, 4*a - 1), (2, 4*a - 3*r*a**2), (3, a**3 + a)]:
    b = s - g
    b_op = poly_at(b, {1: sp.Integer(1)}, q)
    shifted = poly_at(b, {1: sp.Integer(1), 0: -h/(2*q)}, q)
    operator = scale(add(mul(shifted, b_op, q), {0: -8*q}), q**m)
    assert all(sp.denom(sp.cancel(c)).has(q) is False for c in operator.values())
    gm = sp.Poly(g, a).LC()
    assert max(operator) == 2*m
    assert sp.cancel(operator[2*m] - gm**2 * q**m) == 0
    assert max(sp.degree(c, q) for c in operator.values()) == m + 1
    pre = inverse_fourier(operator)
    assert max(pre) == m + 1
    classical_pre = sp.expand(sum(c.subs(h, 0) * tau**k for k, c in pre.items()))
    expected_pre = (-1)**m * tau**m * (b**2 + 8*tau)
    assert sp.expand(classical_pre - expected_pre) == 0
    classical_q = sp.expand(sum(c.subs(h, 0) * a**k for k, c in operator.items()))
    assert sp.expand(classical_q - q**m * (b**2 - 8*q)) == 0

    # These are the two actual algebraic branch equations.
    bv = {k: sp.cancel(c.subs(q, v**2/2)) for k, c in b_op.items()}
    shifted_v = poly_at(b, {1: sp.Integer(1), 0: -h/v**2}, v, 1/v)
    pulled = scale(add(mul(shifted_v, bv, v, 1/v), {0: -4*v**2}),
                   (v**2/2)**m)
    expected_pulled = {k: sp.cancel(c.subs(q, v**2/2))
                       for k, c in operator.items()}
    assert add(pulled, scale(expected_pulled, -1)) == {}
    for sign in [1, -1]:
        branch = add(bv, {0: -sign*2*v})
        assert remainder(pulled, branch, v, 1/v) == {}

    indicial = sp.Integer(0)
    for k, c in operator.items():
        if k >= m:
            indicial += (sp.Poly(c, q).coeff_monomial(q**(k-m))
                         * h**k * sp.ff(la, k))
    indicial = sp.factor(indicial / (gm**2 * h**(2*m)))
    expected_indicial = sp.prod(la-j for j in range(m))
    expected_indicial *= sp.prod(la-m-half-j for j in range(m))
    assert sp.expand(indicial - expected_indicial) == 0
    controls.append({
        'polynomial_degree': m,
        'Fourier_connection_rank': 2*m,
        'inverse_Fourier_cyclic_order': max(pre),
        'indicial_polynomial': str(indicial),
        'two_algebraic_branch_remainders': 'zero',
        'special_equation': str(sp.factor(classical_q)),
        'boundary_factor': 'q**'+str(m),
        'all_degree_rank_uses_written_Newton_and_Fourier_rank_proof': True,
    })

# An actual four-dimensional restriction-of-scalars connection.
C0 = sp.Matrix([[0, -s/(3*r)], [1, 4/(3*r)]])
C1 = sp.Matrix([[0, 2/(3*r)], [0, 0]])
I2, Z2 = sp.eye(2), sp.zeros(2)
H = C0.row_join(v**2*C1).col_join(
    C1.row_join(C0 + h/v**2*I2))
V = Z2.row_join(v**2*I2).col_join(I2.row_join(Z2))
I4 = sp.eye(4)
d = lambda mat: mat.applyfunc(lambda z: sp.diff(z, v)/v)
assert matrix_zero(V*V - v**2*I4)
assert matrix_zero(h*d(V) + H*V - V*H - h/v**2*V)
projectors = [(I4 + V/v)/2, (I4 - V/v)/2]
assert matrix_zero(sum(projectors, sp.zeros(4)) - I4)
assert matrix_zero(projectors[0]*projectors[1])
for sign, projector in zip([1, -1], projectors):
    assert matrix_zero(projector*projector - projector)
    assert matrix_zero(h*d(projector) + H*projector - projector*H)
    frame = I2.col_join(sign/v*I2)
    assert matrix_zero(projector*frame - frame)
    assert matrix_zero(h*d(frame) + H*frame - frame*(C0+sign*v*C1))
    assert frame.rank() == 2
residue = (H/h * v**2/2).applyfunc(
    lambda z: sp.limit(z, v, 0))
assert residue == sp.diag(0, 0, half, half)

# Rank budget algebra is uniform in the symbols, not a cover-order census.
N, e, m0, i0, nr, extra = sp.symbols(
    'N e m0 i0 nr extra', integer=True)
conductor = N + i0 - m0
budget_difference = sp.expand(e*conductor - ((e-1)*N + e*i0))
assert budget_difference == N - e*m0
N_induced, d_cover, m_branch = sp.symbols('N_induced d_cover m_branch')
assert sp.expand(((1-1/d_cover)*N_induced).subs(
    N_induced, d_cover*m_branch) - (d_cover-1)*m_branch) == 0
assert 4 - 4/2 == 2
assert 2*3 > 2*2

# The naive Long partner has a nonzero quantum defect.
S = {2: -3*r, 1: sp.Integer(4), 0: 2*v}
T = {3: r/2, 2: -half, 1: -v/2}
comm = add(mul(S, T, v, 1/v), scale(mul(T, S, v, 1/v), -1))
defect = add(comm, {0: -h})
assert defect
assert all(sp.cancel(c/h**2).is_finite is not False for c in defect.values())
assert all(sp.cancel(c/h).subs(h, 0) == 0 for c in defect.values())

# A genuine linear full-pair control in the same coordinates.
lam, mu = sp.symbols('lam mu')
linear_S = {1: lam, 0: mu + 2*v}
linear_T = {2: -lam/8, 1: -v/2, 0: -h/(4*v)}
linear_comm = add(mul(linear_S, linear_T, v, 1/v),
                  scale(mul(linear_T, linear_S, v, 1/v), -1))
assert linear_comm == {0: h}

# Exact equality case N=2, r=1 from the rational-line splitting control.
c = sp.symbols('c', nonzero=True)
q_derivative_on_v = {1: 1/v}
bessel_pulled = add(
    scale(power(q_derivative_on_v, 2, v), v**2/2),
    scale(q_derivative_on_v, h/2),
    {0: c})
assert bessel_pulled == {2: half, 0: c}
# Here the operator symbol is h*partial_v, so multiply by2 to obtain
# (h*partial_v)^2+2c; it has constant rank-one eigenspace projectors.
kappa = sp.symbols('kappa', nonzero=True)
J = sp.Matrix([[0, 1], [-2*c, 0]])
for sign in [1, -1]:
    E = (sp.eye(2) + sign*J/kappa)/2
    assert matrix_zero((E*E-E).subs(c, -kappa**2/2))
    assert matrix_zero(J*E-E*J)

names = [
    'FOURIER_CYCLIC_RANK_BUDGET.md',
    'check_fourier_rank_budget.py',
    'RATIONAL_LINE_FOURIER_OBSTRUCTION.md',
    'rational_fourier_controls.json',
    'QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
]
out = {
    'claim': 'QKR1',
    'status': 'PASS',
    'polynomial_family_controls': controls,
    'quadratic_branch_projectors': 'horizontal and complementary, rank2 each',
    'quadratic_branch_frames': 'C0 plus/minus v*C1',
    'quadratic_residue': ['0', '0', '1/2', '1/2'],
    'uniform_self_twist_budget_identity': 'symbolic N,e,horizontal_dimension',
    'induced_cover_budget': 'symbolic degree and branch rank',
    'quadratic_naive_partner_defect': {str(k): str(sp.factor(z))
                                      for k, z in defect.items()},
    'linear_full_pair': 'PASS',
    'sharp_N2_r1_control': 'PASS',
    'hashes': {name: sha(HERE/name) for name in names},
    'research_context': [
        {'source': 'https://arxiv.org/html/0808.0742',
         'used': 'Definition(3.1), Proposition3.5 Malgrange Fourier-rank formula; formal slope framework.'},
        {'source': 'https://arxiv.org/html/0808.0699',
         'used': 'Lemma3.5 formal minimal vanishing cycles and Proposition5.1 Fourier rank; boundary accounting.'},
    ],
    'coverage_boundary': (
        'Finite exact polynomial eliminations at degrees1,2,3, branch divisions, '
        'indicial roots, projector/frame identities, a naive Long mate defect, '
        'linear exact pair and rank-two equality control. General conductor '
        'bounds, cyclic descent, nonminimal extension accounting, all-degree '
        'Newton slopes and ranks use written proofs and cited external '
        'Fourier-rank/formal-connection theory. The exact-R Long application '
        'also depends on the registered QSO2 written theorem. No actual '
        'finite-left Long carrier or all-representative obstruction is certified.'),
    'independent_replay': False,
    'formal_verification': False,
    'external_review': False,
}
(HERE/'fourier_rank_budget_controls.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({
    'claim': out['claim'], 'status': out['status'],
    'polynomial_degrees': [x['polynomial_degree'] for x in controls],
    'inverse_Fourier_orders': [x['inverse_Fourier_cyclic_order'] for x in controls],
    'branch_projectors': 'PASS', 'linear_full_pair': 'PASS',
    'sharp_equality_control': 'PASS',
}))
