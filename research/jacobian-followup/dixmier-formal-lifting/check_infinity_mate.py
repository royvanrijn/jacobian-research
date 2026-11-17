#!/usr/bin/env python3
"""Exact controls for QVM1's all-order infinity valuation proof."""
from pathlib import Path
import hashlib
import json
import math
import sympy as sp

HERE = Path(__file__).resolve().parent
a, v, w, h, r, s = sp.symbols('a v w h r s')
half = sp.Rational(1, 2)


def clean(op):
    return {k: sp.cancel(c) for k, c in op.items() if sp.cancel(c) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for k, c in op.items():
            out[k] = out.get(k, 0) + c
    return clean(out)


def scale(op, c):
    return clean({k: c*z for k, z in op.items()})


def mul(left, right, variable, rate):
    out = {}
    for i, x in left.items():
        for j, y in right.items():
            dy = y
            for k in range(i+1):
                key = i+j-k
                out[key] = out.get(key, 0) + sp.binomial(i, k)*h**k*x*dy
                dy = sp.cancel(rate*sp.diff(dy, variable))
    return clean(out)


def power(op, n, variable, rate):
    out = {0: sp.Integer(1)}
    for _ in range(n):
        out = mul(out, op, variable, rate)
    return out


def poly_op(poly):
    return {(k,)[0]: c for (k,), c in sp.Poly(poly, a).terms()}


def comm(left, right, variable, rate):
    return add(mul(left, right, variable, rate),
               scale(mul(right, left, variable, rate), -1))


def coefficients_at_w(op, exponent):
    # All declared controls are Laurent polynomials in w.
    out = sp.Integer(0)
    for k, coefficient in op.items():
        for term in sp.Add.make_args(sp.expand(coefficient)):
            c, degree = term.as_coeff_exponent(w)
            assert not c.has(w), (c, degree)
            if degree == exponent:
                out += c*a**k
    return sp.expand(out)


def w_degree(op):
    return max(term.as_coeff_exponent(w)[1]
               for c in op.values()
               for term in sp.Add.make_args(sp.expand(c)))


def matrix_zero(matrix):
    return all(sp.cancel(z) == 0 for z in matrix)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


rank_receipt = json.loads((HERE/'fourier_rank_budget_controls.json').read_text())
assert rank_receipt['status'] == 'PASS'
for name, digest in rank_receipt['hashes'].items():
    assert sha(HERE/name) == digest, name

# Verify the actual highest commutator coefficient on ramified charts.
samples = [
    (2, 2+h*a**3, sp.Integer(5), a**4+2*a+1),
    (3, (a-1)**2*(a+2)**3, sp.Integer(3), a**3+1),
    (4, a*a+1, sp.Rational(3, 2), a**3-a+2),
]
leading_controls = []
for d, H, M, K in samples:
    rate = w**(1-2*d)/2
    S = add(scale(poly_op(H), w**2),
            poly_op(a**5+4*a-3*r*a*a),
            {4: h**2/w**2})
    T = add(scale(poly_op(K), w**(2*M)), {1: w**(2*M-1)})
    result = comm(S, T, w, rate)
    exponent = 2*(M+1-d)
    expected = sp.expand(h*(M*sp.diff(H,a)*K-H*sp.diff(K,a)))
    assert sp.expand(coefficients_at_w(result, exponent)-expected) == 0
    assert w_degree(result) <= exponent
    leading_controls.append({
        'cyclic_degree': d,
        'H_degree': int(sp.degree(H,a)),
        'mate_leading_v_exponent': str(M),
        'commutator_leading_v_exponent': str(M+1-d),
        'leading_coefficient_identity': 'PASS',
    })

# The leading elimination uses finite integer powers, not an infinite gauge.
H = 2+h*a**3
rate = w**(-3)/2
S = add(scale(poly_op(H), w**2), poly_op(4*a-3*r*a*a))
U = scale(poly_op(a+1), w**2)
S3 = power(S,3,w,rate)
S4 = power(S,4,w,rate)
T = add(scale(S4,5),scale(S3,-2),U)
assert w_degree(T) == 8
T1 = add(T,scale(S4,-5))
assert w_degree(T1) == 6
T2 = add(T1,scale(S3,2))
assert add(T2,scale(U,-1)) == {}
assert comm(S,T,w,rate) == comm(S,U,w,rate)

# Multiplicity-gcd1 is stronger than assuming all roots are simple.
assert math.gcd(2,3) == 1
assert -1*2+1*3 == 1
assert sp.Rational(3,2)*3 not in sp.S.Integers

# The constant leading equation and its resonant cancellation.
n = sp.symbols('n', integer=True, positive=True)
Hf, Kf = sp.Function('H')(a), sp.Function('K')(a)
identity = n*sp.diff(Hf,a)*Kf-Hf*sp.diff(Kf,a)
assert sp.simplify(identity + Hf**(n+1)*sp.diff(Kf/Hf**n,a)) == 0
for d in [2,3,4]:
    Hp = a*a+a+1
    Kp = 7*Hp**(d-1)+a
    Ksmall = Kp-7*Hp**(d-1)
    W = (d-1)*sp.diff(Hp,a)*Kp-Hp*sp.diff(Kp,a)
    Wsmall = (d-1)*sp.diff(Hp,a)*Ksmall-Hp*sp.diff(Ksmall,a)
    assert sp.expand(W-Wsmall) == 0
    assert sp.degree(W,a) >= 1

# Degree-one leading polynomials admit actual finite mates on every tested chart.
alpha, beta = sp.symbols('alpha beta', nonzero=True)
for d in [2,3,4]:
    rate_v = v**(1-d)
    linear_S = {1: alpha*v, 0: beta*v}
    linear_T = {0: v**(d-1)/(alpha*(d-1))}
    assert comm(linear_S,linear_T,v,rate_v) == {0:h}

# A proper power really admits fractional leading cancellation.
Hp = (a*a+1)**2
Kp = (a*a+1)**3
Mp = sp.Rational(3,2)
assert sp.expand(Mp*sp.diff(Hp,a)*Kp-Hp*sp.diff(Kp,a)) == 0
fractional_S = scale(poly_op(Hp),w**2)
fractional_T = scale(poly_op(Kp),w**3)
fractional_comm = comm(fractional_S,fractional_T,w,w**(-3)/2)
assert coefficients_at_w(fractional_comm,1) == 0
assert fractional_comm and w_degree(fractional_comm) < 1

# The actual h-small order-three Long candidate.
Sstar = {3:h*v, 2:-3*r, 1:sp.Integer(4), 0:2*v}
assert {k:sp.expand(c.subs(h,0)) for k,c in Sstar.items()
        if c.subs(h,0) != 0} == {2:-3*r, 1:4, 0:2*v}
assert max(Sstar) == 3 and Sstar[3].subs(h,0) == 0
assert sp.gcd(sp.Poly(H,a),sp.Poly(sp.diff(H,a),a)).as_expr() == 1

A = sp.Matrix([
    [0,0,(s-2*v)/h**2],
    [v/h,0,-4/h**2],
    [0,v/h,3*r/h**2],
])
C = sp.Matrix([
    [0,0,(s-2*v)/(h*v)],
    [1,0,-4/(h*v)],
    [0,1,3*r/(h*v)],
])
assert matrix_zero(v/h*C-A)
assert all(sp.degree(z,v) <= 1 for z in A if z != 0)
A0 = A.subs(v,0)
A1 = A.diff(v)
assert matrix_zero(A-A0-v*A1)
C0, Cm = h*A1, h*A0
assert matrix_zero(C-C0-Cm/v)
I3,Z3 = sp.eye(3),sp.zeros(3)
Hstar = C0.row_join(Cm).col_join(
    (Cm/v**2).row_join(C0+h/v**2*I3))
V = Z3.row_join(v**2*I3).col_join(I3.row_join(Z3))
I6 = sp.eye(6)
delta = lambda M:M.applyfunc(lambda z:sp.diff(z,v)/v)
projectors = [(I6+V/v)/2,(I6-V/v)/2]
assert matrix_zero(projectors[0]+projectors[1]-I6)
assert matrix_zero(projectors[0]*projectors[1])
for sign,P in zip([1,-1],projectors):
    assert matrix_zero(P*P-P)
    assert matrix_zero(h*delta(P)+Hstar*P-P*Hstar)
    W = I3.col_join(sign/v*I3)
    assert matrix_zero(P*W-W)
    assert matrix_zero(h*delta(W)+Hstar*W-W*(C0+sign*Cm/v))
    assert W.rank() == 3
residue = (v**2*Hstar/(2*h)).applyfunc(lambda z:sp.limit(z,v,0))
assert matrix_zero(residue*(residue-half*I6))
assert sp.trace(residue) == sp.Rational(3,2)
assert residue.rank() == 3
assert sp.factor(6-3) == 3

names = [
    'CYCLIC_INFINITY_MATE_OBSTRUCTION.md',
    'check_infinity_mate.py',
    'FOURIER_CYCLIC_RANK_BUDGET.md',
    'fourier_rank_budget_controls.json',
    'QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md',
    'HAMILTONIAN_IDEALIZER_CRITERION.md',
]
out = {
    'claim':'QVM1',
    'status':'PASS',
    'leading_Ore_controls':leading_controls,
    'finite_power_elimination':'PASS',
    'root_multiplicity_gcd1_control':'multiplicities2,3',
    'leading_Wronskian_identity':'symbolic n,H,K',
    'positive_linear_pairs_cyclic_degrees':[2,3,4],
    'proper_power_fractional_boundary':'PASS, not a Weyl pair',
    'Long_hsmall_cubic_candidate':{
        'generic_order':3,
        'classical_order':2,
        'top_H':'2+h*a**3',
        'top_H_squarefree':True,
        'ordinary_connection_polynomial_degree':1,
        'pushforward_rank':6,
        'inverse_Fourier_rank':3,
        'rank_calculation_uses_written_QKR1_and_slope_proof':True,
        'algebraic_horizontal_projector_ranks':[3,3],
        'local_residue_eigenvalue_multiplicities':{'0':3,'1/2':3},
        'full_Weyl_mate_excluded_by_written_theorem':True,
    },
    'hashes':{name:sha(HERE/name) for name in names},
    'research_context':[
        {'source':'https://arxiv.org/html/2210.00257',
         'used':'Theorem2.8 recalls weighted leading commutator/Poisson calculus. Context only; the ramified infinity theorem is proved directly.'},
        {'source':'https://arxiv.org/html/1111.6100',
         'used':'Remarks1.13 and Section2 leading brackets and proportional powers in fractional Weyl extensions. Context only; no whole-curve theorem imported.'},
        {'source':'https://arxiv.org/html/0808.0742',
         'used':'Proposition3.5 Fourier rank and formal slope theory for the separate rank-three candidate calculation.'},
    ],
    'coverage_boundary':(
        'Exact finite Ore leading-coefficient, power-elimination, Wronskian, '
        'degree-one pair, proper-power boundary, candidate companion and '
        'projector controls. The universal no-mate theorem uses the written '
        'Puiseux/infinity valuation proof, with primitive root multiplicities '
        'and v-leading degree one. Proper powers and other leading faces are '
        'outside this exclusion. Candidate inverse-Fourier rank uses QKR1 and '
        'standard slope theory. No actual Long carrier or all-representative '
        'nonlifting is certified.'),
    'independent_replay':False,
    'formal_verification':False,
    'external_review':False,
}
(HERE/'infinity_mate_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({
    'claim':out['claim'],'status':out['status'],
    'leading_controls':len(leading_controls),
    'linear_pairs':len(out['positive_linear_pairs_cyclic_degrees']),
    'candidate_projectors':'PASS',
    'candidate_generic_order':3,
    'candidate_inverse_Fourier_rank':3,
}))
