#!/usr/bin/env python3
"""Exact finite controls for QPC1. Universal coverage uses the written proof."""
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
v, w, h, lam, mu = sp.symbols('v w h lam mu', nonzero=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Keep the exact hypotheses and prior classification inputs bound to this receipt.
for receipt_name in ['zero_mode_controls.json', 'infinity_limits_controls.json']:
    prior = json.loads((HERE / receipt_name).read_text())
    assert prior['status'] == 'PASS', receipt_name
    for name, digest in prior['hashes'].items():
        assert sha(HERE / name) == digest, name


# Sparse normal order: key (a order, Laurent exponent of v).
def ore_mul(left, right, p, H=1):
    out = {}
    for (i, e), x in left.items():
        for (j, f), y in right.items():
            falling = 1
            for k in range(i + 1):
                value = x * y * comb(i, k) * falling * pow(H, k, p)
                key = (i + j - k, e + f - 2 * k)
                out[key] = (out.get(key, 0) + value) % p
                falling = falling * (f - 2 * k) % p
                if falling == 0:
                    break
    return {key: value for key, value in out.items() if value}


def ore_power(op, exponent, p, H=1):
    result = {(0, 0): 1}
    for _ in range(exponent):
        result = ore_mul(result, op, p, H)
    return result


def binom_poly(x, m, p):
    numerator = 1
    for j in range(m):
        numerator = numerator * (x - j) % p
    return numerator * pow(factorial(m) % p, -1, p) % p


def two_v_sum(n, p):
    m = (p + 1) // 2
    cm = 1
    for j in range(m):
        cm = cm * (1 - 2 * j) % p
    total = 0
    inner = 0
    for u in range(p - 1):
        inner = (inner + binom_poly(n * u, m, p)) % p
        total += inner * binom_poly(n * u - m, m, p)
    return cm * cm * total % p


declared_cases = [
    ('n2_p7', 7, [0, 0, 1], 1, 1),
    ('n2_p11', 11, [0, 0, 1], 1, 1),
    ('n3_p11', 11, [0, 0, 0, 1], 1, 1),
    ('n4_p11', 11, [0, 0, 0, 0, 1], 1, 1),
    ('n5_p13', 13, [0, 0, 0, 0, 0, 1], 1, 1),
    ('nonmonic_lower_terms', 11, [5, 7, 1, 3], 2, 4),
]
power_controls = []
for label, p, G, Mu, H in declared_cases:
    n = len(G) - 1
    assert p > 2 * n + 1
    S = {(j, 0): c % p for j, c in enumerate(G) if c % p}
    S[(0, 1)] = Mu
    result = ore_power(S, p, p, H)
    J = (n - 1) * p - (2 * n + 1)
    gamma = (n * n - 1) * pow(24 * n, -1, p) % p
    scale = pow(G[-1], p - 2, p) * Mu**2 * pow(H, p + 1, p) % p
    expected = gamma * scale % p
    actual = result.get((J, -2 * p), 0)
    assert actual == expected != 0, (label, actual, expected)
    assert J > 0 and J % p
    assert two_v_sum(n, p) == gamma
    # Check the coefficient's noncentrality by actual commutation with v.
    V = {(0, 1): 1}
    left, right = ore_mul(result, V, p, H), ore_mul(V, result, p, H)
    assert left != right
    power_controls.append({
        'case': label, 'p': p, 'G_coefficients_ascending': G,
        'mu': Mu, 'h': H, 'target_v_exponent': -2 * p,
        'target_a_order': J, 'coefficient': actual,
        'predicted_coefficient': expected,
        'two_v_word_sum': two_v_sum(n, p),
        'gamma_n': gamma, 'normal_order_terms': len(result),
        'nonzero_commutator_with_v': True,
    })

# At p=2n+1, the target order can be zero; no obstruction is claimed there.
bad_p = 5
boundary = ore_power({(2, 0): 1, (0, 1): 1}, bad_p, bad_p)
assert all(i % bad_p == 0 and e % bad_p == 0 for i, e in boundary)
assert (2 - 1) * bad_p - (2 * 2 + 1) == 0

# Universal rational identity: top coefficients computed before m=1/2 reduction.
n, m = sp.symbols('n m')
ea1 = m * (m - 1) / 2
ea2 = m * (m - 1) * (m - 2) * (3 * m - 1) / 24
ea3 = m**2 * (m - 1)**2 * (m - 2) * (m - 3) / 48
eb1 = m * (3 * m - 1) / 2
pb2 = m * (2 * m - 1) * (7 * m - 1) / 6
pb3 = m**2 * (15 * m**2 - 14 * m + 3) / 4
eb2 = (eb1**2 - pb2) / 2
eb3 = (eb1**3 - 3 * eb1 * pb2 + 2 * pb3) / 6
A = [sp.Integer(1), -ea1 / n, ea2 / n**2, -ea3 / n**3]
B = [sp.Integer(1), -eb1 / n, eb2 / n**2, -eb3 / n**3]
F = [
    A[0] / (m + 1),
    A[0] / 2 + A[1] / m,
    A[0] * m / 12 + A[1] / 2 + A[2] / (m - 1),
    A[1] * (m - 1) / 12 + A[2] / 2 + A[3] / (m - 2),
]
Ah = [sp.cancel(z.subs(m, sp.Rational(1, 2))) for z in A]
Bh = [sp.cancel(z.subs(m, sp.Rational(1, 2))) for z in B]
Fh = [sp.cancel(z.subs(m, sp.Rational(1, 2))) for z in F]
assert Ah == [1, 1/(8*n), 1/(128*n**2), -5/(1024*n**3)]
assert Bh == [1, -1/(8*n), 1/(128*n**2), 5/(1024*n**3)]
D = sp.cancel(sum(Fh[j] * Bh[3 - j] for j in range(4)))
assert sp.cancel(D - (1 - n**2) / (96 * n**3)) == 0
C = sp.cancel(-4 * n**2 * D)
assert sp.cancel(C - (n**2 - 1) / (24 * n)) == 0


# Exact characteristic-zero Ore controls with rational symbolic coefficients.
def clean(op):
    return {j: z for j, x in op.items() if (z := sp.cancel(x)) != 0}


def add(*ops):
    result = {}
    for op in ops:
        for j, x in op.items():
            result[j] = result.get(j, 0) + x
    return clean(result)


def scale(op, x):
    return clean({j: x * y for j, y in op.items()})


def symbolic_mul(left, right, variable=v, rate=None):
    if rate is None:
        rate = 1/v
    result = {}
    for i, x in left.items():
        for j, y in right.items():
            derivative = y
            for k in range(i + 1):
                key = i + j - k
                result[key] = result.get(key, 0) + comb(i, k) * x * h**k * derivative
                derivative = sp.cancel(rate * sp.diff(derivative, variable))
    return clean(result)


def comm(left, right):
    return add(symbolic_mul(left, right), scale(symbolic_mul(right, left), -1))


Saff = {1: lam, 0: mu*v + 7}
Taff = {2: -lam/(2*mu**2), 1: -v/mu, 0: -h/(2*mu*v)}
assert comm(Saff, Taff) == {0: h}
assert comm({1: lam, 0: sp.Integer(7)}, {0: v**2/(2*lam)}) == {0: h}
assert comm({1: lam*v, 0: mu*v+7}, {0: v/lam}) == {0: h}

domain = sp.QQ.frac_field(v, h)


def ntrace(function, cover_degree):
    numerator, denominator = sp.fraction(sp.cancel(function))
    modulus = sp.Poly(w**cover_degree-v, w, domain=domain)
    num = sp.Poly(numerator, w, domain=domain)
    den = sp.Poly(denominator, w, domain=domain)
    reduced = (num * sp.invert(den, modulus)).rem(modulus)
    return sp.cancel(reduced.nth(0))


trace_controls = []
for cover_degree in [2, 3]:
    rate = sp.Rational(1, cover_degree)*w**(1-2*cover_degree)
    T = {
        2: (w**2+2*w+1)/(w+2),
        1: w/(w+1),
        0: w**2+1,
    }
    traced = clean({j: ntrace(x, cover_degree) for j, x in T.items()})
    for label, base_op, cover_op in [
        ('a', {1: sp.Integer(1)}, {1: sp.Integer(1)}),
        ('v', {0: v}, {0: w**cover_degree}),
    ]:
        lhs = symbolic_mul(cover_op, T, w, rate)
        lhs_trace = clean({j: ntrace(x, cover_degree) for j, x in lhs.items()})
        assert lhs_trace == symbolic_mul(base_op, traced)
        rhs = symbolic_mul(T, cover_op, w, rate)
        rhs_trace = clean({j: ntrace(x, cover_degree) for j, x in rhs.items()})
        assert rhs_trace == symbolic_mul(traced, base_op)
    trace_controls.append({'cover_degree': cover_degree, 'left_and_right_a_v': 'PASS'})

projector = sp.Matrix([[1, w], [1/w, 1]]) / 2
assert (projector*projector-projector).applyfunc(sp.cancel) == sp.zeros(2)
traced_projector = projector.applyfunc(lambda x: ntrace(x, 2))
assert traced_projector == sp.eye(2)/2
assert traced_projector*traced_projector != traced_projector


# Actual split fibers of the relative rank-one Weyl algebra in characteristic7.
p = 7
I, N, derivative = sp.eye(p), sp.zeros(p), sp.zeros(p)
for j in range(p-1):
    N[j+1, j] = 1
    derivative[j, j+1] = j+1


def mod(matrix):
    return matrix.applyfunc(lambda x: int(x) % p)


def is_scalar(matrix):
    return matrix == matrix[0, 0]*sp.eye(matrix.rows)


def charpoly_data(matrix):
    coefficients = [int(x) % p for x in matrix.charpoly().all_coeffs()]
    degree = matrix.rows
    derivative_zero = all(
        ((degree-j)*coefficient) % p == 0
        for j, coefficient in enumerate(coefficients[:-1])
    )
    return coefficients, derivative_zero


V = I+N
H = 4
a = mod(H*V.inv_mod(p)*derivative)
assert mod(a*V-V*a) == mod(H*V.inv_mod(p))
Sbad = mod(a*a+V)
cp_bad, derivative_bad = charpoly_data(Sbad)
assert not derivative_bad
assert not is_scalar(mod(Sbad**p))

La, Mu = 3, 2
Sgood = mod(La*a+Mu*V)
Tgood = mod(
    (-La*pow(2*Mu**2, -1, p))*a*a
    - pow(Mu, -1, p)*V*a
    - H*pow(2*Mu, -1, p)*V.inv_mod(p)
)
assert mod(Sgood*Tgood-Tgood*Sgood) == H*I
cp_good, derivative_good = charpoly_data(Sgood)
assert derivative_good and is_scalar(mod(Sgood**p))

# Genuine QPH1 phase fiber: valid matrices need not have scalar p-th powers.
q, aa = I+N, derivative
g = mod(4*aa-3*aa**2)
c0 = mod(8*q-g*g)
inv2, inv4 = pow(2, -1, p), pow(4, -1, p)
t0 = mod(inv2*aa**2-inv4*aa**3)
Sphase = mod(sp.BlockMatrix([[sp.zeros(p), I], [c0, 2*g]]).as_explicit())
Tphase = mod(sp.BlockMatrix([
    [t0, -inv4*aa],
    [I-inv4*aa*c0, t0-inv2*aa*g],
]).as_explicit())
assert mod(Sphase*Tphase-Tphase*Sphase) == sp.eye(2*p)
cp_phase, derivative_phase = charpoly_data(Sphase)
assert derivative_phase
assert cp_phase == [1]+[0]*(2*p-1)+[p-1]
assert not is_scalar(mod(Sphase**p))

bound_names = [
    'UNIVERSAL_P_POWER_OBSTRUCTION.md', 'check_p_power_obstruction.py',
    'QUANTUM_ZERO_MODE_TRANSPORT.md', 'zero_mode_controls.json',
    'ALGEBRAIC_INFINITY_LIMITS.md', 'infinity_limits_controls.json',
    'QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md', 'phase_carrier_controls.json',
]
for name in bound_names:
    assert (HERE/name).is_file(), name
phase_prior = json.loads((HERE/'phase_carrier_controls.json').read_text())
assert phase_prior['status'] == 'PASS'
for name, digest in phase_prior['hashes'].items():
    assert sha(HERE/name) == digest, name

receipt = {
    'claim': 'QPC1', 'status': 'PASS',
    'full_sparse_Ore_power_controls': power_controls,
    'symbolic_top_coefficient_identity': str(sp.factor(D)),
    'universal_gamma_identity': str(sp.factor(C)),
    'positive_classification_families': [
        'constant_H_affine_G_mu_nonzero',
        'zero_H_affine_G_lambda_nonzero',
        'linear_H_constant_G',
    ],
    'linear_trace_bimodule_controls': trace_controls,
    'trace_of_projector_is_not_projector': True,
    'small_prime_boundary': {
        'p': bad_p, 'n': 2, 'target_a_order': 0,
        'p_power_is_central': True,
    },
    'scalar_split_fibers': {
        'p': p, 'dimension': p, 'h': H,
        'bad_S_charpoly_coefficients': cp_bad,
        'bad_S_charpoly_derivative_zero': derivative_bad,
        'bad_S_p_power_scalar': False,
        'positive_pair_commutator': 'h*I',
        'positive_pair_charpoly_coefficients': cp_good,
        'positive_pair_p_power_scalar': True,
    },
    'genuine_phase_matrix_boundary': {
        'source': 'QPH1 equation6', 'p': p, 'dimension': 2*p,
        'commutator': 'I', 'charpoly_coefficients': cp_phase,
        'charpoly_derivative_zero': derivative_phase,
        'p_power_scalar': False,
    },
    'hashes': {name: sha(HERE/name) for name in bound_names},
    'research_context': [
        {
            'source': 'https://arxiv.org/html/1508.02499',
            'used': 'Sections1.2 and1.4: standard Weyl center and good reduction framework.',
        },
        {
            'source': 'https://arxiv.org/html/2601.23110',
            'used': 'Section3.1: rank-one Weyl algebra splitting background. The p^2 lifting theorem is not used.',
        },
    ],
    'coverage_boundary': (
        'Declared finite sparse powers, two-v word sums, symbolic rational top '
        'coefficients, positive mates, trace generator identities on covers2,3, '
        'projector countercontrol and actual p-by-p/2p-by-2p fibers. Universal '
        'word selection, finite-field summation, specialization, arbitrary '
        'finite covers and the complete separable classification use written '
        'proofs, including QVL1 and QZM1. No arbitrary matrix p-centrality, '
        'all-representative Long exclusion or actual finite carrier is certified.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'p_power_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QPC1', 'status': 'PASS', 'full_power_cases': len(power_controls),
    'universal_gamma_identity': str(sp.factor(C)),
    'matrix_boundary': receipt['genuine_phase_matrix_boundary'],
}, indent=2))
