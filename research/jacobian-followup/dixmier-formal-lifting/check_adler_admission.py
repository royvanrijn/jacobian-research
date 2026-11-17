#!/usr/bin/env python3
"""Declared exact controls for QAD1; all-degree coverage is in the written proof."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
v, u, w, h, mu = sp.symbols('v u w h mu')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


prior = json.loads((HERE/'p_power_controls.json').read_text())
assert prior['status'] == 'PASS'
for name, digest in prior['hashes'].items():
    assert sha(HERE/name) == digest, name


def binomial_integer(i, k):
    return Q(prod(i-j for j in range(k)), prod(j+1 for j in range(k)))


# Sparse PDO product in the exact projection modulo (h^3,mu^3).
# Key: a order, v exponent, h degree, potential-amplitude degree.
def sparse_add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, Q(0)) + value
    return {key: value for key, value in out.items() if value}


def sparse_mul(left, right, low):
    out = {}
    for (i, e, hd, md), x in left.items():
        for (j, f, Hd, Md), y in right.items():
            if hd+Hd > 2 or md+Md > 2:
                continue
            dmax = min(i+j-low, 2-hd-Hd)
            falling = 1
            for d in range(dmax+1):
                key = (i+j-d, e+f-2*d, hd+Hd+d, md+Md)
                out[key] = (
                    out.get(key, Q(0))
                    + x*y*binomial_integer(i, d)*falling
                )
                falling *= f-2*d
    return {key: value for key, value in out.items() if value}


def sparse_power(op, exponent, low):
    result = {(0, 0, 0, 0): Q(1)}
    for step in range(exponent):
        result = sparse_mul(result, op, low-(exponent-step-1))
    return result


def sparse_root(n, depth, lower=None):
    L = {(n, 0, 0, 0): Q(1), (0, 1, 0, 1): Q(1)}
    if lower:
        L = sparse_add(L, {(j, 0, 0, 0): Q(x) for j, x in lower.items()})
    R = {(1, 0, 0, 0): Q(1)}
    for k in range(depth+1):
        target = n-1-k
        got = sparse_power(R, n, target)
        difference = sparse_add(L, {key: -x for key, x in got.items()})
        correction = {
            (-k, e, hd, md): x/n
            for (i, e, hd, md), x in difference.items() if i == target
        }
        R = sparse_add(R, correction)
    return R


sparse_controls = []
for n, lower in [(2, None), (3, None), (4, None), (5, None),
                 (3, {2: 2, 1: 3, 0: 5})]:
    j = 2*n+1
    root = sparse_root(n, j, lower)
    result = sparse_power(root, j, -1)
    actual = result.get((-1, -2, 2, 2), Q(0))
    expected = -Q((2*n+1)*(n*n-1), 24*n*n)
    assert actual == expected
    sparse_controls.append({
        'n': n, 'lower_G_coefficients': lower or {},
        'root_power': j, 'projection': 'modulo(h^3,mu^3)',
        'coefficient_of_v^-2_h^2_mu^2': str(actual),
        'predicted_coefficient': str(expected),
    })

# Rational all-n identity, separate from the declared root computations.
n = sp.symbols('n', positive=True, integer=True)
q = 1+1/n
beta = 2+1/n
A = sp.cancel(q*(q-1)*n*(n-1)/2 + q*(q-1)*(q-2)*n**2/3)
assert sp.cancel(A-(n**2-1)/(6*n)) == 0
gamma = (n**2-1)/(24*n)
period_coefficient = sp.cancel(-beta*A/4)
assert sp.cancel(period_coefficient+beta*gamma) == 0


# Rational symbolic PDO calculus; each call explicitly specifies the needed order.
def clean(op):
    return {i: value for i, x in op.items() if (value := sp.cancel(x)) != 0}


def add(*ops):
    result = {}
    for op in ops:
        for i, x in op.items():
            result[i] = result.get(i, 0)+x
    return clean(result)


def scale(op, value):
    return clean({i: value*x for i, x in op.items()})


def mul(left, right, low, variable=u, rate=sp.Integer(1)):
    result = {}
    for i, x in left.items():
        for j, y in right.items():
            derivative = y
            for k in range(max(0, i+j-low+1)):
                if derivative == 0:
                    break
                key = i+j-k
                bc = sp.Rational(binomial_integer(i, k))
                result[key] = result.get(key, 0)+bc*h**k*x*derivative
                derivative = sp.cancel(rate*sp.diff(derivative, variable))
    return clean(result)


def power(op, exponent, low, variable=u, rate=sp.Integer(1)):
    result = {0: sp.Integer(1)}
    upper = max(op, default=0)
    for step in range(exponent):
        result = mul(result, op, low-(exponent-step-1)*upper, variable, rate)
    return result


def root(L, n, depth, variable=u, rate=sp.Integer(1)):
    R = {1: sp.Integer(1)}
    for k in range(depth+1):
        target = n-1-k
        got = power(R, n, target, variable, rate)
        value = sp.cancel((L.get(target, 0)-got.get(target, 0))/n)
        if value != 0:
            R[-k] = value
    return R


def comm(left, right, low=0, variable=u, rate=sp.Integer(1)):
    return add(
        mul(left, right, low, variable, rate),
        scale(mul(right, left, low, variable, rate), -1),
    )


def positive(op):
    return {i: x for i, x in op.items() if i >= 0}


# Full quadratic computation: no amplitude or h-degree projection.
L = {2: sp.Integer(1), 0: mu*v}
R = root(L, 2, 5, v, 1/v)
density = power(R, 5, -1, v, 1/v)[-1]
expected_density = 5*mu*(-3*h**2+2*mu*v**5)*(h**2+mu*v**5)/(32*v**7)
assert sp.cancel(density-expected_density) == 0
assert sp.residue(density*v, v, 0) == -5*h**2*mu**2/32

# The explicit all-order monomial commutator primitive, on declared +/- orders.
monomial_controls = []
f, g = v**3+v**-2, v**2+v**-1
for i, j in [(3, -2), (-1, 2), (2, -4), (0, -1), (4, -1)]:
    N = i+j+1
    residue = comm({i: f}, {j: g}, -1, v, 1/v).get(-1, 0)
    if N >= 1:
        derivatives_f, derivatives_g = [f], [g]
        for _ in range(N):
            derivatives_f.append(sp.cancel(sp.diff(derivatives_f[-1], v)/v))
            derivatives_g.append(sp.cancel(sp.diff(derivatives_g[-1], v)/v))
        primitive = (
            h**N*sp.Rational(binomial_integer(i, N))
            * sum((-1)**k*derivatives_f[k]*derivatives_g[N-1-k]
                  for k in range(N))
        )
        assert sp.cancel(residue-sp.diff(primitive, v)/v) == 0
    else:
        assert residue == 0
    assert sp.residue(residue*v, v, 0) == 0
    monomial_controls.append({'left_order': i, 'right_order': j, 'status': 'PASS'})

# Conjugation trace: invert a genuine order-zero series to the needed precision.
W = {0: sp.Integer(1), -1: v}
U = {-1: v}
Winv = {0: sp.Integer(1)}
term = {0: sp.Integer(1)}
for k in range(1, 4):
    term = mul(term, U, -3, v, 1/v)
    Winv = add(Winv, scale(term, (-1)**k))
P = {2: sp.Integer(1), -1: v**-2}
conjugated = mul(mul(W, P, -1, v, 1/v), Winv, -1, v, 1/v)
assert sp.residue(conjugated.get(-1, 0)*v, v, 0) == 1
assert sp.residue(P[-1]*v, v, 0) == 1

# Finite-cover exact-form trace and the persistence of the local logarithmic class.
domain = sp.QQ.frac_field(v, h)


def ntrace(function, degree):
    numerator, denominator = sp.fraction(sp.cancel(function))
    modulus = sp.Poly(w**degree-v, w, domain=domain)
    numerator = sp.Poly(numerator, w, domain=domain)
    denominator = sp.Poly(denominator, w, domain=domain)
    return sp.cancel((numerator*sp.invert(denominator, modulus)).rem(modulus).nth(0))


cover_controls = []
for degree in [2, 3]:
    function = (w**2+2*w+1)/(w+2)
    delta_cover = sp.Rational(1, degree)*w**(1-2*degree)
    assert sp.cancel(
        ntrace(delta_cover*sp.diff(function, w), degree)
        - sp.diff(ntrace(function, degree), v)/v
    ) == 0
    pulled_log_form = degree/w
    assert sp.residue(pulled_log_form, w, 0) == degree
    cover_controls.append({'degree': degree, 'pulled_dv_over_v_residue': degree})


# A genuinely mixed monic cubic and its all-order basis characterization.
S = {3: sp.Integer(1), 1: u*u+1, 0: u}
R = root(S, 3, 4)
basis = {j: positive(power(R, j, 0)) for j in range(6)}
defects = {j: comm(S, basis[j]) for j in basis}
assert all(max(C, default=-1) <= 1 for C in defects.values())
T = add(scale(basis[5], 2), scale(basis[2], -3),
        scale(basis[1], 7), basis[0])
remainder = dict(T)
coefficients = {}
while remainder:
    order = max(remainder)
    coefficient = remainder[order]
    assert sp.diff(coefficient, u) == 0
    coefficients[order] = coefficient
    remainder = add(remainder, scale(basis[order], -coefficient))
assert coefficients == {5: 2, 2: -3, 1: 7, 0: 1}
assert comm(S, T) == add(*[scale(defects[j], c) for j, c in coefficients.items()])

# Positive all-n mixed family in the original v,a chart: E=v^3*a, u=-1/v.
E = {1: v**3}
Upos = {0: -1/v}
positive_controls = []
for order in [2, 3, 4]:
    Sp = add(power(E, order, 0, v, 1/v), Upos)
    Tp = scale(E, -1)
    assert comm(Sp, Tp, 0, v, 1/v) == {0: h}
    leading = Sp[order]
    assert sp.cancel(leading-v**(3*order)) == 0
    assert sp.cancel(sp.diff(-1/v, v)-v**-2) == 0
    positive_controls.append({'order': order, 'leading_coefficient': str(leading)})

# Harmonic oscillator: all root coefficients lie in K[h,u] by recursion.
Harm = {2: sp.Integer(1), 0: u*u}
Rh = root(Harm, 2, 5)
harmonic_densities = []
for j in [1, 2, 3, 4, 5]:
    res = power(Rh, j, -1).get(-1, sp.Integer(0))
    sp.Poly(res, u, h)  # Polynomiality, hence an exact rational differential.
    primitive = sp.integrate(res, u)
    assert sp.cancel(sp.diff(primitive, u)-res) == 0
    harmonic_densities.append({'root_power': j, 'residue_density': str(sp.factor(res))})
b1, c0 = sp.symbols('b1 c0')
first_order_defect = comm(Harm, {1: b1, 0: c0})
assert first_order_defect == {0: -2*h*u*b1}

# Construct the first three dressing coefficients by the intertwining equation.
Wa = {0: sp.Integer(1)}
Aop = {1: sp.Integer(1)}
for m in range(1, 4):
    error = add(mul(Rh, Wa, -m), scale(mul(Wa, Aop, -m), -1))
    correction = sp.integrate(-error.get(-m, 0)/h, u)
    sp.Poly(correction, u, domain=sp.QQ.frac_field(h))
    Wa[-m] = sp.cancel(correction)
    check = add(mul(Rh, Wa, -m), scale(mul(Wa, Aop, -m), -1))
    assert all(i < -m for i in check)
Ua = {i: x for i, x in Wa.items() if i < 0}
inverse = {0: sp.Integer(1)}
term = {0: sp.Integer(1)}
for k in range(1, 4):
    term = mul(term, Ua, -3)
    inverse = add(inverse, scale(term, (-1)**k))
Tpdo = scale(mul(mul(Wa, {-1: u}, -4), inverse, -4), sp.Rational(1, 2))
pdo_check = comm(Harm, Tpdo, -2)
assert pdo_check == {0: h}, pdo_check

# Universal leading-differential gate, independent of every lower coefficient.
Alead, Blead, Clead = sp.symbols('Alead Blead Clead', nonzero=True)
Hlinear = Alead*v+Blead
Ulinear = n*v/(Alead*(2*n-1))-n**2*Blead/(Alead**2*(n-1)*(2*n-1))
assert sp.cancel(Hlinear*sp.diff(Ulinear, v)+(n-1)/n*sp.diff(Hlinear, v)*Ulinear-v) == 0
Hcentered = Alead*v*v+Clead
Ucentered = n/(2*Alead*(n-1))
assert sp.cancel(Hcentered*sp.diff(Ucentered, v)+(n-1)/n*sp.diff(Hcentered, v)*Ucentered-v) == 0
Hshifted = Alead*v*v+Blead*v+Clead
shifted_defect = sp.cancel(
    Hshifted*sp.diff(Ucentered, v)+(n-1)/n*sp.diff(Hshifted, v)*Ucentered-v
)
assert shifted_defect == Blead/(2*Alead)
degree_U = sp.symbols('degree_U', integer=True, nonnegative=True)
degree_H = sp.symbols('degree_H', integer=True, positive=True)
n_minus2 = sp.symbols('n_minus2', integer=True, nonnegative=True)
highest_multiplier = degree_U+(n_minus2+1)*degree_H/(n_minus2+2)
assert highest_multiplier.is_positive is True
pole_order = sp.symbols('pole_order', integer=True, positive=True)
assert (-pole_order-1/n).is_negative is True
# The nonsquarefree high-degree leading coefficient of the positive family
# passes the rational equation with c=v^3, u=-1/v, Q=c*u=-v^2.
Qpositive = -v*v
assert sp.cancel(sp.diff(Qpositive, v)-3*Qpositive/v-v) == 0
leading_controls = {
    'all_n_linear_H_primitive_U': str(Ulinear),
    'all_n_centered_quadratic_H_primitive_U': str(Ucentered),
    'noncentered_quadratic_incompatibility': str(shifted_defect),
    'universal_highest_multiplier_positive': True,
    'universal_negative_pole_multiplier_nonzero': True,
    'squarefree_degree_at_least3_quadratic_chart': 'written pole/divisibility/degree exclusion',
    'general_delta_v^(1-d)': 'written bound degree_H<=d for squarefree polynomial H',
    'nonsquarefree_positive_Q': str(Qpositive),
}

# A finite differential mate has an additional divisor-order obstruction.
from math import gcd
kk = sp.symbols('kk', positive=True, integer=True)
Hfun = sp.Function('Hfun')(v)
leading_relation = (
    n*Hfun*sp.diff(Hfun**kk, v)
    - n*kk*Hfun**kk*sp.diff(Hfun, v)
)
assert sp.simplify(leading_relation) == 0
Smixed = {3: v**3-v, 1: v*v+1, 0: v**-2}
low_remainder = {2: v+1, 0: v**-3}
large_partner = add(scale(power(Smixed, 2, 0, v, 1/v), 2), low_remainder)
assert max(large_partner) == 6
reduced_partner = add(large_partner, scale(power(Smixed, 2, 0, v, 1/v), -2))
assert reduced_partner == low_remainder
assert comm(Smixed, large_partner, 0, v, 1/v) == comm(Smixed, reduced_partner, 0, v, 1/v)
assert max(comm(Smixed, reduced_partner, 0, v, 1/v)) > 0
assert gcd(3, 1, 1, 1, 3) == 1
assert gcd(4, 2, -3, 1) == 1
assert gcd(6, 2, 4, -6) == 2
divisor_controls = {
    'universal_power_subtraction_leading_relation': 'PASS',
    'mixed_order3_squarefree_leading_coefficient': 'v^3-v',
    'declared_partner_order_reduction': [6, 2],
    'commutator_preserved_under_power_subtraction': True,
    'coprime_zero_pole_order_example': {'operator_order': 4, 'orders': [2, -3, 1], 'gcd': 1},
    'nonprimitive_divisor_boundary': {'operator_order': 6, 'orders': [2, 4, -6], 'gcd': 2},
    'every_nonconstant_squarefree_polynomial_leading_coefficient': 'written all-order finite-mate exclusion',
}

hash_names = [
    'ADLER_PERIODS_AND_FINITE_ADMISSION.md', 'check_adler_admission.py',
    'UNIVERSAL_P_POWER_OBSTRUCTION.md', 'p_power_controls.json',
    'QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md', 'QUANTUM_ZERO_MODE_TRANSPORT.md',
]
receipt = {
    'claim': 'QAD1', 'status': 'PASS',
    'graded_fractional_root_controls': sparse_controls,
    'universal_resolvent_coefficient': str(sp.factor(A)),
    'universal_period_coefficient': str(sp.factor(period_coefficient)),
    'full_quadratic_density': str(sp.factor(density)),
    'full_quadratic_period': str(sp.factor(sp.residue(density*v, v, 0))),
    'monomial_commutator_primitive_controls': monomial_controls,
    'pseudodifferential_conjugation_nonzero_period_preserved': True,
    'finite_cover_controls': cover_controls,
    'leading_differential_gate': leading_controls,
    'finite_divisor_gate': divisor_controls,
    'mixed_cubic_basis_reconstruction': {
        'basis_indices': list(basis),
        'max_defect_order': max(max(C, default=-1) for C in defects.values()),
        'reconstructed_constants': {str(j): str(c) for j, c in coefficients.items()},
    },
    'positive_mixed_families': positive_controls,
    'harmonic_period_insufficiency': {
        'declared_densities': harmonic_densities,
        'all_positive_classes_vanish': 'written polynomial-root recursion proof',
        'no_finite_mate': 'QSO2 plus written first-order contradiction',
        'first_three_dressing_coefficients': {str(j): str(x) for j, x in Wa.items() if j < 0},
        'PDO_mate_commutator_through_order_minus2': 'h',
        'PDO_mate_leading_order': max(Tpdo),
    },
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'research_context': [
        {
            'source': 'https://arxiv.org/html/2405.05421',
            'used': 'Established almost-commuting basis theorem3.5; the direct proof and finite-mate span criterion are written in the note.',
        },
        {
            'source': 'https://eprints.whiterose.ac.uk/178183/7/rnab288.pdf',
            'used': 'Lemma1.1 Adler commutator trace; the monomial primitive is proved and checked directly.',
        },
    ],
    'coverage_boundary': (
        'Declared sparse fractional-root coefficients modulo(h^3,mu^3), '
        'a full quadratic root density, rational all-n coefficient identity, '
        'leading primitive, divisor and pole/degree identities, monomial trace primitives, finite-cover trace, conjugation, a mixed '
        'cubic basis reconstruction, three positive mixed mates and initial '
        'harmonic dressing/PDO coefficients. All-order trace necessity, '
        'PDO sufficiency, constant-linear finite admission, universal '
        'two-potential selection and harmonic hierarchy polynomiality use '
        'written proofs. No finite prefix decides the infinite defect span. '
        'No arbitrary Fourier completion comparison, QRF2 identification, '
        'general matrix criterion or actual finite Long carrier is certified.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'adler_admission_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QAD1', 'status': 'PASS',
    'graded_root_cases': len(sparse_controls),
    'universal_period_coefficient': str(sp.factor(period_coefficient)),
    'full_quadratic_period': str(sp.factor(sp.residue(density*v, v, 0))),
    'harmonic_PDO_mate_leading_order': max(Tpdo),
    'finite_differential_admission': 'constant-linear span of all basis defects',
}, indent=2))
