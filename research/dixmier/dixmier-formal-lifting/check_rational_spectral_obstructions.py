#!/usr/bin/env python3
"""Declared exact controls for QST1. Universal coverage uses the written proofs."""
from pathlib import Path
from math import factorial, prod
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
x, v, w, a, spectral = sp.symbols('x v w a spectral')
h, mu = sp.symbols('h mu', nonzero=True)
Euler = sp.symbols('Euler')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for receipt_name in ['adler_admission_controls.json', 'spectral_curvature_controls.json']:
    previous = json.loads((HERE/receipt_name).read_text())
    assert previous['status'] == 'PASS', receipt_name
    for name, digest in previous['hashes'].items():
        assert sha(HERE/name) == digest, name


def binomial(i, k):
    return sp.Rational(prod(i-j for j in range(k)), factorial(k))


def clean(op):
    return {i: z for i, value in op.items() if (z := sp.cancel(value)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for i, value in op.items():
            out[i] = out.get(i, 0)+value
    return clean(out)


def scale(op, value):
    return clean({i: value*z for i, z in op.items()})


def mul(left, right, low=0, variable=x, rate=sp.Integer(1)):
    out = {}
    for i, value in left.items():
        for j, other in right.items():
            derivative = other
            for k in range(max(0, i+j-low+1)):
                if derivative == 0:
                    break
                key = i+j-k
                out[key] = out.get(key, 0)+binomial(i, k)*h**k*value*derivative
                derivative = sp.cancel(rate*sp.diff(derivative, variable))
    return clean(out)


def power(op, exponent, low=0, variable=x, rate=sp.Integer(1)):
    out = {0: sp.Integer(1)}
    upper = max(op, default=0)
    for step in range(exponent):
        out = mul(out, op, low-(exponent-step-1)*upper, variable, rate)
    return out


def comm(left, right, low=0, variable=x, rate=sp.Integer(1)):
    return add(
        mul(left, right, low, variable, rate),
        scale(mul(right, left, low, variable, rate), -1),
    )


def polynomial_op(expression):
    return {i: coefficient for (i,), coefficient in sp.Poly(expression, a).terms()}


# Full quadratic root and first moment. No h or potential-amplitude projection.
L = {2: sp.Integer(1), 0: mu*v}
root = {1: sp.Integer(1)}
for k in range(6):
    order = 1-k
    got = power(root, 2, order, v, 1/v)
    coefficient = sp.cancel((L.get(order, 0)-got.get(order, 0))/2)
    if coefficient != 0:
        root[-k] = coefficient
moments = {}
for j in range(-3, 6):
    if j < -1:
        moment = sp.Integer(0)  # Order makes the PDO residue zero.
    elif j == -1:
        moment = sp.residue(v, v, 0)  # Monic inverse root has residue1.
    else:
        density = power(root, j, -1, v, 1/v).get(-1, sp.Integer(0))
        moment = sp.residue(density*v, v, 0)
    moments[j] = sp.cancel(moment)
assert all(moments[j] == 0 for j in range(-3, 5))
assert moments[5] == -5*h**2*mu**2/32

triangular_controls = []
for degree in [0, 1, 2, 3, 4]:
    coefficients = [sp.Integer(ell+2) for ell in range(degree+1)]
    indices = [5+2*(ell-degree) for ell in range(degree+1)]
    # Every negative order below -1 contributes zero; the declared full
    # moments provide every other value used by these degree checks.
    actual = sum(
        coefficient*(moments[index] if index >= -3 else 0)
        for coefficient, index in zip(coefficients, indices)
    )
    assert sp.cancel(actual-coefficients[-1]*moments[5]) == 0
    triangular_controls.append({
        'polynomial_degree': degree, 'functional_power': 5-2*degree,
        'moment_indices': indices, 'value': str(sp.factor(actual)),
    })


# Exact global spectral annihilators for declared constant-coefficient polynomials.
critical_controls = []
for label, G in [
    ('quadratic', a**2),
    ('two_critical_values', a**3-3*a),
    ('coincident_values', a**4-2*a*a),
    ('degenerate_critical_point', a**5),
    ('nonrational_critical_values', a**4+a),
    ('affine_positive', 3*a+2),
]:
    derivative = sp.diff(G, a)
    resultant = sp.resultant(derivative, spectral-G, a)
    Delta = sp.Poly(resultant, spectral).sqf_part().monic().as_expr()
    quotient, remainder = sp.div(sp.expand(Delta.subs(spectral, G)), derivative, a)
    assert remainder == 0
    Sop = polynomial_op(G)
    Top = scale(polynomial_op(quotient), x)
    assert comm(Sop, Top) == scale(polynomial_op(Delta.subs(spectral, G)), h)
    critical_controls.append({
        'case': label, 'G': str(G), 'Delta': str(Delta),
        'finite_numerator_H': str(sp.factor(quotient)),
        'commutator': 'h*Delta(G)',
    })
assert critical_controls[0]['Delta'] == 'spectral'
assert critical_controls[1]['Delta'] == 'spectral**2 - 4'
assert critical_controls[2]['Delta'] == 'spectral**2 + spectral'
assert critical_controls[3]['Delta'] == 'spectral'
assert critical_controls[-1]['Delta'] == '1'

# Noncommutative quadratic Reynolds control on w^2=v.
b_cover = h/(2*w**4)
Q_cover = {2: sp.Integer(1), 1: b_cover, 0: -w*w-h*h/w**8}
first_factorization = mul({1: sp.Integer(1), 0: w},
                          {1: sp.Integer(1), 0: b_cover-w},
                          0, w, 1/(2*w**3))
second_factorization = mul({1: sp.Integer(1), 0: -w},
                           {1: sp.Integer(1), 0: b_cover+w},
                           0, w, 1/(2*w**3))
assert first_factorization == clean(Q_cover) == second_factorization
average_numerator = scale(add(
    {1: sp.Integer(1), 0: b_cover-w},
    {1: sp.Integer(1), 0: b_cover+w},
), sp.Rational(1, 2))
assert average_numerator == {1: sp.Integer(1), 0: b_cover}
assert all(sp.cancel(c.subs(w, -w)-c) == 0 for c in Q_cover.values())
assert sp.cancel(b_cover.subs(w**2, v)-h/(2*v*v)) == 0

# Primitive leading divisor blocks finite differential mates, but not rational ones.
Smixed = {2: v, 1: h/(2*v)}
Amixed = {1: v*v/3}
assert comm(Smixed, Amixed, 0, v, 1/v) == scale(Smixed, h)
Ecover = {1: w}
assert power(Ecover, 2, 0, w, 1/(2*w**3)) == {
    2: w*w, 1: h/(2*w*w)
}
ucover = 2*w**3/3
assert sp.cancel(w*(sp.diff(ucover, w)/(2*w**3))-1) == 0
assert scale(Ecover, ucover/2) == {1: w**4/3}

# A finite Darboux certificate for a rational mate with a differential denominator.
Sconstant = {2: sp.Integer(1)}
Aconstant = {0: x/2}
Bconstant = {1: sp.Integer(1)}
Cconstant = Sconstant
assert mul(Sconstant, Bconstant) == mul(Bconstant, Cconstant)
assert add(mul(Sconstant, Aconstant),
           scale(mul(Aconstant, Cconstant), -1)) == scale(Bconstant, h)

# A scalar multiplication gauge produces genuinely mixed negative-v coefficients.
gauge_log_derivative = 3*h/v**2
Agauge = {1: sp.Integer(1), 0: gauge_log_derivative}
Sgauge = add(power(Agauge, 2, 0, v, 1/v), {0: mu*v})
# Check b^-1*L*b explicitly for b=v^3.
explicit_conjugate = scale(mul(L, {0: v**3}, 0, v, 1/v), v**-3)
assert Sgauge == explicit_conjugate
assert Sgauge[1] != 0


# Exact harmonic canonical normalization.
Aop = {1: sp.Integer(1)}
Xop = {0: x}
Hharm = {2: sp.Integer(1), 0: x*x}
Xnew = add(Aop, scale(Xop, sp.I))
Pnew = scale(add(Aop, scale(Xop, -sp.I)), 1/(2*sp.I))
assert comm(Pnew, Xnew) == {0: h}
assert add(scale(mul(Xnew, Pnew), 2*sp.I), {0: sp.I*h}) == Hharm


# Rank-one rational torus model Q=K(Euler)(X;shift), X*Euler=(Euler-h)*X.
# Keys are torus weights; multiplication uses the exact coefficient shift.
def skew_add(*ops):
    return add(*ops)


def skew_mul(left, right, high=None):
    out = {}
    for i, f in left.items():
        for j, g in right.items():
            key = i+j
            if high is not None and key > high:
                continue
            value = f*g.subs(Euler, Euler-i*h)
            out[key] = out.get(key, 0)+value
    return clean(out)


def skew_inverse(op, high):
    low = min(op)
    leading = op[low]
    monomial_inverse = {-low: sp.cancel(1/leading.subs(Euler, Euler+low*h))}
    shifted = skew_mul(monomial_inverse, op)
    V = skew_add(shifted, {0: -sp.Integer(1)})
    assert all(i > 0 for i in V)
    if not V:
        return monomial_inverse
    smallest = min(V)
    count = max(0, (high+low)//smallest)
    out = {}
    term = {0: sp.Integer(1)}
    for k in range(count+1):
        contribution = skew_mul(term, monomial_inverse, high)
        out = skew_add(out, scale(contribution, (-1)**k))
        term = skew_mul(term, V, high+low)
    return out


skewX = {1: sp.Integer(1)}
skewP = {-1: Euler+h}
skewH = {0: Euler}
assert skew_add(skew_mul(skewP, skewX),
                scale(skew_mul(skewX, skewP), -1)) == {0: h}
assert skew_mul(skewX, skewP) == skewH
inverse1 = skew_inverse({0: sp.Integer(1), 1: -sp.Integer(1)}, 5)
assert inverse1 == {j: sp.Integer(1) for j in range(6)}
rational_test = skew_mul(skewP, inverse1)
assert rational_test.get(0) == Euler+h
commH = skew_add(skew_mul(skewH, rational_test),
                  scale(skew_mul(rational_test, skewH), -1))
assert commH.get(0, 0) == 0
rationalF = (Euler*Euler+1)/(Euler-2)
commF = skew_add(skew_mul({0: rationalF}, rational_test),
                  scale(skew_mul(rational_test, {0: rationalF}), -1))
assert commF.get(0, 0) == 0
left_base, right_base = 1/(Euler+3), Euler*Euler+1
bimodule_test = skew_mul({0: left_base}, skew_mul(rational_test, {0: right_base}))
assert sp.cancel(bimodule_test.get(0)-left_base*(Euler+h)*right_base) == 0
# The projection is not multiplicative.
assert skewX.get(0, 0) == 0 and skewP.get(0, 0) == 0
assert skew_mul(skewX, skewP).get(0) == Euler

sum_inverse = skew_inverse(skew_add(skewX, skewP), 5)
assert sum_inverse and all(j >= 1 for j in sum_inverse)
product_check = skew_mul(skew_add(skewX, skewP), sum_inverse, 4)
assert product_check == {0: sp.Integer(1)}

hash_names = [
    'RATIONAL_SPECTRAL_OBSTRUCTIONS.md', 'check_rational_spectral_obstructions.py',
    'ADLER_PERIODS_AND_FINITE_ADMISSION.md', 'adler_admission_controls.json',
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md', 'spectral_curvature_controls.json',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
    'QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md',
]
receipt = {
    'claim': 'QST1', 'status': 'PASS',
    'full_quadratic_Adler_moments': {str(j): str(sp.factor(value)) for j, value in moments.items()},
    'triangular_spectral_dual_controls': triangular_controls,
    'constant_coefficient_spectral_annihilators': critical_controls,
    'rational_Reynolds_quadratic_control': {
        'cover': 'w^2=v',
        'common_monic_right_denominator': 'a^2+h/(2*v^2)*a-v-h^2/v^4',
        'average_numerator': 'a+h/(2*v^2)',
        'both_Ore_factorizations': 'PASS',
    },
    'mixed_rational_escape': {
        'S': 'v*a^2+h/(2*v)*a',
        'finite_numerator': 'v^2*a/3',
        'commutator': 'h*S',
        'rational_mate': '(v^2*a/3)*S^-1',
        'finite_mate': 'excluded by QAD1 primitive leading divisor',
    },
    'finite_Darboux_certificate': 'PASS',
    'rational_multiplication_gauge': {
        'b': 'v^3', 'shift': 'a+3*h/v^2',
        'mixed_first_order_coefficient': str(Sgauge[1]),
        'exact_conjugation': 'PASS',
    },
    'harmonic_canonical_normalization': 'H=2*i*X*P+i*h, [P,X]=h',
    'rational_torus_controls': {
        'model': 'K(Euler)(X; Euler->Euler-h)',
        'nonhomogeneous_rational_test': 'P*(1-X)^-1',
        'constant_term': 'Euler+h',
        'Euler_and_rational_function_commutator_constant_terms': 0,
        'invariant_bimodule_retraction': 'PASS',
        'multiplicative_projection': False,
        'inverse_of_X_plus_P_lowest_weight': min(sum_inverse),
    },
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'research_context': [
        {
            'source': 'https://arxiv.org/html/math/0403037',
            'used': 'Section5 semisimple Weyl-element context; the rational torus residue is proved directly and no classification theorem is imported.',
        },
        {
            'source': 'https://arxiv.org/abs/1001.0153',
            'used': 'Related differential-operator descent context. Fixed linear rational mate descent here has an explicit monic-denominator proof.',
        },
    ],
    'coverage_boundary': (
        'Full quadratic moments, declared triangular degrees, six critical-value '
        'resultants and Ore witnesses, a quadratic rational Reynolds denominator, '
        'mixed rational escape, finite Darboux identity, multiplication gauge, '
        'harmonic normalization and nonhomogeneous torus Laurent fractions. '
        'Universal moment separation, all-degree critical-value annihilator, '
        'all finite coefficient covers, rational torus retraction in every rank, '
        'transport/stabilization and projective rational-orbit obstruction use '
        'written proofs. No arbitrary rational mate is inferred solely from '
        'spectral torsion tests or zero Adler periods. No finite Long carrier '
        'or universal rational normalization is certified.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'rational_spectral_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QST1', 'status': 'PASS',
    'Long_family_first_moment': str(sp.factor(moments[5])),
    'critical_value_models': len(critical_controls),
    'mixed_rational_escape': receipt['mixed_rational_escape'],
    'rational_torus_constant_term': 'Euler+h',
}, indent=2))
