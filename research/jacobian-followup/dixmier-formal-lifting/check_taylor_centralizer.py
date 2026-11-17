#!/usr/bin/env python3
"""Declared QTB1 controls; no finite prefix proves the universal Taylor theorem."""
from pathlib import Path
from math import factorial, comb
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h = sp.symbols('h', nonzero=True)
x, y, s, t, lam = sp.symbols('x y s t lam')
u, v, z, w, root = sp.symbols('u v z w root')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['etale_centralizer_controls.json', 'finite_centralizer_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, digest in prior['hashes'].items():
        assert sha(HERE/rel) == digest, rel


def truncate(poly, order):
    return sp.expand(sum(coefficient*s**i*t**j
                         for (i, j), coefficient in sp.Poly(sp.expand(poly), s, t).terms()
                         if i+j <= order))


def taylor(expr, order):
    return sp.expand(sum(
        sp.diff(expr, x, i, y, j).subs({x: 1, y: 2})*s**i*t**j/(factorial(i)*factorial(j))
        for i in range(order+1) for j in range(order+1-i)
    ))


f = x/(x*x+y)
g = (x+2*y)/(1+x*y)
order = 3
tf, tg = taylor(f, order), taylor(g, order)
assert truncate(taylor(f*g, order)-tf*tg, order) == 0
translated = sp.series(f.subs({x: 1+lam*s, y: 2+lam*t}), lam, 0, order+1).removeO()
assert sp.expand(translated.subs(lam, 1)-tf) == 0
assert truncate(sp.diff(tf, s)-taylor(sp.diff(f, x), order-1), order-1) == 0
assert truncate(sp.diff(tf, t)-taylor(sp.diff(f, y), order-1), order-1) == 0
assert taylor(x, order) == 1+s and taylor(y, order) == 2+t


# Coprime multivariate numerator/denominator can vanish at the same point.
q = x*x+y
b = x/q
assert sp.cancel(sp.diff(b, x)-2*x*sp.diff(b, y)-1/q) == 0
assert x.subs({x: 0, y: 0}) == 0 and q.subs({x: 0, y: 0}) == 0
assert sp.gcd(x, q) == 1


# Singular divisor: distinct formal factors, although total order is two.
node = y*y-x*x-x**3
assert sp.gcd(node, sp.diff(node, y)) == 1
sqrt_prefix = sum(sp.binomial(sp.Rational(1, 2), j)*x**j for j in range(9))
plus, minus = y+x*sqrt_prefix, y-x*sqrt_prefix
assert sp.series(plus*minus-node, x, 0, 11).removeO() == 0
assert sp.expand(plus-minus).coeff(x, 1) == 2
assert min(sum(alpha) for alpha in sp.Poly(node, x, y).monoms()) == 2
assert sp.expand(plus.subs(y, x*sqrt_prefix)).coeff(x, 1) == 2
assert sp.expand(minus.subs(y, -x*sqrt_prefix)).coeff(x, 1) == -2


# A non-Galois quartic field with a quadratic Galois factor.
quartic = z**4+u*z*z+v
quadratic = w*w+u*w+v
assert sp.expand(quadratic.subs(w, z*z)-quartic) == 0
assert sp.discriminant(quadratic, w) == u*u-4*v
assert sp.Poly(u*u-4*v, v).degree() == 1
assert sp.Poly(-u-z*z, u).degree() == 1
assert sp.expand(quartic.subs(v, -z**4-u*z*z)) == 0


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    return tuple(a.index(i) for i in range(len(a)))


def generated(gens, size):
    out = {tuple(range(size))}
    while True:
        new = out | {compose(g, p) for g in gens for p in out}
        if new == out:
            return out
        out = new


def normal_closure(stabilizer, group, size):
    conjugates = {compose(compose(g, p), inverse(g)) for g in group for p in stabilizer}
    return generated(conjugates, size)


d4 = generated([(1, 2, 3, 0), (0, 3, 2, 1)], 4)
d4_stab = {p for p in d4 if p[0] == 0}
d4_normal = normal_closure(d4_stab, d4, 4)
assert len(d4) == 8 and len(d4_stab) == 2 and len(d4_normal) == 4
s3 = generated([(1, 2, 0), (1, 0, 2)], 3)
s3_stab = {p for p in s3 if p[0] == 0}
s3_normal = normal_closure(s3_stab, s3, 3)
assert len(s3) == 6 and len(s3_stab) == 2 and len(s3_normal) == 6
c3 = generated([(1, 2, 0)], 3)
c3_stab = {p for p in c3 if p[0] == 0}
assert len(normal_closure(c3_stab, c3, 3)) == 1


# Rank-one actual Weyl multiplication with polynomial coefficients in x,h.
# These operators commute with an independent second canonical position.
ONE = {0: sp.Integer(1)}
ZERO = {}


def clean(op):
    return {key: value for key, old in op.items() if (value := sp.cancel(old)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, 0)+value
    return clean(out)


def scale(op, value):
    return clean({key: value*coefficient for key, coefficient in op.items()})


def mul(left, right):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            for k in range(i+1):
                key = i+j-k
                out[key] = out.get(key, 0)+a*comb(i, k)*h**k*sp.diff(b, x, k)
    return clean(out)


def power(op, exponent):
    out = ONE
    for _ in range(exponent):
        out = mul(out, op)
    return out


def comm(a, b):
    return add(mul(a, b), scale(mul(b, a), -1))


W = {0: x, 2: h*h}
F1 = add(W, scale(power(W, 2), h*h))
radical = add(ONE, scale(W, 2*h*h))
assert comm(F1, radical) == ZERO
assert power(radical, 2) == add(ONE, scale(F1, 4*h*h))
assert clean({j: coefficient.subs(h, 0) for j, coefficient in F1.items()}) == {0: x}
assert sp.Poly(1+4*h*h*u, u).degree() == 1
radical_controls = []
for k in [2, 3, 5]:
    theta = add(ONE, scale(W, h*h))
    cyclic_F = scale(add(power(theta, k), scale(ONE, -1)), 1/(k*h*h))
    assert power(theta, k) == add(ONE, scale(cyclic_F, k*h*h))
    assert comm(theta, cyclic_F) == ZERO
    assert clean({j: coefficient.subs(h, 0) for j, coefficient in cyclic_F.items()}) == {0: x}
    cyclotomic = sp.cyclotomic_poly(k, root)
    assert sp.rem(root**k-1, cyclotomic, root) == 0
    assert sp.rem((root*z)**k-z**k, cyclotomic, root) == 0
    radical_controls.append({'k': k, 'exact_actual_radical_identity': 'PASS',
                             'classical_first_image': 'x', 'Long_pair': False})


# Selected denominators from the infinite differential-localization boundary.
localization_controls = []
for n in [1, 2, 5]:
    denominator = 1-n*h*x
    assert sp.cancel(sp.diff(1/denominator, x)-n*h/denominator**2) == 0
    assert denominator.subs(x, 1/(n*h)) == 0
    assert denominator.subs(h, 0) == 1
    localization_controls.append(n)
for n, m in [(1, 2), (2, 5), (1, 5)]:
    assert sp.Poly(1-n*h*x, x).gcd(sp.Poly(1-m*h*x, x)).degree() == 0


hash_names = [
    'TAYLOR_BRANCHES_AND_CENTRALIZER_RIGIDITY.md', 'check_taylor_centralizer.py',
    'ETALE_CENTRALIZERS_AND_GALOIS_COMPLETION.md', 'etale_centralizer_controls.json',
    'FINITE_CENTRALIZERS_AND_COMPLETION.md', 'finite_centralizer_controls.json',
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md',
    'QUANTUM_DIVISION_AND_COMPLETION.md', 'IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md',
    'BOUNDED_REALIZATION.md',
]
receipt = {
    'claim': 'QTB1', 'status': 'PASS',
    'Taylor_control': {'point': [1, 2], 'total_order': order,
                       'product_rule': 'PASS', 'derivative_intertwining': 'PASS',
                       'direct_translated_series': 'PASS'},
    'birational_denominator_control': {
        'fraction': 'x/(x^2+y)', 'derivative_certificate': '(partial_x-2*x*partial_y)b=1/(x^2+y)',
        'coprime_but_common_zero': [0, 0], 'identity': 'PASS'},
    'singular_divisor_control': {
        'divisor': 'y^2-x^2-x^3', 'total_origin_order': 2,
        'two_distinct_formal_prime_factors': 'written proof',
        'factorization_prefix_through_x10': 'PASS',
        'global_gcd_with_y_derivative': 1},
    'quartic_Galois_factor_control': {
        'field_polynomial': str(quartic), 'quadratic_subfield': str(quadratic),
        'quadratic_discriminant': 'u^2-4*v', 'quadratic_nontrivial': 'simple v-divisor',
        'quartic_non_Galois': 'written other-root nonsquare argument'},
    'normal_stabilizer_controls': {
        'D4': {'group_order': 8, 'stabilizer_order': 2, 'normal_closure_order': 4,
               'forbidden_Galois_factor_degree': 2},
        'S3': {'group_order': 6, 'stabilizer_order': 2, 'normal_closure_order': 6,
               'passes_normal_generation': True},
        'C3_regular': {'normal_stabilizer_closure_order': 1, 'forbidden_Galois_factor_degree': 3}},
    'actual_momentum_dependent_radical_control': {
        'W': 'x+h^2*p^2', 'F1': 'W+h^2*W^2', 'b': '1+2*h^2*W',
        'b_square': '1+4*h^2*F1', 'actual_Weyl_identity': 'PASS',
        'classical_graph': 'identity positions', 'Long_pair': False},
    'declared_cyclic_and_root_of_unity_controls': radical_controls,
    'infinite_localization_declared_samples': localization_controls,
    'external_dependencies': [
        {'source': 'https://arxiv.org/html/math/0401263',
         'used': 'Apply Corollary1.6 to each P[b] to obtain algebraicity without finite generation of all B.'},
        {'source': 'https://stacks.math.columbia.edu/tag/00MB',
         'used': 'Flatness of Noetherian local completion in the coordinate Taylor map.'},
        {'source': 'https://stacks.math.columbia.edu/tag/0AG0',
         'used': 'Factoriality of the regular local formal power series ring.'},
        {'source': 'https://stacks.math.columbia.edu/tag/09E3',
         'used': 'Galois transitivity and ramification indices of primes over a DVR.'},
        {'source': 'https://stacks.math.columbia.edu/tag/0BMB',
         'used': 'Purity of the branch locus for the normalization of each specified finite field.'},
        {'source': 'https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27',
         'used': 'Characteristic-zero finite etale affine-space covers are constant.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Declared Taylor, denominator, singular factor, quartic subfield, finite '
        'permutation, actual Weyl radical and localization controls only. Universal '
        'Taylor injectivity at arbitrary residue fields, flat regular-sequence '
        'arguments, localization without finite generation, all finite-intermediate '
        'unramified valuations, Galois-factor rigidity, normal-stabilizer generation, '
        'radical descent under constant extension and finite full completion use '
        'written proofs and explicit classical geometric inputs. The cyclic '
        'samples do not prove all k and the localization samples do not prove '
        'infinite nonfinite generation. No new Long field or radical certificate, '
        'full finite lift, carrier or all-representative nonlifting is supplied.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'taylor_centralizer_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QTB1', 'status': 'PASS', 'Taylor_product_and_derivatives': True,
    'D4_normal_stabilizer_closure_order': 4, 'S3_normal_stabilizer_closure_order': 6,
    'momentum_dependent_radical_certificate': 'PASS',
    'Long_certificate_constructed': False,
}, indent=2))
