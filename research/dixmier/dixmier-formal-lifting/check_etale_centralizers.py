#!/usr/bin/env python3
"""Declared QET1 geometry controls; universal etale/purity arguments are written."""
from pathlib import Path
from math import comb
from itertools import permutations
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h = sp.symbols('h', nonzero=True)
x, y, u, v, z, t1, t2 = sp.symbols('x y u v z t1 t2')
A, B, C, w = sp.symbols('A B C w')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['finite_centralizer_controls.json', 'division_index_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, digest in prior['hashes'].items():
        assert sha(HERE/rel) == digest, rel


# The retained cubic field: exact degree proof is in QDI1.
cubic = C*w**3-2*w**2+B*w-2*A
disc = sp.discriminant(cubic, w)
expected_disc = 4*B**2-4*C*B**3+72*A*B*C-64*A-108*A**2*C**2
assert sp.expand(disc-expected_disc) == 0
disc_A = sp.discriminant(disc, A)
assert sp.expand(disc_A-64*(4-3*B*C)**3) == 0
pdisc = sp.Poly(disc, A, domain=sp.QQ.frac_field(B, C))
assert pdisc.degree() == 2
assert pdisc.gcd(pdisc.diff()).degree() == 0
A_from_w = (C*w**3-2*w**2+B*w)/2
assert sp.expand(cubic.subs(A, A_from_w)) == 0


# Fixed sheets are local inertia data, not a globally chosen branch.
def compose(left, right):
    return tuple(left[right[i]] for i in range(len(left)))


def fixed(perm):
    return {i for i, j in enumerate(perm) if i == j}


s3 = list(permutations(range(3)))
transpositions = [p for p in s3 if len(fixed(p)) == 1]
three_cycles = [p for p in s3 if len(fixed(p)) == 0]
assert len(transpositions) == 3 and len(three_cycles) == 2
g1, g2 = (1, 0, 2), (0, 2, 1)
assert fixed(g1) and fixed(g2) and not fixed(g1).intersection(fixed(g2))
generated = {(0, 1, 2)}
while True:
    new = generated | {compose(g, p) for g in (g1, g2) for p in generated}
    if new == generated:
        break
    generated = new
assert len(generated) == 6 and not fixed(compose(g1, g2))
assert fixed((1, 0)) == set()  # quadratic nontrivial inertia


# A genuine quadratic Galois Ore boundary, with nonconstant units.
relation = z*z-1-h*t1
def nabla1(expr):
    return sp.cancel(sp.diff(expr, t1)+h/(2*z)*sp.diff(expr, z))


def nabla2(expr):
    return sp.diff(expr, t2)


assert nabla1(relation) == 0
assert nabla2(relation) == 0
for expr in [t1, t2, z, 1/z, 1/(z+1), t1*t2*z/(z+1)]:
    assert sp.cancel(nabla1(nabla2(expr))-nabla2(nabla1(expr))) == 0
assert nabla1(t1) == 1 and nabla2(t2) == 1
assert nabla1(t2) == 0 and nabla2(t1) == 0
assert sp.expand((z-1)*(z+1)-h*t1-relation) == 0
assert sp.cancel(nabla1((z-1)*(z+1))) == h
assert sp.cancel(z*z/(1+h*t1)-1-relation/(1+h*t1)) == 0
assert sp.cancel(nabla1(-z)+h/(2*z)) == 0
assert sp.diff(1+h*t1, t1) == h
# Reduction z^2=1 with z+1 inverted keeps precisely z=1.
assert sp.cancel((z*z-1)/(z+1)-(z-1)) == 0
assert relation.subs({h: 0, z: 1}) == 0
root_prefix = sum(sp.binomial(sp.Rational(1, 2), n)*h**n*t1**n for n in range(9))
assert sp.series(root_prefix**2-(1+h*t1), h, 0, 9).removeO() == 0
root_degrees = []
for n in range(1, 9):
    coefficient = sp.expand(root_prefix).coeff(h, n)
    assert coefficient != 0 and sp.degree(coefficient, t1) == n
    root_degrees.append(n)
quadratic_root_field = sp.Poly(z*z-(1+h*t1), z)
assert sp.expand(sp.discriminant(quadratic_root_field.as_expr(), z)-4*(1+h*t1)) == 0


# Leading birationality can coexist with actual quadratic field degree.
U = x*x*y+y
V = x*y
assert sp.expand((v*x*x-u*x+v).subs({u: U, v: V})) == 0
assert sp.cancel((x*x*y)/(x*y)-x) == 0
assert sp.cancel((x*y)**2/(x*x*y)-y) == 0
assert sp.expand(sp.Matrix([[sp.diff(U, x), sp.diff(U, y)],
                  [sp.diff(V, x), sp.diff(V, y)]]).det()-y*(x*x-1)) == 0
quad_disc = sp.discriminant(v*x*x-u*x+v, x)
assert sp.expand(quad_disc-(u*u-4*v*v)) == 0
assert sp.factor(quad_disc) == (u-2*v)*(u+2*v)


# Scalar quantum quadratic positions: their rational mate really exists.
# Ordinary differential operators over L(x,y), normal ordered p_x,p_y.
ONE = {(0, 0): sp.Integer(1)}
ZERO = {}


def clean(op):
    return {key: value for key, old in op.items() if (value := sp.cancel(old)) != 0}


def add(*ops):
    out = {}
    for op in ops:
        for key, value in op.items():
            out[key] = out.get(key, 0)+value
    return clean(out)


def scale(op, factor):
    return clean({key: factor*value for key, value in op.items()})


def mul(left, right):
    out = {}
    for (ip, iq), a in left.items():
        for (jp, jq), b in right.items():
            for k in range(ip+1):
                for ell in range(iq+1):
                    value = a*comb(ip, k)*comb(iq, ell)*h**(k+ell)*sp.diff(b, x, k, y, ell)
                    key = (ip+jp-k, iq+jq-ell)
                    out[key] = out.get(key, 0)+value
    return clean(out)


def comm(left, right):
    return add(mul(left, right), scale(mul(right, left), -1))


f1 = x+h*h*x*x
positions = [{(0, 0): f1}, {(0, 0): y}]
fprime = sp.diff(f1, x)
mates = [{(1, 0): -1/fprime}, {(0, 1): sp.Integer(-1)}]
assert comm(*positions) == ZERO and comm(*mates) == ZERO
for i in range(2):
    for j in range(2):
        assert comm(positions[j], mates[i]) == (scale(ONE, h) if i == j else ZERO)
sigma_x = -1/h**2-x
assert sp.cancel(f1.subs(x, sigma_x)-f1) == 0
assert sp.cancel(sigma_x.subs(x, sigma_x)-x) == 0
assert sp.discriminant(h*h*x*x+x-u, x) == 1+4*h*h*u
assert sp.Poly(fprime, x).degree() == 1
assert sp.cancel(-fprime*mates[0][(1, 0)]) == 1
rational_mate_prefix = -sum((-2*h*h*x)**n for n in range(6))
assert sp.series(fprime*rational_mate_prefix+1, h, 0, 12).removeO() == 0
assert sp.degree(sp.expand(rational_mate_prefix).coeff(h, 10), x) == 5


# Actual parameter-pole descent with momentum-dependent regular positions.
def power(op, exponent):
    out = ONE
    for _ in range(exponent):
        out = mul(out, op)
    return out


def h_valuation(coefficient):
    numerator, denominator = sp.fraction(sp.cancel(coefficient))
    return (min(k[0] for k in sp.Poly(numerator, h).monoms())
            -min(k[0] for k in sp.Poly(denominator, h).monoms()))


X = {(0, 0): x}
Y = {(0, 0): y}
PX = {(1, 0): sp.Integer(1)}
PY = {(0, 1): sp.Integer(1)}
regular_F = [add(X, scale(power(PX, 2), h*h)),
             add(Y, scale(power(PY, 2), h*h))]


def eval_at_F(poly):
    out = ZERO
    for (i, j), coefficient in sp.Poly(poly, x, y).terms():
        out = add(out, scale(mul(power(regular_F[0], i), power(regular_F[1], j)), coefficient))
    return out


def remove_poles(op):
    polar_part = ZERO
    orders = []
    while True:
        valuation = min(h_valuation(coefficient) for coefficient in op.values())
        if valuation >= 0:
            return op, polar_part, orders
        leading = clean({
            key: sp.cancel(coefficient*h**(-valuation)).subs(h, 0)
            for key, coefficient in op.items()
        })
        assert set(leading) == {(0, 0)}
        correction = scale(eval_at_F(leading[(0, 0)]), h**valuation)
        op = add(op, scale(correction, -1))
        polar_part = add(polar_part, correction)
        orders.append(valuation)


F1, F2 = regular_F
classical_shift = [
    add(scale(power(F1, 2), 3), power(F2, 2)),
    scale(mul(F1, F2), 2),
]
polar_form = [
    add(scale(mul(F1, F2), 2/h**3), scale(power(F2, 2), 1/h)),
    add(scale(power(F1, 2), 1/h**3), scale(mul(F1, F2), 2/h)),
]
base_mates = [add(scale(PX, -1), classical_shift[0]),
              add(scale(PY, -1), classical_shift[1])]
generic_mates = [add(base_mates[i], polar_form[i]) for i in range(2)]
assert comm(*regular_F) == ZERO
assert comm(*base_mates) == ZERO and comm(*generic_mates) == ZERO
normalized = [remove_poles(op) for op in generic_mates]
for i in range(2):
    assert normalized[i][0] == base_mates[i]
    assert normalized[i][1] == polar_form[i]
    assert normalized[i][2] == [-3, -1]
final_mates = [add(normalized[i][0], scale(classical_shift[i], -1)) for i in range(2)]
assert final_mates == [scale(PX, -1), scale(PY, -1)]
assert comm(*final_mates) == ZERO
regular_extra = scale(mul(power(F1, 2), F2), h*h)
curved_generic = [generic_mates[0], add(generic_mates[1], regular_extra)]
curved_normalized = [remove_poles(op)[0] for op in curved_generic]
expected_curvature = scale(mul(F1, F2), -2*h**3)
assert comm(*curved_generic) == expected_curvature
assert comm(*curved_normalized) == expected_curvature
for family in [generic_mates, base_mates, final_mates, curved_generic, curved_normalized]:
    for i in range(2):
        for j in range(2):
            assert comm(regular_F[j], family[i]) == (scale(ONE, h) if i == j else ZERO)
for expr1, expr2 in [
    (2*x*y/h**3+y*y/h, x*x/h**3+2*x*y/h),
    (3*x*x+y*y, 2*x*y),
]:
    assert sp.diff(expr2, x)-sp.diff(expr1, y) == 0


hash_names = [
    'ETALE_CENTRALIZERS_AND_GALOIS_COMPLETION.md', 'check_etale_centralizers.py',
    'FINITE_CENTRALIZERS_AND_COMPLETION.md', 'finite_centralizer_controls.json',
    'LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md',
    'IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md',
    'QUANTUM_DIVISION_AND_COMPLETION.md', 'division_index_controls.json',
    'BOUNDED_REALIZATION.md',
]
receipt = {
    'claim': 'QET1', 'status': 'PASS',
    'retained_cubic_field_control': {
        'minimal_polynomial': str(cubic), 'field_degree': 3,
        'discriminant': str(sp.factor(disc)), 'discriminant_as_quadratic_in_A': str(sp.factor(disc_A)),
        'discriminant_squarefree_in_A': True, 'normal_closure': 'S3',
        'exact_degree_and_irreducibility': 'retained QDI1 written cubic map proof',
        'rank': 3, 'Long_candidate': False,
    },
    'inertia_controls': {
        'transposition_fixed_sheet_counts': [len(fixed(p)) for p in transpositions],
        'three_cycle_fixed_sheet_counts': [len(fixed(p)) for p in three_cycles],
        'two_fixed_sheet_generators_generate_S3': True, 'global_common_fixed_sheet': False,
        'nontrivial_quadratic_inertia_fixed_sheets': 0,
    },
    'actual_parameter_pole_descent': {
        'positions': ['x+h^2*p_x^2', 'y+h^2*p_y^2'],
        'pole_orders_removed': [-3, -1], 'algorithmic_normalization': 'PASS',
        'flatness_preserved': True, 'prescribed_classical_shift_removed_finitely': True,
        'nonzero_integral_curvature_preserved': '-2*h^3*F1*F2',
        'rational_operator_denominators_removed': False, 'Long_candidate': False,
    },
    'Galois_boundary_with_units': {
        'ring': 'R0[t1,t2,z,z^-1,(z+1)^-1]/(z^2-1-h*t1)',
        'field_degree': 2, 'nonconstant_base_unit': '1+h*t1=z^2',
        'relation_preserved_by_connections': True, 'connections_commute': True,
        'classical_branch': 'z=1', 'formal_root_prefix_coefficient_degrees': root_degrees,
        'full_flat_Ore_tuple': 'positions t_i, mates -p_i',
        'scalar_finite_Weyl_algebra': False,
    },
    'leading_field_boundary': {
        'map': ['x^2*y+y', 'x*y'], 'leading_field_degree': 1, 'actual_field_degree': 2,
        'actual_discriminant': str(sp.factor(quad_disc)), 'Jacobian': 'y*(x^2-1)',
        'Keller': False,
    },
    'scalar_quadratic_position_control': {
        'positions': ['x+h^2*x^2', 'y'], 'classical_positions': ['x', 'y'],
        'centralizer': 'L[x,y], written leading-order proof',
        'centralizer_field_degree': 2, 'field_involution': '-h^-2-x',
        'finite_individual_mate': 'excluded by QET1 and direct written leading-order argument',
        'rational_mate': '-(1+2*h^2*x)^-1*p_x', 'all_rational_Weyl_relations': 'PASS',
        'Long_candidate': False,
    },
    'external_dependencies': [
        {'source': 'https://arxiv.org/html/math/0401263', 'used': 'Corollary1.6: commutative Weyl subalgebras have transcendence degree at most the rank; algebraic-closure base change covers L.'},
        {'source': 'https://stacks.math.columbia.edu/tag/02LR', 'used': 'Normalization form of Zariski Main: the actual etale centralizer is open in the finite normalization.'},
        {'source': 'https://stacks.math.columbia.edu/tag/09E3', 'used': 'Galois transitivity and equal ramification indices over a DVR.'},
        {'source': 'https://stacks.math.columbia.edu/tag/0BMB', 'used': 'Purity upgrades absence of codimension-one branch to finite etaleness.'},
        {'source': 'https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27', 'used': 'Characteristic-zero finite etale affine-space covers are constant, as in QFC1.'},
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Declared discriminants, field charts, permutation inertia, Ore connection and '
        'formal branch identities, a leading-versus-actual degree countercontrol and '
        'a scalar quadratic rational mate and actual parameter-pole normalization. Differential simplicity, etaleness, '
        'height-one boundary retention, birational/Galois polynomiality, all-rank '
        'integral completion and the scalar dimension/unit facts use written proofs '
        'and explicit external geometric inputs. No finite-generation theorem for '
        'all Long centralizers or actual Long field certificate is supplied. '
        'Centralizer field degree is not identified with a full quantum division index.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'etale_centralizer_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QET1', 'status': 'PASS', 'retained_cubic_field_normal_closure': 'S3',
    'scalar_quadratic_rational_relations': 'PASS', 'Galois_boundary_violates_unit_hypothesis': True,
    'leading_birational_actual_quadratic_control': True, 'Long_candidate_constructed': False,
}, indent=2))
