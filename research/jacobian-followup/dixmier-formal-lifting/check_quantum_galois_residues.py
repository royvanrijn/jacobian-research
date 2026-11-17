#!/usr/bin/env python3
"""Exact QGD1 controls; group/valuation exhaustion is a written theorem."""
from pathlib import Path
from itertools import permutations
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
h, x, p, r, s, t, a, u = sp.symbols('h x p r s t a u')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name in ['division_index_controls.json', 'centralizer_index_controls.json']:
    prior = json.loads((HERE/name).read_text())
    assert prior['status'] == 'PASS', name
    for rel, expected in prior['hashes'].items():
        assert sha(HERE/rel) == expected, rel

prior_division = json.loads((HERE/'division_index_controls.json').read_text())
assert prior_division['Long_classical_field_degree']['value'] == 3
prior_centralizer = json.loads((HERE/'centralizer_index_controls.json').read_text())
assert [q['m'] for q in prior_centralizer['genuine_rational_residual_family_controls']] == [2, 3, 5]

# The whole classical Long field, including its independent fourth output.
cubic = r*a**3-2*a**2+s*a+4*t
disc = sp.discriminant(cubic, a)
expected_disc = 4*s*s-4*r*s**3+128*t-144*r*s*t-432*r*r*t*t
assert sp.expand(disc-expected_disc) == 0
disc_in_t = sp.discriminant(disc, t)
assert sp.expand(disc_in_t-256*(4-3*r*s)**3) == 0
assert sp.degree(disc, t) == 2
coefficient_field = sp.QQ.frac_field(r, s)
assert sp.gcd(sp.Poly(disc, t, domain=coefficient_field),
              sp.Poly(sp.diff(disc, t), t, domain=coefficient_field)).degree() == 0
assert sp.discriminant(cubic, a) != 0

# First-order Ore operators are a(x)*p+b(x), with [p,x]=h.
def first_order_comm(left, right):
    la, lb = left
    ra, rb = right
    return (
        sp.cancel(h*(la*sp.diff(ra, x)-ra*sp.diff(la, x))),
        sp.cancel(h*(la*sp.diff(rb, x)-ra*sp.diff(lb, x))),
    )


# A cyclic degree-three quotient is allowed: its residue cover is Galois.
omega = (-1+sp.sqrt(-3))/2
assert sp.simplify(omega**3-1) == 0
assert sp.simplify(omega-1) != 0
Fcoef, Tvalue = 1/(3*x*x), x**3
assert first_order_comm((Fcoef, 0), (0, Tvalue)) == (0, h)
assert sp.simplify(Fcoef.subs(x, omega*x)/omega-Fcoef) == 0
assert sp.simplify(Tvalue.subs(x, omega*x)-Tvalue) == 0
for power in (1, 2):
    assert sp.simplify(sum(omega**(j*power)*x**power for j in range(3))) == 0
assert sp.simplify(sum((omega**j*x)**3 for j in range(3))-3*x**3) == 0
cyclic_disc = sp.discriminant(a**3-t, a)
assert sp.expand(cyclic_disc+27*t*t) == 0
assert sp.simplify(cyclic_disc-(3*sp.sqrt(-3)*t)**2) == 0

# A pole-bearing involution moves the valuation; it cannot be silently discarded.
sigma_x = -1/h-x


def sigma(function):
    return sp.cancel(sp.sympify(function).subs(x, sigma_x))


def sigma_operator(operator):
    oa, ob = operator
    return -sigma(oa), sigma(ob)


U = x+h*x*x
Vcoef = 1/(1+2*h*x)
assert sp.cancel(sigma(sigma_x)-x) == 0
assert first_order_comm((Vcoef, 0), (0, U)) == (0, h)
assert sp.cancel(sigma(U)-U) == 0
assert all(sp.cancel(q) == 0 for q in (
    sigma_operator((Vcoef, 0))[0]-Vcoef,
    sigma_operator((Vcoef, 0))[1],
))
assert first_order_comm((-1, 0), (0, sigma_x)) == (0, h)
assert sp.degree(h*a*a+a-u, a) == 2
assert sp.degree((h*a*a+a-u).subs(h, 0), a) == 1


def h_order(polynomial):
    poly = sp.Poly(sp.expand(polynomial), h)
    assert not poly.is_zero
    return min(monomial[0] for monomial, coefficient in poly.terms())


def h_valuation(function):
    function = sp.cancel(function)
    if function == 0:
        return None
    num, den = sp.fraction(function)
    return h_order(num)-h_order(den)


def h_residue(function):
    function = sp.cancel(function)
    if function == 0:
        return sp.S.Zero
    num, den = sp.fraction(function)
    vn, vd = h_order(num), h_order(den)
    assert vn >= vd
    if vn > vd:
        return sp.S.Zero
    numerator = sp.Poly(num, h).coeff_monomial(h**vn)
    denominator = sp.Poly(den, h).coeff_monomial(h**vd)
    return sp.cancel(numerator/denominator)


assert h_valuation(x) == 0 and h_valuation(sigma_x) == -1
branch_separator_argument = h*x*x
assert h_valuation(branch_separator_argument) == 1
assert h_valuation(sigma(branch_separator_argument)) == -1
selector = 1/(1+branch_separator_argument)
assert h_valuation(selector) == 0 and h_residue(selector) == 1
assert h_valuation(sigma(selector)) == 1 and h_residue(sigma(selector)) == 0
selected_lift = selector**2*x
assert h_valuation(selected_lift) == 0 and h_residue(selected_lift) == x
assert h_valuation(sigma(selected_lift)) == 1
averaged = sp.cancel(selected_lift+sigma(selected_lift))
assert sp.cancel(sigma(averaged)-averaged) == 0
assert h_valuation(averaged) == 0 and h_residue(averaged) == x
# This second rational control also fails QCI1 finite admission: V has degree zero.
vnum, vden = sp.fraction(Vcoef)
assert sp.degree(vnum, x)-sp.degree(vden, x)+1 == 0

# Finite normalizer controls distinguish a cubic source from its S3 Galois closure.
def compose(left, right):
    return tuple(left[right[j]] for j in range(len(left)))


def inverse(element):
    return tuple(element.index(j) for j in range(len(element)))


identity = (0, 1, 2)
symmetric = set(permutations(range(3)))
point_stabilizer = {g for g in symmetric if g[0] == 0}


def normalizer(group, subgroup):
    return {g for g in group
            if {compose(compose(g, q), inverse(g)) for q in subgroup} == subgroup}


assert len(symmetric) == 6 and len(point_stabilizer) == 2
assert len(normalizer(symmetric, point_stabilizer)) == 2
assert len(normalizer(symmetric, point_stabilizer))//len(point_stabilizer) == 1
cyclic = {identity, (1, 2, 0), (2, 0, 1)}
cyclic_stabilizer = {g for g in cyclic if g[0] == 0}
assert cyclic_stabilizer == {identity}
assert len(normalizer(cyclic, cyclic_stabilizer)) == 3
assert 6//2 == 3  # Classical cubic source degree, not a Galois source cover.

hash_names = [
    'FINITE_QUANTUM_SYMMETRIES_AND_GALOIS_RESIDUES.md',
    'check_quantum_galois_residues.py',
    'QUANTUM_DIVISION_AND_COMPLETION.md', 'division_index_controls.json',
    'CENTRALIZER_INDEX_AND_RATIONAL_POLARIZATION.md', 'centralizer_index_controls.json',
    'FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md',
    'PROOF.md',
]
receipt = {
    'claim': 'QGD1', 'status': 'PASS',
    'Long_classical_cubic': {
        'polynomial': str(cubic), 'discriminant': str(sp.expand(disc)),
        'discriminant_in_t': str(sp.factor(disc_in_t)),
        'squarefree_in_t': 'PASS',
        'degree': 3,
        'irreducibility': 'written cubic function-field degree proof, retained QDI1 chart',
        'Galois_closure_group': 'S3, written irreducible-cubic criterion',
        'classical_deck_group_order': 1,
        'independent_fourth_output': 'retained QDI1 birational chart',
    },
    'cyclic_rational_positive_control': {
        'F': 'p/(3*x^2)', 'T': 'x^3',
        'automorphism': ['x -> omega*x', 'p -> omega^-1*p'],
        'exact_Weyl_relation': 'PASS', 'fixed_generators': 'PASS',
        'group_order': 3, 'valuation_stabilizer_order': 3,
        'valuation_branches': 1, 'generic_index': 3, 'completed_index': 3,
        'residue_extension_Galois': True,
        'finite_carrier': 'excluded separately by QCI1 bounded powers',
        'Long_lift': False,
    },
    'pole_bearing_rational_involution_control': {
        'U': str(U), 'V_coefficient_of_p': str(Vcoef),
        'sigma_x': str(sigma_x), 'sigma_p': '-p',
        'order_two': 'PASS', 'exact_Weyl_relation': 'PASS', 'fixed_generators': 'PASS',
        'v_x': 0, 'v_sigma_x': -1,
        'branch_selector': str(selector),
        'selector_valuations': [0, 1], 'selector_residues': ['1', '0'],
        'selected_lift': str(selected_lift),
        'averaged_invariant': 'selected_lift + sigma(selected_lift)',
        'averaged_invariance_and_residue_x': 'PASS',
        'group_order': 2, 'valuation_stabilizer_order': 1,
        'valuation_branches': 2, 'generic_index': 2, 'completed_index': 1,
        'residue_extension_Galois': True,
        'Bernstein_degree_V': 0,
        'finite_carrier': 'excluded separately by QCI1 bounded powers',
        'Long_lift': False,
    },
    'finite_monodromy_controls': {
        'S3_order': 6, 'point_stabilizer_order': 2,
        'point_stabilizer_normalizer_order': 2, 'cubic_source_deck_order': 1,
        'cyclic_regular_deck_order': 3,
        'larger_Galois_closure_with_intermediate_source_excluded': False,
        'quantum_realizations': False,
    },
    'research_context': [
        {
            'source': 'https://arxiv.org/html/1708.07923',
            'used': (
                'Tikaradze Theorem1 is comparison context for actual invariant orders; '
                'it is not used to infer division-ring quotient rigidity.'
            ),
        },
        {
            'source': 'https://arxiv.org/html/2610.15469v2',
            'used': (
                'General Weyl index symmetry is context only. The finite-quotient '
                'left/right equality is proved directly here using Artin bounds '
                'and valuation-branch residue independence.'
            ),
        },
    ],
    'hashes': {name: sha(HERE/name) for name in hash_names},
    'coverage_boundary': (
        'Exact discriminants and squarefreeness, retained classical cubic chart, '
        'first-order rational Weyl and finite-group action identities, a cyclic '
        'Galois residue control, a pole-bearing valuation-moving involution, '
        'branch selectors, averaged residues and finite subgroup normalizers. '
        'Universal branch approximation, residue Galois descent, characteristic-zero '
        'kernel rigidity, Artin dimension bounds, generic/completed index accounting, '
        'centralizer and automorphism rigidity, all-degree Long exclusion, stabilization '
        'and rational source/target transport use written proofs. The controls do not '
        'prove group or operator exhaustion. The conclusion excludes direct finite '
        'fixed-division-ring models; a larger Galois closure with an intermediate '
        'subgroup-fixed source, different nonsymmetric embeddings, and general '
        'higher-rank matrix actions remain open. No actual Long carrier is built.'
    ),
    'independent_replay': False, 'formal_verification': False, 'external_review': False,
}
(HERE/'quantum_galois_residue_controls.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({
    'claim': 'QGD1', 'status': 'PASS',
    'Long_discriminant_squarefree': True,
    'cyclic_degree3_positive_control': 'PASS',
    'valuation_moving_involution_and_branch_selector': 'PASS',
    'S3_intermediate_source_boundary': 'PASS',
    'Long_lift_constructed': False,
}, indent=2))
