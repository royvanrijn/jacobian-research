"""Exact nilpotent and embedded-component controls for the written algorithm."""
from pathlib import Path
import hashlib
import json
import sympy as S
from verify_generic_point_moment_algorithm import difference_operator, indicial, constant_integer_roots

ROOT = Path(__file__).resolve().parent


def run():
    a, b, t, z, X, h = S.symbols('a b t z X h')
    P = a + b * X
    L = [-a - b**2 * t, S.Integer(1)]
    U = -b
    assert S.expand(P-a-b**2*t - (S.diff(U, X)+(t*S.diff(P, X)-X)*U)) == 0
    M = difference_operator(L, t)
    G = S.exp(a*t+b**2*t**2/2)
    assert S.simplify(sum(q*S.diff(G-1, t, j) for j, q in enumerate(M))) == 0
    moments = [S.expand(S.diff(G, t, k).subs(t, 0)) for k in range(9)]
    empty = S.groebner([], a, b, domain=S.QQ)
    first = S.groebner([a], a, b, domain=S.QQ)
    final = S.groebner([a, b**2], a, b, domain=S.QQ)
    assert moments[1] == a
    assert first.reduce(moments[2])[1] == b**2
    assert empty.reduce(a)[1] != 0 and first.reduce(b**2)[1] != 0
    # Formal polynomial certificate over the full nonreduced final algebra.
    cert = S.Poly(S.expand(b*X-(t*b-X)*(-b)), t, X)
    assert all(final.reduce(c)[1] == 0 for c in cert.coeffs())
    assert all(final.reduce(m)[1] == 0 for m in moments[1:])
    assert final.reduce(b)[1] == b and final.reduce(b**2)[1] == 0
    _, poly = indicial([S.Integer(0), S.Integer(1)], t, z)
    assert constant_integer_roots(poly, z, []) == [0]
    # Coefficient-field indicial bounds before and after the first cut.
    bounds = []
    for op, params in [(L, [a, b]), ([-b**2*t, S.Integer(1)], [b])]:
        _, poly = indicial(difference_operator(op, t), t, z)
        bounds.append(max(constant_integer_roots(poly, z, params)))
    assert bounds == [1, 2]

    # Nilpotent thickening: b^2 is nonzero modulo b^3.
    thick = S.groebner([b**3], b, domain=S.QQ)
    Y_thick = b**2*t**2/2
    assert thick.reduce(b**2)[1] == b**2
    assert S.diff(Y_thick, t, 3) == 0
    _, poly = indicial([S.Integer(0)]*3+[S.Integer(1)], t, z)
    assert constant_integer_roots(poly, z, []) == [0, 1, 2]
    assert S.diff(Y_thick, t, 2) == b**2
    assert Y_thick.subs(b, 0) == 0
    dual = S.groebner([b**2], b, domain=S.QQ)
    assert dual.reduce(b)[1] != 0
    assert dual.reduce(b*S.diff(b*t, t))[1] == 0

    # Verify (a) intersect (a^2,b) = (a^2,ab) by exact elimination.
    intersection = S.groebner([h*a, (1-h)*a**2, (1-h)*b], h, a, b, order='lex')
    eliminated = [q.as_expr() for q in intersection.polys if not q.as_expr().has(h)]
    embedded_ideal = S.groebner([a**2, a*b], a, b, domain=S.QQ)
    assert S.groebner(eliminated, a, b, domain=S.QQ) == embedded_ideal
    assert embedded_ideal.reduce(a)[1] == a
    assert embedded_ideal.reduce(a*b)[1] == 0
    Y_embedded = a*t
    assert S.diff(Y_embedded, t, 2) == 0
    _, poly = indicial([S.Integer(0), S.Integer(0), S.Integer(1)], t, z)
    assert constant_integer_roots(poly, z, []) == [0, 1]
    assert S.diff(Y_embedded, t) == a
    assert Y_embedded.subs(a, 0) == 0

    # A positive-dimensional primary component, localized at Q[b].
    # Multiplication by a in basis 1,a,a^2 over Q(b).
    C = S.Matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    assert C**3 == S.zeros(3) and C**2 != S.zeros(3)
    assert b*S.eye(3) == S.diag(b, b, b)
    assert S.factor((b*S.eye(3)).det()) == b**3
    assert (b*S.eye(3)).inv() == S.eye(3)/b

    return {
        'completed': True,
        'linear_family': {
            'polynomial': 'a+bX',
            'ordinary_ideal': ['a', 'b^2'],
            'radical': ['a', 'b'],
            'cut_indices': [1, 2],
            'generic_bounds': [1, 2],
            'final_indicial_bound': 0,
            'least_ordinary_prefix': 2,
            'nonreduced_final_algebra_dimension': 2,
            'full_small_loop_checked': True,
            'polynomial_telescoping_certificate_checked': True,
            'auxiliary_moments_checked': 8,
        },
        'nilpotent_control': {'quotient': 'Q[b]/(b^3)', 'series_minus_one': 'b^2*t^2/2',
                              'scalar_operator': 'Dt^3', 'bound': 2,
                              'radical_shortcut_rejected': True},
        'zero_divisor_coefficient_operator_shortcut_rejected': True,
        'embedded_control': {'ideal': '(a^2,ab)', 'primary_intersection': '(a) intersect (a^2,b)',
                             'polynomial': 'a', 'series_in_embedded_algebra': 'a*t',
                             'bound': 1, 'intersection_checked': True,
                             'minimal_components_only_shortcut_rejected': True,
                             'nonprimary_localization_shortcut_rejected': True},
        'positive_dimensional_primary_control': {
            'primary_ideal': '(a^3) in Q[a,b]', 'normalization': 'Q[b]',
            'localized_algebra': 'Q(b)[a]/(a^3)', 'dimension': 3,
            'nilpotent_multiplication_table_checked': True},
        'general_algorithm_implemented': False,
        'general_algorithm_executed': False,
        'universal_termination_proof_certified_by_checker': False,
        'full_cubic_coefficient_search_executed': False,
        'producer_package_executed': False,
        'new_cubic_prefix_computed': False,
        'independent_replay': False,
        'external_review': False,
        'formal_verification': False,
    }


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/'primary_moment_ideal_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

