"""Exact scalar polynomiality and full-source controls; no general module solver."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as S
from verify_generic_point_moment_algorithm import compose, indicial, constant_integer_roots

ROOT = Path(__file__).resolve().parent


def infinity_indicial(operator, t, z):
    parts = []
    for j, q in enumerate(operator):
        if q != 0:
            p = S.Poly(q, t)
            parts.append((p.degree()-j, p.LC(), j))
    if not parts:
        raise ValueError('zero operator')
    h = max(p[0] for p in parts)
    C = S.expand(sum(c*S.ff(z, j) for w, c, j in parts if w == h))
    assert C != 0
    return h, C


def polynomiality(f, L, t, z, params=()):
    assert S.simplify(sum(q*S.diff(f, t, j) for j, q in enumerate(L))) == 0
    _, C = infinity_indicial(L, t, z)
    roots = constant_integer_roots(C, z, params)
    bound = max(roots) if roots else -1
    candidate = S.series(f, t, 0, bound+1).removeO() if bound >= 0 else S.Integer(0)
    rhs = S.expand(sum(q*S.diff(candidate, t, j) for j, q in enumerate(L)))
    M = L if rhs == 0 else compose([-S.diff(rhs, t), rhs], L, t)
    _, origin = indicial(M, t, z)
    origin_roots = constant_integer_roots(origin, z, params)
    jet_bound = max(origin_roots) if origin_roots else -1
    difference = S.series(f-candidate, t, 0, jet_bound+1).removeO() if jet_bound >= 0 else S.Integer(0)
    first = next((k for k in range(jet_bound+1) if S.expand(difference).coeff(t, k) != 0), None)
    accepted = first is None
    if not accepted:
        assert first > bound
    assert accepted == (S.simplify(f-candidate) == 0)
    return {'infinity_indicial': str(S.factor(C)), 'polynomial_degree_bound': bound,
            'origin_jet_bound': jet_bound, 'polynomial': accepted,
            'first_tail_index': first, 'candidate': str(candidate)}


def gaussian(expr, variables):
    Z, W, T = variables
    value = S.Integer(0)
    for (a, b, c), coefficient in S.Poly(S.expand(expr), Z, W, T).terms():
        if a == b and c % 2 == 0:
            value += coefficient*factorial(a)*S.Rational(factorial(c),2**(c//2)*factorial(c//2))
    return S.cancel(value)


def run():
    t, z, a = S.symbols('t z a')
    examples = [
        ('constant', S.Integer(1), [0, 1]),
        ('polynomial', t**5+2*t**2, [0]*6+[1]),
        ('exponential', S.exp(t), [-1, 1]),
        ('marked_quartic', t/(1-t), [-1, t*(1-t)]),
        ('simple_pole', 1/(1-t), [-1, 1-t]),
        ('delayed_rational', t**5/(1-t), [-5+4*t, t*(1-t)]),
        ('zero', S.Integer(0), [t, 1]),
    ]
    results = {}
    for name, f, raw in examples:
        results[name] = polynomiality(f, list(map(S.sympify, raw)), t, z)
    assert [results[n]['polynomial'] for n, _, _ in examples] == [True, True, False, False, False, False, True]
    assert results['delayed_rational']['polynomial_degree_bound'] == 4
    assert results['delayed_rational']['first_tail_index'] == 5
    # The candidate truncation alone would wrongly accept the delayed series.
    assert S.series(t**5/(1-t), t, 0, 5).removeO() == 0
    _, generic_C = infinity_indicial([-a, 1+t], t, z)
    assert constant_integer_roots(generic_C, z, [a]) == []
    specialized = []
    for N in [0, 1, 7]:
        out = polynomiality((1+t)**N, [S.Integer(-N), 1+t], t, z)
        assert out['polynomial'] and out['polynomial_degree_bound'] == N
        specialized.append(N)

    Z, W, T, alpha, beta, gamma = S.symbols('Z W T alpha beta gamma')
    B = W*(1+Z)-(1+3*Z/2+Z**2/2)*T**2
    q = 1-t
    g = (t+beta)/q
    R = 1+t*(3-t)*beta+t*beta**2
    assert S.cancel(q**2*(1+t*(2+3*g+g**2))-R) == 0
    full = R**S.Rational(-1, 2)*S.exp(alpha*(t+beta)/q+gamma**2*q**2/(2*R))
    multipliers = [((0, 0, 0), S.Integer(1)), ((1, 0, 0), Z), ((0, 1, 0), W),
                   ((0, 0, 1), T), ((0, 0, 2), T**2), ((2, 0, 0), Z**2),
                   ((1, 1, 0), Z*W)]
    source_images = []
    for (az, bw, ct), Q in multipliers:
        f = S.diff(full, alpha, az, beta, bw, gamma, ct).subs({alpha: 0, beta: 0, gamma: 0})
        source_images.append(S.cancel(f))
    assert source_images[0] == 1
    assert S.cancel(source_images[1]-t/q) == 0
    assert S.cancel(source_images[2]+t*(3-t)/2) == 0
    assert source_images[3] == 0
    assert S.cancel(source_images[4]-q**2) == 0
    power = S.Integer(1)
    checked = 0
    for m in range(7):
        for ((az, bw, ct), Q), f in zip(multipliers, source_images):
            exact = gaussian(Q*power, (Z, W, T))
            expected = S.diff(f, t, m).subs(t, 0)
            assert S.cancel(exact-expected) == 0, (m, Q, exact, expected)
            checked += 1
        power = S.expand(power*B)
    # Principal parts generate the point-module quotient.
    for k in range(1, 7):
        assert S.cancel(S.diff(1/q, t, k-1)-factorial(k-1)/q**k) == 0
    assert S.cancel((t-1)/q) == -1

    # Entire P=W integration module: exact monomial image formula.
    safe = 0
    for az in range(6):
        for bw in range(6-az):
            f = factorial(az)*t**(az-bw)/factorial(az-bw) if az >= bw else S.Integer(0)
            for m in range(8):
                assert gaussian(Z**az*W**(bw+m), (Z, W, T)) == S.diff(f, t, m).subs(t, 0)
                safe += 1
    return {
        'completed': True, 'scalar_polynomiality_examples': results,
        'delayed_jet_shortcut_rejected': True,
        'generic_specialization_shortcut_rejected': True,
        'specialized_polynomial_degrees': specialized,
        'quartic_source_identity_checked': True, 'quartic_exact_moments_checked': checked,
        'quartic_pure_series': '1', 'quartic_marked_Z_series': 't/(1-t)',
        'quartic_full_gaussian_image_written_proof': 'K[t,(1-t)^(-1)]',
        'principal_parts_checks': 6, 'safe_W_monomial_moments_checked': safe,
        'general_algorithm_implemented': False, 'general_algorithm_executed': False,
        'module_generator_algorithm_executed': False,
        'universal_finite_multiplier_proof_certified_by_checker': False,
        'producer_package_executed': False, 'new_cubic_prefix_computed': False,
        'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/'finite_multiplier_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

