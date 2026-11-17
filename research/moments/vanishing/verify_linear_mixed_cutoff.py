"""Exact b-function-lift and degree-recursion controls, not a general solver."""
from pathlib import Path
from itertools import product
from math import factorial
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def multiindices(n, max_degree):
    return [a for a in product(range(max_degree+1), repeat=n) if sum(a)<=max_degree]


def twist_derivative(A, x, Psi):
    return S.expand(S.diff(A, x)+S.diff(Psi, x)*A)


def theta_action(A, xs, Psi):
    return S.expand(-sum(twist_derivative(x*A, x, Psi) for x in xs))


def apply_certificate(b, terms, xs, t, s, Psi):
    coeffs = S.Poly(b, s)
    amplitudes = [S.Integer(1)]
    for _ in range(coeffs.degree()):
        amplitudes.append(theta_action(amplitudes[-1], xs, Psi))
    total = sum(coeffs.nth(j)*amplitudes[j] for j in range(len(amplitudes)))
    for c, k, j, u, v in terms:
        A = S.Integer(1)
        for x, count in zip(xs, v):
            for _ in range(count):
                A = twist_derivative(A, x, Psi)
        A *= S.prod(x**a for x, a in zip(xs, u))
        for _ in range(j):
            A = twist_derivative(A, t, Psi)
        total += c*t**k*A
    return S.expand(total)


def extract(b, terms, s, base_degrees):
    if b == 0:
        raise ValueError('zero b polynomial')
    roots = sorted(int(x) for x in S.polys.polytools.ground_roots(S.Poly(b, s)).keys() if x.is_Integer and x>=0)
    if not roots or 0 not in roots:
        raise ValueError('Gaussian certificate must include root zero')
    B = max(roots)
    slopes = [S.Rational(0)]
    for c, k, j, u, v in terms:
        if c == 0:
            continue
        drop = sum(v)-sum(u)
        if drop<1:
            raise ValueError('every remainder term must lower amplitude degree')
        slopes.append(S.Rational(k-j, drop))
    slope = max(slopes)
    offset = max([S.Rational(0)]+[S.Rational(N)-slope*sum(alpha) for alpha, N in base_degrees.items()])
    pure = max([slope]+[S.Rational(N, sum(alpha)) for alpha, N in base_degrees.items() if sum(alpha)>0])
    assert base_degrees.get(tuple([0]*len(next(iter(base_degrees))))) == 0
    return B, slope, offset, pure


def image_W(alpha, t):
    a, b = alpha
    return S.Rational(factorial(a), factorial(a-b))*t**(a-b) if a>=b else S.Integer(0)


def recurrence(b, terms, alpha, t, s, image):
    rhs = S.Integer(0)
    for c, k, j, u, v in terms:
        eta = tuple(a+u0 for a, u0 in zip(alpha, u))
        if any(v0>e for e, v0 in zip(eta, v)):
            continue
        factor = S.prod(S.ff(e, v0) for e, v0 in zip(eta, v))
        child = tuple(e-v0 for e, v0 in zip(eta, v))
        rhs -= c*(-1)**sum(v)*factor*t**k*S.diff(image(child, t), t, j)
    lhs = b.subs(s, sum(alpha))*image(alpha, t)
    assert S.expand(lhs-rhs) == 0, (alpha, lhs, rhs)
    return rhs


def run():
    Z, W, X, Y, t, s = S.symbols('Z W X Y t s')
    R0 = [(-2, 0, 0, (0, 0), (1, 1)), (1, 1, 0, (0, 0), (1, 0))]
    R1 = [(2, 0, 0, (1, 0), (2, 1)), (2, 0, 0, (0, 1), (1, 2)),
          (6, 0, 0, (0, 0), (1, 1)), (-1, 1, 0, (1, 0), (2, 0)),
          (-1, 1, 0, (0, 1), (1, 1)), (-3, 1, 0, (0, 0), (1, 0)),
          (1, 0, 1, (0, 0), (1, 0)), (1, 0, 0, (0, 0), (2, 0))]
    zero_terms = [(-1, 0, 0, (0, 0), (2, 0)), (-1, 0, 0, (0, 0), (0, 2))]
    assert apply_certificate(s, R0, [Z, W], t, s, t*W-Z*W) == 0
    assert apply_certificate(s*(s-1), R1, [Z, W], t, s, t*W-Z*W) == 0
    assert apply_certificate(s, zero_terms, [X, Y], t, s, -(X**2+Y**2)/2) == 0
    out0 = extract(s, R0, s, {(0, 0): 0})
    out1 = extract(s*(s-1), R1, s, {(0, 0): 0, (1, 0): 1})
    outzero = extract(s, zero_terms, s, {(0, 0): 0})
    assert out0 == (0, 1, 0, 1) and out1 == (1, 1, 0, 1) and outzero == (0, 0, 0, 0)
    checks = 0
    for alpha in multiindices(2, 7):
        for b, terms in [(s, R0), (s*(s-1), R1)]:
            recurrence(b, terms, alpha, t, s, image_W)
            checks += 1

    def gaussian_number(a):
        return S.Rational(factorial(a), 2**(a//2)*factorial(a//2)) if a%2 == 0 else S.Integer(0)

    def zero_image(alpha, t):
        return S.prod(gaussian_number(a) for a in alpha)

    for alpha in multiindices(2, 7):
        recurrence(s, zero_terms, alpha, t, s, zero_image)
        checks += 1
    # b(-d) would give the wrong scalar in the nonminimal certificate.
    alpha = (2, 0)
    rhs = recurrence(s*(s-1), R1, alpha, t, s, image_W)
    mutated = (s*(s-1)).subs(s, -sum(alpha))*image_W(alpha, t)
    assert S.expand(mutated-rhs) == 4*t**2
    rejected = 0
    for b, terms in [(S.Integer(0), R0), (s, [(1, 0, 0, (1, 0), (1, 0))])]:
        try:
            extract(b, terms, s, {(0, 0): 0})
        except ValueError:
            rejected += 1
    assert rejected == 2
    # Rational slope extraction includes both degree drops and derivative order.
    rational_terms = [(1, 3, 1, (0, 0), (3, 0))]
    abstract = extract(s, rational_terms, s, {(0, 0): 0})
    assert abstract == (0, S.Rational(2, 3), 0, S.Rational(2, 3))
    # This arithmetic-only control is NOT an annihilator certificate.
    return {
        'completed': True,
        'exponential_amplitude_identities_checked': 3,
        'normal_order_recurrences_checked': checks,
        'W_degree_zero_base_certificate': {'b': 's', 'base_bound': 0, 'slope': '1', 'pure_cutoff': '1'},
        'W_nonminimal_certificate': {'b': 's*(s-1)', 'base_bound': 1, 'slope': '1', 'pure_cutoff': '1'},
        'zero_polynomial_certificate': {'b': 's', 'base_bound': 0, 'slope': '0', 'pure_cutoff': '0'},
        'wrong_b_sign_mutation_detected': True,
        'zero_b_rejected': True, 'nonlowering_remainder_rejected': True,
        'rational_net_shift_arithmetic_only_control': '2/3,not an annihilator',
        'general_b_function_search_executed': False, 'general_algorithm_implemented': False,
        'general_algorithm_executed': False, 'universal_linear_cutoff_proof_certified_by_checker': False,
        'producer_package_executed': False, 'new_cubic_prefix_computed': False,
        'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/'linear_mixed_cutoff_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

