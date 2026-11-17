"""Exact differential/indicial components and one full small descent."""
from pathlib import Path
import hashlib
import json
from math import comb
import sympy as S

ROOT = Path(__file__).resolve().parent


def compose(left, right, t):
    out = [S.Integer(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            for k in range(i + 1):
                out[i - k + j] += comb(i, k) * a * S.diff(b, t, k)
    return [S.cancel(a) for a in out]


def difference_operator(operator, t):
    p0 = operator[0]
    return operator if p0 == 0 else compose([-S.diff(p0, t), p0], operator, t)


def indicial(operator, t, z):
    pieces = []
    for j, coef in enumerate(operator):
        if coef == 0:
            continue
        poly = S.Poly(coef, t)
        valuation = min(m[0] for m in poly.monoms())
        pieces.append((valuation - j, poly.nth(valuation), j))
    if not pieces:
        raise ValueError('The annihilator must be nonzero after reduction.')
    mu = min(p[0] for p in pieces)
    result = S.expand(sum(c * S.ff(z, j) for v, c, j in pieces if v == mu))
    assert result != 0
    return mu, result


def constant_integer_roots(poly, z, params):
    numerator = S.fraction(S.cancel(poly))[0]
    components = S.Poly(numerator, *params).coeffs() if params else [numerator]
    gcd = S.Poly(0, z, domain=S.QQ)
    for c in components:
        gcd = S.gcd(gcd, S.Poly(c, z, domain=S.QQ))
    assert not gcd.is_zero
    return sorted(int(r) for r in S.polys.polytools.ground_roots(gcd).keys() if r.is_Integer and r >= 0)


def run():
    a, b, c, t, z, X, Y = S.symbols('a b c t z X Y')
    params = [a, b, c]
    v = b ** 2 + c ** 2
    P = a + b * X + c * Y
    # Full polynomial telescoping certificate, including the Gaussian terms.
    rhs = (t * b - X) * (-b) + (t * c - Y) * (-c)
    assert S.expand(P - a - v * t - rhs) == 0
    original = [-a - v * t, S.Integer(1)]
    full = difference_operator(original, t)
    series = S.exp(a * t + v * t ** 2 / 2) - 1
    assert S.simplify(sum(q * S.diff(series, t, j) for j, q in enumerate(full))) == 0

    tree = []
    moments = [S.Integer(0), a, a ** 2 + v]
    for constraints, expected_bound, first in [([], 1, 1), ([a], 2, 2), ([a, v], 0, None)]:
        ideal = S.groebner(constraints, *params, domain=S.QQ)

        def reduce_poly(expr):
            expr = S.Poly(S.expand(expr), t)
            return S.expand(sum(ideal.reduce(coef)[1] * t ** mon[0] for mon, coef in expr.terms()))

        # Recompute at this prime node, before valuations or integer-root tests.
        reduced = [reduce_poly(q) for q in original]
        operator = difference_operator(reduced, t)
        mu, poly = indicial(operator, t, z)
        roots = constant_integer_roots(poly, z, params)
        bound = max(roots) if roots else -1
        assert bound == expected_bound
        nonzero = [n for n in range(1, bound + 1) if ideal.reduce(moments[n])[1] != 0]
        assert (nonzero[0] if nonzero else None) == first
        tree.append({'prime_generators': [str(x) for x in constraints],
                     'indicial': str(S.factor(poly)), 'nonnegative_integer_roots': roots,
                     'bound': bound, 'first_nonzero_moment': first,
                     'accepted': first is None})

    # A generic bound transported to the exceptional child would be wrong.
    assert moments[1].subs({a: 0, b: 1, c: 0}) == 0
    assert moments[2].subs({a: 0, b: 1, c: 0}) == 1
    rejected_zero_operator = False
    try:
        indicial([S.Integer(0)], t, z)
    except ValueError:
        rejected_zero_operator = True
    assert rejected_zero_operator
    return {'completed': True, 'telescoping_certificate': True,
            'difference_from_one_operator': True, 'linear_family_descent': tree,
            'null_radical': ['a', 'b**2+c**2'], 'least_linear_prefix': 2,
            'generic_bound_specialization_mutation_detected': True,
            'zero_operator_rejected': True,
            'general_algorithm_implemented': False, 'general_algorithm_executed': False,
            'universal_termination_proof_certified_by_checker': False,
            'full_cubic_coefficient_search_executed': False,
            'producer_package_executed': False,
            'independent_replay': False, 'external_review': False, 'formal_verification': False}


if __name__ == '__main__':
    result = run()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'generic_point_moment_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
