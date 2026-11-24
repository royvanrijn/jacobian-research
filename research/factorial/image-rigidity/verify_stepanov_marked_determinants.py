#!/usr/bin/env python3
"""Exact controls for factorial marked-rank and determinant mechanisms.

Standard library only. The all-parameter arguments are in
STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md;
these finite controls verify formulas, inverse markers and explicit countercontrols.
"""
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial, prod
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("Assertions are required; do not run with python -O")


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p or [F(0)]


def add(p, q):
    return trim([(p[i] if i < len(p) else 0) +
                 (q[i] if i < len(q) else 0)
                 for i in range(max(len(p), len(q)))])


def scale(p, c):
    return trim([F(c) * x for x in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def evaluate(p, x):
    value = F(0)
    for c in reversed(p):
        value = value*x+c
    return value


def rising_polynomial(start, count):
    p = [F(1)]
    for j in range(count):
        p = mul(p, [start+j, 1])
    return p


def determinant(matrix):
    """Fraction elimination, with an independent permutation control below."""
    a = [list(map(F, row)) for row in matrix]
    value = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            value = -value
        value *= a[k][k]
        for i in range(k+1, len(a)):
            ratio = a[i][k]/a[k][k]
            for j in range(k+1, len(a)):
                a[i][j] -= ratio*a[k][j]
    return value


def polynomial_determinant(matrix):
    n = len(matrix)
    result = [F(0)]
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        term = [F((-1)**inversions)]
        for i in range(n):
            term = mul(term, matrix[i][p[i]])
        result = add(result, term)
    return result


def alpha_factorial(a):
    return prod(factorial(x) for x in a)


def exponent(alpha, beta, n, k):
    return tuple((n-k)*a+k*b for a, b in zip(alpha, beta))


def moment(alpha, beta, n, marker=None):
    marker = marker or (0,)*len(alpha)
    return [F(comb(n, k)*alpha_factorial(tuple(x+y for x, y in
            zip(exponent(alpha, beta, n, k), marker)))) for k in range(n+1)]


def direct_moment(alpha, beta, n, marker=None):
    """Independent repeated multiplication in the original monomial ring."""
    d = len(alpha)
    marker = marker or (0,)*d
    phase = {tuple(alpha)+(0,): 1, tuple(beta)+(1,): 1}
    power = {(0,)*(d+1): 1}
    for _ in range(n):
        out = {}
        for a, ac in power.items():
            for b, bc in phase.items():
                e = tuple(x+y for x, y in zip(a, b))
                out[e] = out.get(e, 0)+ac*bc
        power = out
    result = [F(0)]*(n+1)
    for e, c in power.items():
        result[e[-1]] += c*alpha_factorial(tuple(x+y for x, y in zip(e[:-1], marker)))
    return trim(result)


def dual_marker(nodes, k):
    """Convert the k-th Lagrange polynomial to the rising basis (z+1)_r."""
    p = [F(1)]
    for j, node in enumerate(nodes):
        if j != k:
            p = scale(mul(p, [-node, 1]), F(1, nodes[k]-node))
    remainder = p[:]
    q = [F(0)]*len(nodes)
    for r in reversed(range(len(nodes))):
        coefficient = remainder[r] if r < len(remainder) else F(0)
        q[r] = coefficient
        remainder = add(remainder, scale(rising_polynomial(1, r), -coefficient))
    assert remainder == [0]
    assert all(evaluate(p, z) == int(j == k) for j, z in enumerate(nodes))
    return trim(q)


def marked_operator(p, alpha, beta, n, marker):
    return trim([c*prod(prod(n*a+1+(b-a)*k+j for j in range(eta))
                        for a, b, eta in zip(alpha, beta, marker))
                 for k, c in enumerate(p)])


def gaussian_add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def gaussian_mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gaussian_power_moments(poly, top):
    power = [(F(1), F(0))]
    moments = []
    for _ in range(top):
        out = [(F(0), F(0)) for _ in range(len(power)+len(poly)-1)]
        for i, a in enumerate(power):
            for j, b in enumerate(poly):
                out[i+j] = gaussian_add(out[i+j], gaussian_mul(a, b))
        power = out
        moments.append(tuple(sum(c[j]*factorial(i) for i, c in enumerate(power))
                             for j in range(2)))
    return moments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('stepanov_marked_determinants_verification.json'))
    args = parser.parse_args()
    cases = [((0,), (1,)), ((0,), (2,)), ((1,), (2,)), ((2,), (5,)),
             ((3,), (7,)), ((1, 0), (0, 1)), ((0, 2), (3, 0)),
             ((2, 2), (2, 5)), ((1, 0, 2), (0, 2, 1)),
             ((0, 0, 0), (1, 1, 1))]
    checks = dict(vandermonde=0, direct_contractions=0, dual_markers=0,
                  marked_operator=0, raising_recurrence=0, factorial_gram=0)
    sample_markers = []
    for alpha, beta in cases:
        i = next(j for j in range(len(alpha)) if alpha[j] != beta[j])
        delta = beta[i]-alpha[i]
        for n in range(1, 6):
            exponents = [exponent(alpha, beta, n, k) for k in range(n+1)]
            nodes = [e[i] for e in exponents]
            matrix = [[alpha_factorial(tuple(x+(r if j == i else 0)
                                            for j, x in enumerate(e)))
                       for e in exponents] for r in range(n+1)]
            expected = (prod(alpha_factorial(e) for e in exponents)
                        * delta**(n*(n+1)//2)*prod(factorial(k) for k in range(1, n+1)))
            assert determinant(matrix) == expected != 0
            checks['vandermonde'] += 1
            p = moment(alpha, beta, n)
            for r in range(n+1):
                marker = tuple(r if j == i else 0 for j in range(len(alpha)))
                actual = moment(alpha, beta, n, marker)
                assert actual == direct_moment(alpha, beta, n, marker)
                checks['direct_contractions'] += 1
                assert actual == marked_operator(p, alpha, beta, n, marker)
                checks['marked_operator'] += 1
            for k in (0, n):
                q = dual_marker(nodes, k)
                contraction = [sum(q[r]*matrix[r][j] for r in range(len(q)))
                               for j in range(n+1)]
                assert contraction == [alpha_factorial(e) if j == k else 0
                                       for j, e in enumerate(exponents)]
                checks['dual_markers'] += 1
                if alpha == (1,) and beta == (2,) and n == 2:
                    sample_markers.append({'k': k, 'coefficients': list(map(str, q))})
            raising = add(marked_operator(p, alpha, beta, n, alpha),
                          [0]+marked_operator(p, alpha, beta, n, beta))
            assert raising == moment(alpha, beta, n+1)
            checks['raising_recurrence'] += 1
    sets = [[(i,) for i in range(k+1)] for k in range(7)]
    sets += [[a for a in product(range(k+1), repeat=2) if sum(a) <= k]
             for k in range(4)]
    for S in sets:
        S = sorted(S, key=lambda a: (sum(a), a))
        H = [[alpha_factorial(tuple(x+y for x, y in zip(a, b))) for b in S] for a in S]
        expected = prod(alpha_factorial(a)**2 for a in S)
        assert determinant(H) == expected
        if len(S) <= 4:
            assert polynomial_determinant([[[x] for x in row] for row in H]) == [expected]
        checks['factorial_gram'] += 1
    p = [moment((1,), (2,), n) for n in range(5)]
    H0 = add(mul(p[0], p[2]), scale(mul(p[1], p[1]), -1))
    H1 = add(mul(p[1], p[3]), scale(mul(p[2], p[2]), -1))
    assert H0 == [1, 8, 20]
    assert H1 == scale(mul([1, 6], [1, 12, 60, 72]), 2)
    assert [evaluate(q, F(-1, 6)) for q in p[1:4]] == [F(2, 3)]*3
    spatial_H1 = add(mul(moment((1,), (2,), 1, (0,)), moment((1,), (2,), 1, (2,))),
                      scale(mul(moment((1,), (2,), 1, (1,)), moment((1,), (2,), 1, (1,))), -1))
    assert spatial_H1 == [2, 12, 12]
    f = [(F(2), F(-2)), (F(-4), F(2)), (F(1), F(0))]
    pure = gaussian_power_moments(f, 3)
    assert pure == [(0, 0), (0, 0), (32, 80)]
    marks = [tuple(sum(c[j]*factorial(i+r) for i, c in enumerate(f))
                   for j in range(2)) for r in range(3)]
    assert marks == [(0, 0), (0, 2), (4, 8)]
    marked_det = gaussian_add(gaussian_mul(marks[0], marks[2]),
                              tuple(-x for x in gaussian_mul(marks[1], marks[1])))
    assert marked_det == (4, 0)
    # An operator/degree constraint does not supply actual moment seeds.
    abstract_p = [11, -20, 9]
    abstract_next = add(marked_operator(abstract_p, (1,), (2,), 2, (1,)),
                        [0]+marked_operator(abstract_p, (1,), (2,), 2, (2,)))
    assert evaluate(abstract_p, 1) == evaluate(abstract_next, 1) == 0
    assert -20+2*9 == -2 != 0
    # Corruption controls: both the determinant factor and the inverse marker matter.
    spaced_matrix = [[factorial(e+r) for e in (2, 4, 6)] for r in range(3)]
    omitted_spacing_factor = prod(factorial(e) for e in (2, 4, 6))*factorial(2)
    assert determinant(spaced_matrix) == 8*omitted_spacing_factor
    assert determinant(spaced_matrix) != omitted_spacing_factor
    bad_marker = dual_marker([2, 3, 4], 0)
    bad_marker[0] += 1
    assert sum(bad_marker[r]*factorial(3+r) for r in range(len(bad_marker))) != 0
    receipt = {
        'passed': True,
        'scope': 'Exact bounded controls plus explicit counterexamples only; universal algebraic proofs are in STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md.',
        'coefficient_domains': ['Q', 'Q[t]', 'Q(i), i^2=-1'],
        'supports': [{'alpha': a, 'beta': b} for a, b in cases],
        'power_range': [1, 5], 'checks': checks,
        'sample_dual_markers_for_alpha1_beta2_n2': sample_markers,
        'countercontrols': {
            'pure_H0_X_plus_tX2': list(map(str, H0)),
            'pure_H1_X_plus_tX2': list(map(str, H1)),
            'rational_H1_zero': {'t': '-1/6', 'moments_1_2_3': ['2/3']*3},
            'spatial_H1_X_plus_tX2': list(map(str, spatial_H1)),
            'three_term_two_zero_phase': 'X^2+(-4+2i)X+2-2i',
            'three_term_pure_moments_1_2_3': [[str(x) for x in q] for q in pure],
            'three_term_spatial_markers_0_1_2': [[str(x) for x in q] for q in marks],
            'three_term_spatial_2x2_determinant': '4',
            'abstract_operator_seed': abstract_p,
            'abstract_operator_next': list(map(str, abstract_next)),
            'abstract_common_root': '1', 'abstract_seed_derivative_at_root': '-2'
        },
        'limitations': ['No SFC counterexample is claimed.',
                        'The three-term control lies outside SFC(2).',
                        'The abstract operator seed is not the factorial-moment seed.',
                        'No propagation from two scalar zero moments to n+1 marked zero moments is proved.',
                        'No positivity, Riemann hypothesis or analytic asymptotic input is used.']
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'passed': True, 'checks': checks, 'output': str(args.output)}))


if __name__ == '__main__':
    main()
