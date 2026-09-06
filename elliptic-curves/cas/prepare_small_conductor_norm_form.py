#!/usr/bin/env python3
"""Exact maximal-order binary cubic and Hessian reduction for a norm sieve.

Only four coefficients and exact generator maps are computed. No factor base,
class relations, Selmer classes, or point searches are enumerated.
"""
import argparse
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path

import certify_small_conductor_class_rank_target as target

ROOT, ART = target.ROOT, target.ART
D = ROOT / 'artifacts/local/elliptic-curves/small-conductor-descent-profile-v1'
OUT = ART / 'small_conductor_reduced_norm_form_v1.json'


def inverse(matrix):
    n = len(matrix)
    rows = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        k = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[k] = rows[k], rows[j]
        pivot = rows[j][j]
        rows[j] = [x/pivot for x in rows[j]]
        for i in range(n):
            if i != j:
                q = rows[i][j]
                rows[i] = [x-q*y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def determinant(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def discriminant(c):
    a, b, d, e = c
    return b*b*d*d - 4*a*d**3 - 4*b**3*e - 27*a*a*e*e + 18*a*b*d*e


def hessian(c):
    a, b, d, e = c
    return [b*b-3*a*d, b*d-9*a*e, d*d-3*b*e]


def multiply_matrix(m, n):
    return [m[0]*n[0]+m[1]*n[2], m[0]*n[1]+m[1]*n[3],
            m[2]*n[0]+m[3]*n[2], m[2]*n[1]+m[3]*n[3]]


def transform(c, m):
    out = [0]*4
    for i, coefficient in enumerate(c):
        for j in range(4-i):
            for k in range(i+1):
                out[j+k] += (coefficient*comb(3-i, j)*comb(i, k)
                             *m[0]**(3-i-j)*m[1]**j*m[2]**(i-k)*m[3]**k)
    return out


def reduce_form(c):
    matrix = [1, 0, 0, 1]
    current = c
    history = []
    for _ in range(1000):
        a, b, d = hessian(current)
        if a <= 0 or 4*a*d-b*b <= 0:
            raise ArithmeticError('positive definite Hessian required')
        shift = (a-b)//(2*a)
        if shift:
            m = [1, shift, 0, 1]
            current = transform(current, m)
            matrix = multiply_matrix(matrix, m)
            history.append(m)
        a, b, d = hessian(current)
        if a <= d:
            if not -a < b <= a:
                raise ArithmeticError('quadratic reduction interval failed')
            return current, matrix, history
        m = [0, -1, 1, 0]
        current = transform(current, m)
        matrix = multiply_matrix(matrix, m)
        history.append(m)
    raise ArithmeticError('finite reduction limit exceeded')


def expected():
    proof = target.expected()
    polynomial = list(map(F, proof['two_division_cubic_ascending']))
    basis = [list(map(F, r)) for r in target.original.cert.read(D / 'basis_export.json')]
    one = [F(1), F(0), F(0)]
    def multiply(u, v):
        c = [sum((u[i]*v[j] for i in range(3) for j in range(3) if i+j == k), F(0)) for k in range(5)]
        for k in (4, 3):
            q = c[k]
            for i in range(3):
                c[k-3+i] -= q*polynomial[i]
        return c[:3]
    def coords(u, bs):
        inv = inverse([list(row) for row in zip(*bs)])
        return [sum(x*y for x, y in zip(row, u)) for row in inv]
    def product(u, v, bs):
        return coords(multiply(u, v), bs)
    pair = product(basis[1], basis[2], basis)
    w = [basis[1][i]-pair[2]*one[i] for i in range(3)]
    t = [basis[2][i]-pair[1]*one[i] for i in range(3)]
    normal = [one, w, t]
    table = [[product(u, v, normal) for v in normal] for u in normal]
    if any(x.denominator != 1 for rows in table for row in rows for x in row):
        raise ArithmeticError('ring basis multiplication not integral')
    a, b, c, d = -table[1][1][2], table[1][1][1], -table[2][2][2], table[2][2][1]
    if table[1][2] != [-a*d, 0, 0] or table[1][1] != [-a*c, b, -a] or table[2][2] != [-b*d, d, -c]:
        raise ArithmeticError('binary cubic ring identities failed')
    form = list(map(int, [a, b, c, d]))
    reduced, matrix, history = reduce_form(form)
    if matrix[0]*matrix[3]-matrix[1]*matrix[2] != 1 or transform(form, matrix) != reduced:
        raise ArithmeticError('exact determinant-one transport failed')
    index = 1/determinant([list(row) for row in zip(*basis)])
    field_disc = F(proof['cubic_discriminant'])/index**2
    if discriminant(form) != field_disc or discriminant(reduced) != field_disc:
        raise ArithmeticError('ring/form discriminants differ')
    # Verify the complete homogeneous norm identity through four coefficients.
    # Norm(a*m+n*w)=a^2 F(m,n), using determinant of multiplication by a*m+n*w.
    for m, n in [(0, 1), (1, 1), (-1, 1), (2, 1), (1, 0)]:
        u, v = matrix[0]*m+matrix[1]*n, matrix[2]*m+matrix[3]*n
        element = [a*u*one[i]+v*w[i] for i in range(3)]
        norm = determinant([list(row) for row in zip(*(multiply(element, e) for e in [one, [0, 1, 0], [0, 0, 1]]))])
        value = sum(coef*m**(3-i)*n**i for i, coef in enumerate(reduced))
        if norm != a*a*value:
            raise ArithmeticError('norm-polynomial identity failed')
    original = list(map(int, reversed(polynomial)))
    return {
        'schema': 'elliptic-curves.small-conductor-reduced-norm-form.v1', 'status': 'PASS',
        'sources': {str(p.relative_to(ROOT)): target.original.cert.hashed(p)
                    for p in [Path(__file__).resolve(), target.OUTPUT, D / 'basis_export.json', D / 'field.log']},
        'original_monic_cubic_descending': list(map(str, original)),
        'maximal_order_basis_ascending': [list(map(str, row)) for row in basis],
        'maximality_scope': 'Maximality was checked separately by PARI nfcertify in field.log using complete proved discriminant support. This checker independently verifies the ring multiplication, index and discriminant identities; those identities alone do not establish maximality.',
        'normal_basis_ascending': [list(map(str, row)) for row in normal],
        'multiplication_table': [[[str(x) for x in row] for row in rows] for rows in table],
        'initial_binary_cubic_descending': list(map(str, form)),
        'reduced_binary_cubic_descending': list(map(str, reduced)),
        'sl2_matrix': [[matrix[0], matrix[1]], [matrix[2], matrix[3]]],
        'reduction_steps': history, 'reduced_hessian': list(map(str, hessian(reduced))),
        'defining_order_index': str(index), 'field_discriminant': str(field_disc),
        'integral_norm_generator': {'fixed_a': str(a), 'w_power_basis': list(map(str, w)),
                                    'formula': 'u=M00*m+M01*n; v=M10*m+M11*n; beta=a*u+v*w; Norm(beta)=a^2*F_reduced(m,n).'},
        'maximum_coefficient_bits': {'original_monic': max(abs(v).bit_length() for v in original),
                                     'reduced_binary': max(abs(v).bit_length() for v in reduced)},
        'box10_polynomial_value_upper_bound_bits': {'original_monic': (1000*sum(map(abs, original))).bit_length(),
                                                   'reduced_binary': (1000*sum(map(abs, reduced))).bit_length()},
        'claim_boundary': 'Exact norm-sieve coordinate preparation on one field. The smaller varying polynomial is accompanied by the exact fixed norm square factor a^2. No smoothness rate, nontrivial class relation, certified factor-base generation, class-rank upper bound, Selmer bound, new point or end-to-end speedup has been proved.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if json.loads(json.dumps(data)) != target.original.cert.read(OUT):
            raise ArithmeticError('norm-form certificate differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve norm-form preparation')
        target.original.cert.write(OUT, data)
    print('EXACT NORM FORM', data['reduced_binary_cubic_descending'])
    print('COEFFICIENT BITS', data['maximum_coefficient_bits'], 'BOX10 BITS', data['box10_polynomial_value_upper_bound_bits'])
