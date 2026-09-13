#!/usr/bin/env python3
"""Verify the literal NS0031 stable-group counter-witness; standard library only.

No lattice search or period-curve construction. The original finite modular
arithmetic is replayed separately by verify_ns0031_marking_arithmetic.py.
This checks an integral discriminant-kernel reflection whose projective spin
representative has nonsquare determinant, invalidating the old containment in
rational norm-one units. It produces no rational K3 or nonexistence theorem.
"""
from fractions import Fraction as F
from math import isqrt
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WITNESS = ROOT/'artifacts/generated-results/elkies-k3-ns0031-period-group-counterwitness-v1.json'
ORIGINAL = ROOT/'artifacts/generated-results/elkies-k3-ns0031-qq-marking-obstruction-v1.json'


def require(test, message):
    if not test:
        raise ValueError(message)


def transpose(a):
    return list(map(list, zip(*a)))


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def determinant(a):
    if len(a) == 2:
        return a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return sum((-1)**j*a[0][j]*determinant([row[:j]+row[j+1:] for row in a[1:]])
               for j in range(len(a)))


def combination(basis, coefficients):
    return [[sum(c*a[i][j] for c, a in zip(coefficients, basis))
             for j in range(len(basis[0][0]))] for i in range(len(basis[0]))]


def verify(witness, original):
    require(witness['schema'] == 'elkies-k3.ns0031-period-group-counterwitness.v1', 'wrong witness schema')
    require(witness['status'] == 'PASS_STABLE_REFLECTION_OUTSIDE_RATIONAL_NORM_ONE_IMAGE', 'wrong witness status')
    require(witness['surface_id'] == original['surface_id'] == 'K3-d1b1381f87d69f1c', 'wrong surface')
    g = witness['transcendental_gram']
    require(g == original['transcendental_lattice']['gram'] == [[0,0,4],[0,74,1],[4,1,-2]],
            'wrong literal transcendental lattice')
    v = witness['negative_root']
    require(v == [0,0,1], 'wrong negative root')
    gv = mul(g, [[x] for x in v])
    norm = sum(x*y[0] for x, y in zip(v, gv))
    require(norm == -2, 'root does not have norm minus two')
    identity = [[int(i == j) for j in range(3)] for i in range(3)]
    reflection = [[identity[i][j]-F(2*v[i]*gv[j][0], norm)
                   for j in range(3)] for i in range(3)]
    require(reflection == witness['reflection'], 'reflection formula mismatch')
    require(all(x.denominator == 1 for row in reflection for x in map(F, row)), 'nonintegral reflection')
    require(determinant(reflection) == -1 and mul(reflection, reflection) == identity, 'wrong reflection determinant/order')
    require(mul(mul(transpose(reflection), g), reflection) == g, 'not an isometry')
    # Columns of G^-1 generate T^dual. (R-I)G^-1 integral means identity on A_T.
    difference = [[reflection[i][j]-identity[i][j] for j in range(3)] for i in range(3)]
    kernel = witness['discriminant_kernel_integral_matrix']
    require(all(type(x) is int for row in kernel for x in row), 'nonintegral discriminant action witness')
    require(mul(kernel, g) == difference, 'discriminant-kernel identity failed')
    plane = witness['fixed_positive_plane_basis']
    require(mul(plane, transpose(reflection)) == plane, 'positive plane is not fixed pointwise')
    plane_gram = mul(mul(plane, g), transpose(plane))
    require(plane_gram == witness['fixed_positive_plane_gram'], 'positive-plane Gram mismatch')
    require(plane_gram[0][0] > 0 and determinant(plane_gram) > 0, 'plane is not positive definite')

    order = original['even_clifford_order']['integral_basis_after_diag_4_1_conjugation']
    require(order == [[[1,0],[0,1]],[[0,4],[0,0]],[[4,0],[0,0]],[[1,1],[37,0]]], 'wrong Clifford embedding')
    # z=e0 wedge e1 wedge e2 has z^2=148. Compute ei*z in the old even basis.
    b = [combination(order, c) for c in [[0,2,0,0],[74,F(1,2),-37,0],[0,-1,F(-1,2),2]]]
    require(b == witness['trace_zero_clifford_images'], 'Clifford-vector identification mismatch')
    for i in range(3):
        for j in range(3):
            left, right = mul(b[i], b[j]), mul(b[j], b[i])
            require([[left[a][c]+right[a][c] for c in range(2)] for a in range(2)]
                    == [[148*g[i][j]*int(a==c) for c in range(2)] for a in range(2)],
                    'Clifford anticommutator mismatch')
    require(determinant([[x[0][0],x[0][1],x[1][0]] for x in b]) != 0,
            'Clifford vectors do not span the trace-zero algebra')
    a = witness['projective_spin_matrix']
    det = determinant(a)
    require(det == witness['projective_spin_determinant'] == 37, 'wrong spin determinant')
    inverse = [[F(a[1][1],det),F(-a[0][1],det)],[F(-a[1][0],det),F(a[0][0],det)]]
    for j in range(3):
        require(mul(mul(a,b[j]),inverse) == combination(b,[-reflection[i][j] for i in range(3)]),
                'adjoint action differs from the projective reflection')
    # Spanning all trace-zero matrices makes any other rational representative
    # a scalar multiple of A. Its determinant is 37*c^2, never 1 for rational c.
    require(det > 0 and isqrt(det)**2 != det, 'determinant is a rational square')
    coordinates = witness['order_conjugation_coordinates']
    require(all(type(x) is int for row in coordinates for x in row), 'nonintegral order transport')
    for j in range(4):
        require(mul(mul(a,order[j]),inverse) == combination(order,[row[j] for row in coordinates]),
                'order normalizer identity failed')
    require(mul(coordinates,coordinates) == [[int(i==j) for j in range(4)] for i in range(4)],
            'order transport is not an integral involution')
    return 'PASS NS0031 stable reflection outside rational norm-one image; rational K3 existence remains UNKNOWN.'


def main():
    witness = json.loads(WITNESS.read_text())
    raw = ORIGINAL.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == witness['inputs'][str(ORIGINAL.relative_to(ROOT))],
            'original arithmetic certificate changed')
    print(verify(witness, json.loads(raw)))


if __name__ == '__main__':
    main()
