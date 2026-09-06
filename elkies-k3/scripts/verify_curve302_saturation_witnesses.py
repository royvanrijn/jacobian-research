#!/usr/bin/env python3
"""Independent standard-library replay of the finite saturation witnesses.

This verifies the lattice argument, not the cubic geometry or section maps;
those are separate dependencies checked by the Sage construction replay.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CERTIFICATE = ROOT / 'artifacts/generated-results/elkies-k3-curve302-section-saturation-v1.json'


def determinant(a):
    a = [[Fraction(c) for c in r] for r in a]
    result = Fraction(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return Fraction(0)
        if i != j:
            a[i], a[j] = a[j], a[i]
            result = -result
        pivot = a[i][i]
        result *= pivot
        for j in range(i+1, len(a)):
            scale = a[j][i]/pivot
            a[j] = [x-scale*y for x, y in zip(a[j], a[i])]
    return result


def pairing(a, v, w):
    return sum(v[i]*a[i][j]*w[j] for i in range(len(a)) for j in range(len(a)))


def verify(certificate):
    for path, expected in certificate['input_sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == expected
    baseline = json.loads((ROOT/next(iter(certificate['input_sha256']))).read_text())
    G = baseline['generic_section_lattice']['height_gram']
    H = certificate['K3']['height_gram']
    T = [[Fraction(c) for c in r] for r in certificate['K3']['new_in_old_rational_matrix']]
    assert all(pairing(G, T[i], T[j]) == H[i][j] for i in range(9) for j in range(9))
    assert determinant(T) == Fraction(1, 3) and determinant(H) == 512
    assert all(determinant([r[:i] for r in H[:i]]) > 0 for i in range(1, 10))
    A = certificate['source_basis']['height_gram']
    assert all(2*A[i][j] == H[i][j] for i in range(8) for j in range(8))
    assert determinant(A) == 1 and all(A[i][i] % 2 == 0 for i in range(8))
    old_image = baseline['specialization_overlap']['matrix_9_by_31_rows']
    image = certificate['K3']['specialization_matrix']
    assert all(sum(T[i][j]*old_image[j][k] for j in range(9)) == image[i][k]
               for i in range(9) for k in range(31))
    assert abs(determinant([r[:9] for r in image])) == 3
    assert all(not any(r[9:]) for r in image)
    # This visible unit minor plus determinant 3 proves the stated Smith form.
    assert abs(determinant([r[1:9] for r in image[:7]+[image[8]]])) == 1
    possibilities = set()
    for v in product((0, 1), repeat=9):
        if not any(v):
            continue
        if all(sum(H[i][j]*v[j] for j in range(9)) % 2 == 0 for i in range(9)):
            if pairing(H, v, v) % 8 == 0:
                possibilities.add(v)
    supplied = certificate['saturation_certificate']['height_two_witnesses']
    assert len(supplied) == len(possibilities) == 71
    assert {tuple(r['parity']) for r in supplied} == possibilities
    for row in supplied:
        w = row['twice_vector']
        assert len(w) == 9 and all(type(x) is int for x in w)
        assert tuple(x % 2 for x in w) == tuple(row['parity'])
        assert pairing(H, w, w) == 8


if __name__ == '__main__':
    verify(json.loads(CERTIFICATE.read_text()))
    print('PASS: index 3, source determinant 1, K3 determinant 512, all 71 root-forcing cosets, specialization Smith factors 1^8,3')
