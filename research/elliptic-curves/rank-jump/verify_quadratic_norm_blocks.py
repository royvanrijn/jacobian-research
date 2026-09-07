#!/usr/bin/env python3
"""Independent polynomial, permutation and symplectic checks for norm blocks."""
import argparse
from pathlib import Path
import retrospective as r
import quadratic_norm_blocks as source

OUTPUT = r.OUT/'rank_jump_quadratic_norm_block_verification_v1.json'


def compute():
    from sage.all import QQ, GF, PolynomialRing, NumberField, matrix, identity_matrix, GL
    data, inp = r.read(source.OUTPUT), r.read(source.INPUT)
    for name, digest in data['bindings'].items():
        assert r.digest((r.ROOT/name).read_bytes()) == digest
    R = PolynomialRing(QQ, 't')
    t = R.gen()
    K = NumberField(t*t+1, 'ii')
    ii = K.gen()
    P = PolynomialRing(K, names=('A', 'B', 'a', 'b', 'c', 'd', 'x', 'y'))
    A, B, a, b, c, d, x, y = P.gens()
    E = lambda X, Y: Y*Y-X**3-A*X-B
    Et = y*y-x**3-A*x+B
    assert E(-x, ii*y) == -Et
    real = c*c-d*d-a**3+3*a*b*b-A*a-B
    imag = 2*c*d-3*a*a*b+b**3-A*b
    assert E(a+ii*b, c+ii*d) == real+ii*imag

    M = matrix(QQ, [[1, 1], [1, -1]])
    swap = matrix(QQ, [[0, 1], [1, 0]])
    descent = matrix(QQ, [[1, 0], [0, -1]])
    assert M*M == 2*identity_matrix(QQ, 2)
    assert M.transpose()*M == 2*identity_matrix(QQ, 2)
    assert swap*M == M*descent and M.det()**2 == 4

    F = GF(2)
    G = GL(2, F)
    cycle, transposition = G([0, 1, 1, 1]), G([0, 1, 1, 0])
    assert cycle.order() == 3 and transposition.order() == 2
    centralizers = {}
    for label, generators in [('S3', [cycle, transposition]), ('C3', [cycle])]:
        C = sorted([list(map(int, g.matrix().list())) for g in G
                    if all(g*h == h*g for h in generators)])
        assert C == sorted(data['GL2_centralizers'][label])
        centralizers[label] = len(C)

    rows = []
    for old, out in zip(inp['rows'], data['rows'], strict=True):
        f = R(list(map(QQ, old['cubic_ascending'])))
        cert = out['cubic_certificate']
        assert f.discriminant() == QQ(cert['discriminant'])
        assert not f.discriminant().is_square()
        assert f.change_ring(GF(cert['irreducible_reduction_prime'])).is_irreducible()
        mat = matrix(F, old['CT_matrix'])
        flat = out['nondegenerate_block_masks'] + [v['strict_basis_mask'] for v in old['radical_basis']]
        T = matrix(F, [[v >> i & 1 for v in flat] for i in range(mat.nrows())])
        assert T.is_invertible() and mat == mat.transpose()
        assert all(mat[i, i] == 0 for i in range(mat.nrows()))
        rank = out['CT_rank']
        expected = matrix(F, [[int(i < rank and j < rank and i//2 == j//2 and i != j)
                              for j in range(mat.nrows())] for i in range(mat.nrows())])
        assert T.transpose()*mat*T == expected and mat.rank() == rank
        assert out['local_global_elliptic_norm_defect_order_at_least'] == 2**rank
        rows.append({'case_index': out['case_index'], 'independent_CT_rank': int(mat.rank()),
                     'S3_verified': True, 'all_block_masks_verified': True})
    small = R(inp['small_control_cubic_ascending'])
    assert small.discriminant() == 163**2
    assert small.change_ring(GF(data['small_control']['irreducible_reduction_prime'])).is_irreducible()
    from sage.all import EllipticCurve
    assert EllipticCurve([0, -11, 0, -14, -1]).j_invariant() not in [0, 1728]
    return {
        'schema': 'rank-jump.quadratic-norm-block-verification.v1',
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (Path(__file__), source.OUTPUT)},
        'status': 'PASS', 'twist_isomorphism_identity': True, 'Weil_restriction_equations': True,
        'isogeny_degree': 4, 'descent_equivariance': True, 'polarization_pullback_multiplier': 2,
        'centralizer_sizes': centralizers, 'rows': rows,
        'boundary': 'No norm preimages, new CT entries or full Sha groups computed. The norm-defect identity follows from the written cohomology argument.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['build', 'check'])
    args = p.parse_args()
    data = compute()
    if args.mode == 'check':
        assert r.read(OUTPUT) == data
        print('PASS independent quadratic norm-block verification')
    else:
        r.write_new(OUTPUT, data)
        print(data['status'])
