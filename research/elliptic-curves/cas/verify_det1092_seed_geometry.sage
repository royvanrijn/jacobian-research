#!/usr/bin/env sage-python
"""Generic-only identities for the seed-density and limiting-lattice theorem.

No specialization, exceptional point, V3 landscape, CVP enumeration or search.
Bound: one 18-dimensional Gram, one 4-dimensional resultant matrix, and
polynomial identities of degree at most four (apart from the existing j map).
Run under timeout 25s. --write creates an immutable arithmetic certificate;
the default rederives it and checks equality. Analytic deductions are proved
in the accompanying note, not asserted to be proved by this script.
"""
import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, EllipticCurve, PolynomialRing, block_diagonal_matrix, gcd, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_seed_geometry_v1.json'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
COVER = ART / 'det1092_orbit8044_rank18_base_change_v2.json'
REPLAY = ART / 'det1092_orbit8044_rank18_base_change_replay_v2.json'
CHART = ART / 'orbit8044_seed_factory_v1/parameter-chart.json'


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(A):
    return [list(map(str, row)) for row in A.rows()]


def verify():
    parent, cover, replay, chart = map(read, [PARENT, COVER, REPLAY, CHART])
    require(cover['status'] == 'PASS_EXPLICIT_Q_RATIONAL_BASE_CHANGE_RANK_AT_LEAST_18', 'cover status')
    require(replay['status'] == 'PASS_INDEPENDENT_COEFFICIENT_AND_INTERSECTION_HEIGHT_REPLAY', 'cover replay status')
    require(replay['inputs'][str(COVER.relative_to(ROOT))] == sha(COVER), 'cover replay binding')
    require(chart['status'] == 'PASS_EXACT_ORBIT8044_PARAMETER_CHART', 'factory chart status')
    require(chart['inputs'][str(COVER.relative_to(ROOT))] == sha(COVER), 'factory cover binding')

    G = matrix(QQ, parent['generic_height_gram'])
    w = vector(QQ, cover['lift']['trace_word'])
    require(G.nrows() == G.ncols() == 17 and G.is_symmetric(), 'generic Gram size/symmetry')
    require(G.det() == 1092 and w * G * w == 10, 'generic determinant/trace norm')
    g_minors = [G[:i, :i].det() for i in range(1, 18)]
    require(all(d > 0 for d in g_minors), 'generic positivity')
    cross = G*w
    H = (2*G).augment(matrix(QQ, 17, 1, list(cross)))
    H = H.stack(matrix(QQ, 1, 18, list(cross) + [8]))
    require(H == matrix(QQ, cover['independence']['gram']), 'saved cover Gram')
    T = matrix.identity(QQ, 18)
    for i in range(17):
        T[i, 17] = w[i]/2
    D = block_diagonal_matrix(2*G, matrix(QQ, [[3]]))
    require(H == T.transpose()*D*T, 'complete-square identity')
    schur = H[17, 17] - cross*(2*G).inverse()*cross
    median = sorted(H[i, i] for i in range(17))[8]
    require(schur == 3 and median == 8 and H.det() == 429391872, 'Schur/scale/determinant')

    P = PolynomialRing(ZZ, names=('a', 'b'))
    a, b = P.gens()
    def form(raw, aa=a, bb=b):
        c0, c1, c2 = map(ZZ, raw)
        return c0*bb**2 + c1*aa*bb + c2*aa**2
    F, Y, Z = [form(chart[key]) for key in ['s_numerator', 's_denominator', 'z_numerator']]
    q0, q1, q2 = map(ZZ, chart['primitive_conic'])
    require(Z**2 == q0*Y**2 + q1*F*Y + q2*F**2, 'conic identity')
    require(q1**2 - 4*q0*q2 != 0, 'smooth conic')
    # Column convention: [aF,bF,aY,bY] = [a^3,a^2b,ab^2,b^3] * M.
    monomials = [a**3, a*a*b, a*b*b, b**3]
    columns = [a*F, b*F, a*Y, b*Y]
    M = matrix(ZZ, [[v.monomial_coefficient(m) for v in columns] for m in monomials])
    resultant = M.det()
    require(resultant != 0, 'degree-two morphism must be basepoint free')
    adj = M.adjugate()
    require(M*adj == resultant*matrix.identity(ZZ, 4), 'adjugate identity')
    for j in [0, 3]:
        require(sum(columns[i]*adj[i,j] for i in range(4)) == resultant*monomials[j], 'homogeneous Bezout identity')
    K = max(sum(abs(adj[i,j]) for i in range(4)) for j in [0,3])
    U = max(sum(abs(ZZ(x)) for x in chart[k]) for k in ['s_numerator', 's_denominator'])

    # Exact deck involution from F(u)Y(v)-F(v)Y(u).
    V = PolynomialRing(QQ, names=('u', 'v'))
    u, v = V.gens()
    def affine(raw, x):
        return sum(QQ(c)*x**i for i, c in enumerate(raw))
    fu, yu = affine(chart['s_numerator'], u), affine(chart['s_denominator'], u)
    fv, yv = affine(chart['s_numerator'], v), affine(chart['s_denominator'], v)
    quotient, remainder = (fu*yv-fv*yu).quo_rem(u-v)
    require(remainder == 0, 'deck division')
    aa, bb, cc = [quotient.monomial_coefficient(m) for m in [u*v, u, V.one()]]
    require(quotient == aa*u*v + bb*(u+v) + cc, 'deck bilinear form')
    deck = matrix(QQ, [[-bb,-cc],[aa,bb]])
    deck = matrix(ZZ, deck / gcd(deck.list()))
    require(deck.det() != 0 and deck*deck == -deck.det()*matrix.identity(ZZ, 2), 'deck involution')
    da, db = deck*vector(P, [a,b])
    fd, yd, zd = [form(chart[key], da, db) for key in ['s_numerator','s_denominator','z_numerator']]
    require(fd*Y == F*yd and zd*Y == -Z*yd, 'deck fixes s and changes conic sign')

    R = PolynomialRing(QQ, 't')
    Kt = R.fraction_field()
    ainv = [Kt(R(x['numerator']))/R(x['denominator']) for x in parent['a_invariants']]
    j = EllipticCurve(Kt, ainv).j_invariant()
    j_degree = max(j.numerator().degree(), j.denominator().degree())
    require(j_degree > 0, 'nonisotriviality')
    return dict(
        schema='det1092-seed-geometry-v1',
        status='PASS_GENERIC_SEED_GEOMETRY_IDENTITIES',
        classification='verified application and elementary new deductions; analytic proof in canonical note',
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,COVER,REPLAY,CHART]},
        checker_sha256=sha(Path(__file__)),
        lattice=dict(generic_G=rows(G), trace_word=list(map(str,w)), cover_H=rows(H),
            generic_positive_leading_minors=list(map(str,g_minors)),
            complete_square_transform=rows(T), schur=str(schur),
            inherited_median=str(median), normalized_schur=str(schur/median),
            determinant=str(H.det()), parent_j_degree=int(j_degree), pullback_j_degree=2*int(j_degree)),
        parameter_map=dict(numerator=chart['s_numerator'], denominator=chart['s_denominator'],
            conic_ordinate=chart['z_numerator'], conic=chart['primitive_conic'],
            degree=2, sylvester_columns=rows(M), resultant=str(resultant), adjugate=rows(adj),
            height_lower_denominator=str(K), height_upper_multiplier=str(U),
            height_inequality='H(u)^2 / K <= H(S(u)) <= U H(u)^2 for every rational projective u',
            deck_matrix=rows(deck)),
        bounds=dict(gram_dimension=18, resultant_dimension=4, conic_identity_degree=4,
            specializations=0, cvps=0, point_searches=0, subgroup_states=0),
        boundary='No effective specialization threshold, rank upper bound, new curve rank, or finite V3 success guarantee. No exceptional point or V3 artifact is a computational input.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = verify()
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.write and not OUT.exists():
        with OUT.open('x') as output:
            output.write(text)
    else:
        require(OUT.read_text() == text, 'immutable arithmetic certificate mismatch')
    print(result['status'], flush=True)
    print('normalized Schur:', result['lattice']['normalized_schur'], flush=True)
    print('degree:', result['parameter_map']['degree'], 'CVPs: 0; point searches: 0', flush=True)
