#!/usr/bin/env sage-python
"""Two retrospective cubic squareclasses; no class-group or norm factoring.

Run with timeout 25s. Protocol is checkpointed before arithmetic. A three-
dimensional exact trace Gram is LLL-reduced, then only its three basis vectors
are examined. No shortest-vector or principal-ideal assertion is made.
"""
import hashlib, json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, NumberField, matrix, pari, lcm

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_seed_half_ideal_v1'
ARITH = ART / 'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
SEED = ART / 'det1092_marked_kummer_transport_v1.json'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name, data):
    p = OUT / name
    if p.exists(): assert json.loads(p.read_text()) == data
    else:
        with p.open('x') as stream:
            json.dump(data, stream, indent=2, sort_keys=True); stream.write('\n')
def entries(A): return [[str(x) for x in row] for row in A.rows()]
def coefficients(x): return [str(x[i]) for i in range(3)]
def stripped(q, S):
    q = abs(QQ(q)); powers = []
    for p in S:
        e = q.valuation(p)
        if e: powers.append([str(p), int(e)]); q /= QQ(p)**e
    return {'S_valuations': powers, 'remaining': str(q)}

OUT.mkdir(exist_ok=True)
save('protocol.json', {
    'classification': 'retrospective diagnostic, not prospective selection',
    'inputs': {str(p.relative_to(ROOT)): sha(p) for p in [ARITH, SEED, Path(__file__)]},
    'limits': {'seconds': 25, 'classes': 2, 'LLL_vectors_per_class': 3,
               'point_searches': 0, 'class_group_computations': 0,
               'unrestricted_factorizations': 0},
    'method': 'I=(alpha,sqrt(N(alpha))); reduce I inverse using Tr(alpha*u*v) if positive definite, otherwise Tr(alpha^2*u*v); test only three reduced basis vectors.',
    'boundary': 'No global minimality, nonprincipality, or prospective source follows from reduction.'})
d = json.loads(ARITH.read_text()); seed = json.loads(SEED.read_text())['cases'][0]
assert seed['label'] == '302-first-unlock'
R = PolynomialRing(QQ, 'z'); z = R.gen(); f = R(d['cubic_ascending'])
K = NumberField(f, 'theta'); theta = K.gen()
S = list(map(ZZ, d['S_finite']))
pari.allocatemem(64000000, 268435456, silent=True)
nf = pari.nfinit([pari(f), list(map(int, S))])
assert ZZ(nf[3]) == ZZ(d['defining_order_index'])
assert ZZ(nf.disc()) == ZZ(d['field_discriminant'])
assert f.discriminant() == ZZ(nf[3])**2 * ZZ(nf.disc())
zk = [K(R(str(v))) for v in nf.nf_get_zk()]
assert zk == [K(R(s)) for s in d['maximal_order_basis']]
seed_poly = R(seed['beta_branch']); den = lcm([q.denominator() for q in seed_poly.list()])
assert den.is_square()
cases = [('302-first-unlock', K(seed_poly * den)),
         ('generic-section-0', K(R(d['generic_classes'][0]['beta_ascending'])))]
def to_pari(a): return pari.Mod(pari(R(a.list())), pari(f))
def ideal_matrix(a): return matrix(QQ, a.sage())
results = []
for label, alpha in cases:
    norm = ZZ(alpha.norm()); assert norm.is_square() and norm > 0
    y = norm.sqrt(); pa = to_pari(alpha)
    I = pari.idealhnf(nf, pa, pari(y))
    Iinv = pari.idealinv(nf, I)
    defect = pari.idealdiv(nf, pari.idealpow(nf, I, 2), pa)
    defect_norm = QQ(pari.idealnorm(nf, defect))
    support = stripped(defect_norm, S)
    assert support['remaining'] == '1'
    B = ideal_matrix(Iinv)
    basis = [sum((B[i,j] * zk[i] for i in range(3)), K.zero()) for j in range(3)]
    weight = alpha
    G = matrix(QQ, 3, 3, lambda i,j: (weight*basis[i]*basis[j]).trace())
    positive = all(G[:i,:i].det() > 0 for i in [1,2,3])
    if not positive:
        weight = alpha**2
        G = matrix(QQ, 3, 3, lambda i,j: (weight*basis[i]*basis[j]).trace())
        assert all(G[:i,:i].det() > 0 for i in [1,2,3])
    scale = lcm([v.denominator() for v in G.list()])
    U = matrix(ZZ, pari.qflllgram(pari(matrix(ZZ, scale*G))).sage())
    assert abs(U.det()) == 1
    reduced = U.transpose()*G*U
    candidates = []
    for j in range(3):
        gamma = sum((U[i,j]*basis[i] for i in range(3)), K.zero())
        beta = alpha*gamma**2
        assert gamma and beta.norm() == (QQ(y)*gamma.norm())**2
        coords = matrix(QQ, 3, 3, lambda i,j: zk[j][i]).solve_right(beta.vector())
        candidates.append({'gamma': coefficients(gamma), 'beta': coefficients(beta),
                           'norm_square_root': str(QQ(y)*gamma.norm()),
                           'trace': str(beta.trace()),
                           'integral_basis_coordinates': list(map(str,coords)),
                           'is_integral': all(c.denominator() == 1 for c in coords),
                           'is_original_up_to_rational_square': gamma in QQ,
                           'weighted_norm': str(reduced[j,j])})
    item = {'label': label, 'alpha': coefficients(alpha), 'norm_square_root': str(y),
            'half_ideal_basis': entries(ideal_matrix(I)),
            'inverse_half_ideal_basis': entries(B),
            'half_ideal_norm': str(pari.idealnorm(nf,I)),
            'square_defect_basis': entries(ideal_matrix(defect)),
            'square_defect_norm': str(defect_norm), 'square_defect_support': support,
            'exact_square_root_ideal': ideal_matrix(defect) == matrix.identity(QQ,3),
            'weight': 'alpha' if positive else 'alpha^2',
            'trace_Gram': entries(G), 'LLL_transform': entries(U),
            'reduced_Gram': entries(reduced), 'candidates': candidates}
    save(label+'.json', item); results.append(item)
    print(label, 'square ideal', item['exact_square_root_ideal'],
          'defect norm', defect_norm, 'original representatives',
          [v['is_original_up_to_rational_square'] for v in candidates], flush=True)
save('summary.json', {'status':'PASS_TWO_RETROSPECTIVE_HALF_IDEAL_REDUCTIONS',
    'classification': 'verified application', 'checker_sha256': sha(Path(__file__)),
    'cases': [x['label'] for x in results],
    'boundary': 'Squareclass equivalence only. No class-group computation, no proof of shortest representative, and no prospective nongeneric class source.'})
