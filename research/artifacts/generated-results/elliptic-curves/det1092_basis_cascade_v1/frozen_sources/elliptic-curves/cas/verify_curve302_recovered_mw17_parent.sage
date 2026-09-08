#!/usr/bin/env sage-python
"""Replay the explicit MW17 parent, saturation and exact302 specialization.

One worker, 120 seconds. No search and no backend rerun. The retained full
Frobenius polynomial is external ToricControlledReduction evidence, bound to
this exact model; two moments are independently recomputed from fibre counts.
Use --build once to freeze the compact proof output, then run without flags.
"""
import argparse
import json
import runpy
import signal
import sys
from hashlib import sha256
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/generated-results/elliptic-curves'
INPUT = BASE / 'curve302_recovered_mw17_parent_v1.json'
OUTPUT = BASE / 'curve302_recovered_mw17_parent_proof_v1.json'
LOADER = ROOT / 'elliptic-curves/cas/load_curve302_recovered_parent.sage'
PUBLIC = ROOT / 'elliptic-curves/cas/icarm_curve302.py'
CORE = BASE / 'curve302_integral_shell_inputs_v1.json'
RECOVERY = BASE / 'curve302_det1092_recovery_v1.json'
PARSER = ROOT / 'elkies-k3/scripts/parse_toric_controlled_reduction_output.py'


def build():
    data = json.loads(INPUT.read_text())
    E, points, t0 = runpy.run_path(str(LOADER))['load_curve302_recovered_parent'](INPUT)
    K = E.base_ring()
    R = K.ring()
    assert t0 == 0 and len(points) == 17
    assert E.a_invariants()[:3] == (K(1), K(1), K(1))
    u, r, s, tau = map(K, data['raw_short_to_normalized_u_r_s_t'])
    assert [u, r, s, tau] == list(map(K, [-6, 15, -3, -108]))
    # Inverse of raw_short.change_weierstrass_model([u,r,s,tau]).
    J = E.change_weierstrass_model([1/u, -r/u**2, -s/u, (s*r-tau)/u**3])
    assert J.a_invariants()[:3] == (K(0), K(0), K(0))
    A, B = map(R, [J.a4(), J.a6()])
    delta = R(J.discriminant())
    assert (A.degree(), B.degree(), delta.degree()) == (8, 12, 24)
    assert delta.is_squarefree() and A.gcd(delta).degree() == 0
    # Thus the minimal elliptic surface is K3, with 24 I1 fibres and smooth
    # infinity. There are no geometric correction terms or nonzero torsion.
    raw_points = [J(u**2*P[0]+r, u**3*P[1]+s*u**2*P[0]+tau) for P in points]

    def height(x):
        return max(x.numerator().degree(), x.denominator().degree()+4)

    H = matrix(QQ, 17)
    for i, P in enumerate(raw_points):
        H[i, i] = height(P[0])
        for j in range(i):
            Q = raw_points[j]
            assert P[0] != Q[0]
            x_difference = ((P[1]+Q[1])/(P[0]-Q[0]))**2-P[0]-Q[0]
            H[i, j] = H[j, i] = (H[i, i]+H[j, j]-height(x_difference))/2
    assert H == matrix(QQ, data['generic_height_gram'])
    assert H.is_positive_definite() and H.det() == 1092
    assert all(a in ZZ for a in H.list()) and all(a in 2*ZZ for a in H.diagonal())
    H = matrix(ZZ, H)
    assert list(H.smith_form()[0].diagonal()) == [1]*16+[1092]
    # det=4*3*7*13: the only possible first enlargement has index two.
    parity = H.change_ring(GF(2)).right_kernel_matrix()
    assert parity.nrows() == 1
    half = vector(QQ, list(map(ZZ, parity.row(0))))/2
    obstruction_norm = half*H*half
    assert obstruction_norm in ZZ and obstruction_norm % 2 == 1

    pub = runpy.run_path(str(PUBLIC))
    E0 = EllipticCurve(QQ, list(map(QQ, pub['GENERAL_WEIERSTRASS_COEFFICIENTS'])))
    assert [a(t0) for a in E.a_invariants()] == list(E0.a_invariants())
    public = [E0(QQ(x), QQ(y)) for x, y in pub['POINTS']]
    embedding = matrix(ZZ, data['basis_embedding_in_public_D'])
    core = matrix(ZZ, json.loads(CORE.read_text())['target_core']['public_embedding'])
    words = matrix(ZZ, data['basis_words_in_recovered_core'])
    assert words.nrows() == words.ncols() == 17 and abs(words.det()) == 1
    assert embedding == core*words.transpose()
    assert embedding.nrows() == 31 and embedding.ncols() == 17
    assert list(embedding.smith_form()[0].diagonal()) == [1]*17
    recovered = matrix(ZZ, json.loads(RECOVERY.read_text())['recovery']['gram'])
    assert H == words*recovered*words.transpose()
    for P, col in zip(points, embedding.columns()):
        image = E0([c(t0) for c in P.xy()])
        expected = sum((n*Q for n, Q in zip(col, public) if n), E0(0))
        assert image == expected
    print('PASS exact equation, all 17 heights, saturation obstruction and 302 images', flush=True)

    f = data['frobenius']
    external_paths = []
    for name in ['input.txt', 'output.txt', 'independent-moments.json', 'backend.log']:
        path = ROOT / f[name+'_path']
        assert sha256(path.read_bytes()).hexdigest() == f[name+'_sha256']
        external_paths.append(path)
    p = ZZ(f['prime'])
    assert p == 149 and f['local_scaling_exponent'] == 0
    assert f['local_base_map'] == {'r0': 0, 'rinf': 1}
    F = GF(p)
    Rp = PolynomialRing(F, 'v')
    v = Rp.gen()
    a0, b0 = Rp(A), Rp(B)
    a = sum(c*v**i*(v+1)**(8-i) for i, c in enumerate(a0))
    b = sum(c*v**i*(v+1)**(12-i) for i, c in enumerate(b0))
    d = -16*(4*a**3+27*b**2)
    assert (a.degree(), b.degree(), d.degree()) == (8, 12, 24)
    assert d.is_squarefree() and b.is_squarefree() and a.gcd(d).degree() == 0
    assert b[0] and b[12] and d[0] and d[24]
    assert list(map(int, a.list())) == f['reduced_A']
    assert list(map(int, b.list())) == f['reduced_B']
    # The four-facet reflexive simplex has unique interior point (1,1,1).
    # Face nondegeneracy: y=0 uses squarefree delta, x=0 squarefree B;
    # endpoint faces use nonzero delta; remaining faces are binomial with
    # exponents invertible at149. These checks also give good K3 reduction.
    facets = [(1,0,0), (0,1,0), (0,0,1), (-1,-4,-6)]
    offsets = [0,0,0,12]
    assert [sum(row)+c for row, c in zip(facets, offsets)] == [1]*4
    terms = {(0,0,2): 1, (0,3,0): int(p-1)}
    terms.update({(i,1,0): int(-c) for i,c in enumerate(a) if c})
    terms.update({(i,0,0): int(-c) for i,c in enumerate(b) if c})
    mons = sorted(terms)
    coefs = [terms[m] for m in mons]
    def nv(xs): return '['+' '.join(map(str, xs))+']'
    def nm(xs): return '['+''.join(nv(row) for row in xs)+']'
    label = 'curve302-det1092-recovered-p149'
    expected_input = ':'.join([label,nm(mons),nv(coefs),nm(facets),nv(offsets),str(p)])+'\n'
    assert (ROOT/f['input.txt_path']).read_text() == expected_input
    parsed = runpy.run_path(str(PARSER))['parse_readfile_output']((ROOT/f['output.txt_path']).read_text())
    for key, value in [('label',label),('monomials',mons),('coefficients',coefs),
                       ('halfspace_A',facets),('halfspace_b',offsets),('prime',p),('hodge_numbers',[1,18,1])]:
        assert parsed[key] == value
    S = PolynomialRing(QQ, 'T')
    T = S.gen()
    P = S(parsed['frobenius_coefficients'])
    assert P == (T-p)**17*(T+p)*(T**2+248*T+p**2)
    assert list(map(str, P.list())) == f['primitive_polynomial']
    assert P.degree() == 20 and 248**2 < 4*p**2
    # Remaining normalized quadratic is z^2+(248/149)z+1. Its nonintegral
    # trace excludes roots of unity. The -p factor adds one geometric class
    # but no F149-rational class; two ambient classes are rational.
    assert QQ(248)/p not in ZZ
    assert f['arithmetic_NS_upper_bound'] == 19 and f['geometric_NS_upper_bound'] == 20
    assert f['arithmetic_MW_upper_bound'] == 17
    assert f['backend_commit'] == '74cda9e8148cd8e9a3928fc15a558c9a70b67cc1'

    moments = []
    expected_traces = [-P[19], P[19]**2-2*P[18]]
    for n in [1,2]:
        Fn = GF(p**n, 'v')
        q = Fn.order()
        Rn = PolynomialRing(Fn, 't')
        an, bn = Rn(a), Rn(b)
        total = ZZ(0)
        singular = 0
        for parameter in list(Fn)+[None]:
            aa, bb = (an(parameter),bn(parameter)) if parameter is not None else (an[8],bn[12])
            if 4*aa**3+27*bb**2:
                total += q+1-EllipticCurve(Fn,[aa,bb]).cardinality()
            else:
                assert aa
                node = -3*bb/(2*aa)
                total += 1 if (3*node).is_square() else -1
                singular += 1
        assert -total == expected_traces[n-1]
        moments.append(dict(extension_degree=n, primitive_trace=int(-total),
                            surface_point_count=int(q*q+1+2*q-total), singular_fibres=singular,
                            field_modulus=list(map(int,Fn.modulus().list()))))
    assert moments == json.loads((ROOT/f['independent-moments.json_path']).read_text())['records']
    print('PASS Frobenius model binding, rank17 capacity gate and two independent moments', flush=True)
    paths = [Path(__file__), INPUT, LOADER, PUBLIC, CORE, RECOVERY, PARSER,
             ROOT/'elliptic-curves/cas/reconstruct_curve302_recovered_parent.sage',
             BASE/'curve302_recovered_mw17_construction_inputs_v1.json']+external_paths
    return dict(schema='curve302.recovered-mw17-parent-proof.v1', status='PASS_FULL_ARITHMETIC_MW17_PARENT',
                input_sha256={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in paths},
                generic_arithmetic_MW_rank=17, generic_torsion_order=1, full_basis_size=17,
                generic_geometric_MW_rank_interval=[17,18], rational_NS_rank=19,
                geometric_NS_rank_interval=[19,20], absolute_rational_NS_determinant=1092,
                height_determinant=1092, discriminant_group_invariants=[1092],
                only_possible_enlargement_prime=2, unique_order_two_coset_norm=str(obstruction_norm),
                full_generic_saturation=True, geometric_fibres='24 I1',
                specialization_parameter='0', specialization_model='literal public302',
                specialization_Smith_factors=[1]*17, displayed_D_quotient='Z^14',
                primitive_rank17_core_recovered=True, frobenius_prime=149,
                primitive_frobenius_factorization='(T-149)^17*(T+149)*(T^2+248*T+22201)',
                independent_fibre_count_moments=moments,
                independent_full_frobenius_implementation=False,
                original_discoverer_provenance='UNKNOWN', E302_exact_rank='UNKNOWN')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', action='store_true')
    args = parser.parse_args()
    signal.alarm(120)
    result = build()
    if args.build:
        assert not OUTPUT.exists()
        OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    assert result == json.loads(OUTPUT.read_text())
    print(result['status'])
