#!/usr/bin/env sage-python
"""Complete the constructed curve-302 K3's arithmetic MW9 basis.

Build with --toric-output PATH; --check replays the retained exact data.
Add --recompute-frobenius to repeat the pinned external calculation (300s
limit). Two independent PARI fibre-count moments are always recomputed.
The complete Frobenius calculation is not independently implemented here.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from sage.all import (AA, EllipticCurve, GF, PolynomialRing, QQ, ZZ,
                      cyclotomic_polynomial, euler_phi, matrix, vector)
from sage.env import SAGE_VERSION
from parse_toric_controlled_reduction_output import parse_readfile_output

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json'
SAT = ROOT / 'artifacts/generated-results/elkies-k3-curve302-section-saturation-v1.json'
OUT = ROOT / 'artifacts/generated-results/elkies-k3-curve302-full-mw-basis-v1.json'
COMMIT = '74cda9e8148cd8e9a3928fc15a558c9a70b67cc1'
BACKEND = ROOT / ('artifacts/local/tools/ToricControlledReduction-' + COMMIT)
EXE = BACKEND / 'build/examples/readfile.exe'
PINS = {BASE: '51c9ef1f87dcb31effb83ab71f5b50a18ee56c351b8d4764f973143e43b64a77',
        SAT: 'a572da9b528454d51e77fcaf96e44ce95fc217c454a5bfa92c62518dfd0dd886'}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def pack(f):
    return {'numerator': list(map(str, f.numerator().list())),
            'denominator': list(map(str, f.denominator().list()))}


def certify(raw, recompute=False):
    for path, pinned in PINS.items():
        assert digest(path) == pinned, ('changed dependency', path)
    base, sat = [json.loads(path.read_text()) for path in (BASE, SAT)]
    R = PolynomialRing(QQ, 't'); t = R.gen(); K = R.fraction_field()
    model = base['weierstrass_family']
    A, B = [R(model[key + '_coefficients_low_to_high']) for key in ('A', 'B')]
    delta = -16*(4*A**3 + 27*B**2)
    assert A.degree() == 8 and B.degree() == 12
    assert delta.degree() == 24 and delta.is_squarefree()

    # Reflexive tetrahedron: its unique interior lattice point is (1,1,1).
    # The shifted primitive facet constants are all one. Nondegeneracy on
    # every face follows from squarefree B and delta, nonzero endpoints,
    # and nonsingular elliptic curves at t=0 and infinity.
    p = 47; F = GF(p); Rp = PolynomialRing(F, 't')
    a, b = [Rp([F(c) for c in f]) for f in (A, B)]
    d = -16*(4*a**3 + 27*b**2)
    assert (a.degree(), b.degree(), d.degree()) == (8, 12, 24)
    assert d.is_squarefree() and b.is_squarefree() and a.gcd(d).degree() == 0
    assert b[0] and b[12] and d[0] and d[24]
    facets = [(1,0,0), (0,1,0), (0,0,1), (-1,-4,-6)]
    offsets = [0,0,0,12]
    assert [sum(row) + c for row, c in zip(facets, offsets)] == [1]*4
    interior = [(i,j,k) for i in range(1,12) for j in range(1,3)
                for k in range(1,2) if i+4*j+6*k < 12]
    assert interior == [(1,1,1)]
    terms = {(0,0,2): 1, (0,3,0): p-1}
    terms.update({(i,1,0): int(-c) for i,c in enumerate(a) if c})
    terms.update({(i,0,0): int(-c) for i,c in enumerate(b) if c})
    monomials = sorted(terms); coefficients = [terms[m] for m in monomials]
    label = 'curve302-k3-p47'
    def ntlvec(v): return '[' + ' '.join(map(str,v)) + ']'
    def ntlmat(v): return '[' + ''.join(ntlvec(row) for row in v) + ']'
    raw_input = ':'.join([label, ntlmat(monomials), ntlvec(coefficients),
                          ntlmat(facets), ntlvec(offsets), str(p)]) + '\n'
    parsed = parse_readfile_output(raw)
    for key, expected in [('label', label), ('monomials', monomials),
                          ('coefficients', coefficients), ('halfspace_A', facets),
                          ('halfspace_b', offsets), ('prime', p),
                          ('hodge_numbers', [1,18,1])]:
        assert parsed[key] == expected, ('toric input binding', key)
    if recompute:
        assert subprocess.check_output(['git','rev-parse','HEAD'], cwd=BACKEND,
                                       text=True).strip() == COMMIT
        assert not subprocess.check_output(['git','diff','HEAD','--'], cwd=BACKEND)
        with tempfile.TemporaryDirectory(prefix='curve302-frobenius-') as tmp:
            inp, out = Path(tmp)/'input', Path(tmp)/'output'
            inp.write_text(raw_input)
            # Use the executing Sage runtime, without assuming a user prefix.
            env = os.environ.copy()
            lib = str(Path(sys.executable).resolve().parent.parent / 'lib')
            env['LD_LIBRARY_PATH'] = lib + ':' + env.get('LD_LIBRARY_PATH', '')
            run = subprocess.run([str(EXE),str(inp),str(out)], env=env,
                                 capture_output=True, text=True, timeout=300)
            assert run.returncode == 0, run.stderr
            assert out.exists() and out.stat().st_size, 'backend produced no result'
            assert parse_readfile_output(out.read_text()) == parsed

    ZT = PolynomialRing(ZZ, 'T'); T = ZT.gen()
    P = ZT(parsed['frobenius_coefficients'])
    assert P.degree() == 20 and P.is_monic()
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
    normalized = Rz(P(p*z)/p**20)
    assert normalized[0] in (-1,1)
    assert normalized == normalized[0]*z**20*normalized(1/z)
    hits = []; remainder = normalized
    # phi(m)>=sqrt(m/2): exhaustive for factors of degree at most 20.
    for order in range(1, 801):
        if euler_phi(order) > 20: continue
        cyclo = Rz(cyclotomic_polynomial(order)); exponent = 0
        while remainder % cyclo == 0:
            remainder //= cyclo; exponent += 1
        if exponent: hits.append({'order': order, 'multiplicity': exponent,
                                  'degree': int(cyclo.degree())})
    assert hits == [{'order':1,'multiplicity':9,'degree':1},
                    {'order':2,'multiplicity':1,'degree':1}]
    assert remainder.degree() == 10
    assert remainder == z**10*remainder(1/z)
    RW = PolynomialRing(QQ, 'W'); W = RW.gen()
    cheb = [RW(2), W]
    for i in range(2,6): cheb.append(W*cheb[-1]-cheb[-2])
    trace_poly = RW(remainder[5]) + sum(remainder[5+i]*cheb[i] for i in range(1,6))
    assert remainder == z**5*trace_poly(z+1/z)
    roots = trace_poly.roots(AA)
    assert sum(m for _,m in roots) == 5 and all(-2 < r < 2 for r,_ in roots)
    residual = P // ((T-p)**9*(T+p))
    assert P == (T-p)**9*(T+p)*residual

    # Independent point counting uses PARI cardinalities on each smooth
    # fibre, with the elementary nodal formula at any singular fibre.
    expected_moments = [-P[19], P[19]**2-2*P[18]]
    moments = []
    for n in (1,2):
        field = GF(p**n,'v'); q = field.order()
        RF = PolynomialRing(field,'t')
        af, bf = [RF([field(c) for c in f]) for f in (A,B)]
        trace = ZZ(0); singular = 0
        for aa, bb in [(af(v),bf(v)) for v in field] + [(af[8],bf[12])]:
            if 4*aa**3+27*bb**2:
                trace += q+1-EllipticCurve(field,[aa,bb]).cardinality()
            else:
                assert aa
                node = -3*bb/(2*aa)
                trace += 1 if (3*node).is_square() else -1
                singular += 1
        assert -trace == expected_moments[n-1]
        moments.append({'extension_degree': n, 'primitive_trace': int(-trace),
                        'surface_point_count': int(q*q+1+2*q-trace),
                        'singular_rational_fibres': singular,
                        'field_modulus': list(map(int,field.modulus().list())) if n>1 else None})

    # Expand all nine points in the pinned compact K3 gauge; no local probe
    # file supplies any coordinate. Reconstruct the ninth from the cubic map.
    c = base['pointed_cubic_family']
    N,D = [R(c[key+'_coefficients_low_to_high']) for key in ('N','D')]
    u = K(N/D); U = PolynomialRing(QQ,'u')
    def evaluate(record):
        return K(U(record['numerator'])(u)/U(record['denominator'])(u))
    points = [(D**2*evaluate(x), D**3*evaluate(y))
              for x,y in sat['source_basis']['Weierstrass_basis_xy']]
    mp = sat['cubic_to_Weierstrass']
    M = matrix(K, [[evaluate(v) for v in row] for row in mp['inverse_frame_matrix']])
    L = vector(K, [R(row) for row in c['moving_section_coordinates_low_to_high']])
    v,w,z0 = M*L; aa,bb = evaluate(mp['a']),evaluate(mp['b'])
    k,r,s,q0 = map(evaluate,mp['intermediate_to_compact_u_r_s_t'])
    xx,yy = -aa*bb*v/z0, aa*bb**2*v*w/z0**2
    points.append((K(D**2*(xx-r)/k**2), K(D**3*(yy-s*(xx-r)-q0)/k**3)))
    assert len(points) == 9
    assert all(y*y == x**3+A*x+B for x,y in points)
    assert N(1) == 0 and D(1) != 0 and delta(1) != 0
    # The compact source fibre is exactly the pinned short model, with its
    # sign fixed by eight specialized points in the saturation dependency.
    assert list(map(QQ,sat['source_basis']['fibre_to_short_u_r_s_t'])) == [1,0,0,0]
    sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
    from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS
    E302 = EllipticCurve(QQ, list(map(QQ, GENERAL_WEIERSTRASS_COEFFICIENTS)))
    public = [E302(list(map(QQ,pt))) for pt in POINTS]
    expected = [public[i]-public[0] for i in range(1,8)]
    expected += [-2*public[0]-public[1], public[8]-public[0]]
    image = []
    for x,y in points:
        X,Y = x(1)/D(1)**2, y(1)/D(1)**3
        xp = (X-15)/36; yp = (Y/108-xp-1)/2
        image.append(E302([xp,yp]))
    assert image == expected
    H = matrix(ZZ,sat['K3']['height_gram'])
    assert H.det() == 512 and H.rank() == 9
    assert sat['K3']['saturated_in_full_geometric_MW'] is True
    return {
        'schema': 'elkies-k3.curve302-full-mw-basis.v1',
        'status': 'PASS_EXACT_CONSTRUCTED_K3_FULL_ARITHMETIC_MW9_BASIS',
        'source_sha256': digest(Path(__file__)),
        'input_sha256': {str(path.relative_to(ROOT)): pin for path,pin in PINS.items()},
        'software': {'sage': SAGE_VERSION, 'toric_commit': COMMIT,
                     'toric_repository': 'https://github.com/edgarcosta/ToricControlledReduction',
                     'original_backend_binary_sha256': '0207af5145b4e9e85e5790ae0f91d91787b8d170d6d88e5d46a55c46c6627289',
                     'strict_parser_sha256': digest(Path(__file__).with_name('parse_toric_controlled_reduction_output.py'))},
        'protocol': {'fixed_object': 'pinned curve302 nine-direction K3',
                     'prime_budget': [47,53,127], 'used_primes': [47],
                     'maximum_seconds_per_prime': 300, 'workers': 1,
                     'stop': 'arithmetic upper bound matches saturated rank9',
                     'failure': 'missing, invalid or timed-out output leaves upper bound UNKNOWN'},
        'frobenius': {'prime': p, 'raw_input': raw_input, 'raw_output': raw.strip()+'\n',
                      'newton_vertices': [[0,0,0],[12,0,0],[0,3,0],[0,0,2]],
                      'interior_point': [1,1,1], 'nondegenerate_face_gates': True,
                      'hodge_vector': [1,18,1], 'primitive_degree': 20,
                      'full_H2_degree': 22, 'complement_dimension': 2,
                      'characteristic_coefficients': list(map(str,P.list())),
                      'residual_degree10_coefficients': list(map(str,residual.list())),
                      'cyclotomic_factors_normalized': hits,
                      'weil_circle_certified': True, 'independent_point_counts': moments,
                      'complete_polynomial_independently_implemented': False},
        'surface': {'equation': 'y^2=x^3+A(t)*x+B(t)',
                    'A_coefficients_low_to_high': list(map(str,A.list())),
                    'B_coefficients_low_to_high': list(map(str,B.list())),
                    'singular_fibres_geometric': '24 I1; smooth at infinity'},
        'basis': {'labels': ['P2-P1','P3-P1','P4-P1','P5-P1','P6-P1','P7-P1','P8-P1','Q12-P1','S'],
                  'Weierstrass_xy': [[pack(x),pack(y)] for x,y in points],
                  'height_gram': [list(map(int,row)) for row in H.rows()],
                  'determinant': 512, 'torsion_order': 1,
                  'full_arithmetic_MW_basis': True, 'generic_rank_Qt': 9,
                  'generic_geometric_rank_bounds': [9,10]},
        'specialization': {'t0': '1', 'source_parameter': '0',
                           'D_at_1': str(D(1)),
                           'to_public': 'X=x/D(1)^2; Y=y/D(1)^3; xp=(X-15)/36; yp=(Y/108-xp-1)/2',
                           'public_a_invariants': list(map(str,E302.a_invariants())),
                           'all_nine_transports_checked': True,
                           'public_coordinate_matrix': sat['K3']['specialization_matrix'],
                           'displayed_public_quotient': 'Z^22 + Z/3'},
        'conclusions': {'rational_NS_rank': 11, 'geometric_NS_rank_bounds': [11,12],
                        'quadratic_source_twist_rank_Qu': 1,
                        'any_Q_elliptic_fibration_on_this_K3_MW_rank_upper_bound': 9,
                        'original_parent_provenance': 'UNKNOWN',
                        'rank_E302_unconditional_upper_bound': 'UNKNOWN'},
        'proof': 'Good-reduction specialization injects NS into H2. The 20-dimensional toric factor has nine p-eigenvalues and ten p-times-root-of-unity eigenvalues. The complementary space has dimension two. Thus rational NS rank is at most11 and geometric NS rank at most12. The fibre, zero section and nine rational independent sections give rational NS rank at least11. Shioda-Tate gives arithmetic MW rank9; the pinned saturation and trivial torsion make the nine sections a full basis. No Tate converse, BSD or GRH is used.',
        'boundary': 'Explicit constructed parent with full arithmetic MW9 basis and t=1 fibre E302; this does not identify the discoverers original construction or certify the full geometric MW basis.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--toric-output', type=Path)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--recompute-frobenius', action='store_true')
    args = parser.parse_args()
    previous = json.loads(args.output.read_text()) if args.check else None
    if args.toric_output:
        raw = args.toric_output.read_text()
    elif previous:
        raw = previous['frobenius']['raw_output']
    else:
        parser.error('build requires --toric-output')
    result = certify(raw, args.recompute_frobenius)
    if args.check:
        assert result == previous, 'certificate mismatch'
    else:
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'])
    print('rank Q(t)=9; rational NS rank=11; geometric MW rank in [9,10]; t0=1')


if __name__ == '__main__':
    main()
