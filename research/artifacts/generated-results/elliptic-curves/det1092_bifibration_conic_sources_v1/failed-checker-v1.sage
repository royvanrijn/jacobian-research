#!/usr/bin/env sage-python
"""Independent three-source map, no-incidence and one-address seed replay.

No producer imports, polynomial factor/root discovery, point search, or
later exceptional input. Sage finite groups replace the producer's manual
group engine. All-prime statements are proved in the companion note;
the polynomial identities and finite-cycle tests here are not that proof.
One action, hard 25-second limit. Immutable input and output checkpoints.
"""
import argparse
import hashlib
import json
import signal
from pathlib import Path
from sage.all import (QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix,
                      vector, gcd, lcm, prime_range)

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_bifibration_conic_sources_v1'
CONICS = ART/'det1092_rational_bisection_index_v1'
MASKS = [8044, 47755, 103186]
R = PolynomialRing(QQ, 'u')
K = R.fraction_field()
u = R.gen()


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(v): return K(R(v['numerator']))/R(v['denominator'])
def degree(v): return max(v.numerator().degree(), v.denominator().degree())
def key(P): return tuple(int(v) for v in P)
def point(E, v): return E(0) if v is None else E(list(map(QQ, v)))
def hashes(row):
    for name, digest in row['inputs'].items(): assert sha(ROOT/name) == digest


def save(name, value):
    path = OUT/name
    if path.exists(): assert read(path) == value, 'immutable checkpoint changed'
    else:
        with path.open('x') as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write('\n')


def check_root_certificate(f, proof, pool):
    """Check supplied lifts/lattices, never run the producer's Gauss solver."""
    assert f.degree() > 0 and f[0] and f.leading_coefficient() > 0
    assert all(c in ZZ for c in f.list()) and gcd(list(f)) == 1
    H = max(abs(ZZ(f[0])), abs(ZZ(f.leading_coefficient())))
    assert ZZ(proof['numerator_bound']) == abs(f[0])
    assert ZZ(proof['denominator_bound']) == f.leading_coefficient()
    assert ZZ(proof['common_bound']) == H
    attempts = []
    for p in pool:
        assert ZZ(p).is_prime(proof=True)
        if f.leading_coefficient() % p == 0:
            attempts.append(dict(prime=p, gate='leading_coefficient_zero'))
            continue
        roots = [i for i in range(p) if f(i) % p == 0]
        simple = all(f.derivative()(i) % p for i in roots)
        attempts.append(dict(prime=p, gate='PASS' if simple else 'multiple_root', roots=roots))
        if simple: break
    else: raise AssertionError('no eligible prime in frozen pool')
    assert attempts == proof['prime_attempts'] and p == proof['prime']
    assert roots == proof['all_roots_mod_p']
    assert [row['initial_root'] for row in proof['lifts']] == roots
    for row in proof['lifts']:
        M, r = ZZ(row['modulus']), ZZ(row['lift'])
        exponent = M.valuation(p)
        assert exponent > 0 and ZZ(p)**exponent == M and M > 8*H*H
        assert 0 <= r < M and r % p == row['initial_root'] and f(r) % M == 0
        a, b = [vector(ZZ, v) for v in row['reduced_basis']]
        assert abs(matrix(ZZ, [a, b]).det()) == M
        assert all((v[0]-r*v[1]) % M == 0 for v in [a, b])
        assert 0 < a*a <= b*b and 2*abs(a*b) <= a*a
        if a*a > 2*H*H:
            assert row['decision'] == 'NO_BOUNDED_NUMERATOR_DENOMINATOR_VECTOR'
        else:
            assert a[1] and row['decision'] == 'ONLY_POSSIBLE_RATIO_IS_NOT_A_ROOT'
            ratio = QQ(a[0])/a[1]
            assert ratio == QQ(row['candidate']) and f(ratio) != 0
    assert proof['rational_roots'] == []
    return dict(prime=int(p), modular_roots=len(roots), rational_roots=0)


def geometry(generic):
    old = read(ART/'det1092_norm8_singular_members_v1/independent-replay.json')
    hashes(old)
    assert old['status'] == 'PASS_INDEPENDENT_NORM8_RATIONAL_SINGULAR_MEMBER_OBSTRUCTION'
    assert old['Kodaira_configuration'] == dict(I1=14, I2=5)
    data = read(ART/'det1092_norm8_singular_members_v1/geometry.json')
    e, d, c, b, a = map(R, generic['quartic_t_coefficients_in_z'])
    I = 12*a*e-3*b*d+c*c
    J = 72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    A, B = -27*I, -27*J
    delta = -16*(4*A**3+27*B*B)
    repeated, simple = R(data['repeated_support']), R(data['simple_support'])
    assert A == R(generic['Jacobian_A']) == R(data['A']) and A.degree() == 8
    assert B == R(generic['Jacobian_B']) == R(data['B']) and B.degree() == 12
    assert delta == QQ(data['delta_scalar'])*repeated**2*simple and delta.degree() == 22
    assert repeated.degree() == 4 and simple.degree() == 14
    assert gcd(repeated, repeated.derivative()) == gcd(simple, simple.derivative()) == 1
    assert gcd(repeated, simple) == gcd(A, delta) == 1
    assert A[8] != 0  # scaled c4 is a unit at infinity; discriminant order2.
    assert 4-QQ(5)/2 == QQ(3)/2 and 19-2-5 == old['geometric_generic_MW_rank'] == 12
    ell = R.gen()
    bound = 1-ell**2+QQ(19)/2*(ell-1)**2
    assert bound == (ell-1)*(17*ell-21)/2
    assert bound-30 == (ell-3)*(17*ell+13)/2
    assert 1-4+QQ(14)/2 == 4 and 1-3+QQ(14)/2 == 5
    # Regression only: every affine lift of a transvection has at least
    # (ell-1)^2 permutation defect for odd primes (>=1 at ell2).
    tests = []
    for prime in [2, 3, 5, 7]:
        defects = []
        for aa in range(prime):
            for bb in range(prime):
                unseen = {(x, y) for x in range(prime) for y in range(prime)}
                cycles = 0
                while unseen:
                    v = next(iter(unseen)); cycles += 1
                    while v in unseen:
                        unseen.remove(v)
                        v = ((v[0]+v[1]+aa) % prime, (v[1]+bb) % prime)
                defects.append(prime*prime-cycles)
        assert min(defects) >= (1 if prime == 2 else (prime-1)**2)
        tests.append(dict(prime=prime, minimum_defect=min(defects), maximum_defect=max(defects)))
    return dict(Kodaira_configuration=dict(I1=14, I2=5), generic_rank=12,
                nonzero_section_height_lower_bound='3/2',
                primitive_division_genus_lower_bound_ell2=4,
                nonzero_two_torsion_genus=5,
                odd_prime_division_genus_lower_bound='(ell-1)*(17*ell-21)/2',
                all_prime_claim_source='written argument, not finite regression', cycle_regressions=tests)


def maps_and_incidence(parent, pencil, protocol):
    old_replay = read(CONICS/'independent-replay.json'); hashes(old_replay)
    assert old_replay['status'] == 'PASS_INDEPENDENT_RATIONAL_BISECTION_INDEX_AND_THREE_MAPS'
    ip = read(OUT/'incidence-protocol.json'); hashes(ip)
    old_ip = read(ART/'det1092_signed_source_orbits_v1/protocol.json')
    assert ip['carrier_label'] == old_ip['carrier_label']
    assert ip['primes'] == old_ip['limits']['prime_pool']
    G = matrix(ZZ, parent['generic_height_gram']); w = vector(ZZ, [0]*14+[1, -1, 0])
    assert G.det() == 1092 and w*G*w == 8
    old_curve = EllipticCurve(K, [dec(v) for v in parent['a_invariants']])
    # A displayed original section is a degree1 source/zero for the alternate pencil.
    origin = vector(ZZ, [0]*14+[1, 0, 0]); assert origin*G*origin-w*G*origin == 1
    rows = []
    for mask in MASKS:
        source = read(CONICS/('orbit-%d.json' % mask)); hashes(source)
        mp = read(OUT/('map-%d.json' % mask)); inc = read(OUT/('incidence-%d.json' % mask))
        assert mp['source_sha256'] == sha(CONICS/('orbit-%d.json' % mask))
        assert mp['generic_protocol_sha256'] == sha(OUT/'generic-protocol.json')
        assert inc['map_sha256'] == sha(OUT/('map-%d.json' % mask))
        assert inc['incidence_protocol_sha256'] == sha(OUT/'incidence-protocol.json')
        a = vector(ZZ, source['word']); assert a*G*a == 10
        expected = int(8-w*G*a)
        assert expected == next(v['alternate_degree'] for v in protocol['sources'] if v['mask'] == mask)
        T, W, Z = [dec(mp[v]) for v in ['T', 'conic_ordinate', 'Z']]
        assert T == dec(source['base_map']) and W == dec(source['conic_ordinate'])
        q = R(source['q']); assert q.degree() == 2 and q.discriminant() != 0
        assert W*W == q(T) and degree(T) == mp['original_degree'] == 2
        # Degree1 common fibre over (T,W) verifies birationality to the conic.
        V = PolynomialRing(K, 'v'); v = V.gen()
        fT = V(T.numerator())-T*V(T.denominator())
        fW = V(W.numerator())-W*V(W.denominator())
        assert gcd(fT, fW).degree() == 1
        x0, x1, y0, y1 = [R(coeff)(T) for coeff in source['elliptic_quadratic_maps']]
        x, y = x0+x1*W, y0+y1*W
        assert x1 and [dec(v) for v in mp['original_point']] == [x, y]
        values = [v(T) for v in old_curve.a_invariants()]
        a1, a2, a3, a4, a6 = values
        assert y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6
        X, Y = x+old_curve.b2()(T)/12, y+(a1*x+a3)/2
        assert [dec(v) for v in mp['short_point']] == [X, Y]
        h, shift = R(pencil['pole_h'])(T), R(pencil['shift'])(T)
        cx, cy = R(pencil['nx'])(T)/h**2, R(pencil['ny'])(T)/h**3
        slope = (Y+cy)/(X-cx)
        assert Z == (slope+shift/h)/h and degree(Z) == mp['alternate_degree'] == expected
        # Reconstruct the alternate quartic coordinate and check its equation.
        alt_W = (2*X+cx-slope**2)/h
        alt_F = sum(R(coeff)(Z)*T**i for i, coeff in enumerate(pencil['quartic_t_coefficients_in_z']))
        assert alt_W**2 == alt_F
        f = Z.numerator()-QQ(ip['carrier_label'])*Z.denominator()
        f *= lcm([v.denominator() for v in f]); f /= gcd(list(f))
        f = R(f)
        if f.leading_coefficient() < 0: f = -f
        assert f == R(inc['polynomial']) and f.degree() == expected
        assert inc['rational_preimages'] == [] and not inc['infinity_preimage']
        assert inc['status'] == 'EXACT_NO_RATIONAL_SOURCE_POINT'
        proof = check_root_certificate(f, inc['root_proof'], ip['primes'])
        rows.append(dict(mask=mask, bidegree=[2, expected], conic_birational=True,
                         first_carrier_rational_incidence=False, no_root_certificate=proof))
    selected = min(rows, key=lambda row: row['bidegree'][1])
    assert selected['mask'] == protocol['selected']['mask'] == 47755
    assert selected['bidegree'] == [2, 2] and protocol['parameter'] == '0'
    return rows


def specialized_seed(parent, protocol):
    result = read(OUT/'seed-result.json'); frame = read(OUT/'seed-generic-frame.json')
    mp = read(OUT/'map-47755.json')
    assert result['frame_sha256'] == sha(OUT/'seed-generic-frame.json')
    assert result['map_sha256'] == sha(OUT/'map-47755.json')
    assert result['generic_protocol_sha256'] == sha(OUT/'generic-protocol.json')
    parameter = QQ(protocol['parameter']); tau = dec(mp['T'])(parameter)
    assert result['mask'] == 47755 and QQ(result['parameter']) == parameter == 0
    assert tau == QQ(result['original_parameter'])
    assert dec(mp['Z'])(parameter) == QQ(result['alternate_parameter'])
    old = EllipticCurve(QQ, [dec(v)(tau) for v in parent['a_invariants']])
    E = EllipticCurve(QQ, [-old.c4()/48, -old.c6()/864]); assert E.discriminant()
    basis = []
    for row in parent['basis_weierstrass_coordinates']:
        x, y = [dec(v)(tau) for v in row]
        basis.append(E([x+old.b2()/12, y+(old.a1()*x+old.a3())/2]))
    P = E([dec(v)(parameter) for v in mp['short_point']])
    assert point(E, result['original']) == point(E, result['terminal']) == P
    assert result['status'] == 'NEW_INDEPENDENT_DIRECTION' and result['reason'] == 'FINITE_FOOTPRINT_ESCAPE'
    assert result['steps'] == 0 and result['history'] == [] and not any(map(ZZ, result['chain_word']))
    assert list(map(QQ, frame['curve'])) == list(E.a_invariants())
    assert [point(E, v) for v in frame['basis']] == basis and frame['candidate_inputs'] == 0
    assert frame['prime_cap'] == protocol['prime_cap'] == 1009
    groups = []; seen_odd = False; boundaries = []
    for row in frame['records']:
        p = row['prime']; assert ZZ(p).is_prime(proof=True) and 2 < p <= 1009
        assert all(a.denominator() % p for a in E.a_invariants()) and E.discriminant() % p
        ep = EllipticCurve(GF(p), E.a_invariants()); points = sorted(ep.points(), key=key)
        doubles = {key(2*Q): 2*Q for Q in points}
        labels = {v: 0 for v in doubles}; reps = [ep(0)]; dim = 0
        for Q in points:
            if key(Q) in labels: continue
            new = []
            for i, V in enumerate(reps):
                N = Q+V; new.append(N)
                for D in doubles.values():
                    v = key(N+D); assert v not in labels
                    labels[v] = i | (1 << dim)
            reps += new; dim += 1
        assert len(labels) == len(points) == len(doubles)*2**dim == row['order']
        assert dim == row['dimension']
        seen_odd |= bool(len(points) % 2)
        if p == frame['no_two_torsion_prime']: assert len(points) == 13 and p == 7
        groups.append((p, ep, labels, dim)); boundaries.append((sum(v[3] for v in groups), seen_odd))
    def code(Q):
        bits = []
        for p, ep, labels, dim in groups:
            den = lcm([v.denominator() for v in Q]); coords = [ZZ(v*den) for v in Q]
            common = gcd(coords); red = ep([int(v/common % p) for v in coords])
            label = labels[key(red)]
            bits.extend((label >> i) & 1 for i in range(dim))
        return vector(GF(2), bits)
    M = matrix(GF(2), [code(Q) for Q in basis]).transpose()
    N = matrix(GF(2), [code(Q) for Q in basis+[P]]).transpose()
    assert M.rank() == 17 and N.rank() == 18 and seen_odd
    assert frame['no_two_torsion_prime'] == 7
    for nrows, odd in boundaries[:-1]: assert not (odd and M[:nrows].rank() == 17)
    last = groups[-1][0]; assert last == 109
    good = [int(p) for p in prime_range(3, last+1)
            if all(a.denominator() % p for a in E.a_invariants()) and E.discriminant() % p]
    assert good == [row['prime'] for row in frame['records']]
    return dict(mask=47755, parameter='0', generic_rank=17, augmented_rank=18,
                generic_only_earliest_complete_prime_prefix=good,
                odd_order_certificate=dict(prime=7, order=13),
                quotient_matrix=[list(map(int, row)) for row in N.rows()],
                original_parameter=str(tau), alternate_specialized_rank='NOT_TESTED')


def replay():
    protocol = read(OUT/'generic-protocol.json'); hashes(protocol)
    assert protocol['point_searches'] == protocol['parameter_replacements'] == 0
    dependencies = [OUT/'generic-protocol.json', OUT/'incidence-protocol.json',
        OUT/'seed-generic-frame.json', OUT/'seed-result.json', CONICS/'independent-replay.json',
        ART/'det1092_norm8_singular_members_v1/independent-replay.json',
        *[OUT/('%s-%d.json' % (kind, mask)) for mask in MASKS for kind in ['map', 'incidence']]]
    save('replay-protocol.json', dict(checker_sha256=sha(Path(__file__)), seconds=25,
        inputs={str(p.relative_to(ROOT)): sha(p) for p in dependencies},
        method='Exact polynomial maps, supplied rational-root lattice proof checks, independent Sage finite quotients.'))
    parent = read(ART/'curve302_recovered_mw17_parent_v1.json')
    pencil = read(ART/'det1092_norm8_seed_cover_v2/generic.json')
    return dict(status='PASS_INDEPENDENT_BIFIBRATION_CONICS_AND_FIXED_M18_SEED',
        classification='Verified application and new deduction; not a prospective302 seed',
        geometry=geometry(pencil), sources=maps_and_incidence(parent, pencil, protocol),
        seed=specialized_seed(parent, protocol),
        replay_protocol_sha256=sha(OUT/'replay-protocol.json'),
        point_searches=0, full_atlas_runs=0, later_exceptional_inputs=0,
        boundary='Three fixed original bisections, not all rational curves. No alternate specialized rank or rank upper bound. First302 carrier is retrospective incidence input only.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    signal.alarm(25)
    result = replay()
    if args.write: save('independent-replay.json', result)
    else: assert result == read(OUT/'independent-replay.json')
    print(result['status'], flush=True)
    print('bidegrees', [v['bidegree'] for v in result['sources']], 'fixed seed rank', result['seed']['augmented_rank'], flush=True)
