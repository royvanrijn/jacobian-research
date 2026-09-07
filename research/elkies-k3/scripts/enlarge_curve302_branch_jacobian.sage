#!/usr/bin/env sage-python
"""Explicit half of the natural divisor, index2 enlargement and2-saturation.

Uses Stoll's linear halving equations, followed by an exact Cantor doubling
check. The known rational torsion and the2-primary group at53 prove that
the enlarged rank-one subgroup is2-saturated, without assuming full rank.
The |nH|<=1024 sieve retains every old |nD|<=512 word since D=2H.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import HyperellipticCurve, PolynomialRing, QQ, ZZ, matrix, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
SOURCE = ART/'elkies-k3-curve302-branch-jacobian-generator-v1.json'
HELPER = ROOT/'elkies-k3/scripts/certify_curve302_branch_jacobian_generator.sage'
OUT = ART/'elkies-k3-curve302-branch-jacobian-enlargement-v1.json'
LIMIT = 1024


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def coeff(poly):
    return list(map(str, poly.list()))


def compute():
    data = json.loads(SOURCE.read_text())
    for path, expected in data['input_sha256'].items():
        assert digest(ROOT/path) == expected
    helper = runpy.run_path(str(HELPER))
    ring = PolynomialRing(QQ, 'z')
    z = ring.gen()
    F = ring(data['odd_degree_chart']['polynomial_low_to_high'])
    A, B = [ring(p) for p in data['reduced_generator_u_v']]
    roots = list(map(QQ, data['torsion_finite_branch_roots']))
    c = F.leading_coefficient()
    assert A.degree() == 3 and A.gcd(A.derivative()).degree() == A.gcd(F).degree() == 0
    values = [(-c)**3*A(root) for root in roots]
    assert all(value != 0 and value.is_square() for value in values)
    square_roots = [value.sqrt() for value in values]
    S = ring.lagrange_polynomial(list(zip(roots, square_roots)))
    assert (S*S-(-c)**3*A) % F == 0
    # Stoll Prop5.1, coprime degree3 case: v=ws modF, v=ub modA,
    # deg(u)<=1, deg(v)<=4, deg(w)<=3. These are ten linear equations.
    columns = []
    for i in range(2):
        columns.append([QQ(0)]*7+list((-(z**i*B) % A).padded_list(3)))
    for i in range(5):
        columns.append(list((z**i % F).padded_list(7))+list((z**i % A).padded_list(3)))
    for i in range(4):
        columns.append(list((-(z**i*S) % F).padded_list(7))+[QQ(0)]*3)
    system = matrix(QQ, columns).transpose()
    kernel = system.right_kernel()
    assert kernel.dimension() == 1
    solution = kernel.basis()[0]
    u, v, w = ring(list(solution[:2])), ring(list(solution[2:7])), ring(list(solution[7:]))
    normalizer = w.leading_coefficient()
    u, v, w = u/normalizer, v/normalizer, w/normalizer
    assert w.degree() == 3 and u.gcd(w) == 1
    assert u*u*F == v*v-(-c)**3*A*w*w
    r = (-v*u.inverse_mod(w)) % w
    group = HyperellipticCurve(F).jacobian()(QQ)
    D = group([A, B])
    raw_half = group([w, r])
    assert 2*raw_half == D
    masks = [n for n in range(128) if n.bit_count() <= 3]
    candidates = []
    halves = []
    for mask in masks:
        T = group([prod((z-roots[i] for i in range(7) if mask >> i & 1), ring(1)), ring(0)])
        H = raw_half+T
        assert 2*H == D
        bits = max(max(abs(q.numerator()).nbits(), q.denominator().nbits()) for q in H[0].list()+H[1].list())
        candidates.append({'torsion_mask': mask, 'coefficient_max_bits': int(bits), 'mumford_degree': int(H[0].degree())})
        halves.append(H)
    selected = min(range(64), key=lambda i: (candidates[i]['coefficient_max_bits'], candidates[i]['torsion_mask']))
    H = halves[selected]
    point_candidates = [{'torsion_mask': row['torsion_mask'], 'mumford': [coeff(point[0]), coeff(point[1])]}
                        for row, point in zip(candidates, halves) if point[0].degree() <= 1]
    finite = []
    for old in data['finite_tables']:
        table = helper['finite_table'](F, H[0], H[1], roots, old['prime'], ZZ(old['jacobian_order_bound']))
        previous = {(r['n_residue'], r['torsion_mask']) for r in old['complete_word_image']}
        current = {(r['n_residue'], r['torsion_mask']) for r in table['complete_word_image']}
        for n in range(old['D_order']):
            for mask in masks:
                assert ((n, mask) in previous) == ((2*n % table['D_order'], mask) in current)
        finite.append(table)
        print('HALF_FINITE_IMAGE', table['prime'], table['D_order'], len(current), flush=True)
    at53 = next(t for t in finite if t['prime'] == 53)
    two_valuation = ZZ(at53['jacobian_order_bound']).valuation(2)
    h_valuation = ZZ(at53['D_order']).valuation(2)
    assert two_valuation == 9 and h_valuation == 4
    # Six independent rational two-torsion directions remain independent
    # at53. The sum of the six2-primary invariant-factor exponents is9,
    # so their maximum is at most9-5=4. H attains this maximum.
    assert two_valuation-(6-1) == h_valuation
    allowed = [{(r['n_residue'], r['torsion_mask']) for r in t['complete_word_image']} for t in finite]
    excluded = [0]*len(finite)
    survivors = []
    for n in range(-LIMIT, LIMIT+1):
        for mask in masks:
            for i, (table, allowed_here) in enumerate(zip(finite, allowed)):
                if (n % table['D_order'], mask) not in allowed_here:
                    excluded[i] += 1
                    break
            else:
                survivors.append({'n': n, 'torsion_mask': mask})
    known = [r for r in survivors if r['n'] == 0 and r['torsion_mask'].bit_count() <= 1]
    pending = [r for r in survivors if r not in known]
    assert sum(excluded)+len(survivors) == (2*LIMIT+1)*64
    return {'schema': 'curve302.branch-jacobian-enlargement.v1',
            'status': 'INDEX_TWO_ENLARGEMENT_PROVED_TWO_SATURATED' if not pending and not point_candidates else 'INDEX_TWO_ENLARGEMENT_WITH_POINT_CANDIDATES',
            'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), SOURCE, HELPER]},
            'halving_reference': 'Michael Stoll, Chabauty without the Mordell-Weil group, Proposition5.1',
            'halving_reference_url': 'https://www.mathe2.uni-bayreuth.de/stoll/papers/ratpts-selmer-2016-12-07.pdf',
            'kummer_values': list(map(str, values)), 'kummer_square_roots': list(map(str, square_roots)),
            'interpolated_square_root_S': coeff(S), 'linear_system_rank': int(system.rank()),
            'linear_solution_u_v_w': [coeff(u), coeff(v), coeff(w)],
            'raw_half_mumford': [coeff(raw_half[0]), coeff(raw_half[1])],
            'all_half_height_records': candidates, 'selected_half': candidates[selected],
            'selected_half_mumford': [coeff(H[0]), coeff(H[1])],
            'all64_halves_double_to_D': True, 'index_of_old_group_in_new': 2,
            'rank_lower_bound': 1,
            'two_saturation': {'prime': 53, 'jacobian_order': at53['jacobian_order_bound'],
                               'two_torsion_rank': 6, 'two_primary_order_valuation': int(two_valuation),
                               'max_two_primary_exponent': int(h_valuation),
                               'H_reduction_order': at53['D_order'],
                               'proof': 'Every odd multiple ofH plus rational two-torsion has2-primary order16 at53, whereas twice any element has2-primary order at most8. Thus none can have a rational half. If2X=2nH+T, nonzeroT would give rational4-torsion, already excluded; whenT=0, X-nH is rational two-torsion and lies in the subgroup. Therefore <H,J(Q)[2]> is2-saturated inJ(Q).'},
            'finite_tables': finite,
            'old_word_embedding': 'nD+T=(2n)H+T; every old prime-table membership is replayed under this embedding.',
            'limits': {'abs_n_max': LIMIT, 'torsion_classes': 64, 'workers': 1, 'seconds': 300},
            'word_count': (2*LIMIT+1)*64, 'exclusions_by_prime': excluded,
            'known_Weierstrass_survivors': known, 'pending_nonWeierstrass_words': pending,
            'curve_points_among_halves': point_candidates,
            'boundary': 'This is an explicit2-saturated rank-one subgroup, not a full Mordell-Weil basis or an exact rank determination. Other independent directions, odd-index enlargements and global nonWeierstrass points remain UNKNOWN. No new cover or302 parent is constructed.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    signal.alarm(300)
    result = compute()
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.check:
        assert OUT.read_text() == text
    else:
        assert not OUT.exists(), 'Use --check; preserve the certificate'
        OUT.write_text(text)
    print('PASS_HALF_AND_TWO_SATURATION', result['exclusions_by_prime'], result['pending_nonWeierstrass_words'], flush=True)
