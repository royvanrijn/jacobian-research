#!/usr/bin/env sage-python
"""Natural nontorsion divisor and a bounded Jacobian-word point sieve.

The quartic identity already supplies a rational degree-four divisor. Reduce
it on the frozen odd-degree chart and prove it is not2-torsion. Since all
rational torsion has exponent2, this proves positive rank and rules out the
previously proposed rank-zero gate. Test nD+T, |n|<=512 and all64 rational
two-torsion classes, at47 and53. No full MW-group coverage is asserted.
"""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import runpy
import signal
from sage.all import GF, HyperellipticCurve, PolynomialRing, QQ, ZZ, prod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
CARRIER = ART/'elkies-k3-curve302-rational-branch-carrier-v1.json'
PROTOCOL = ART/'elkies-k3-curve302-branch-carrier-protocol-v1.json'
PROBE = ART/'elkies-k3-curve302-branch-carrier-probe-v1.json'
TORSION = ART/'elkies-k3-curve302-branch-jacobian-simple-v1.json'
OUT = ART/'elkies-k3-curve302-branch-jacobian-generator-v1.json'
LIMIT = 512


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def coeff(poly):
    return list(map(str, poly.list()))


def point_key(point):
    return tuple(map(int, point[0].list())), tuple(map(int, point[1].list()))


def finite_table(F, U, V, roots, p, order_bound):
    ring = PolynomialRing(GF(p), 'z')
    z = ring.gen()
    f, u, v = ring(F), ring(U), ring(V)
    assert f.degree() == 7 and f.gcd(f.derivative()).degree() == 0
    assert u.degree() == U.degree()
    group = HyperellipticCurve(f).jacobian()(GF(p))
    D = group([u, v])
    zero = group(0)
    assert order_bound*D == zero
    order = ZZ(order_bound)
    for prime, _ in ZZ(order_bound).factor():
        while order % prime == 0 and (order//prime)*D == zero:
            order //= prime
    assert order <= 10000, 'Declared cyclic-table limit'
    multiples = {}
    P = zero
    for n in range(int(order)):
        assert point_key(P) not in multiples
        multiples[point_key(P)] = n
        P += D
    assert P == zero
    finite_roots = list(map(GF(p), roots))
    assert len(set(finite_roots)) == 7 and all(f(r) == 0 for r in finite_roots)
    torsion = {}
    for size in range(4):
        for subset in combinations(range(7), size):
            mask = sum(1 << i for i in subset)
            T = group([prod((z-finite_roots[i] for i in subset), ring(1)), ring(0)])
            assert 2*T == zero
            torsion[mask] = T
    assert len(torsion) == len({point_key(T) for T in torsion.values()}) == 64
    curve_points = [(zero, {'infinity': True})]
    for x in GF(p):
        value = f(x)
        if value.is_square():
            for y in value.sqrt(all=True):
                curve_points.append((group([z-x, ring(y)]), {'x': int(x), 'y': int(y)}))
    records = []
    for mask, T in sorted(torsion.items()):
        for point, coordinates in curve_points:
            residue = multiples.get(point_key(point-T))
            if residue is not None:
                records.append({'n_residue': residue, 'torsion_mask': mask, 'curve_point': coordinates})
    assert len(records) == len({(r['n_residue'], r['torsion_mask']) for r in records})
    return {'prime': p, 'jacobian_order_bound': int(order_bound), 'D_order': int(order),
            'D_mumford': [coeff(u), coeff(v)],
            'twice_D_mumford': [coeff((2*D)[0]), coeff((2*D)[1])],
            'twice_D_is_nonzero': bool(2*D != zero),
            'curve_point_count': len(curve_points), 'complete_word_image': records}


def control():
    helper = runpy.run_path(str(ROOT/'elkies-k3/scripts/probe_curve302_branch_carrier.sage'))
    ring = helper['R']
    z = ring.gen()
    f = 40320*prod(z-i for i in range(8))
    row = helper['chart_pool'](f, list(map(QQ, range(8))))['selected'][0]
    a, b, c, d = map(QQ, row['matrix_a_b_c_d'])
    scale = QQ(row['square_scale'])
    F = ring(row['polynomial_low_to_high'])
    x = (b-8*d)/(8*c-a)
    y = QQ(40320)*(c*x+d)**4/scale
    assert F(x) == y*y and y != 0
    roots = [(b-m*d)/(m*c-a) for m in map(QQ, range(8)) if m*c != a]
    for p in [11, 13, 19, 23, 31]:
        if any(q.denominator() % p == 0 for q in [x, y]+roots):
            continue
        fp = F.change_ring(GF(p))
        if fp.degree() != 7 or fp.gcd(fp.derivative()).degree() != 0 or GF(p)(y) == 0:
            continue
        order = ZZ(helper['pari'](fp).hyperellcharpoly().subst('x', 1))
        result = finite_table(F, z-x, ring(y), roots, p, order)
        assert any(r['n_residue'] == 1 and r['torsion_mask'] == 0 and r['curve_point'].get('y', 0) != 0
                   for r in result['complete_word_image'])
        return {'curve': 'y²=40320*product(x-i,i=0..7)', 'point': ['8', '40320'],
                'positive_nontorsion_label_not_required': True,
                'nonWeierstrass_n1_T0_retained': True, 'finite_table': result}
    raise ArithmeticError('No good positive-control prime in declared list')


def compute():
    files = [Path(__file__), CARRIER, PROTOCOL, PROBE, TORSION,
             ROOT/'elkies-k3/scripts/probe_curve302_branch_carrier.sage']
    stamps = {str(p.relative_to(ROOT)): digest(p) for p in files}
    data = {}
    for p in [CARRIER, PROTOCOL, PROBE, TORSION]:
        data[p] = json.loads(p.read_text())
        for path, expected in data[p]['input_sha256'].items():
            assert digest(ROOT/path) == expected
    carrier = data[CARRIER]
    row = data[PROTOCOL]['chart_pool']['selected'][0]
    ring = PolynomialRing(QQ, 'z')
    z = ring.gen()
    a, b, c, d = map(QQ, row['matrix_a_b_c_d'])
    scale = QQ(row['square_scale'])
    F = ring(row['polynomial_low_to_high'])
    dc = ring(carrier['D4_low_to_high'])
    db = ring(carrier['branch_discriminant_in_u']['linear'])
    transform = lambda g: sum(g[i]*(a*z+b)**i*(c*z+d)**(4-i) for i in range(5))
    initial_u = transform(dc).monic()
    initial_v = (transform(db)/(4*scale)) % initial_u
    quotient, remainder = (F-initial_v**2).quo_rem(initial_u)
    assert remainder == 0
    U = quotient.monic()
    V = (-initial_v) % U
    assert U.degree() == 3 and V.degree() == 2 and (F-V*V) % U == 0
    group = HyperellipticCurve(F).jacobian()(QQ)
    D = group([U, V])
    assert D != -D
    roots = []
    indices = []
    for idx, m in enumerate(map(QQ, carrier['R8_rational_roots'])):
        if m*c == a:
            assert idx == row['root_triple_infinity_zero_one'][0]
            continue
        root = (b-m*d)/(m*c-a)
        assert F(root) == 0
        roots.append(root)
        indices.append(idx)
    assert len(roots) == 7
    positive = control()
    print('PASS_NONWEIERSTRASS_CONTROL', positive['finite_table']['prime'], flush=True)
    tables = []
    for arithmetic in data[PROBE]['jacobian_arithmetic']['records']:
        p = arithmetic['prime']
        table = finite_table(F, U, V, roots, p, ZZ(arithmetic['jacobian_order']))
        assert table['curve_point_count'] == arithmetic['extension_point_counts_n1_n2_n3'][0]
        assert table['twice_D_is_nonzero']
        tables.append(table)
        print('FINITE_WORD_IMAGE', p, table['D_order'], len(table['complete_word_image']), flush=True)
    masks = sorted({r['torsion_mask'] for r in tables[0]['complete_word_image']})
    # All64 torsion labels, including those never occurring in the Abel image.
    masks = [n for n in range(128) if n.bit_count() <= 3]
    assert len(masks) == 64
    allowed = [{(r['n_residue'], r['torsion_mask']) for r in t['complete_word_image']} for t in tables]
    survivors = []
    excluded = [0]*len(tables)
    for n in range(-LIMIT, LIMIT+1):
        for mask in masks:
            for idx, (table, allowed_here) in enumerate(zip(tables, allowed)):
                if (n % table['D_order'], mask) not in allowed_here:
                    excluded[idx] += 1
                    break
            else:
                survivors.append({'n': n, 'torsion_mask': mask})
    known = []
    pending = []
    for record in survivors:
        if record['n'] == 0 and record['torsion_mask'].bit_count() <= 1:
            known.append(record)
        else:
            pending.append(record)
    total = (2*LIMIT+1)*64
    assert sum(excluded)+len(survivors) == total
    return {'schema': 'curve302.branch-jacobian-generator.v1',
            'status': 'POSITIVE_JACOBIAN_RANK_AND_COMPLETE_BOUNDED_WORD_SIEVE' if not pending else 'POSITIVE_JACOBIAN_RANK_WITH_PENDING_WORDS',
            'input_sha256': stamps,
            'limits': {'abs_n_max': LIMIT, 'torsion_classes': 64, 'cyclic_table_order_cap': 10000,
                       'workers': 1, 'seconds': 300, 'rational_large_multiple_computation': False},
            'odd_degree_chart': row,
            'quartic_divisor_u_v': [coeff(initial_u), coeff(initial_v)],
            'cantor_reduction_quotient': coeff(quotient),
            'reduced_generator_u_v': [coeff(U), coeff(V)],
            'generator_not_equal_to_negative': True,
            'rank_lower_bound': 1,
            'nontorsion_proof': 'Reduced Mumford V is nonzero, so D differs from-D. Its finite reductions also have nonzero double. The complete rational torsion group has exponent2; D is therefore nontorsion.',
            'torsion_finite_branch_roots': list(map(str, roots)),
            'torsion_root_indices_in_carrier': indices,
            'torsion_labels': 'Masks on the seven listed finite roots with at most3 set bits. Larger subsets are replaced by their seven-root complement; the sum of all seven Weierstrass classes is zero.',
            'positive_control': positive, 'finite_tables': tables,
            'word_count': total, 'exclusions_by_prime': excluded,
            'known_Weierstrass_survivors': known, 'pending_nonWeierstrass_words': pending,
            'boundary': 'The rank-zero premise for this genus-three Jacobian is false. Its exact rank, saturation of the subgroup generated byD and torsion, and global curve rational points remain UNKNOWN. The finite sieve only concerns |n|<=512 in that explicit subgroup. No nonWeierstrass point, new cover or alternative parent is claimed.'}


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
    print('PASS_JACOBIAN_GENERATOR', result['exclusions_by_prime'], result['pending_nonWeierstrass_words'], flush=True)
