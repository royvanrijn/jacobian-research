#!/usr/bin/env python3
"""Independent integer and finite-field replay; no Sage or point search."""
import argparse
from pathlib import Path
from math import prod
import retrospective as r

HERE = Path(__file__).resolve().parent
INPUT = r.OUT / 'rank_jump_soluble_quartet_compression_inputs_v1.json'
SUPPORT = r.OUT / 'rank_jump_collision_prime_lift_v1.json'
SOURCE = r.OUT / 'rank_jump_collision_defect_v2.json'
FAILED = r.OUT / 'rank_jump_collision_defect_v1.json'
OUTPUT = r.OUT / 'rank_jump_collision_defect_verification_v1.json'


def trim(a):
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a


def remainder(a, b, p):
    a = a.copy()
    while a != [0] and len(a) >= len(b):
        shift = len(a)-len(b); c = a[-1]*pow(b[-1], -1, p) % p
        for i in range(len(b)): a[i+shift] = (a[i+shift]-c*b[i]) % p
        trim(a)
    return a


def gcd(a, b, p):
    a = trim([x % p for x in a]); b = trim([x % p for x in b])
    while b != [0]: a, b = b, remainder(a, b, p)
    c = pow(a[-1], -1, p)
    return [x*c % p for x in a]


def evaluate(q, x): return sum(a*x**i for i, a in enumerate(q))


def legendre(a, p):
    z = pow(a % p, (p-1)//2, p)
    return 0 if z == 0 else 1 if z == 1 else -1


def root_count(g, p):
    if p == 2: return sum(evaluate(g, x) % p == 0 for x in range(p))
    if len(g) == 1: return 0
    if len(g) == 2: return 1
    assert len(g) == 3
    return 1 + legendre(g[1]**2-4*g[0]*g[2], p)


def compute():
    source = r.read(SOURCE); support = r.read(SUPPORT); inputs = r.read(INPUT)
    for doc in (source, support, r.read(FAILED)):
        for path, sha in doc['bindings'].items():
            assert r.digest((r.ROOT/path).read_bytes()) == sha
    assert all(x['execution']['status'] == 'FAILED' for x in r.read(FAILED)['rows'])
    summaries = []
    for case, record, old in zip(inputs['cases'], source['rows'], support['rows'], strict=True):
        assert record['execution']['status'] == 'COMPLETE'
        result = record['result']; old = old['result']
        assert result['id'] == case['id'] == old['id']
        qs = [list(map(int, c['form'])) for c in case['covers']]
        assert [x['prime'] for x in result['rows']] == old['collision_primes']
        masks = []; toggles = []; pruned = []; complex_rows = []
        for row in result['rows']:
            p = int(row['prime']); roots = set()
            assert [x['indices'] for x in row['pairs']] == [[j, i] for i in range(4) for j in range(i)]
            for pair in row['pairs']:
                i, j = pair['indices']; g = gcd(qs[i], qs[j], p)
                assert g == pair['gcd']
                rr = list(map(int, pair['finite_roots']))
                assert rr == sorted(set(rr)) and all(0 <= x < p for x in rr)
                assert len(rr) == root_count(g, p)
                assert all(evaluate(g, x) % p == 0 for x in rr)
                assert pair['infinity'] == (qs[i][2] % p == qs[j][2] % p == 0)
                roots.update(rr)
            clusters = [{'chart': 'finite', 'root': str(x),
                         'mask': sum(1 << i for i, q in enumerate(qs) if evaluate(q, x) % p == 0)}
                        for x in sorted(roots)]
            infinity = sum(1 << i for i, q in enumerate(qs) if q[2] % p == 0)
            if infinity.bit_count() >= 2: clusters.append({'chart': 'infinity', 'root': '0', 'mask': infinity})
            assert clusters == row['collision_clusters']
            upper = sorted({0} | {m for c in clusters for m in range(1, 16)
                                  if m & c['mask'] == m and m.bit_count() % 2 == 0})
            assert upper == row['necessary_parity_masks']
            assert row['prunable'] == (not clusters)
            if not clusters: pruned.append(str(p))
            incident = [(pair['indices'], int(e)) for pair in old['pair_resultants']
                        for pp, e in pair['prime_factorization'] if int(pp) == p]
            ordinary = p >= 5 and len(incident) == 1 and incident[0][1] == 1
            ordinary = ordinary and len(clusters) == 1 and clusters[0]['chart'] == 'finite'
            if ordinary:
                x = int(clusters[0]['root']); pair = incident[0][0]
                ordinary = all((qs[i][1]+2*qs[i][2]*x) % p for i in pair)
            assert bool(ordinary) == row['ordinary_simple_pair_collision']
            if not ordinary:
                complex_rows.append({'prime': str(p), 'necessary_parity_masks': upper})
                continue
            a = [(evaluate(qs[i], x)//p) % p for i in pair]
            b = [(qs[i][1]+2*qs[i][2]*x) % p for i in pair]
            assert (a[0]*b[1]-a[1]*b[0]) % p
            leading = b[0]*b[1]*prod(evaluate(qs[i], x) for i in range(4) if i not in pair) % p
            assert int(row['odd_pair_residue_count']) == (p-2-legendre(leading, p))//2 > 0
            odd = row['odd_pair_witness']; zero = row['zero_parity_witness']
            assert odd['status'] == zero['status'] == 'CERTIFIED'
            xx = int(odd['residue_mod_p_squared']); assert xx == x+p*odd['s']
            assert 0 <= odd['s'] < min(p, 64)
            values = [evaluate(q, xx) for q in qs]
            mask = sum(1 << i for i in pair); assert odd['mask'] == mask
            assert all(v % p == 0 and v % (p*p) != 0 if i in pair else v % p != 0
                       for i, v in enumerate(values))
            product_unit = prod(values)//(p*p) % p
            assert str(product_unit) == odd['product_unit'] and legendre(product_unit, p) == 1
            xx = int(zero['residue_mod_p']); assert 0 <= xx < min(p, 64)
            product_unit = prod(evaluate(q, xx) for q in qs) % p
            assert zero['mask'] == 0 and str(product_unit) == zero['product_unit']
            assert legendre(product_unit, p) == 1
            assert upper == [0, mask]
            masks.append(mask); toggles.append({'prime': str(p), 'exact_local_parity_image': upper})
        assert pruned == result['pruned_primes']
        assert result['retained_support'] == [x for x in old['collision_primes'] if x not in pruned]
        assert sorted(set(masks)) == result['certified_pair_masks'] == [3, 5, 6, 9, 10, 12]
        assert r.rank(masks) == result['native_mask_span_dimension'] == 3
        assert len(masks) == result['separate_prime_defect_witnesses']
        summaries.append({'id': case['id'], 'pruned_primes': pruned, 'exact_ordinary_local_images': toggles,
                          'other_prime_upper_bounds_only': complex_rows,
                          'native_parity_span_dimension': 3, 'independent_place_coordinates': len(toggles)})
    assert [s['independent_place_coordinates'] for s in summaries] == [13, 17, 18]
    return {'schema': 'rank-jump.collision-defect-verification.v1', 'status': 'PASS', 'rows': summaries,
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                         for p in (INPUT, SUPPORT, SOURCE, FAILED, Path(__file__), HERE/'retrospective.py')},
            'boundary': 'Exact local valuation images at ordinary primes. Separate local choices give no global rational-point existence or original Mordell-Weil rank.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check'])
    args = parser.parse_args(); result = compute()
    if args.mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS: all six pair defects; native span 3; exact local toggles at 13,17,18 primes')
