#!/usr/bin/env python3
"""Bounded exact geometry and square-character audit of three frozen quartets."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb, prod
from pathlib import Path
import retrospective as r
from cover_experiment import evaluate, mul, sqrtq, trim

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'SOLUBLE_QUARTET_COMPRESSION_PROTOCOL.json'
SOURCE = r.OUT / 'rank_jump_local_solubility_block_inputs_v1.json'
COHORT = r.OUT / 'rank_jump_solubility_first_v1.json'
INPUT = r.OUT / 'rank_jump_soluble_quartet_compression_inputs_v1.json'
OUTPUT = r.OUT / 'rank_jump_soluble_quartet_compression_v1.json'


def capture():
    source = r.read(SOURCE)
    cases = []
    for group in source['observed_groups']:
        if group['mask'].bit_count() == 4:
            cases.append({'id': group['source_id'], 'observed_parameter': group['published_parameter'],
                          'covers': [q for i, q in enumerate(source['covers']) if group['mask'] >> i & 1]})
    labels = ['orbit-030cb', 'orbit-03da0', 'orbit-07086', 'orbit-11278']
    cases.append({'id': 'obstructed_ABCD', 'observed_parameter': None,
                  'covers': [next(q for q in source['covers'] if q['label'] == label) for label in labels]})
    fields = ['source_id', 'published_parameter', 'observed_quotient_rank', 'generic_rank']
    rows = [{key: row[key] for key in fields} for row in r.read(COHORT)['rows']]
    assert len(cases) == 3 and len(rows) == 32
    r.write_new(INPUT, {'schema': 'rank-jump.soluble-quartet-compression-inputs.v1',
                       'cases': cases, 'cohort': rows,
                       'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                                    for p in (PROTOCOL, SOURCE, COHORT)}})


def remainder(f, g):
    f = trim(list(map(F, f))); g = trim(list(map(F, g)))
    while f != [0] and len(f) >= len(g):
        c = f[-1] / g[-1]; k = len(f) - len(g)
        for i, a in enumerate(g):
            f[k + i] -= c * a
        f = trim(f)
    return f


def gcd_degree(f, g):
    f, g = list(map(F, f)), list(map(F, g))
    while g != [0]:
        f, g = g, remainder(f, g)
    return len(f) - 1


def analyse():
    inp = r.read(INPUT)
    for path, sha in inp['bindings'].items():
        assert r.digest((r.ROOT / path).read_bytes()) == sha
    rows = []
    for case in inp['cases']:
        forms = [q['form'] for q in case['covers']]
        assert len(forms) == 4 and all(len(f) == 3 and f[-1] for f in forms)
        discriminants = [f[1] ** 2 - 4 * f[0] * f[2] for f in forms]
        assert all(discriminants)
        pair_gcds = [gcd_degree(f, g) for f, g in combinations(forms, 2)]
        assert pair_gcds == [0] * 6
        characters = []
        for mask in range(1, 16):
            polynomial = [F(1)]
            for i, f in enumerate(forms):
                if mask >> i & 1:
                    polynomial = mul(polynomial, f)
            derivative = [i * polynomial[i] for i in range(1, len(polynomial))]
            assert gcd_degree(polynomial, derivative) == 0
            size = mask.bit_count()
            assert len(polynomial) == 2 * size + 1
            characters.append({'mask': mask, 'weight': size, 'genus': size - 1,
                               'coefficients': list(map(str, polynomial))})
        assert sum(c['genus'] for c in characters) == 17
        fibres = []
        individual_witnesses = [None] * 4
        for old in inp['cohort']:
            t = F(old['published_parameter']); values = [evaluate(f, t) for f in forms]
            assert all(values), 'Branch specialization must be analysed separately'
            kernel = [0]; roots = {}
            for mask in range(1, 16):
                v = prod(values[i] for i in range(4) if mask >> i & 1)
                root = sqrtq(v)
                if root is not None:
                    kernel.append(mask); roots[str(mask)] = str(root)
            dim = r.rank(kernel)
            assert len(kernel) == 2 ** dim and all(a ^ b in kernel for a in kernel for b in kernel)
            for i in range(4):
                if 1 << i in kernel and individual_witnesses[i] is None:
                    individual_witnesses[i] = {'parameter': str(t), 'root': roots[str(1 << i)]}
            fibres.append({**old, 'values': list(map(str, values)), 'square_character_masks': kernel,
                           'product_roots': roots, 'radicand_squareclass_rank': 4 - dim,
                           'full_lift': dim == 4, 'all_four_product_square': 15 in kernel})
        assert all(individual_witnesses)
        observed = case['observed_parameter']
        if observed is not None:
            assert next(f for f in fibres if f['published_parameter'] == observed)['full_lift']
        rows.append({'id': case['id'], 'labels': [q['label'] for q in case['covers']],
                     'quadratic_discriminants': list(map(str, discriminants)),
                     'pairwise_gcd_degrees': pair_gcds, 'generic_squareclass_rank': 4,
                     'carrier_degree_over_parameter_line': 16, 'carrier_genus': 17,
                     'rational_gonality': 8, 'geometric_gonality': 8,
                     'minimum_geometric_degree_to_genus_one': 4,
                     'individual_conic_witnesses': individual_witnesses,
                     'character_quotients': characters, 'fibres': fibres,
                     'specialized_squareclass_rank_counts': dict(sorted(Counter(str(f['radicand_squareclass_rank']) for f in fibres).items())),
                     'full_lift_fibres': [f['source_id'] for f in fibres if f['full_lift']],
                     'proper_character_only_fibres': [f['source_id'] for f in fibres if len(f['square_character_masks']) > 1 and not f['full_lift']],
                     'all_four_product_false_positives': [f['source_id'] for f in fibres if f['all_four_product_square'] and not f['full_lift']]})
    genus_table = []
    for n in range(2, 5):
        genus = 1 + 2 ** (n - 1) * (n - 2)
        previous = 0 if n == 2 else 1 + 2 ** (n - 2) * (n - 3)
        assert genus - (2 * previous - 1) == 2 ** (n - 1)
        assert sum(comb(n, k) * (k - 1) for k in range(1, n + 1)) == genus
        genus_table.append({'radicands': n, 'genus': genus, 'forget_one_genus': previous,
                            'gonality': 2 ** (n - 1), 'minimum_degree_to_genus_one': 2 ** (n - 2)})
    return {'schema': 'rank-jump.soluble-quartet-compression.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                         for p in (INPUT, PROTOCOL, Path(__file__), HERE / 'retrospective.py', HERE / 'cover_experiment.py')},
            'layer': 'solubility', 'product_square_tests': 1440, 'rows': rows, 'tower_genus_table': genus_table,
            'theorem_boundary': 'Genus and map-degree conclusions use Riemann-Hurwitz and Castelnuovo-Severi as proved in the companion note. Exact branch conditions and rational conic witnesses are replayed here. The radicand squareclass rank is not an elliptic Selmer or MW rank.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['capture', 'build', 'check'])
    mode = parser.parse_args().mode
    if mode == 'capture':
        capture()
    else:
        result = analyse()
        if mode == 'build': r.write_new(OUTPUT, result)
        else: assert result == r.read(OUTPUT)
        for row in result['rows']:
            print(row['id'], 'genus', row['carrier_genus'], 'squareclass ranks', row['specialized_squareclass_rank_counts'],
                  'product false positives', row['all_four_product_false_positives'])
