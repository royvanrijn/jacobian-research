#!/usr/bin/env python3
"""Small exact synthesis replay. No search, factorization, Sage or heights.

General theorems are proved in the note. This verifies their retained numerical
applications and binds inherited geometric evidence without overstating replay.
"""
import argparse
import hashlib
import json
from math import prod
from pathlib import Path

import mechanism_theory as t
import verify_collision_defect as collision
import verify_paired_quartet_relations as quartet

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / 'artifacts/generated-results/elliptic-curves'
OUTPUT = OUT / 'rank_jump_mechanism_theory_v1.json'


def read(path):
    return json.loads(path.read_text())


def require(value, message):
    if not value:
        raise ValueError(message)


def compute():
    # These independent legacy verifiers use assertions as proof checks.
    require(__debug__, 'Replay must not run with Python -O')
    require(quartet.verify() == read(quartet.OUTPUT), 'quartet exact replay differs')
    require(collision.compute() == read(collision.OUTPUT), 'collision exact replay differs')
    panel_path = OUT / 'rank_jump_degree_one_relation_panel_v1.json'
    panel_verification = OUT / 'rank_jump_degree_one_relation_panel_verification_v1.json'
    input_path = OUT / 'rank_jump_soluble_quartet_compression_inputs_v1.json'
    support_path = OUT / 'rank_jump_collision_prime_lift_v1.json'
    source_path = OUT / 'rank_jump_collision_defect_v2.json'
    relations_path = OUT / 'rank_jump_paired_quartet_relations_v1.json'
    panel = read(panel_path)
    inherited = read(panel_verification)
    # Hash validation is not a rerun of the Sage lattice/geometric proof.
    for doc in (panel, inherited):
        for path, expected in doc['bindings'].items():
            require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,
                    'panel input binding differs: ' + path)
    graph_rows = []
    for row in panel['rows']:
        edges = [(*e['indices'], *e['signs']) for e in row['signed_quotient_edges']]
        result = t.signed_relation_bound(len(row['labels']), edges)
        require(result['relation_rank'] == row['signed_relation_rank_over_Q'], 'panel relation rank')
        require(result['quotient_upper_bound'] == row['native_quotient_rank_upper_bound_from_pair_relations'], 'panel bound')
        graph_rows.append({'parameter': row['parameter'], 'edge_count': len(edges), **result})
    require(len(graph_rows) == 165, 'frozen panel length')
    require(sum(row['edge_count'] for row in graph_rows) == 18, 'frozen panel edge count')
    losses = []
    for record in read(relations_path)['rows']:
        r = record['result']
        dimension = r['basis_rank']
        generic = [[int(i == j) for j in range(dimension)] for i in range(r['generic_rank'])]
        q = t.quotient_dimension(generic, r['coordinates_in_witness_basis'])
        require(q == r['exact_quotient_rank'] == 3, 'quartet quotient rank')
        k = len(r['labels'])  # independent singleton characters, written J3 + pinned maps
        losses.append({'id': r['id'], 'generic_marked_increment': k,
                       'specialized_quotient_dimension': q, 'relative_kernel_dimension': k-q,
                       'witness_directions_outside_marked_mechanism': dimension-r['generic_rank']-q})
    collisions = []
    for case, support, source in zip(read(input_path)['cases'], read(support_path)['rows'],
                                     read(source_path)['rows'], strict=True):
        forms = [c['form'] for c in case['covers']]
        support = support['result']
        source = source['result']
        require(case['id'] == support['id'] == source['id'], 'case alignment')
        exact_resultants = []
        for pair in support['pair_resultants']:
            i, j = pair['indices']
            actual = t.resultant(forms[i], forms[j])
            require(actual == int(pair['resultant']), 'pair resultant')
            exact_resultants.append(actual)
        require(abs(prod(exact_resultants)) == t.quadratic_support(forms), 'total resultant support')
        clusters = [cluster['mask'] for row in source['rows'] for cluster in row['collision_clusters']]
        minimum = t.minimum_square_tests(4, clusters)
        require(minimum['minimum_tests'] == 3 and minimum['pair_masks'] == [3, 5, 6, 9, 10, 12], 'complete collision graph')
        count = prod(len(row['necessary_parity_masks']) for row in source['rows'])
        result = {'id': case['id'], 'collision_prime_count': len(support['collision_primes']),
                  'valuation_tuple_upper_bound': count,
                  'positivity_on_rational_product_points': support['real_lift_map_surjective'],
                  'minimum_tests_from_collision_bounds_given_positivity': minimum}
        if case['observed_parameter'] is not None:
            parameter = t.F(case['observed_parameter'])
            lift = t.lift_at_parameter(forms, parameter.numerator, parameter.denominator)
            require(lift['status'] == 'NATIVE_LIFT', 'retained rational native lift')
            result['retained_parameter'] = str(parameter)
            result['retained_lift'] = lift
        else:
            result['retained_lift'] = 'No supplied rational product point; UNKNOWN by this gate.'
        collisions.append(result)
    require([x['valuation_tuple_upper_bound'] for x in collisions] == [262144, 18874368, 301989888], 'mask product counts')
    paths = [Path(__file__), HERE/'mechanism_theory.py', HERE/'test_mechanism_theory.py',
             HERE/'RANK_JUMP_MECHANISM_THEOREMS.md', HERE/'SEARCH_THEOREM_GATES.json',
             panel_path, panel_verification, input_path, support_path, source_path, relations_path,
             quartet.OUTPUT, collision.OUTPUT, Path(quartet.__file__), Path(collision.__file__)]
    return {'schema': 'rank-jump.mechanism-theory.v1', 'status': 'PASS',
            'quartet_specialization_accounting': losses,
            'signed_relation_panel': graph_rows, 'collision_test_compression': collisions,
            'bindings': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'replayed': ['quartet rational point identities and finite Kummer independence',
                         'all ordinary collision parity witnesses and necessary masks',
                         'signed graph and rational row rank on all 165 panel addresses',
                         'pair resultants and retained homogeneous native roots',
                         'minimum fixed square tests from all collision clusters'],
            'inherited_not_recomputed_here': ['generic native maps, lattice minimum and panel point identities: run verify_degree_one_relation_panel.py check',
                                             'real sign surjectivity: run verify_collision_prime_lift.py check',
                                             'generic geometry: written J3 and existing carrier certificates'],
            'boundary': 'No full specialized rank, global product-curve image, new search, prospective enrichment claim, formal proof-assistant verification, or external review.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['build', 'check'])
    args = parser.parse_args()
    result = compute()
    if args.mode == 'build':
        with OUTPUT.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    else:
        require(read(OUTPUT) == result, 'certificate differs; preserve it and investigate')
    print('PASS: 165 signed graphs, two quotient losses, three collision test bounds; native roots and independent point certificates replayed')
