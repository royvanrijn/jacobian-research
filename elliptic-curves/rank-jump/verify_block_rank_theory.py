#!/usr/bin/env python3
"""Exact finite applications of J7-J10; no new points, classes, or pairings."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path

import block_rank_theory as b
import retrospective as r
import verify_paired_quartet_relations as old

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT/'artifacts/generated-results/elliptic-curves'
OUTPUT = OUT/'rank_jump_block_rank_theory_v1.json'


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compute():
    require(__debug__, 'retained verifiers require Python without -O')
    require(old.verify() == read(old.OUTPUT), 'independent ambient witness replay differs')
    paths = [Path(__file__), HERE/'block_rank_theory.py', HERE/'mechanism_theory.py',
             HERE/'test_block_rank_theory.py', HERE/'INDEPENDENT_BLOCKS_AND_RANK_OBSTRUCTIONS.md',
             HERE/'SEARCH_THEOREM_GATES_V2.json', Path(old.__file__), old.INPUT, old.SOURCE, old.OUTPUT,
             HERE/'retrospective.py']
    quartet_rows = []
    for case, source in zip(read(old.INPUT)['cases'], read(old.SOURCE)['rows'], strict=True):
        result = source['result']
        require(case['id'] == result['id'], 'case mismatch')
        model, generic = r.short(case['model'], case['basis'][:17])
        _, points = r.short(case['model'], [x['point'] for x in case['lifts']])
        primes = [x['prime'] for x in case['rank_certificate']['signatures']]
        tp = case['rank_certificate']['no_rational_2_torsion_prime']
        valid_primes = set(r.primes(max(primes+[tp])))
        require(all(p in valid_primes for p in primes+[tp]), 'nonprime reduction modulus')
        require(r.roots_at(model[3], model[4], tp) == (), 'rational 2-torsion not excluded')
        blocks = [(p, r.roots_at(model[3], model[4], p)) for p in primes]
        require(all(roots is not None for _, roots in blocks), 'bad reduction in finite signatures')
        width = sum(len(roots) for _, roots in blocks)
        sigs = [r.point_signature(model, P, blocks) for P in generic+points]
        columns = [[v >> j & 1 for j in range(width)] for v in sigs]
        # Replay the relation directly on just these 21 points.
        aa = F(model[3])
        for coeffs, word in zip(result['kernel_integer_vectors'], result['kernel_generic_coordinates'], strict=True):
            total = None
            for n, P in zip(coeffs, points, strict=True):
                total = old.add(aa, total, old.mul(aa, int(n), tuple(map(F, P))))
            for n, P in zip(word, generic, strict=True):
                require(F(n).denominator == 1, 'nonintegral relation word')
                total = old.add(aa, total, old.mul(aa, -int(F(n)), tuple(map(F, P))))
            require(total is None, 'native relation failed exact group law')
        sandwich = b.signature_sandwich(17, 0, columns[:17], columns[17:], result['kernel_integer_vectors'])
        require(sandwich['exact_quotient_rank'] == 3, 'J7 does not close quartet rank')
        # All three unordered partitions into two pairs; oracle coordinates
        # are explicitly retrospective and used only after ambient replay.
        coordinates = result['coordinates_in_witness_basis']
        dim = result['basis_rank']
        base = [[int(i == j) for j in range(dim)] for i in range(17)]
        partitions = []
        for j in (1, 2, 3):
            first = [0, j]
            second = [i for i in range(4) if i not in first]
            overlap = b.block_overlap(base, [coordinates[i] for i in first], [coordinates[i] for i in second])
            require(overlap['union_rank'] == 3, 'union rank mismatch')
            partitions.append({'first_indices': first, 'second_indices': second, **overlap})
        quartet_rows.append({'id': case['id'], 'point_count_in_direct_certificate': 21,
                             'no_rational_two_torsion_prime': tp, 'finite_signature_width': width,
                             'sandwich': sandwich, 'pair_partitions': partitions})

    block_path = OUT/'rank_jump_block_inputs_v1.json'
    block_input = read(block_path)
    paths.append(block_path)
    ct_rows = []
    for case in block_input['ct']:
        source_path = ROOT/case['source']
        require(digest(source_path) == block_input['bindings'][case['source']], 'CT source hash differs')
        raw = source_path.read_bytes()
        original = json.loads(gzip.decompress(raw) if source_path.suffix == '.gz' else raw)
        full = original['arithmetic']['matrix'] if case['u'] == -1 else original['ct']['matrix']
        require(full == case['matrix'], 'matrix projection differs from arithmetic source')
        b.alternating(full)
        n = len(full)
        full_rank = r.rank([r.pack(row) for row in full])  # independent packed elimination
        require(full_rank == b.mod_rank(full), 'independent matrix ranks disagree')
        prefixes = []
        for d in range(n+1):
            record = b.radical_partner_bound([row[:d] for row in full[:d]],
                                             [row[d:] for row in full[:d]], n-d)
            require(record['certified_pairing_rank_lower_bound'] <= full_rank, 'J9 exceeds actual rank')
            if n-d == 1:
                require(record['certified_pairing_rank_lower_bound'] == full_rank, 'one-partner exact rank identity')
            prefixes.append(record)
        ct_rows.append({'u': case['u'], 'retained_selmer_subspace_dimension': n,
                        'retained_pairing_rank': full_rank, 'prefixes': prefixes,
                        'target_gate_without_absolute_selmer_upper': b.rank_exclusion(None, 0, full_rank, 20)})
        paths.append(source_path)

    pair_path = OUT/'rank_jump_degree_one_relation_panel_v1.json'
    triple_path = OUT/'rank_jump_low_degree_triple_panel_v1.json'
    triple_verified = OUT/'rank_jump_low_degree_triple_panel_verification_v1.json'
    pairs = read(pair_path)
    triples = read(triple_path)
    require(read(triple_verified)['status'] == 'PASS', 'missing verified triple result')
    for doc in (pairs, triples, read(triple_verified)):
        for path, sha in doc['bindings'].items():
            require(digest(ROOT/path) == sha, 'relation source binding differs: '+path)
    pair = next(x for x in pairs['rows'] if x['parameter'] == '3/8')
    triple = next(x for x in triples['rows'] if x['parameter'] == '3/8')
    require(pair['labels'] == triple['labels'], 'pair/triple label mismatch')
    n = len(pair['labels'])
    def rows(relations):
        out = []
        for relation in relations:
            row = [0]*n
            for i, sign in zip(relation['indices'], relation['signs'], strict=True):
                row[i] = sign
            out.append(row)
        return out
    rp = rows(pair['signed_quotient_edges'])
    rt = rows(triple['canonical_signed_relations'])
    pair_rank, triple_rank, combined_rank = b.qrank(rp), b.qrank(rt), b.qrank(rp+rt)
    require((pair_rank, triple_rank, combined_rank) == (7, 11, 13), 'combined relation ranks')
    require(triple['combined_relation_rank'] == combined_rank, 'triple combined bound')
    paths += [pair_path, triple_path, triple_verified]
    return {'schema': 'rank-jump.block-rank-theory.v1', 'status': 'PASS',
            'quartet_independence_and_overlap': quartet_rows,
            'radical_partner_replays': ct_rows,
            'combined_relation_control': {'parameter': '3/8', 'native_labels': n,
               'pair_rank': pair_rank, 'triple_rank': triple_rank, 'combined_rank': combined_rank,
               'relation_span_intersection_dimension': pair_rank+triple_rank-combined_rank,
               'new_constraints_beyond_pairs': combined_rank-pair_rank,
               'native_quotient_upper_bound': n-combined_rank},
            'bindings': {str(p.relative_to(ROOT)): digest(p) for p in paths},
            'scope': 'Exact finite applications of written J7-J10. Seven retained arithmetic matrices, not seven fresh CT computations. Absolute Selmer bounds remain UNKNOWN; no new curve rank or class construction.',
            'separate_arithmetic_replay': ['sage -python elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage --check',
                                          'sage -python elliptic-curves/rank-jump/verify_low_degree_triple_panel.py check']}


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
        require(read(OUTPUT) == result, 'certificate differs; investigate rather than overwrite')
    print('PASS: two direct 21-point quotient proofs; six block partitions; seven CT matrices; combined relation ranks 7,11,13')
