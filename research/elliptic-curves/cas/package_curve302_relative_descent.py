#!/usr/bin/env python3
"""Deterministic manifest and binary-rank audit for the relative descent result."""
import argparse
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUTPUT = ART / 'curve302_relative_descent_manifest_v1.json'
SOURCE_NAMES = [
    'curve302_relative_ideal_anatomy.sage',
    'curve302_even_ideal_extension.sage',
    'complete_curve302_even_ideal_extension.sage',
    'probe_curve302_even_artin_kernel.sage',
    'verify_curve302_relative_ideal_anatomy.sage',
    'curve302_seed_norm_artin_obstruction.sage',
    'verify_det1092_seed_artin.sage',
    'rank_jump_matched_descent_panel.py',
    'complete_rank_jump_matched_dyadic.sage',
    'verify_rank_jump_matched_descent_panel.sage',
    'package_curve302_relative_descent.py',
    'r17_60_arithmetic.py',
    'compact_atlas_specialization.py',
    'certify_compact_r17_candidates.py',
    'memory_rank_certificate.py',
    'select_r17_60_panel.py',
    'v3_warm_support.py',
    'research_runtime/local_kummer.py',
]
ARTIFACT_NAMES = [
    'curve302_descent_anatomy_manifest_v1.json',
    'curve302_relative_ideal_anatomy_v1.json',
    'curve302_even_ideal_extension_v1.json',
    'curve302_even_ideal_extension_complete_v1.json',
    'curve302_seed_norm_artin_obstruction_v1.json',
    'curve302_recovered_mw17_parent_v1.json',
    'det1092_seed_artin_v2/independent-replay.json',
    'rank_jump_matched_descent_panel_protocol_v1.json',
    'rank_jump_matched_descent_panel_v1.json',
    'rank_jump_matched_descent_dyadic_completion_v1.json',
    'rank_jump_matched_descent_panel_verification_v1.json',
    'compact_six_r17_atlas_v1.json',
]


def read(p):
    return json.loads(p.read_text())


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def rank(rows):
    pivots = {}
    for row in rows:
        word = sum(int(bit) << i for i, bit in enumerate(row))
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                break
    return len(pivots)


def compute():
    paths = [ROOT / 'elliptic-curves/cas' / n for n in SOURCE_NAMES]
    paths += [ART / n for n in ARTIFACT_NAMES]
    paths += [ROOT / 'elliptic-curves/notes/CURVE302_RELATIVE_DESCENT_AND_MATCHED_FIBRES_2026-09-10.md']
    bindings = {}
    # Preserve and validate transitive explicit file bindings, including the
    # earlier independent seed replay's complete ideal/norm-form inputs.
    while paths:
        p = paths.pop()
        key = str(p.relative_to(ROOT))
        if key in bindings:
            continue
        bindings[key] = sha(p)
        if p.suffix != '.json':
            continue
        d = read(p)
        if not isinstance(d, dict):
            continue
        for field in ['bindings', 'inputs']:
            for name, digest in d.get(field, {}).items():
                if not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest):
                    continue
                target = ROOT / name
                assert sha(target) == digest, (key, name)
                paths.append(target)
    b = read(ART / 'curve302_relative_ideal_anatomy_v1.json')
    f = read(ART / 'curve302_even_ideal_extension_complete_v1.json')
    s = read(ART / 'curve302_seed_norm_artin_obstruction_v1.json')
    panel = read(ART / 'rank_jump_matched_descent_panel_verification_v1.json')
    A = b['artin_matrix']
    assert [rank([row[:8] for row in A]), rank([row[8:] for row in A]), rank(A)] == [8, 10, 18]
    assert rank(f['artin_matrix']) == 19 and len(f['kernel_words']) == 3
    relative = read(ART / 'det1092_seed_artin_v2/relative-summary.json')
    seed_vectors = []
    for c, correction in zip(s['cases'], relative['compatible_words']):
        assert c['integral_norm_equation'] == 'NO_INTEGER_SOLUTION'
        assert all(e['status'] == 'PASS' for e in c['entries'])
        bits = [e['bit'] for e in c['entries']]
        assert any(bits)
        seed_vectors.append([(bit + sum(r[i] for i, v in enumerate(correction) if v)) % 2
                             for bit, r in zip(bits, A)])
    assert seed_vectors[0] == seed_vectors[1]
    assert rank([r[:8] + [bit] for r, bit in zip(A, seed_vectors[0])]) == 9
    assert panel['completed_cases'] == 7 and len(panel['cases']) == 10
    done = [r for r in panel['cases'] if r['status'].startswith('PASS')]
    assert all(r['full_Selmer_boundary_closed'] for r in done)
    assert all(r['ordinary_unramified_relative_dimension'] == r['rank_lower_bound'] - 17 for r in done)
    return {
        'schema': 'curve302.relative-descent.manifest.v1',
        'bindings': dict(sorted(bindings.items())),
        'software': {'Sage': '10.9', 'PARI': '2.17.3'},
        'results': {
            'generic_even_dimension': 8, 'generic_half_ideal_dimension': 8,
            'strict_dimension': 10, 'strict_unit_kernel': 0,
            'strict_relative_half_ideal_dimension': 10,
            'generic_plus_strict_ordinary_split_dimension': 18,
            'known_even_dimension': 22, 'known_ordinary_unramified_dimension': 20,
            'generic_ordinary_unramified_dimension': 6,
            'ordinary_class_two_rank_lower_bound': 20,
            'ordinary_split_elementary_dimension_lower_bound': 19,
            'known_even_unit_kernel_interval': [0, 2],
            'known_even_half_ideal_dimension_interval': [20, 22],
            'relative_half_ideal_dimension_interval': [12, 14],
            'first_seed_ideal_outside_full_inherited_image': True,
            'generic_plus_first_seed_split_dimension': 9,
            'historical_norm_equations_integer_solubility': [False, False],
            'matched_roster': 10, 'matched_local_anatomy_completed': 7,
            'matched_incomplete_support': 3,
            'matched_full_half_ideal_anatomy': 'UNKNOWN',
            'full_class_groups_and_selmer_upper_bounds': 'UNKNOWN',
        },
        'replay_from_repository_root': [
            'python3 research/elliptic-curves/cas/package_curve302_relative_descent.py --check',
            'timeout 120 sage -python research/elliptic-curves/cas/verify_curve302_relative_ideal_anatomy.sage --negative-controls',
            'timeout 60 sage -python research/elliptic-curves/cas/verify_det1092_seed_artin.sage',
            'timeout 60 sage -python research/elliptic-curves/cas/curve302_seed_norm_artin_obstruction.sage --check',
            'timeout 180 sage -python research/elliptic-curves/cas/verify_rank_jump_matched_descent_panel.sage --check',
        ],
        'limits': {'workers': 1, 'seconds_per_control': 60,
                   'independent_replay_seconds': 180, 'initial_panel_replay_seconds': 120, 'new_factorizations': 0,
                   'point_searches': 0, 'new_parameters': 0, 'bnf_calls': 0},
        'boundary': 'Retrospective exact arithmetic. No full class-group computation, exact elliptic rank, prospective predictor, shared carrier, or propagation theorem.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.check:
        assert read(OUTPUT) == result
    else:
        with OUTPUT.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print('PASS relative descent manifest;', len(result['bindings']), 'file bindings')
