#!/usr/bin/env python3
"""Deterministic portable manifest for the exact strict302 descent block."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUTPUT = ART/'curve302_descent_anatomy_manifest_v1.json'
SOURCES = [
    'elliptic-curves/cas/curve302_descent_anatomy.sage',
    'elliptic-curves/cas/complete_curve302_descent_anatomy.sage',
    'elliptic-curves/cas/curve302_descent_remaining_class.sage',
    'elliptic-curves/cas/verify_curve302_descent_anatomy.sage',
    'elliptic-curves/cas/package_curve302_descent_anatomy.py',
    'elliptic-curves/cas/icarm_curve302.py',
    'elliptic-curves/cas/verify_curve302_recovered_quotient_local_filtration.sage',
    'elliptic-curves/cas/research_runtime/local_kummer.py',
    'elliptic-curves/rank-jump/verify_strict_half_ideals.py',
    'elliptic-curves/rank-jump/verify_half_ideal_artin.py',
    'elliptic-curves/rank-jump/verify_generic_only_class_anchors.py',
    'elliptic-curves/rank-jump/STRICT_SELMER_AND_ARTIN_BLOCKS.md',
    'elliptic-curves/rank-jump/DERIVATIVE_RECIPROCITY_AND_COMPLETE_BOUNDARY.md',
    'elliptic-curves/notes/CURVE302_DESCENT_ANATOMY_2026-09-10.md',
]
ARTIFACTS = [
    'curve302_descent_anatomy_v1.json',
    'curve302_descent_remaining_class_v1.json',
    'curve302_recovered_quotient_local_filtration_v1.json',
    'rank_jump_curve302_strict_constructor_arithmetic_v1.json',
    'rank_jump_generic_only_class_anchors_v1.json',
    'rank_jump_generic_only_class_anchors_verification_v1.json',
]


def compute():
    paths = [ROOT/name for name in SOURCES] + [ART/name for name in ARTIFACTS]
    d = json.loads((ART/ARTIFACTS[0]).read_text())
    last = json.loads((ART/ARTIFACTS[1]).read_text())
    assert d['artin_rank'] == 9 and last['remaining_class_detected']
    assert last['unit_kernel_dimension_interval'] == [0, 0]
    return {
        'schema': 'curve302.descent-anatomy.manifest.v1',
        'bindings': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'software': d['half_ideal_packet']['software'],
        'scope': 'Retrospective exact anatomy of the displayed strict ten-dimensional rational Kummer block.',
        'results': {'strict_dimension': 10, 'unit_kernel_dimension': 0,
                    'ordinary_half_ideal_dimension': 10, 'narrow_half_ideal_dimension': 10,
                    'ordinary_elementary_direct_factor_dimension': 10,
                    'S_elementary_direct_factor_dimension_lower_bound': 9,
                    'ordinary_class_2_rank_lower_bound': 16, 'S_class_2_rank_lower_bound': 10,
                    'full_selmer_dimension': '21+c_S', 'c_S': 'UNKNOWN; at least 10',
                    'exact_rank': 'UNKNOWN; at least 31',
                    'image_modulo_full_generic_even_ideal_image': 'UNKNOWN'},
        'replay_from_repository_root': [
            'python3 research/elliptic-curves/cas/package_curve302_descent_anatomy.py --check',
            'timeout 120 sage -python research/elliptic-curves/cas/verify_curve302_descent_anatomy.sage --negative-controls',
        ],
        'producer_from_repository_root': [
            'timeout 120 sage -python research/elliptic-curves/cas/curve302_descent_anatomy.sage',
            'timeout 60 sage -python research/elliptic-curves/cas/complete_curve302_descent_anatomy.sage',
            'timeout 60 sage -python research/elliptic-curves/cas/curve302_descent_remaining_class.sage',
        ],
        'limits': {'workers': 1, 'half_ideal_and_initial_artin_seconds': 120,
                   'completion_seconds': 60, 'remaining_class_seconds': 60,
                   'independent_replay_seconds': 120, 'point_searches': 0,
                   'new_parameters': 0, 'bnf_calls': 0},
        'boundary': 'No complete class group, full Selmer upper bound, matched-fibre predictor, shared carrier, or propagation theorem is claimed.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    result = compute()
    if args.check:
        assert json.loads(OUTPUT.read_text()) == result
    else:
        with OUTPUT.open('x') as out:
            json.dump(result, out, indent=2, sort_keys=True)
            out.write('\n')
    print('PASS curve302 descent anatomy manifest')
