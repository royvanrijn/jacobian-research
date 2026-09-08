#!/usr/bin/env python3
"""Certify the translation-class corollary, using exact finite reduction columns."""
import argparse
import json
from pathlib import Path

import certify_compact_r17_candidates as cert
import study_mw16_rank27_visibility as study
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

ROOT, ART = study.ROOT, study.ART
OUT = ART / 'mw16_rank27_translation_classes_v1.json'


def expected():
    data = cert.read(study.INPUT)
    proof = data['point_proof']
    points = [tuple(map(cert.F, p)) for p in proof['discovery_points']]
    model = tuple(map(cert.F, proof['discovery_curve']))
    old = proof['rank_certificate']
    checked = checked_rank(model, points, [s['prime'] for s in old['signatures']],
                           old['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(checked)) != old:
        raise ArithmeticError('exact mod2 injectivity proof differs')
    # Recompute every centre and both-sign chart geometry as part of this proof.
    visibility = study.expected(data)
    if visibility != cert.read(study.OUTPUT):
        raise ArithmeticError('pointwise geometry certificate differs')
    initial = [r['centre']['mask'] for r in data['arms']['initial']]
    adaptive = [r['centre']['parity'] for r in data['arms']['adaptive']]
    if len(set(initial + adaptive)) != 344 or any(m >> 16 for m in initial):
        raise ArithmeticError('distinct initial/adaptive parity roster required')
    if any(not m >> 16 or m >> 26 for m in adaptive):
        raise ArithmeticError('nonzero ten-bit adaptive quotient required')
    winner = data['arms']['adaptive'][85]['centre']
    return {'schema': 'elliptic-curves.mw16-rank27-translation-classes.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(ROOT)): cert.hashed(p) for p in
                         (study.INPUT, study.OUTPUT, Path(study.__file__), Path(__file__).resolve())},
            'model': proof['discovery_curve'], 'certified_mod2_column_rank': 27,
            'initial_subgroup_rank': 26, 'generic_subgroup_rank': 16,
            'available_classes_in_initial_subgroup_mod2': 2**26,
            'lifts_per_generic_parity': 2**10,
            'certified_distinct_translation_classes': 344,
            'generic_classes_tested': 43, 'adaptive_classes_tested': 301,
            'winning_generic_mask': winner['generic_mask'],
            'winning_quotient_word': winner['quotient_word'],
            'winning_parity': winner['parity'],
            'winning_quotient_support_one_based': [j+1 for j in range(10)
                                                 if (winner['quotient_word'] >> j) & 1],
            'winning_generic_mask_was_in_original_roster': winner['generic_mask'] in initial,
            'lemma': 'For i_C(R)=C-R and tau_T(R)=R+T, tau_T i_C tau_-T = i_(C+2T). '
                     'Hence i_C and i_D are conjugate by a Q-rational translation iff D-C is '
                     'in 2E(Q). Independent finite quotient columns inject H/2H into '
                     'E(Q)/2E(Q), so distinct coefficient parities in this certified H '
                     'give distinct rational-translation classes. For C in H and P outside '
                     'H tensor Q, the two points P and C-P give opposite nonzero classes '
                     'modulo H tensor Q and have the same pointed-map coordinate.',
            'corollary': 'Every adaptive centre here is rational-translation inequivalent '
                         'to every original generic centre. Each fixed generic parity '
                         'has 1024 distinct lifts in the rank26 subgroup. Only 344 of '
                         '67108864 classes in that subgroup were tested by these two rosters.',
            'claim_boundary': 'Translation equivalence of pointed involutions, not abstract '
                              'isomorphism of genus-one curves or insolubility of a 2-cover. '
                              'All pointed curves already have known endpoints. The successful '
                              'generic mask was not originally searched, and horizontal maps '
                              'also differ: this is not a controlled causal ablation of quotient '
                              'bits alone. No new rank, upper bound, saturation at all primes, '
                              'or global quotient-coset visibility follows.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = expected()
    if args.check:
        if cert.read(OUT) != result:
            raise ArithmeticError('translation class certificate differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve translation class proof')
        checkpoint(OUT, result)
    print('EXACTLY 344 RATIONAL-TRANSLATION CLASSES; 1024 LIFTS PER GENERIC PARITY')
