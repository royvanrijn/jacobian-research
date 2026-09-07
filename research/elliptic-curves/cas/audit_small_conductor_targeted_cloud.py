#!/usr/bin/env python3
"""Independently audit the twelve-attempt retained cloud modulo3 and5."""
import argparse
import json

import study_small_conductor_target as study
import audit_retained_cloud_modl as odd

INPUT = study.D / 'point_only_cloud.json'
OUTPUT = study.ART / 'small_conductor_targeted_cloud_modl_v2.json'


def expected():
    study.protocol()
    result = study.cert.read(study.D / 'result.json')
    proof = study.cert.read(study.PROOF)
    if study.cert.read(study.D / 'ledger.json')['status'] != 'PASS':
        raise ArithmeticError('terminal exact admission replay required')
    seed = [tuple(map(study.cert.F, p)) for p in proof['short_points']]
    seen = {(x, abs(y)) for x, y in seed}
    points = list(seed)
    for row in result['charts']:
        for raw in row['search']['finite_curve_points']:
            x, y = study.cert.F(raw['x']), study.cert.F(raw['y'])
            key = x, abs(y)
            if key not in seen:
                seen.add(key)
                points.append((x, y))
    return {'schema': 'elliptic-curves.small-conductor-point-only-cloud.v2',
            'status': 'COMPLETE_DECLARED_FINITE_AUDIT',
            'family': proof['family'], 'parameter': '3/17',
            'curve': proof['short_model'], 'points': [list(map(str, p)) for p in points],
            'rank_lower_bound': result['rank_lower_bound'],
            'bindings': {str(p.relative_to(study.ROOT)): study.cert.hashed(p)
                         for p in [study.PROOF, study.D / 'result.json', study.D / 'ledger.json']},
            'scope': 'Point-only union of the certified22 seed and every exact retained witness from the twelve attempts. Mod2 admission and orbit-compression histories have completed replay. This is not a chronological search transcript or a claim that timed-out boxes completed.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if data != study.cert.read(INPUT):
            raise ArithmeticError('retained union differs')
        audit = study.cert.read(OUTPUT)
        if audit['input_sha256'] != study.cert.hashed(INPUT) or audit['points'] != data['points'] or audit['curve'] != data['curve']:
            raise ArithmeticError('odd audit input differs')
        odd.check(OUTPUT)
    else:
        if INPUT.exists() or OUTPUT.exists():
            raise FileExistsError('preserve complete-cloud audit')
        study.checkpoint(INPUT, data)
        odd.build(INPUT, OUTPUT)
