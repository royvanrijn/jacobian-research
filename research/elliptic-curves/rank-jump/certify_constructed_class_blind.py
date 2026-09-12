#!/usr/bin/env python3
"""Fresh finite independence proof for the blind lifts over generic sixteen."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
import verify_constructed_class_blind as exact

sys.path.insert(0, str(exact.ROOT/'elliptic-curves/cas'))


def certify(output):
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import find_two_torsion_certificate_prime
    start = time.monotonic()
    paths = [exact.EXPERIMENT/'generic.json', exact.EXPERIMENT/'independent-cover-replay.json',
             Path(__file__).resolve()]
    generic, covers = [json.loads(p.read_text()) for p in paths[:2]]
    assert covers['status'] == 'PASS'
    for name, h in covers['bindings'].items():
        assert exact.sha(exact.ROOT/name) == h
    model = tuple(map(F, generic['curve']))
    points = tuple(tuple(map(F, p)) for p in generic['points'])
    points += tuple(tuple(map(F, c['point_on_short_model'])) for c in covers['cases'])
    assert len(points) == 18
    admission = FinitePointAdmission(model, points, prime_bound=1000)
    proof = checked_rank(model, points, admission.primes, find_two_torsion_certificate_prime(model))
    result = {'status': 'PASS_GENERIC16_PLUS_TWO_BLIND_LIFTS',
        'curve': list(map(str, model)), 'points': [list(map(str, p)) for p in points],
        'generic_rank': 16, 'combined_rank_lower_bound': 18, 'proof': proof,
        'wall_seconds': time.monotonic()-start,
        'bindings': {str(p.relative_to(exact.ROOT)): exact.sha(p) for p in paths},
        'boundary': 'Fresh exact finite quotient rank18 and no rational2-torsion prove independence over the supplied generic sixteen. This re-establishes two known directions, not a new curve rank record or exact rank.'}
    with output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write('\n')
    print('PASS independent generic16 plus two blind lifts: rank18')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    certify(p.parse_args().output.resolve())
