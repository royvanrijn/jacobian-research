#!/usr/bin/env python3
"""Whole residue-ring replay of all ten exclusions of a second scale."""
import argparse
from pathlib import Path
import audit_mw16_postscale_reduction as audit
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

OUT = audit.first.support.ART/'mw16_postscale_reduction_replay_v1.json'


def evaluate(a, r, m):
    total = 0
    for c in a[::-1]:
        total = (r*total+c) % m
    return total


def expected():
    p = cert.read(audit.D/'protocol.json')
    d = cert.read(audit.OUT)
    if (p['sources'] != audit.sources() or d['sources'] != audit.sources()
            or d['protocol_sha256'] != cert.hashed(audit.D/'protocol.json')
            or cert.read(audit.D/'ledger.json')['status'] != 'PASS'
            or d['status'] != 'PASS_SINGLE_SCALE_CLASSIFICATION'
            or len(p['balls']) != 10 or len(d['rows']) != 10):
        raise ArithmeticError('complete ten-ball first replay required')
    count = 0
    for ball, row in zip(p['balls'], d['rows']):
        if any(row[k] != ball[k] for k in ('family', 'prime', 'chart')):
            raise ArithmeticError('ball identity differs')
        q = ball['prime']
        a = list(map(int, ball['A_divided_coefficients']))
        b = list(map(int, ball['B_divided_coefficients']))
        tree = row['second_scale_tree']
        if (not row['no_second_scale'] or tree['scale_balls']
                or tree['status'] != 'COMPLETE_RESIDUE_CLASSIFICATION'
                or not 1 <= len(tree['levels']) <= 3):
            raise ArithmeticError('complete second-scale exclusions through depth3 required')
        previous = None
        for depth, level in enumerate(tree['levels'], 1):
            modulus = q**depth
            if modulus > 2197:
                raise ArithmeticError('fixed whole-ring bound exceeded')
            # Entire residue ring, independent of the producer's branch tree.
            live = [r for r in range(modulus)
                    if evaluate(a, r, q**min(depth, 4)) == 0
                    and evaluate(b, r, q**min(depth, 6)) == 0]
            count += modulus
            candidates = list(range(q)) if previous is None else sorted(r+q**(depth-1)*j for r in previous for j in range(q))
            actual = {'depth': depth, 'modulus': modulus,
                      'excluded_residues': [r for r in candidates if r not in set(live)],
                      'admitted_residues': [], 'unresolved_residues': live}
            if level != actual:
                raise ArithmeticError('whole-ring necessary congruences differ')
            previous = live
        if previous:
            raise ArithmeticError('second-scale possibility remains')
        cells = []
        for z in range(q):
            aa, bb = evaluate(a, z, q), evaluate(b, z, q)
            cells.append({'next_digit': z, 'reduced_A': aa, 'reduced_B': bb,
                          'good_reduction_after_one_scale': (4*aa**3+27*bb**2) % q != 0,
                          'parameter_residue': ball['residue']+ball['modulus']*z,
                          'parameter_modulus': ball['modulus']*q})
        if row['cells'] != cells or row['first_residue'] != ball['residue'] or row['first_modulus'] != ball['modulus']:
            raise ArithmeticError('complete local reduction cells differ')
    paths = [Path(__file__).resolve(), Path(audit.__file__), audit.OUT,
             audit.D/'protocol.json', audit.D/'ledger.json']
    return {'schema': 'elliptic-curves.mw16-postscale-reduction-replay.v1',
            'status': 'PASS', 'sources': {str(p.relative_to(audit.ROOT)): cert.hashed(p) for p in paths},
            'first_scale_balls': 10, 'whole_ring_residues_tested': count,
            'no_second_scale_in_any_ball': True,
            'claim_boundary': 'Independent whole residue-ring necessary congruences exclude a second '
                'p-scale everywhere in all ten admitted first-scale balls. Every good/bad next-digit '
                'cell replays. This closes local scaling at the declared5/13 pairs; it is not a '
                'global minimal-model, point, score-optimality, conductor or rank certificate.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if cert.read(OUT) != data:
            raise ArithmeticError('postscale whole-ring replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve independent postscale replay')
        checkpoint(OUT, data)
    print('ALL TEN MW16 SECOND-SCALE EXCLUSIONS PASS', data['whole_ring_residues_tested'], flush=True)
