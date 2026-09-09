#!/usr/bin/env python3
"""Cheap root-number diagnostic on exact-conductor inventory rows, without factoring.

Parity of an existing lower bound is not a rank theorem. This table is only
for choosing bounded follow-up work; no new point or rank claim is made.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from local_conductor_database import load_conductor_inventory, CURRENT, ROOT


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    from sage.all import QQ, EllipticCurve, pari
    if output.exists():
        raise FileExistsError('preserve root-number diagnostic')
    start = time.monotonic()
    rows = []
    bindings = {str(CURRENT.relative_to(ROOT)): sha(CURRENT),
                str(Path(__file__).resolve().relative_to(ROOT)): sha(Path(__file__))}
    for row in load_conductor_inventory():
        info = row['conductor_information']
        if info['status'] != 'EXACT':
            continue
        curve = pari(EllipticCurve(QQ, row['curve']))
        local, sign = [], -1
        for p in map(int, info['bad_primes']):
            w = int(curve.ellrootno(p))
            if w not in (-1, 1):
                raise ArithmeticError('invalid local root number')
            local.append({'prime': str(p), 'root_number': w})
            sign *= w
        certificate = ROOT/info['certificate']
        bindings[info['certificate']] = sha(certificate)
        rank = row['rank_lower_bound']
        rows.append({'id': row['id'], 'family': row['family'], 'parameter': row['parameter'],
            'rank_lower_bound': rank, 'curve': row['curve'], 'root_number': sign,
            'differs_from_lower_bound_parity': sign != (-1)**rank,
            'exact_conductor': info['conductor'], 'conductor_certificate': info['certificate'],
            'local': local, 'infinity': -1})
    if any(sha(ROOT/p) != h for p,h in bindings.items()):
        raise ArithmeticError('root diagnostic input changed')
    result = {'status': 'COMPLETE_EXACT_CONDUCTOR_ROOT_DIAGNOSTIC', 'rows': rows,
        'bindings': bindings, 'wall_seconds': time.monotonic()-start, 'point_searches': 0,
        'claim_boundary': 'PARI local root numbers using certified complete bad-prime lists. '
            'A sign different from the parity of a lower bound is a scheduling diagnostic, '
            'not an algebraic rank lower bound, rank upper bound or proof of an extra point.'}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('ROOT_PRIOR', len(rows), 'rows', result['wall_seconds'], 'seconds')
    print([(r['id'], r['rank_lower_bound'], r['root_number']) for r in rows if r['differs_from_lower_bound_parity']])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output)
