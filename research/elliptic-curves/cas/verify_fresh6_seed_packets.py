#!/usr/bin/env python3
"""Replay self-contained native M17/M18 packets without a point search or Sage."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import compact_atlas_specialization as spec
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank


def verify(path):
    data=json.loads(path.read_text())
    families={r['family']:r for r in json.loads(spec.ATLAS.read_text())['families']}
    if len(data['records'])!=6:raise ArithmeticError('fixed six-case cohort required')
    models=[];successes=0
    for row in data['records']:
        model,initial=spec.specialize(families[row['family']],row['parameter'])
        packet=row['packet'];points=tuple(tuple(map(F,p)) for p in packet['points'])
        if model!=tuple(map(F,packet['curve'])) or points[:17]!=initial:
            raise ArithmeticError('native generic basis differs')
        old=packet['proof'];rank=packet['rank_lower_bound']
        if rank not in (17,18) or len(points)!=rank:raise ArithmeticError('first-seed packet rank differs')
        proof=checked_rank(model,points,[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.loads(json.dumps(proof))!=old:raise ArithmeticError('finite independence proof differs')
        if any(cert.isomorphic(model,other) for other in models):raise ArithmeticError('duplicate curve in six-case cohort')
        models.append(model);successes+=rank==18
        print('EXACT_NATIVE_SEED',row['id'],rank,flush=True)
    if successes!=data['certified_M18_count']:raise ArithmeticError('M18 count differs')
    print('PASS_FRESH6_SEED_PACKETS',successes,'of6',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result',required=True,type=Path)
    verify(parser.parse_args().result.resolve())
