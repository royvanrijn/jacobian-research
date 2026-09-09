#!/usr/bin/env python3
"""Portable exact rank replay. No raw search files, point enumeration or Sage needed."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import json
from r17_60_arithmetic import native_check
import certify_compact_r17_candidates as independent


def verify(path):
    result=json.loads(Path(path).read_text())
    if result['status']!='PASS_CERTIFIED_SEARCH':
        raise ArithmeticError('not a certified rank result')
    packet=result['packet']
    native_check(packet,result)
    if len(packet['points'])!=packet['rank_lower_bound'] or result['rank_lower_bound']!=packet['rank_lower_bound']:
        raise ArithmeticError('reported lower bound differs')
    proof=packet['proof']
    second=independent.checked_rank(tuple(map(F,packet['curve'])),
            tuple(tuple(map(F,p)) for p in packet['points']),[s['prime'] for s in proof['signatures']],
            proof['no_rational_2_torsion_prime'])
    if second['rank_lower_bound']!=result['rank_lower_bound']:
        raise ArithmeticError('second independent finite certificate failed')
    print('FOUNDRY_PORTABLE_RANK_PASS',result['id'],result['rank_lower_bound'])
    return result['rank_lower_bound']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('certificates',type=Path,nargs='+')
    for path in p.parse_args().certificates:verify(path)
