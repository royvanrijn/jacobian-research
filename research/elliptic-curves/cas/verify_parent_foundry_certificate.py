#!/usr/bin/env python3
"""Portable rank replay with two finite implementations; no Sage or point search.

Optional --parent also checks the labelled generic specialization prefix.
This verifies rank lower bounds, not a catalogue or world novelty assertion.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from memory_rank_certificate import checked_rank
import certify_compact_r17_candidates as independent


def verify(path,parent_path=None):
    data=json.loads(Path(path).read_text());packet=data.get('packet',data)
    if not isinstance(packet,dict):raise ValueError('self-contained packet required')
    model=tuple(map(F,packet['curve']));points=tuple(tuple(map(F,p)) for p in packet['points'])
    rank=packet['rank_lower_bound'];proof=packet['proof']
    if len(points)!=rank or data.get('rank_lower_bound',rank)!=rank:
        raise ArithmeticError('rank header differs from point count')
    primes=[s['prime'] for s in proof['signatures']];torsion=proof['no_rational_2_torsion_prime']
    first=checked_rank(model,points,primes,torsion)
    second=independent.checked_rank(model,points,primes,torsion)
    if first['rank_lower_bound']!=rank or second['rank_lower_bound']!=rank:
        raise ArithmeticError('independent finite rank replay failed')
    if parent_path:
        from parent_foundry_worker import specialize
        parent=json.loads(Path(parent_path).read_text());parent=parent.get('parent',parent)
        expected,generic=specialize(parent,data.get('parameter',packet.get('parameter')))
        if expected!=model or generic!=points[:len(generic)]:
            raise ArithmeticError('generic specialization binding differs')
    print('PARENT_FOUNDRY_RANK_PASS',rank,str(path),flush=True)
    return rank


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('certificate',type=Path)
    p.add_argument('--parent',type=Path);args=p.parse_args();verify(args.certificate,args.parent)
