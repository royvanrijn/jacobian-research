#!/usr/bin/env python3
"""Apply the independent private-atom proof to the corrected producer output."""
import argparse
from math import gcd
from pathlib import Path
import verify_generic_chord_ramification as old
import retrospective as r

INPUT=r.OUT/'rank_jump_generic_chord_ramification_v2.json'
OUTPUT=r.OUT/'rank_jump_generic_chord_ramification_verification_v2.json'


def verify():
    old.INPUT=INPUT
    result=old.verify()
    for row in r.read(INPUT)['rows']:
        atoms=list(map(int,row['atoms']))
        assert all(gcd(a,b)==1 for i,a in enumerate(atoms) for b in atoms[i+1:])
    prior=r.read(old.OUTPUT)
    assert result['rows']==prior['rows']
    result['schema']='rank-jump.generic-chord-ramification-verification.v2'
    result['pairwise_coprime_atoms_independently_verified']=True
    result['rows_identical_to_v1_verification']=True
    for p in [Path(__file__),old.OUTPUT]:result['bindings'][str(p.relative_to(r.ROOT))]=r.digest(p.read_bytes())
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','check']);args=p.parse_args()
    result=verify()
    if args.mode=='capture':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],[x['forms_verified'] for x in result['rows']])
