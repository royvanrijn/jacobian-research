#!/usr/bin/env python3
"""Standalone algebra checker; Python standard library only.

This checks finite nonabelian images of the supplied presentations.
Named-knot/cover identification and bridge upper bounds have separate replay.
"""
import json
from pathlib import Path
import sys


def check(case):
    n=case['degree'];identity=list(range(n));images=[case['A'],case['B']]
    if any(sorted(p)!=identity for p in images):raise ValueError('not permutations')
    inverses=[[p.index(i) for i in identity] for p in images]
    def word_value(word):
        result=identity[:]
        for g in word:
            if g not in [1,2,-1,-2]:raise ValueError('invalid generator')
            p=images[g-1] if g>0 else inverses[-g-1]
            result=[result[p[i]] for i in identity]
        return result
    for index,word in enumerate(case['relators']):
        if word_value(word)!=identity:raise ValueError(f'relator {index} fails')
    commutator=word_value([1,2,-1,-2])
    if commutator==identity:raise ValueError('generators commute')
    witness=case['commutator_witness'];point=witness['point']
    if not 0<=point<n or commutator[point]!=witness['image'] or commutator[point]==point:
        raise ValueError('commutator witness fails')
    return {'name':case['name'],'degree':n,'relators':len(case['relators']),'moved_point':point,'image':commutator[point]}


def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('finite_cover_certificates.json')
    data=json.loads(path.read_text());assert data['schema']==1
    print(json.dumps({'status':'PASS finite nonabelian quotient certificates','cases':[check(c) for c in data['cases']]},indent=2))


if __name__=='__main__':main()
