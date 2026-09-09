#!/usr/bin/env python3
"""Replay diagnostic dependence identities with Fraction arithmetic, without Sage."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from half_lattice_pointed_sieve import linear_combination_python
from alternate_quartic_covers import point_on_short_curve
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]


def run(folder):
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    result,frame,protocol=[read(folder/n) for n in ('result.json','frame.json','protocol.json')]
    assert result['protocol_sha256']==sha(folder/'protocol.json')
    assert result['frame_sha256']==sha(folder/'frame.json')
    assert all(sha(ROOT/n)==h for n,h in protocol['inputs'].items())
    model=tuple(map(F,frame['curve']))
    basis=[tuple(map(F,p)) for p in frame['basis']]
    assert all(point_on_short_curve(model,p) for p in basis)
    assert [r['point'] for r in result['rows']]==protocol['points']
    for row in result['rows']:
        p=tuple(map(F,row['point']));r=row['result']
        assert point_on_short_curve(model,p)
        assert r['status']=='INHERITED_RATIONAL_SPAN'
        multiplier=int(r['relation_multiplier'])
        assert multiplier != 0 and multiplier%2==1
        assert linear_combination_python(model,[p],[multiplier])==linear_combination_python(
            model,basis,list(map(int,r['relation_word'])))
    paths=[folder/n for n in ('result.json','frame.json','protocol.json')]
    paths += [Path(__file__),Path(__file__).with_name('half_lattice_pointed_sieve.py'),
              Path(__file__).with_name('alternate_quartic_covers.py')]
    output={'status':'PASS_INDEPENDENT_RATIONAL_DEPENDENCE_IDENTITIES',
        'relations':len(result['rows']), 'bindings':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'claim_boundary':'Exact rational-span membership of these sampled points only. '
            'The halving algorithm itself is not replayed; the final identities suffice.'}
    target=folder/'relations-verified.json'
    if target.exists():assert read(target)==output
    else:checkpoint(target,output)
    print(output['status'],output['relations'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True)
    run(p.parse_args().folder.resolve())
