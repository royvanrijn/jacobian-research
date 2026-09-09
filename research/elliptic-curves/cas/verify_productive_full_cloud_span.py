#!/usr/bin/env python3
"""Exact rank of the retained point subgroup, not an upper bound for E(Q)."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from verify_productive_cloud_relations import run as verify_relations
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]


def run(folder):
    verify_relations(folder)
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    p=read(folder/'protocol.json');r=read(folder/'result.json');f=read(folder/'frame.json')
    terminals=[ROOT/n for n in p['inputs'] if n.endswith('/terminal.json')]
    assert len(terminals)==1
    t=read(terminals[0]);assert t['curve']==f['curve'] and t['points']==f['basis']
    model=tuple(map(F,t['curve']));basis=[tuple(map(F,x)) for x in t['points']]
    proof=t['proof'];fresh=checked_rank(model,basis,[x['prime'] for x in proof['signatures']],proof['no_rational_2_torsion_prime'])
    assert json.loads(json.dumps(fresh))==proof
    # Reconstruct all distinct signed points directly from the retained receipts.
    clouds=[ROOT/n for n in p['inputs'] if n.endswith('/cloud.json')]
    assert set(clouds)==set(terminals[0].parent.glob('epoch-*/cloud.json'))
    points={(F(x['x']),abs(F(x['y']))) for cloud in clouds for chart in read(cloud)['charts']
            for x in chart['search']['finite_curve_points']}
    native={(x,abs(y)) for x,y in basis};extra=points-native
    selected=[tuple(map(F,x)) for x in p['points']]
    assert len(selected)==len(set(selected)) and set(selected)==extra
    assert len(r['rows'])==len(extra)
    result={'status':'PASS_EXACT_RETAINED_POINT_SUBGROUP_RANK','rank':len(basis),
        'nonbasis_points':len(extra),'sign_normalized_cloud_with_basis':len(points|native),
        'curve':t['curve'],'basis':t['points'],'rank_certificate':proof,
        'relations':[{'point':x['point'],'multiplier':x['result']['relation_multiplier'],
                      'word':x['result']['relation_word']} for x in r['rows']],
        'bindings':{str(q.relative_to(ROOT)):sha(q) for q in
            [folder/'relations-verified.json',folder/'result.json',folder/'protocol.json',Path(__file__),*clouds]},
        'claim_boundary':'Every retained point lies in this certified rational span. Exact '
            'rank of the retained subgroup only; no upper bound for the elliptic curve.'}
    q=folder/'full-span-verified.json'
    if q.exists():assert read(q)==result
    else:checkpoint(q,result)
    print('FULL_RETAINED_SPAN_VERIFIED',len(basis),len(extra))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,required=True)
    run(p.parse_args().folder.resolve())
