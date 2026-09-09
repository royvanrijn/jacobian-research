#!/usr/bin/env python3
"""Prepare a new sealed seed after cloud reconciliation, retaining exact parents."""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
from memory_rank_certificate import checked_rank

ROOT=Path(__file__).resolve().parents[2]
def run(preparation,packet,output):
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    old=read(preparation/'prepared.json')
    if not all(sha(preparation/n)==h for n,h in old['files'].items()):raise ArithmeticError('old preparation changed')
    prior=read(next(preparation.glob('seed-M*.json')));seed=read(packet)
    if seed['status']!='PASS_REPLAYED_CLOUD_RECONCILIATION' or seed['curve']!=prior['curve'] or seed['points'][:len(prior['points'])]!=prior['points']:raise ArithmeticError('native seed prefix differs')
    proof=seed['proof'];points=tuple(tuple(map(F,p)) for p in seed['points'])
    fresh=checked_rank(tuple(map(F,seed['curve'])),points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(fresh))!=proof or len(points)!=seed['rank_lower_bound']:raise ArithmeticError('reconciled certificate differs')
    protocol=read(preparation/'protocol.json');protocol['initial_rank']=len(points)
    protocol['continuation_inputs']={str(p.relative_to(ROOT)):sha(p) for p in [packet,preparation/'prepared.json',Path(__file__)]}
    protocol['continuation_rule']='Rebuild from reconciled exact subgroup, retaining the same frozen generic parents; no point search in preparation.'
    output.mkdir(exist_ok=False)
    def save(name,data):(output/name).write_text(json.dumps(data,indent=2)+'\n')
    for n in old['files']:
        if n in ('protocol.json','anchor-bank.json') or n.startswith('seed-M'):continue
        (output/n).write_bytes((preparation/n).read_bytes())
    save('protocol.json',protocol)
    bank=read(preparation/'anchor-bank.json');bank['protocol_sha256']=sha(output/'protocol.json');save('anchor-bank.json',bank)
    save(f'seed-M{len(points)}.json',{k:seed[k] for k in ['curve','points','proof','rank_lower_bound','generic_rank']})
    save('prepared.json',{'status':'PASS_RECONCILED_V3_PREPARATION','point_searches':0,'files':{p.name:sha(p) for p in sorted(output.glob('*.json'))}})
    print('PREPARED_RECONCILED',len(points),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['preparation','packet','output']:p.add_argument('--'+n,required=True,type=Path)
    a=p.parse_args();run(a.preparation.resolve(),a.packet.resolve(),a.output.resolve())
