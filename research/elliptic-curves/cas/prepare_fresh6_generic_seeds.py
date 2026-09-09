#!/usr/bin/env python3
"""Exact native generic17 packets for the frozen fresh-six cohort; no searches."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import compact_atlas_specialization as spec
import certify_compact_r17_candidates as cert
from future_point_admission import FinitePointAdmission
from memory_rank_certificate import checked_rank
from v3_warm_engine import certified_state
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
CAS=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def run(selection,folder,replay=False):
    selected=read(selection/'result.json');verified=read(selection/'verified.json')
    assert verified['status']=='PASS_EXACT_SELECTION_REPLAY' and verified['result_sha256']==sha(selection/'result.json')
    p=read(selection/'protocol.json')
    assert all(sha(ROOT/n)==h for category in ('inputs','sources') for n,h in p[category].items())
    if replay:protocol=read(folder/'protocol.json')
    else:
        folder.mkdir(exist_ok=False)
        paths=[selection/n for n in ('protocol.json','result.json','verified.json')]+[spec.ATLAS]
        protocol=dict(schema='fresh6-generic17-packets.v1',point_searches=0,
            inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},
            sources={str((CAS/n).relative_to(ROOT)):sha(CAS/n) for n in
                ('prepare_fresh6_generic_seeds.py','compact_atlas_specialization.py',
                 'future_point_admission.py','memory_rank_certificate.py','v3_warm_engine.py')},
            rule='Specialize exactly the17 atlas sections on each selected equation, then '
                 'certify finite quotient independence and absence of rational2-torsion. '
                 'No extra point or point search is included.')
        checkpoint(folder/'protocol.json',protocol)
    assert all(sha(ROOT/n)==h for category in ('inputs','sources') for n,h in protocol[category].items())
    families={r['family']:r for r in read(spec.ATLAS)['families']}
    records=[]
    for row in selected['rows']:
        model,points=spec.specialize(families[row['family']],row['parameter'])
        assert model==tuple(map(F,row['model']))
        path=folder/row['id']/'seed-M17.json'
        if replay:
            packet=read(path)
            assert packet['curve']==list(map(str,model)) and packet['points']==[list(map(str,p)) for p in points]
            old=packet['proof']
            proof=checked_rank(model,points,[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
            assert json.loads(json.dumps(proof))==old
        else:
            admission=FinitePointAdmission(model,points,prime_bound=1000)
            torsion=cert.find_two_torsion_certificate_prime(model,prime_bound=200)
            proof=checked_rank(model,points,admission.primes,torsion)
            packet=dict(id=row['id'],family=row['family'],parameter=row['parameter'],
                curve=list(map(str,model)),points=[list(map(str,p)) for p in points],
                generic_rank=17,rank_lower_bound=17,proof=proof,
                claim_boundary='The17 generic sections only; no extra direction or record claimed.')
            path.parent.mkdir();checkpoint(path,packet)
        certified_state(model,points,proof)
        records.append(dict(id=row['id'],family=row['family'],parameter=row['parameter'],
            rank_lower_bound=17,path=str(path.relative_to(ROOT)),sha256=sha(path)))
        print('GENERIC17_VERIFIED',row['id'],flush=True)
    assert all(sha(ROOT/n)==h for category in ('inputs','sources') for n,h in protocol[category].items())
    result=dict(status='PASS_SIX_GENERIC17_PACKETS',point_searches=0,records=records,
                protocol_sha256=sha(folder/'protocol.json'))
    if replay:
        assert result==read(folder/'prepared.json')
        checkpoint(folder/'verified.json',dict(status='PASS_INDEPENDENT_GENERIC17_REPLAY',
            prepared_sha256=sha(folder/'prepared.json'),source_sha256=sha(Path(__file__))))
    else:checkpoint(folder/'prepared.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True);p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--replay',action='store_true');a=p.parse_args()
    run(a.selection.resolve(),a.folder.resolve(),a.replay)
