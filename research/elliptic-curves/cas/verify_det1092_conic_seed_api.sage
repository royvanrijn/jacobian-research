#!/usr/bin/env sage-python
"""Bounded producer interface tests; the separate finite-group proof is primary.

One new address (n=0), its parameter roundtrip, and nine old noncoverage
addresses. Validator-only comparison follows formula-only production.
"""
import argparse,hashlib,json,runpy,signal
from pathlib import Path
from sage.all import QQ,EllipticCurve
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_conic_seed_progression_v1'
PRODUCER=ROOT/'elliptic-curves/cas/det1092_conic_seed_v2.sage'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def run():
    produce=runpy.run_path(str(PRODUCER))['produce']
    seed=produce(argparse.Namespace(n='0',parameter=None))
    assert seed['status']=='PROVED_RANK_AT_LEAST_18_WITH_NON_GENERIC_SEED'
    assert seed==read(DIR/'seed-at-n0-v2.json')
    # Validation data are not an input to production.
    inp=read(DIR/'input.json')
    assert seed['elliptic_a_invariants']==inp['equation_at_zero']
    assert seed['generic_points']+[seed['seed']]==inp['points_at_zero']
    E=EllipticCurve(QQ,list(map(QQ,seed['elliptic_a_invariants'])))
    assert E.discriminant()
    for P in seed['generic_points']+[seed['seed']]:assert E(list(map(QQ,P)))
    roundtrip=produce(argparse.Namespace(n=None,parameter=seed['parameter']))
    assert roundtrip['status']==seed['status'] and roundtrip['parameter']==seed['parameter']
    assert roundtrip['elliptic_a_invariants']==seed['elliptic_a_invariants']
    assert QQ(roundtrip['conic_parameter_u'])==QQ(roundtrip['B'])*QQ(roundtrip['n'])
    outcomes=[]
    roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
    for row in roster['cases']:
        result=produce(argparse.Namespace(n=None,parameter=row['parameter']))
        assert result['status']=='UNKNOWN_OUTSIDE_PROVED_PROGRESSION_IMAGE'
        assert not result['rational_preimages']
        outcomes.append({'label':row['label'],'parameter':row['parameter'],'status':result['status']})
    data={'status':'PASS_PRODUCER_POSITIVE_ROUNDTRIP_AND_NINE_NONCOVERAGE_CASES',
        'classification':'verified interface application; mathematical proof is the independent progression replay',
        'new_addresses':1,'new_parameter':seed['parameter'],
        'positive_roundtrip_passed':True,'controls':outcomes,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PRODUCER,DIR/'replay.json',DIR/'seed-at-n0-v2.json',Path(__file__)]}}
    payload=json.dumps(data,indent=2,sort_keys=True)+'\n';path=DIR/'api-replay.json'
    if path.exists():assert path.read_text()==payload
    else:path.write_text(payload)
    print(data['status'],flush=True)
if __name__=='__main__':signal.alarm(25);run()
