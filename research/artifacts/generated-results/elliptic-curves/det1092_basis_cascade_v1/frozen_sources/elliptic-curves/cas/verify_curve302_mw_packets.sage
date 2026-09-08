#!/usr/bin/env sage-python
"""Replay every finite and rational exclusion from frozen integer packets.

No numerical height enumeration is required to verify the bounded exclusion.
This shares the chart/detector implementation; it is not independent code.
"""
from hashlib import sha256
import argparse
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import signal
import sys
import numpy as np
from sage.all import EllipticCurve,QQ,ZZ,matrix

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
ART=ROOT/'artifacts/generated-results/elliptic-curves'
SEARCH=ART/'curve302_mw_packet_rank17_h190_v1.json'
PACKETS=ART/'curve302_mw_packet_rank17_h190_inputs_v1.npz'
OUT=ART/'curve302_mw_packet_rank17_h190_replay_v1.json'
HELPER=CAS/'select_mw_quartic_packet.sage'
sys.path.insert(0,str(CAS))
from icarm_curve302 import POINTS,GENERAL_WEIERSTRASS_COEFFICIENTS


def digest(p):return sha256(p.read_bytes()).hexdigest()


def build():
    loader=importlib.machinery.SourceFileLoader('packet_replay_helpers',str(HELPER))
    spec=importlib.util.spec_from_loader(loader.name,loader);helper=importlib.util.module_from_spec(spec);loader.exec_module(helper)
    d=json.loads(SEARCH.read_text())
    for path,expected in d['input_sha256'].items():assert digest(ROOT/path)==expected
    arrays=np.load(PACKETS,allow_pickle=False);rows=arrays['rows'].astype(np.int64);keys=arrays['keys'].astype(np.int64)
    assert sha256(rows.astype('<i8').tobytes()).hexdigest()==d['enumeration']['packet_array_sha256']
    assert sha256(keys.astype('<i8').tobytes()).hexdigest()==d['enumeration']['key_array_sha256']
    assert rows.shape==(130706,13,17) and len(set(map(int,keys)))==len(keys)
    masks=(rows%2)@(1<<np.arange(17,dtype=np.int64));assert np.all(masks==keys[:,None])
    B=matrix(ZZ,d['basis_rows']);assert B.rank()==17
    assert all(abs(n)==1 for n in B.smith_form()[0].diagonal())
    E0=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)));E=E0.short_weierstrass_model();iso=E0.isomorphism_to(E)
    public=[iso(E0(list(map(QQ,P)))) for P in POINTS]
    exclusions=d['modular_exclusions'];lookup={int(k):i for i,k in enumerate(keys)}
    checked=0
    for p in d['policy']['primes']:
        chart=helper.modular_chart(E,public,B,p);assert chart is not None
        for k,prime in exclusions.items():
            if prime!=p:continue
            assert helper.modular_pairs(chart(rows[lookup[int(k)]]),p) is False
            checked+=1
            if checked%10000==0:print('MODULAR_REPLAY',checked,flush=True)
    assert checked==130685
    basispoints=[sum((n*P for n,P in zip(v,public)),E(0)) for v in B.rows()]
    def point(v):return sum((ZZ(n)*P for n,P in zip(v,basispoints)),E(0))
    exact_keys=[];detector=helper.load_detector()
    for record in d['exact']:
        k=record['class'];packet=rows[lookup[k]];assert packet.tolist()==record['rows'];C=-point(packet[0]);assert C and C[1]
        z=[]
        for v in packet:
            delta=v-packet[0];assert np.all(delta%2==0);R=point(delta//2)
            if R.is_zero() or R==C:z.append(None)
            elif R[0]==C[0]:
                assert R==-C;z.append(-(3*C[0]**2+E.a4())/(2*C[1]))
            else:z.append((R[1]+C[1])/(R[0]-C[0]))
        result=detector(z);assert result==record['detector'] and not result['hits'];exact_keys.append(k)
    assert len(exact_keys)==21 and not d['unresolved_classes']
    assert set(map(int,exclusions)).isdisjoint(exact_keys)
    assert set(map(int,exclusions))|set(exact_keys)==set(map(int,keys))
    return {'schema':'curve302.mw-packet-replay.v1','status':'PASS_ALL_130706_FROZEN_PACKET_EXCLUSIONS',
        'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),SEARCH,PACKETS,HELPER,helper.PAIR,CAS/'icarm_curve302.py']},
        'modular_exclusions_replayed':checked,'rational_exclusions_replayed':len(exact_keys),
        'primitive_rank17_embedding_verified':True,'all_packet_parities_verified':True,
        'boundary':'Exact exclusion of a six-disjoint-equal-gap configuration for each frozen13-point packet. No numerical completeness claim, no exclusion of other points or classes or subspaces, no exclusion of all302 parents. Replay shares the chart and detector code.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();signal.alarm(300)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)
