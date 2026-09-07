#!/usr/bin/env sage-python
"""Calibrate variable-size pair collisions and replay saved302 exclusions.

Compare all4092 control packets at1019, recognize a24-point supplied-image
control at10007, and check all107094 recorded302 modular exclusions by the
new pole-enumeration algorithm. One worker,300 seconds. No new302 search.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
import numpy as np
from sage.all import ZZ,matrix,pari

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
DRIVER=CAS/'search_curve302_lower_rank_packets.sage'
LIB=CAS/'mw_pair_collision_filter.py'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_mw_pair_collision_control_v1.json'


def build():
    m=runpy.run_path(str(DRIVER))
    from mw_pair_collision_filter import PairCollisionFilter
    helper=m['helper'];B,G,bound,cap,paths=m['old'].inputs(245)
    keys,packets,counts=m['old'].select(B,G,bound,cap)
    control=m['ART']/'curve302_mw_packet_245_control_v1.json';c=json.loads(control.read_text())
    selected={str(int(k)):v.tolist() for k,v in zip(keys,packets)}
    assert sha256(json.dumps(selected,sort_keys=True).encode()).hexdigest()==c['packet_sha256']
    E,public=m['target'](245);p=1019;chart=helper.modular_chart(E,public,B,p);new=PairCollisionFilter(p)
    counts={'true':0,'false':0,'unresolved':0}
    for packet in packets:
        z=chart(packet);old=helper.modular_pairs(z,p);answer=new(z);assert old==answer
        counts['true' if answer is True else 'false' if answer is False else 'unresolved']+=1
    assert counts=={'true':1,'false':3732,'unresolved':359}
    positive=next(x for x in c['exact'] if x['detector']['hits']);chosen=sorted({i for pair in positive['detector']['hits'][0]['pairs'] for i in pair})
    packet=[positive['rows'][i] for i in chosen]
    def ray(v):
        v=tuple(map(int,v));return v if next(x for x in v if x)>0 else tuple(-x for x in v)
    seen=set(map(ray,packet))
    for multiple in range(1,11):
        for j in range(12):
            if len(packet)==24:break
            v=np.asarray(packet[0],dtype=np.int64).copy();v[j]+=2*multiple
            if ray(v) not in seen:packet.append(v.tolist());seen.add(ray(v))
        if len(packet)==24:break
    assert len(packet)==24
    packet.sort(key=lambda v:sha256(str(v).encode()).hexdigest());packet=np.asarray(packet,dtype=np.int64)
    p=10007;new=PairCollisionFilter(p);chart=helper.modular_chart(E,public,B,p);hits=new(chart(packet),witness=True);assert hits
    known={ray(positive['rows'][i]) for i in chosen}
    first13_known=sum(ray(v) in known for v in packet[:13]);assert first13_known==7
    assert helper.modular_pairs(chart(packet[:13]),p) is False
    verified=[]
    for h in hits:
        ids=sorted({i for pair in h['pairs'] for i in pair});ids.append(next(i for i in range(24) if i not in ids))
        part=packet[ids];z,point=m['rational_chart'](E,public,B,part);detected=helper.load_detector()(z)
        verified.extend(m['verify_hits'](E,part,point,detected))
    assert verified and all(v['roots']==[0,106,344,475,594,731] and v['T']=='5801/10' for v in verified)
    print('PASS_4092_EQUIVALENCE_AND_24_POINT_RECOVERY',flush=True)
    protocol,inputs=m['load_inputs']();result=json.loads(m['OUT'].read_text());assert result['protocol_sha256']==m['digest'](m['PROTOCOL'])
    E,public=m['target'](302);replayed=[]
    for entry,B,keys,rows in inputs:
        rec=next(r for r in result['candidates'] if r['candidate_id']==entry['candidate_id']);lookup={int(k):i for i,k in enumerate(keys)};count=0
        for p in m['POLICY']['primes']:
            chart=helper.modular_chart(E,public,B,p);new=PairCollisionFilter(p)
            for key,prime in rec['modular_exclusions'].items():
                if prime!=p:continue
                assert new(chart(rows[lookup[int(key)]])) is False;count+=1
        replayed.append({'candidate_id':entry['candidate_id'],'modular_exclusions_verified':count})
        print('REPLAYED',replayed[-1],flush=True)
    assert sum(r['modular_exclusions_verified'] for r in replayed)==107094
    paths=[Path(__file__),LIB,DRIVER,m['HELPER'],helper.PAIR,control,m['PROTOCOL'],m['OUT'],*paths]
    return {'schema':'curve302.pair-collision-control.v1','status':'PASS_VARIABLE_SIZE_RECOGNITION_AND_107094_RECHECKS',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        'control_4092_comparison_at1019':counts,'control_24_point_rows':packet.tolist(),
        'control_first13_known_images':first13_known,'control_first13_modular_excluded':True,
        'control_24_point_prime':10007,'control_24_point_modular_hits':hits,'control_24_point_verified_families':verified,
        'target_modular_rechecks':replayed,
        'boundary':'The new finite-field pole enumeration agrees with the old algorithm on4092 packets and on every107094 saved302 modular exclusion. The24-point example supplies the known12 images mixed with12 distractors; it calibrates recognition, not blind point selection. Generic identity and exact12-image transport are reverified. No larger302 point packet or new parent has been searched. Elliptic reduction and rational family verification share existing code.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();signal.alarm(300);pari.allocatemem(1000000000,8000000000)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)
