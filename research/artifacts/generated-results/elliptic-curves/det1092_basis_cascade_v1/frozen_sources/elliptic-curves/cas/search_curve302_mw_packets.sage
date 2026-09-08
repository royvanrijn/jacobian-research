#!/usr/bin/env sage-python
"""Bounded six-pair reconstruction in the retained302 rank17 candidate.

Calibrate optimized packet extraction on245 before admitting302. Numerical
height selection only; exact exclusions apply to the frozen13-point packets.
"""
import argparse
from decimal import Decimal
from hashlib import sha256
import gzip
import importlib
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import signal
import sys
import time
import numpy as np
from sage.all import EllipticCurve, QQ, RealField, ZZ, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]; CAS=Path(__file__).resolve().parent
sys.path.insert(0,str(CAS))
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves/curve302-mw-packets'
HELPER=CAS/'select_mw_quartic_packet.sage'
POLICY={'target':302,'rank':17,'height_bound':190,'ray_cap':4000000,
        'packet_size':13,'primes':[1019,1031,1033,1039], 'exact_packet_cap':32,
        'seconds':900,'workers':1,
        'selection':'minimum determinant among retained rank17 first-pass finalists; first13 height-ordered rays per parity class'}


def digest(p):return sha256(p.read_bytes()).hexdigest()
def module(path,name):
    loader=importlib.machinery.SourceFileLoader(name,str(path))
    spec=importlib.util.spec_from_loader(name,loader);m=importlib.util.module_from_spec(spec);loader.exec_module(m);return m


def inputs(target):
    cloud=ART/f'low_height_mw_sublattices_v1_{target}_cloud.json.gz'
    if target==245:
        selector=ART/'latent_lattice_graph_walk_calibration_v1.json'
        entry=next(e for e in json.loads(selector.read_text())['default_controls'] if e['family']=='Fermigier_rank12')
        B=matrix(ZZ,entry['selected_primitive_embedding_matrix_rows']); bound=84; cap=200000
    else:
        selector=ART/'low_height_mw_sublattices_v1_302_selection.json'
        choices=[e for e in json.loads(selector.read_text())['finalists'] if e['rank']==17]
        entry=min(choices,key=lambda e:Decimal(e['determinant']))
        B=matrix(ZZ,entry['primitive_basis_rows']);bound=POLICY['height_bound'];cap=POLICY['ray_cap']
    H=matrix(RealField(280),json.loads(gzip.decompress(cloud.read_bytes()))['height_gram'])
    return B,B*H*B.transpose(),bound,cap,[cloud,selector]


def select(B,G,bound,cap):
    U=matrix(ZZ,pari(G).qflllgram());R=U.transpose()*G*U
    q=pari(R).qfminim(bound,cap,2)
    print('ENUMERATED',int(q[0]),flush=True)
    raw=np.asarray(matrix(ZZ,q[2]),dtype=np.int64).T
    assert int(q[0])==2*len(raw) and len(raw)<cap
    change=np.asarray(U,dtype=np.int64)
    assert int(np.max(np.abs(raw)))*int(np.max(np.abs(change)))*B.nrows()<2**62
    a=raw@change.T
    del raw,q
    leading=a[np.arange(len(a)),np.argmax(a!=0,axis=1)]
    a[leading<0]*=-1
    heights=np.einsum('ij,jk,ik->i',a,np.asarray(G,dtype=float),a)
    masks=(a%2)@(1<<np.arange(B.nrows(),dtype=np.int64))
    order=np.lexsort(tuple(a[:,j] for j in reversed(range(B.nrows())))+(heights,masks))
    mm=masks[order];starts=np.r_[0,np.flatnonzero(mm[1:]!=mm[:-1])+1];ends=np.r_[starts[1:],len(order)]
    eligible=starts[ends-starts>=13];keys=mm[eligible]
    picked=a[order[eligible[:,None]+np.arange(13)[None,:]]]
    return keys,picked,{'ray_count':len(a),'occupied_parity_classes':len(starts),'eligible_packets':len(keys),
                       'packet_array_sha256':sha256(picked.astype('<i8').tobytes()).hexdigest(),
                       'key_array_sha256':sha256(keys.astype('<i8').tobytes()).hexdigest()}


def main(target,check=False):
    start=time.monotonic(); helper=module(HELPER,'packet_helpers');B,G,bound,cap,paths=inputs(target)
    LOCAL.mkdir(parents=True,exist_ok=True)
    stamp={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),HELPER,helper.PAIR,*paths,CAS/f'icarm_curve{target}.py']}
    control_path=ART/'curve302_mw_packet_245_control_v1.json'
    if target==302:
        control=json.loads((LOCAL/'optimized245.json').read_text())
        assert control['status']=='PASS_PACKET_EXTRACTION_EQUALS_245_CONTROL'
        assert control['input_sha256'][str(Path(__file__).relative_to(ROOT))]==digest(Path(__file__))
        assert control['input_sha256'][str(HELPER.relative_to(ROOT))]==digest(HELPER)
    keys,rows,enumeration=select(B,G,bound,cap)
    np.savez_compressed(LOCAL/f'packets{target}.npz',keys=keys,rows=rows)
    if target==245:
        control=json.loads(control_path.read_text());old=json.loads((helper.LOCAL/'packets.json').read_text())
        assert helper.digest(HELPER)==control['input_sha256'][str(HELPER.relative_to(ROOT))]
        assert {str(int(k)):v.tolist() for k,v in zip(keys,rows)}==old['packets']
        result={'status':'PASS_PACKET_EXTRACTION_EQUALS_245_CONTROL','input_sha256':stamp,
                'control_sha256':digest(control_path),'enumeration':enumeration}
        (LOCAL/'optimized245.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        print(result['status'],flush=True);return
    target_module=importlib.import_module(f'icarm_curve{target}')
    E0=EllipticCurve(QQ,list(map(QQ,target_module.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    E=E0.short_weierstrass_model();iso=E0.isomorphism_to(E)
    public=[iso(E0(list(map(QQ,P)))) for P in target_module.POINTS]
    remaining=list(range(len(keys)));exclusions={};stages=[]
    protocol={'policy':POLICY,'input_sha256':stamp,'enumeration':enumeration,'basis_rows':[list(map(int,v)) for v in B.rows()]}
    (LOCAL/'protocol302.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n')
    for p in POLICY['primes']:
        chart=helper.modular_chart(E,public,B,p); rejected=[];unresolved=[]
        for position,index in enumerate(remaining):
            zs=None if chart is None else chart(rows[index])
            decision=helper.modular_pairs(zs,p)
            if decision is False:rejected.append(index)
            elif decision is None:unresolved.append(index)
            if (position+1)%5000==0:
                print('PRIME_PROGRESS',p,position+1,len(remaining),round(time.monotonic()-start,1),flush=True)
                (LOCAL/f'progress_{p}.json').write_text(json.dumps({'processed':position+1,'remaining_indices':remaining,'rejected_indices':rejected,'unresolved_indices':unresolved}))
        rejected_set=set(rejected)
        for index in rejected:exclusions[str(int(keys[index]))]=p
        remaining=[i for i in remaining if i not in rejected_set]
        stages.append({'prime':p,'excluded':len(rejected),'unresolved':len(unresolved),'remaining':len(remaining)})
        (LOCAL/f'after_{p}.json').write_text(json.dumps({'stages':stages,'exclusions':exclusions,'remaining_indices':remaining},sort_keys=True)+'\n')
        print('MODULAR',stages[-1],flush=True)
    detector=helper.load_detector();exact=[];pending=[]
    basispoints=[sum((n*P for n,P in zip(v,public)),E(0)) for v in B.rows()]
    def point(v):return sum((ZZ(n)*P for n,P in zip(v,basispoints)),E(0))
    for pos,index in enumerate(remaining):
        k=int(keys[index]);packet=rows[index]
        if pos>=POLICY['exact_packet_cap']:pending.append(k);continue
        C=-point(packet[0]);z=[]
        if C.is_zero() or C[1]==0:pending.append(k);continue
        for v in packet:
            delta=v-packet[0];assert np.all(delta%2==0);R=point(delta//2)
            if R.is_zero() or R==C:z.append(None)
            elif R[0]==C[0]:
                assert R==-C;z.append(-(3*C[0]**2+E.a4())/(2*C[1]))
            else:z.append((R[1]+C[1])/(R[0]-C[0]))
        if len(set(z))!=13:pending.append(k);continue
        result=detector(z)
        exact.append({'class':k,'rows':packet.tolist(),'detector':result})
        print('EXACT',k,'hits',len(result['hits']),flush=True)
        (LOCAL/'exact302.json').write_text(json.dumps(exact,sort_keys=True)+'\n')
    result={**protocol,'schema':'curve302.mw-packet-search.v1',
        'status':'COMPLETE_FROZEN_PACKET_SCREEN' if not pending else 'INCOMPLETE_EXACT_PACKETS',
        'optimized_control_sha256':digest(LOCAL/'optimized245.json'),
        'stages':stages,'modular_exclusions':exclusions,'exact':exact,'unresolved_classes':pending,
        'boundary':'Only the first13 numerical-height-ordered rays in each eligible parity class of this one retained rank17 candidate are screened. A negative result excludes the six-pair pattern in those packets, not other vectors, higher heights, other MW subspaces, or all parents. Any hit still requires generic construction and exact fibre/section transport.'}
    out=ART/'curve302_mw_packet_rank17_h190_v1.json'
    if check:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],'seconds',round(time.monotonic()-start,1),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--target',type=int,choices=[245,302],required=True);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    signal.alarm(POLICY['seconds']);pari.allocatemem(1000000000,8000000000)
    main(args.target,args.check)
