#!/usr/bin/env sage-python
"""Recognize Mestre incidence in five frozen rank14--16 candidate subspaces.

One numerical height190 selection per subspace, <=4m rays each, first13
rays per parity class. Four exact modular stages, <=32 rational packets per
subspace. One worker, 900 seconds per invocation, resumable by subspace.
--prepare calibrates extraction and transport on245 and freezes integer
inputs. --search executes; --replay verifies from inputs without enumeration.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import gzip
import importlib
import json
from pathlib import Path
import signal
import sys
import time
import numpy as np
from sage.all import EllipticCurve,PolynomialRing,QQ,RealField,ZZ,matrix,pari,prod

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
sys.path.insert(0,str(CAS))
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves/curve302-mw14-16-packets'
OLD=CAS/'search_curve302_mw_packets.sage';HELPER=CAS/'select_mw_quartic_packet.sage'
INTAKE=ART/'curve302_inverse_kihara_and_rank14_16_intake_v1.json'
PROTOCOL=ART/'curve302_mw14_16_packet_protocol_v1.json'
OUT=ART/'curve302_mw14_16_packets_v1.json'
REPLAY=ART/'curve302_mw14_16_packets_replay_v1.json'
POLICY={'candidate_ids':[32,58,59,60,86],'height_bound':190,'ray_cap_per_subspace':4000000,
        'packet_size':13,'primes':[1019,1031,1033,1039],'exact_packet_cap_per_subspace':32,
        'seconds_per_invocation':900,'workers':1,
        'selection':'All five retained primitive rank14--16 subspaces outside the old rank17 core; first13 numerical-height-ordered rays per parity class.'}


def digest(p):return sha256(p.read_bytes()).hexdigest()
def write(path,d):path.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
def module(path,name):
    import importlib.machinery,importlib.util
    loader=importlib.machinery.SourceFileLoader(name,str(path));spec=importlib.util.spec_from_loader(name,loader)
    m=importlib.util.module_from_spec(spec);loader.exec_module(m);return m

old=module(OLD,'old_packet_driver');helper=module(HELPER,'lower_rank_packet_helper')


def target(n):
    public=importlib.import_module(f'icarm_curve{n}')
    E0=EllipticCurve(QQ,list(map(QQ,public.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    E=E0.short_weierstrass_model();iso=E0.isomorphism_to(E)
    return E,[iso(E0(list(map(QQ,p)))) for p in public.POINTS]


def rational_chart(E,public,B,packet):
    basis=[sum((n*p for n,p in zip(v,public)),E(0)) for v in B.rows()]
    def point(v):return sum((ZZ(n)*p for n,p in zip(v,basis)),E(0))
    C=-point(packet[0]);z=[]
    if C.is_zero() or not C[1]:return None,point
    for v in packet:
        delta=v-packet[0];assert np.all(delta%2==0);R=point(delta//2)
        if R.is_zero() or R==C:z.append(None)
        elif R[0]==C[0]:
            assert R==-C;z.append(-(3*C[0]**2+E.a4())/(2*C[1]))
        else:z.append((R[1]+C[1])/(R[0]-C[0]))
    return z,point


def verify_hits(E,packet,point,result):
    from mestre_root_tuples import SixRootMestreConstruction
    from nagao_1994 import primitive_visible_points,quartic_point_to_short_jacobian
    verified=[]
    for hit in result['hits']:
        cons=SixRootMestreConstruction(tuple(Fraction(str(a)) for a in hit['roots']))
        RT=PolynomialRing(QQ,'T');T=RT.gen();RX=PolynomialRing(RT,'X');X=RX.gen()
        samples=[cons.primitive_quartic_coefficients(Fraction(i)) for i in range(1,8)]
        coeff=[RT.lagrange_polynomial([(QQ(i+1),QQ(samples[i][j])) for i in range(7)]) for j in range(5)]
        quartic=sum(coeff[i]*X**i for i in range(5));product=prod((X-QQ(r)-T)*(X-QQ(r)+T) for r in hit['roots']);g=X**6
        for j in range(5,-1,-1):g+=(product[6+j]-(g*g)[6+j])/2*X**j
        assert g*g-product==T*T*QQ(cons.quartic_content)*quartic
        t0=Fraction(str(hit['T']));J=EllipticCurve(QQ,list(map(QQ,cons.primitive_jacobian_coefficients(t0))))
        points=[J(list(map(QQ,quartic_point_to_short_jacobian(cons,t0,p)))) for p in primitive_visible_points(cons,t0)]
        chosen=sorted({i for pair in hit['pairs'] for i in pair});images={i:point(packet[i]) for i in chosen};matches=[]
        for transport in J.isomorphisms(E):
            ids=[[i for i in chosen if transport(p)==images[i] or transport(p)==-images[i]] for p in points]
            if all(len(v)==1 for v in ids) and len({v[0] for v in ids})==12:
                matches.append({'u_r_s_t':list(map(str,transport.tuple())),'indices':[v[0] for v in ids]})
        assert matches
        verified.append({'roots':hit['roots'],'T':hit['T'],'generic_identity_checked':True,'twelve_images_match_up_to_sign':matches})
    return verified


def prepare():
    assert not PROTOCOL.exists(),'Frozen protocol already exists; search or replay it.'
    B,G,bound,cap,paths=old.inputs(245);keys,rows,counts=old.select(B,G,bound,cap)
    snapshot=json.loads((helper.LOCAL/'packets.json').read_text())
    assert {str(int(k)):v.tolist() for k,v in zip(keys,rows)}==snapshot['packets']
    control=ART/'curve302_mw_packet_245_control_v1.json';c=json.loads(control.read_text())
    for path,h in c['input_sha256'].items():assert digest(ROOT/path)==h
    hit=next(e for e in c['exact'] if e['detector']['hits']);packet=rows[list(map(int,keys)).index(hit['class'])]
    E,public=target(245);zs,point=rational_chart(E,public,B,packet)
    recovered=helper.load_detector()(zs);assert recovered==hit['detector']
    assert verify_hits(E,packet,point,recovered)==hit['verified_families']
    print('PASS_245_EXTRACTION_AND_GENERIC_TRANSPORT',flush=True)
    source=json.loads(INTAKE.read_text());candidates={e['candidate_index']:e for e in source['rank14_16_intake'] if not e['contained_in_tested_rank17_core']}
    assert sorted(candidates)==POLICY['candidate_ids']
    cloud=ART/'low_height_mw_sublattices_v1_302_cloud.json.gz'
    H=matrix(RealField(280),json.loads(gzip.decompress(cloud.read_bytes()))['height_gram'])
    entries=[]
    for idx in POLICY['candidate_ids']:
        entry=candidates[idx];B=matrix(ZZ,entry['primitive_basis_rows']);assert B.rank()==entry['rank']
        assert all(abs(v)==1 for v in B.smith_form()[0].diagonal())
        keys,rows,counts=old.select(B,B*H*B.transpose(),POLICY['height_bound'],POLICY['ray_cap_per_subspace'])
        assert max(abs(int(rows.min())),abs(int(rows.max())))<128
        path=ART/f'curve302_mw14_16_packets_candidate{idx}_v1.npz'
        assert not path.exists();np.savez_compressed(path,keys=keys.astype('<u4'),rows=rows.astype('i1'))
        entries.append({'candidate_id':idx,'rank':entry['rank'],'basis_rows':entry['primitive_basis_rows'],
                        'enumeration':counts,'packet_path':str(path.relative_to(ROOT)),'packet_sha256':digest(path)})
        write(LOCAL/'preparation_checkpoint.json',{'completed':entries})
        print('PREPARED',idx,entry['rank'],counts,flush=True)
    paths=[Path(__file__),OLD,HELPER,helper.PAIR,INTAKE,cloud,control,CAS/'icarm_curve302.py',CAS/'icarm_curve245.py',CAS/'mestre_root_tuples.py',CAS/'nagao_1994.py']
    result={'schema':'curve302.lower-rank-packets.protocol.v1','status':'FROZEN_AFTER_245_CONTROL','policy':POLICY,
            'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in paths},'candidate_inputs':entries,
            'control':{'all_4092_packets_equal':True,'recovered_class':hit['class'],'family_and_twelve_images_reverified':True},
            'boundary':'Numerical selection defines frozen packets only; no complete canonical-height-ball claim. The five subspace ranks are candidate dimensions, not proved generic ranks. Any hit requires explicit family and generic-basis work.'}
    write(PROTOCOL,result)


def load_inputs():
    d=json.loads(PROTOCOL.read_text());assert d['policy']==POLICY
    for p,h in d['input_sha256'].items():assert digest(ROOT/p)==h
    result=[]
    for e in d['candidate_inputs']:
        path=ROOT/e['packet_path'];assert digest(path)==e['packet_sha256'];a=np.load(path,allow_pickle=False)
        keys=a['keys'].astype(np.int64);rows=a['rows'].astype(np.int64);B=matrix(ZZ,e['basis_rows'])
        assert rows.shape==(len(keys),13,e['rank']) and B.rank()==e['rank']
        assert all(abs(v)==1 for v in B.smith_form()[0].diagonal())
        assert sha256(rows.astype('<i8').tobytes()).hexdigest()==e['enumeration']['packet_array_sha256']
        assert sha256(keys.astype('<i8').tobytes()).hexdigest()==e['enumeration']['key_array_sha256']
        assert len(set(map(int,keys)))==len(keys)
        assert np.all((rows%2)@(1<<np.arange(e['rank'],dtype=np.int64))==keys[:,None])
        result.append((e,B,keys,rows))
    return d,result


def search():
    protocol,inputs=load_inputs();E,public=target(302);results=[];detector=helper.load_detector()
    for entry,B,keys,rows in inputs:
        idx=entry['candidate_id'];checkpoint=LOCAL/f'completed_{idx}.json'
        if checkpoint.exists():
            saved=json.loads(checkpoint.read_text());assert saved['protocol_sha256']==digest(PROTOCOL)
            results.append(saved);continue
        remaining=list(range(len(keys)));exclusions={};stages=[]
        for p in POLICY['primes']:
            chart=helper.modular_chart(E,public,B,p);rejected=[];unresolved=[]
            for pos,i in enumerate(remaining):
                decision=helper.modular_pairs(None if chart is None else chart(rows[i]),p)
                if decision is False:rejected.append(i)
                elif decision is None:unresolved.append(i)
                if (pos+1)%5000==0:
                    write(LOCAL/f'progress_{idx}.json',{'prime':p,'processed':pos+1,'total':len(remaining)})
                    print('PROGRESS',idx,p,pos+1,len(remaining),flush=True)
            rejected_set=set(rejected)
            for i in rejected:exclusions[str(int(keys[i]))]=p
            remaining=[i for i in remaining if i not in rejected_set]
            stages.append({'prime':p,'excluded':len(rejected),'unresolved':len(unresolved),'remaining':len(remaining)})
            write(LOCAL/f'stages_{idx}.json',{'stages':stages,'modular_exclusions':exclusions,'remaining_indices':remaining})
            print('MODULAR',idx,stages[-1],flush=True)
        exact=[];pending=[]
        for pos,i in enumerate(remaining):
            key=int(keys[i]);packet=rows[i]
            if pos>=POLICY['exact_packet_cap_per_subspace']:pending.append(key);continue
            zs,point=rational_chart(E,public,B,packet)
            if zs is None or len(set(zs))<13:pending.append(key);continue
            result=detector(zs);verified=verify_hits(E,packet,point,result)
            exact.append({'class':key,'detector':result,'verified_families':verified})
            write(LOCAL/f'exact_{idx}.json',exact);print('EXACT',idx,key,'hits',len(verified),flush=True)
        result={'candidate_id':idx,'rank':entry['rank'],'packet_count':len(keys),'protocol_sha256':digest(PROTOCOL),
                'stages':stages,'modular_exclusions':exclusions,'exact':exact,'unresolved_classes':pending}
        write(checkpoint,result);results.append(result)
    hits=sum(len(e['verified_families']) for r in results for e in r['exact']);pending=sum(len(r['unresolved_classes']) for r in results)
    result={'schema':'curve302.lower-rank-packets.search.v1','status':'VERIFIED_FAMILY_HIT' if hits else ('INCOMPLETE_PACKETS' if pending else 'COMPLETE_FROZEN_PACKET_EXCLUSION'),
            'protocol_sha256':digest(PROTOCOL),'candidates':results,'family_hits':hits,'unresolved_count':pending,
            'boundary':'Exact results for five sets of frozen thirteen-point packets only. No exclusion of other point choices, heights, classes without a packet, or all parents. No generic rank follows from a candidate dimension.'}
    write(OUT,result);print(result['status'],flush=True)


def replay():
    protocol,inputs=load_inputs();result=json.loads(OUT.read_text());assert result['protocol_sha256']==digest(PROTOCOL)
    E,public=target(302);detector=helper.load_detector();records=[]
    for entry,B,keys,rows in inputs:
        record=next(r for r in result['candidates'] if r['candidate_id']==entry['candidate_id'])
        lookup={int(k):i for i,k in enumerate(keys)};checked=[]
        for p in POLICY['primes']:
            chart=helper.modular_chart(E,public,B,p)
            for k,prime in record['modular_exclusions'].items():
                if prime!=p:continue
                assert chart is not None and helper.modular_pairs(chart(rows[lookup[int(k)]]),p) is False
                checked.append(int(k))
        exact=[]
        for rec in record['exact']:
            packet=rows[lookup[rec['class']]];zs,point=rational_chart(E,public,B,packet);answer=detector(zs)
            assert answer==rec['detector'] and verify_hits(E,packet,point,answer)==rec['verified_families'];exact.append(rec['class'])
        pending=record['unresolved_classes'];allkeys=checked+exact+pending
        assert len(allkeys)==len(set(allkeys))==len(keys) and set(allkeys)==set(lookup)
        records.append({'candidate_id':entry['candidate_id'],'modular_checked':len(checked),'rational_checked':len(exact),'unresolved':len(pending)})
        print('REPLAYED',records[-1],flush=True)
    output={'schema':'curve302.lower-rank-packets.replay.v1','status':'PASS_COMPLETE_FROZEN_PACKET_REPLAY',
            'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [PROTOCOL,OUT,Path(__file__)]},'candidates':records,
            'boundary':'Replays every claimed finite/rational result from the frozen integer words; shares chart and detector implementation. No numerical height enumeration is needed or certified.'}
    write(REPLAY,output);print(output['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','search','replay']);args=p.parse_args()
    LOCAL.mkdir(parents=True,exist_ok=True);signal.alarm(POLICY['seconds_per_invocation']);pari.allocatemem(1000000000,8000000000)
    globals()[args.mode]()
