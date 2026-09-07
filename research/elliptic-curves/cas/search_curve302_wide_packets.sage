#!/usr/bin/env sage-python
"""Frozen 12--24-point recognition on five lower-rank302 spaces and core17.

Height190 unchanged; <=4m rays per space; primes10007,10009,10037;
<=64 rational packets per space; one worker,900s per invocation.
Preparation checks containment of every previous13-point packet. Modular
checkpoints permit restart without repeating certified exclusions.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import gzip
import itertools
import json
from pathlib import Path
import runpy
import signal
import sys
import numpy as np
from sage.all import QQ,ZZ,RealField,matrix,pari,gcd,lcm

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
sys.path.insert(0,str(CAS))
LOWER=CAS/'search_curve302_lower_rank_packets.sage';m=runpy.run_path(str(LOWER))
from mw_pair_collision_filter import PairCollisionFilter
ART=m['ART'];LOCAL=ROOT/'artifacts/local/elliptic-curves/curve302-wide-packets'
LIB=CAS/'mw_pair_collision_filter.py';CONTROL=ART/'curve302_mw_pair_collision_control_v1.json'
PROTOCOL=ART/'curve302_wide_packet_protocol_v1.json';OUT=ART/'curve302_wide_packets_v1.json';REPLAY=ART/'curve302_wide_packets_replay_v1.json'
POLICY={'candidate_ids':['32','58','59','60','86','core17'],'height_bound':190,'ray_cap_per_space':4000000,
        'minimum_packet_size':12,'maximum_packet_size':24,'primes':[10007,10009,10037],
        'exact_packet_cap_per_space':64,'seconds_per_invocation':900,'workers':1,
        'selection':'First up to24 numerical-height-ordered rays in every class with at least12, on the previous five spaces and rank17 core; height bound unchanged.'}


def digest(p):return sha256(p.read_bytes()).hexdigest()
def write(p,d):
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n');tmp.replace(p)


def detect_rational(z):
    """Complete rational six-disjoint-pair detector for 12--24 points."""
    n=len(z);assert 12<=n<=24 and len(set(z))==n
    offset=QQ(0)
    while offset in z:offset+=1
    w=[QQ(0) if v is None else 1/(v-offset) for v in z];poles={None};equations=0
    for i,j,k,l in itertools.combinations(range(n),4):
        for ids in [(i,j,k,l),(i,k,j,l),(i,l,j,k)]:
            a,b,c,d=[w[h] for h in ids]
            for sign in [-1,1]:
                ab=a-b;cd=sign*(c-d);A=ab-cd;B=-ab*(c+d)+cd*(a+b);C=ab*c*d-cd*a*b;equations+=1
                if not A:
                    assert B or C
                    if B:poles.add(-C/B)
                else:
                    disc=B*B-4*A*C
                    if disc>=0 and disc.is_square():
                        ds=disc.sqrt();poles.update([(-B+ds)/(2*A),(-B-ds)/(2*A)])
    hits=[]
    pairmodule=m['module'](m['helper'].PAIR,'wide_rational_matchings')
    for pole in sorted(poles,key=lambda k:(k is not None,k or 0)):
        xs=[v if pole is None else (None if v==pole else 1/(v-pole)) for v in w];pairs={}
        for i,j in itertools.combinations(range(n),2):
            if xs[i] is None or xs[j] is None:continue
            diff=abs(xs[i]-xs[j])
            if diff:pairs.setdefault(diff,[]).append((i,j))
        for diff,edges in sorted(pairs.items()):
            if len(edges)<6:continue
            for matching in pairmodule.disjoint_matchings(edges,6):
                centers=sorted((xs[i]+xs[j])/2 for i,j in matching);assert len(set(centers))==6
                den=lcm([a.denominator() for a in centers]+[diff.denominator()]);nums=[ZZ((a-centers[0])*den) for a in centers];div=gcd(nums)
                roots=tuple(a//div for a in nums);roots=min(roots,tuple(roots[-1]-a for a in roots[::-1]))
                hits.append({'pole':None if pole is None else str(pole),'pairs':[list(p) for p in matching],
                             'roots':list(map(int,roots)),'T':str(diff*den/(2*div))})
    return {'finite_chart_offset':str(offset),'quadratic_equations':equations,'candidate_poles_including_infinity':len(poles),'hits':hits}


def select(B,G,bound,cap):
    U=matrix(ZZ,pari(G).qflllgram());R=U.transpose()*G*U;q=pari(R).qfminim(bound,cap,2)
    raw=np.asarray(matrix(ZZ,q[2]),dtype=np.int64).T
    assert int(q[0])==2*len(raw) and len(raw)<cap
    change=np.asarray(U,dtype=np.int64);assert int(np.max(np.abs(raw)))*int(np.max(np.abs(change)))*B.nrows()<2**62
    a=raw@change.T;del raw,q
    leading=a[np.arange(len(a)),np.argmax(a!=0,axis=1)];a[leading<0]*=-1
    heights=np.einsum('ij,jk,ik->i',a,np.asarray(G,dtype=float),a);masks=(a%2)@(1<<np.arange(B.nrows(),dtype=np.int64))
    order=np.lexsort(tuple(a[:,j] for j in reversed(range(B.nrows())))+(heights,masks))
    mm=masks[order];starts=np.r_[0,np.flatnonzero(mm[1:]!=mm[:-1])+1];ends=np.r_[starts[1:],len(order)]
    eligible=starts[ends-starts>=12];keys=mm[eligible];lengths=np.minimum(24,(ends-starts)[ends-starts>=12])
    rows=np.zeros((len(keys),24,B.nrows()),dtype=np.int64)
    for j in range(24):
        good=lengths>j;rows[good,j]=a[order[eligible[good]+j]]
    counts={'ray_count':len(a),'occupied_classes':len(starts),'eligible_packets':len(keys),
            'packet_size_histogram':{str(n):int(np.sum(lengths==n)) for n in range(12,25)},
            'rows_sha256':sha256(rows.astype('<i8').tobytes()).hexdigest(),
            'keys_sha256':sha256(keys.astype('<i8').tobytes()).hexdigest(),
            'lengths_sha256':sha256(lengths.astype('<i8').tobytes()).hexdigest()}
    return keys,lengths,rows,counts


def calibrate():
    c=json.loads(CONTROL.read_text())
    for p,h in c['input_sha256'].items():assert digest(ROOT/p)==h
    B,G,bound,cap,paths=m['old'].inputs(245);E,public=m['target'](245)
    c13=json.loads((ART/'curve302_mw_packet_245_control_v1.json').read_text());positive=next(x for x in c13['exact'] if x['detector']['hits'])
    packet=np.asarray(positive['rows'],dtype=np.int64);zs,point=m['rational_chart'](E,public,B,packet)
    assert detect_rational(zs)==positive['detector']
    packet=np.asarray(c['control_24_point_rows'],dtype=np.int64);zs,point=m['rational_chart'](E,public,B,packet)
    result=detect_rational(zs);verified=m['verify_hits'](E,packet,point,result)
    assert verified and all(v['roots']==[0,106,344,475,594,731] and v['T']=='5801/10' for v in verified)
    # Twelve-point endpoint: remove all distractors and use precisely one recovered matching.
    chosen=sorted({i for pair in result['hits'][0]['pairs'] for i in pair});small=packet[chosen]
    z12,point12=m['rational_chart'](E,public,B,small);result12=detect_rational(z12)
    assert m['verify_hits'](E,small,point12,result12)
    return {'size13_equals_previous':True,'size24_complete_rational_detector':result,'size24_verified':verified,'size12_verified':True}


def prepare():
    assert not PROTOCOL.exists(),'Frozen protocol exists; use search or replay.'
    calibration=calibrate();print('PASS_12_13_24_RATIONAL_CONTROLS',flush=True)
    lower_protocol,oldinputs=m['load_inputs']();old17=json.loads((ART/'curve302_mw_packet_rank17_h190_v1.json').read_text())
    entries=[]
    for entry,B,keys,rows in oldinputs:entries.append((str(entry['candidate_id']),B,keys,rows,entry['enumeration']['ray_count']))
    a=np.load(ART/'curve302_mw_packet_rank17_h190_inputs_v1.npz',allow_pickle=False)
    entries.append(('core17',matrix(ZZ,old17['basis_rows']),a['keys'].astype(np.int64),a['rows'].astype(np.int64),old17['enumeration']['ray_count']))
    cloud=ART/'low_height_mw_sublattices_v1_302_cloud.json.gz';H=matrix(RealField(280),json.loads(gzip.decompress(cloud.read_bytes()))['height_gram'])
    sources=[Path(__file__),LOWER,LIB,CONTROL,m['PROTOCOL'],ART/'curve302_mw_packet_rank17_h190_v1.json',ART/'curve302_mw_packet_rank17_h190_inputs_v1.npz',cloud]
    stamps={str(p.relative_to(ROOT)):digest(p) for p in sources}
    draft={'policy':POLICY,'input_sha256':stamps,'calibration':calibration};write(LOCAL/'declared_protocol.json',draft)
    results=[]
    for label,B,oldkeys,oldrows,oldcount in entries:
        path=ART/f'curve302_wide_packets_{label}_v1.npz';saved=LOCAL/f'prepared_{label}.json'
        if saved.exists():
            item=json.loads(saved.read_text());assert item['input_sha256']==stamps and digest(path)==item['packet_sha256'];results.append(item);continue
        assert not path.exists()
        keys,lengths,rows,counts=select(B,B*H*B.transpose(),190,4000000)
        eligible=lengths>=13
        assert np.array_equal(keys[eligible],oldkeys) and np.array_equal(rows[eligible,:13],oldrows)
        assert counts['ray_count']==oldcount and np.max(np.abs(rows))<128
        np.savez_compressed(path,keys=keys.astype('<u4'),lengths=lengths.astype('u1'),rows=rows.astype('i1'))
        item={'label':label,'rank':B.nrows(),'basis_rows':[list(map(int,v)) for v in B.rows()], 'enumeration':counts,
              'previous13_packets_preserved':len(oldkeys),'packet_path':str(path.relative_to(ROOT)),
              'packet_sha256':digest(path),'input_sha256':stamps}
        write(saved,item);results.append(item);print('PREPARED',label,counts,flush=True)
    write(PROTOCOL,{'schema':'curve302.wide-packets.protocol.v1','status':'FROZEN_AFTER_CONTROLS_AND_OLD_PACKET_CONTAINMENT',**draft,'inputs':results})


def load_inputs():
    d=json.loads(PROTOCOL.read_text());assert d['policy']==POLICY
    for p,h in d['input_sha256'].items():assert digest(ROOT/p)==h
    # Transitively verify the shared reduction, model and exact-transport inputs.
    m['load_inputs']()
    out=[]
    for e in d['inputs']:
        path=ROOT/e['packet_path'];assert digest(path)==e['packet_sha256'];a=np.load(path,allow_pickle=False)
        keys=a['keys'].astype(np.int64);lengths=a['lengths'].astype(np.int64);rows=a['rows'].astype(np.int64);B=matrix(ZZ,e['basis_rows'])
        assert rows.shape==(len(keys),24,e['rank']) and np.all((lengths>=12)&(lengths<=24))
        assert len(set(map(int,keys)))==len(keys) and B.rank()==e['rank'] and all(abs(x)==1 for x in B.smith_form()[0].diagonal())
        for name,array in [('rows',rows),('keys',keys),('lengths',lengths)]:assert sha256(array.astype('<i8').tobytes()).hexdigest()==e['enumeration'][name+'_sha256']
        masks=(rows%2)@(1<<np.arange(e['rank'],dtype=np.int64));valid=np.arange(24)[None,:]<lengths[:,None]
        assert np.all(masks[valid]==np.broadcast_to(keys[:,None],masks.shape)[valid]) and np.all(rows[~valid]==0)
        out.append((e,B,keys,lengths,rows))
    return d,out


def search():
    protocol,inputs=load_inputs();E,public=m['target'](302);results=[]
    for entry,B,keys,lengths,rows in inputs:
        label=entry['label'];path=LOCAL/f'state_{label}.json'
        state=json.loads(path.read_text()) if path.exists() else {'label':label,'protocol_sha256':digest(PROTOCOL),'completed_primes':[],'exclusions':{},'exact':[],'pending':[],'stages':[]}
        assert state['protocol_sha256']==digest(PROTOCOL)
        for p in POLICY['primes']:
            if p in state['completed_primes']:continue
            chart=m['helper'].modular_chart(E,public,B,p);test=PairCollisionFilter(p);positive=unresolved=processed=0
            for i,k in enumerate(keys):
                if str(int(k)) in state['exclusions']:continue
                packet=rows[i,:lengths[i]];decision=test(None if chart is None else chart(packet))
                if decision is False:state['exclusions'][str(int(k))]=p
                elif decision is None:unresolved+=1
                else:positive+=1
                processed+=1
                if processed%1000==0:
                    write(path,state);print('PROGRESS',label,p,processed,'excluded_total',len(state['exclusions']),flush=True)
            stage={'prime':p,'excluded':sum(v==p for v in state['exclusions'].values()),'positive':positive,'unresolved':unresolved,'remaining':len(keys)-len(state['exclusions'])}
            state['stages'].append(stage);state['completed_primes'].append(p);write(path,state);print('STAGE',label,stage,flush=True)
        done={e['class'] for e in state['exact']}|set(state['pending'])
        for i,k in enumerate(keys):
            key=int(k)
            if str(key) in state['exclusions'] or key in done:continue
            if len(state['exact'])>=POLICY['exact_packet_cap_per_space']:state['pending'].append(key);write(path,state);continue
            packet=rows[i,:lengths[i]];zs,point=m['rational_chart'](E,public,B,packet)
            if zs is None or len(set(zs))!=len(zs):state['pending'].append(key);write(path,state);continue
            answer=detect_rational(zs);verified=[];failures=[]
            for hit in answer['hits']:
                try:verified.extend(m['verify_hits'](E,packet,point,{'hits':[hit]}))
                except (AssertionError,ValueError,ArithmeticError) as exc:failures.append({'hit':hit,'failure_type':type(exc).__name__,'message':str(exc)})
            state['exact'].append({'class':key,'size':int(lengths[i]),'detector':answer,'verified_families':verified,'unverified_patterns':failures});write(path,state)
            print('EXACT',label,key,'size',int(lengths[i]),'patterns',len(answer['hits']),'families',len(verified),flush=True)
        state['packet_count']=len(keys);write(path,state);results.append(state)
        if any(e['verified_families'] for e in state['exact']):break
    hits=sum(len(e['verified_families']) for r in results for e in r['exact']);pending=sum(len(r['pending'])+sum(bool(e['unverified_patterns']) for e in r['exact']) for r in results)
    status='VERIFIED_FAMILY_HIT' if hits else ('COMPLETE_FROZEN_PACKET_EXCLUSION' if not pending and len(results)==6 else 'INCOMPLETE')
    write(OUT,{'schema':'curve302.wide-packets.search.v1','status':status,'protocol_sha256':digest(PROTOCOL),'results':results,'family_hits':hits,'unresolved_count':pending,
               'boundary':'Only the frozen12--24-point subsets in six specified spaces. No exclusion of other points, entire spaces, generic ranks or all parents. Numerical selection does not certify a complete canonical-height ball.'})
    print(status,flush=True)


def replay():
    protocol,inputs=load_inputs();d=json.loads(OUT.read_text());assert d['protocol_sha256']==digest(PROTOCOL);E,public=m['target'](302);records=[]
    for entry,B,keys,lengths,rows in inputs:
        matches=[s for s in d['results'] if s['label']==entry['label']]
        if not matches:continue
        r=matches[0];lookup={int(k):i for i,k in enumerate(keys)};checked=[]
        for p in POLICY['primes']:
            chart=m['helper'].modular_chart(E,public,B,p);test=PairCollisionFilter(p)
            for k,prime in r['exclusions'].items():
                if prime!=p:continue
                i=lookup[int(k)];assert chart is not None and test(chart(rows[i,:lengths[i]])) is False;checked.append(int(k))
            print('REPLAY_PRIME',entry['label'],p,len(checked),flush=True)
        exact=[]
        for rec in r['exact']:
            i=lookup[rec['class']];packet=rows[i,:lengths[i]];zs,point=m['rational_chart'](E,public,B,packet);result=detect_rational(zs)
            assert result==rec['detector']
            if rec['verified_families']:assert m['verify_hits'](E,packet,point,result)==rec['verified_families']
            exact.append(rec['class'])
        allkeys=checked+exact+r['pending'];assert len(allkeys)==len(set(allkeys))==len(keys) and set(allkeys)==set(lookup)
        records.append({'label':entry['label'],'modular_checked':len(checked),'rational_checked':len(exact),'pending':len(r['pending'])})
    write(REPLAY,{'schema':'curve302.wide-packets.replay.v1','status':'PASS_REPLAY_OF_RECORDED_RESULTS','input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),PROTOCOL,OUT,LIB]},'records':records,
                  'boundary':'Frozen integer input and exact finite/rational claims replay. Shared code, not an independent implementation or numerical enumeration proof.'})
    print('PASS_REPLAY_OF_RECORDED_RESULTS',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','search','replay']);args=p.parse_args();LOCAL.mkdir(parents=True,exist_ok=True)
    signal.alarm(POLICY['seconds_per_invocation']);pari.allocatemem(1000000000,8000000000);globals()[args.mode]()
