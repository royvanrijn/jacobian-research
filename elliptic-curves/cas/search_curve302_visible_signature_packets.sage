#!/usr/bin/env sage-python
"""Frozen12--24-point recognition on new signature spaces106 and130.

Reuse the calibrated wide-packet engine with a new two-space protocol.
Numerical height190, <=4m rays per space, three primes, <=64 rational
packets per space, one worker and900 seconds per invocation.
"""
import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path
import runpy
import signal
import numpy as np
from sage.all import RealField,ZZ,matrix,pari

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent
ENGINE=CAS/'search_curve302_wide_packets.sage'
wide=runpy.run_path(str(ENGINE));g=wide['search'].__globals__
ART=g['ART'];INTAKE=ART/'curve302_visible_subspace_intake_v1.json'
PROTOCOL=ART/'curve302_signature_packet_protocol_v1.json';OUT=ART/'curve302_signature_packets_v1.json';REPLAY=ART/'curve302_signature_packets_replay_v1.json'
LOCAL=ROOT/'artifacts/local/elliptic-curves/curve302-signature-packets'
POLICY={**g['POLICY'],'candidate_ids':['106','130'],
        'selection':'First up to24 numerical-height-ordered rays in each class with at least12, on the two maximal untested visible-signature spaces106/rank13 and130/rank11. Candidate dimension does not prescribe generic parent rank.'}
g.update(PROTOCOL=PROTOCOL,OUT=OUT,REPLAY=REPLAY,LOCAL=LOCAL,POLICY=POLICY,__file__=str(Path(__file__).resolve()))


def digest(p):return sha256(p.read_bytes()).hexdigest()


def prepare():
    assert not PROTOCOL.exists(),'Frozen protocol exists; use search or replay.'
    calibration=g['calibrate']();print('PASS_12_13_24_RATIONAL_CONTROLS',flush=True)
    d=json.loads(INTAKE.read_text())
    for p,h in d['input_sha256'].items():assert digest(ROOT/p)==h
    assert d['maximal_new_candidate_ids']==[106,130]
    entries=[next(e for e in d['records'] if e['candidate_id']==idx) for idx in [106,130]]
    cloud=ART/'low_height_mw_sublattices_v1_302_cloud.json.gz';H=matrix(RealField(280),json.loads(gzip.decompress(cloud.read_bytes()))['height_gram'])
    paths=[Path(__file__),ENGINE,g['LOWER'],g['LIB'],g['CONTROL'],INTAKE,cloud]
    stamps={str(p.relative_to(ROOT)):digest(p) for p in paths};results=[]
    g['write'](LOCAL/'declared_protocol.json',{'policy':POLICY,'input_sha256':stamps})
    for e in entries:
        B=matrix(ZZ,e['primitive_basis_rows']);assert B.rank()==e['rank'] and all(abs(x)==1 for x in B.smith_form()[0].diagonal())
        keys,lengths,rows,counts=g['select'](B,B*H*B.transpose(),190,4000000)
        assert not len(rows) or np.max(np.abs(rows))<128
        path=ART/f"curve302_signature_packets_{e['candidate_id']}_v1.npz";assert not path.exists()
        np.savez_compressed(path,keys=keys.astype('<u4'),lengths=lengths.astype('u1'),rows=rows.astype('i1'))
        result={'label':str(e['candidate_id']),'rank':e['rank'],'basis_rows':e['primitive_basis_rows'],'enumeration':counts,
                'packet_path':str(path.relative_to(ROOT)),'packet_sha256':digest(path),'input_sha256':stamps}
        results.append(result);g['write'](LOCAL/'preparation_checkpoint.json',results);print('PREPARED',e['candidate_id'],counts,flush=True)
    g['write'](PROTOCOL,{'schema':'curve302.signature-packets.protocol.v1','status':'FROZEN_AFTER_RATIONAL_CONTROLS','policy':POLICY,'input_sha256':stamps,'calibration':calibration,'inputs':results})


def search():
    # The inherited search's legacy finalizer assumes six inputs. Recompute
    # completion below against this protocol's exact two-label domain.
    g['search']();d=json.loads(OUT.read_text());protocol,inputs=g['load_inputs']()
    seen=[]
    for e,B,keys,lengths,rows in inputs:
        records=[r for r in d['results'] if r['label']==e['label']]
        if not records:continue
        assert len(records)==1;r=records[0];accounted=list(map(int,r['exclusions']))+[x['class'] for x in r['exact']]+r['pending']
        assert len(accounted)==len(set(accounted))==len(keys) and set(accounted)==set(map(int,keys));seen.append(e['label'])
    complete=set(seen)==set(POLICY['candidate_ids']) and not d['unresolved_count']
    d['schema']='curve302.signature-packets.search.v1'
    d['status']='VERIFIED_FAMILY_HIT' if d['family_hits'] else ('COMPLETE_FROZEN_PACKET_EXCLUSION' if complete else 'INCOMPLETE')
    g['write'](OUT,d);print('FINAL',d['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','search','replay']);args=p.parse_args()
    LOCAL.mkdir(parents=True,exist_ok=True);signal.alarm(900);pari.allocatemem(1000000000,8000000000)
    if args.mode=='prepare':prepare()
    elif args.mode=='search':search()
    else:g['replay']()
