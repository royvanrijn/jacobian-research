#!/usr/bin/env sage-python
"""Replay the245 recognizer after enlarging its primitive subspace to14--16.

Frozen packets are checked exactly by default. Optional numerical selection
replay uses the original height84 bound, <=4m rays per rank and300 seconds.
This is a development robustness control, not a held-out statistical test.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
import numpy as np
from sage.all import QQ,ZZ,matrix,pari

ROOT=Path(__file__).resolve().parents[2]
DRIVER=ROOT/'elliptic-curves/cas/search_curve302_lower_rank_packets.sage'
INPUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_mw14_16_enlarged_control_inputs_v1.json'
OUT=INPUT.with_name('curve302_mw14_16_enlarged_control_v1.json')


def build(recompute=False):
    m=runpy.run_path(str(DRIVER));data=json.loads(INPUT.read_text());assert data['height_bound']==84 and data['ray_cap']==4000000
    B,G,bound,cap,paths=m['old'].inputs(245)
    control=m['ART']/'curve302_mw_packet_245_control_v1.json';c=json.loads(control.read_text())
    for p,h in c['input_sha256'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    positive=next(e for e in c['exact'] if e['detector']['hits'])
    chosen=sorted({i for pair in positive['detector']['hits'][0]['pairs'] for i in pair})
    W=matrix(ZZ,[positive['rows'][i] for i in chosen])*B
    E,public=m['target'](245);records=[];basis=matrix.identity(ZZ,20)
    def ray(v):
        v=tuple(map(int,v));return v if next(x for x in v if x)>0 else tuple(-x for x in v)
    for r in data['records']:
        rank=r['rank']
        for v in basis.rows():
            if B.nrows()==rank:break
            if B.stack(matrix(ZZ,[v])).rank()>B.nrows():B=B.stack(matrix(ZZ,[v])).row_module().saturation().basis_matrix()
        assert [list(map(int,v)) for v in B.rows()]==r['basis_rows']
        assert B.rank()==rank and all(abs(x)==1 for x in B.smith_form()[0].diagonal())
        C=matrix(ZZ,B.solve_left(W));assert C*B==W and [list(map(int,v)) for v in C.rows()]==r['oracle_rows']
        packet=np.asarray(r['selected_rows'],dtype=np.int64);assert packet.shape==(13,rank)
        key=sum((int(x)%2)<<i for i,x in enumerate(C[0]))
        assert key==r['control_class'] and np.all((packet%2)@(1<<np.arange(rank,dtype=np.int64))==key)
        assert len(set(map(ray,C.rows())) & set(map(ray,packet)))==12
        if recompute:
            import gzip
            from sage.all import RealField
            H=matrix(RealField(280),json.loads(gzip.decompress(paths[0].read_bytes()))['height_gram'])
            keys,packets,counts=m['old'].select(B,B*H*B.transpose(),84,4000000)
            assert counts==r['counts'] and packets[list(map(int,keys)).index(key)].tolist()==r['selected_rows']
        zs,point=m['rational_chart'](E,public,B,packet);detected=m['helper'].load_detector()(zs)
        verified=m['verify_hits'](E,packet,point,detected)
        assert len(verified)==1 and verified[0]['roots']==positive['verified_families'][0]['roots'] and verified[0]['T']==positive['verified_families'][0]['T']
        records.append({'rank':rank,'parity_class':key,'visible_images_retained':12,'selection_counts':r['counts'],
                        'detector':detected,'verified_families':verified})
    return {'schema':'curve302.enlarged-control.v1','status':'PASS_245_RECOVERY_IN_RANK14_15_16_CONTAINERS',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),DRIVER,INPUT,control,*paths]},
        'records':records,'boundary':'The known245 parent subspace is enlarged deterministically with public-basis directions. All twelve visible construction images survive the first13 selection, and exact family/section recovery succeeds in these three examples. Selection uses no image list or parity class; the retained oracle is opened afterwards to audit visibility. This does not establish uniform sensitivity on302 or every larger container. Numerical enumeration is not an interval-certified height-ball census.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');p.add_argument('--recompute-selection',action='store_true');args=p.parse_args();signal.alarm(300);pari.allocatemem(1000000000,8000000000)
    result=build(args.recompute_selection)
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)
