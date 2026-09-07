#!/usr/bin/env python3
"""Exact F2[S3] blocks and H1/H2 of the standard module for all subgroups."""
import argparse
from itertools import permutations,product
from pathlib import Path
import retrospective as r

OUTPUT=r.OUT/'rank_jump_standard_s3_class_module_v1.json'


def compute():
    G=list(permutations(range(3)));unit=G.index((0,1,2));assert unit==0
    table=[[G.index(tuple(g[h[i]] for i in range(3))) for h in G] for g in G]
    def act(g,v):
        bits=[v&1,(v>>1)&1,(v&1)^((v>>1)&1)];out=[0]*3
        for i in range(3):out[G[g][i]]=bits[i]
        assert out[2]==out[0]^out[1];return out[0]|out[1]<<1
    assert all(act(table[g][h],v)==act(g,act(h,v)) for g,h,v in product(range(6),range(6),range(4)))
    subgroups=[]
    for mask in range(1,64,2):
        H=[g for g in range(6) if mask>>g&1]
        if all(table[g][h] in H for g,h in product(H,repeat=2)):subgroups.append(H)
    assert sorted(map(len,subgroups))==[1,2,2,2,3,6]
    records=[]
    for H in subgroups:
        non=[g for g in H if g!=unit];tuples={d:list(product(non,repeat=d)) for d in (1,2,3)}
        pos={d:{v:i for i,v in enumerate(tuples[d])} for d in tuples}
        def val(word,d,key):return (word>>(2*pos[d][key]))&3 if key in pos[d] else 0
        def d1(word):
            out=0
            for i,(g,h) in enumerate(tuples[2]):out|=(act(g,val(word,1,(h,)))^val(word,1,(table[g][h],))^val(word,1,(g,)))<<(2*i)
            return out
        def d2(word):
            out=0
            for i,(g,h,k) in enumerate(tuples[3]):
                v=act(g,val(word,2,(h,k)))^val(word,2,(table[g][h],k))^val(word,2,(g,table[h][k]))^val(word,2,(g,h))
                out|=v<<(2*i)
            return out
        d0=[sum((act(g,v)^v)<<(2*i) for i,g in enumerate(non)) for v in (1,2)]
        a=[d1(1<<i) for i in range(2*len(non))];b=[d2(1<<i) for i in range(2*len(tuples[2]))]
        assert all(d1(v)==0 for v in d0) and all(d2(v)==0 for v in a)
        rank0,rank1,rank2=r.rank(d0),r.rank(a),r.rank(b)
        h1=2*len(non)-rank1-rank0;h2=2*len(tuples[2])-rank2-rank1
        assert h1==h2==0
        records.append({'elements':[list(G[g]) for g in H],'order':len(H),'cochain_dimensions':[2,2*len(non),2*len(tuples[2]),2*len(tuples[3])],
            'differential_ranks':[rank0,rank1,rank2],'H0_dimension':2-rank0,'H1_dimension':h1,'H2_dimension':h2})
    def algebra_mul(a,b):
        v=0
        for g in range(6):
            if a>>g&1:
                for h in range(6):
                    if b>>h&1:v^=1<<table[g][h]
        return v
    tau=G.index((1,2,0));tau2=table[tau][tau];e0=(1<<unit)|(1<<tau)|(1<<tau2);e1=(1<<tau)|(1<<tau2)
    assert algebra_mul(e0,e0)==e0 and algebra_mul(e1,e1)==e1 and algebra_mul(e0,e1)==0 and e0^e1==1
    assert all(algebra_mul(e,g)==algebra_mul(g,e) for e in (e0,e1) for g in range(64))
    def matrix_image(word):
        cols=[]
        for v in (1,2):
            y=0
            for g in range(6):
                if word>>g&1:y^=act(g,v)
            cols.append(y)
        return cols[0]|cols[1]<<2
    def matrix_act(m,v):return (m&3 if v&1 else 0)^((m>>2)&3 if v&2 else 0)
    def matrix_mul(a,b):return matrix_act(a,b&3)|(matrix_act(a,(b>>2)&3)<<2)
    ideal0={algebra_mul(e0,w) for w in range(64)};ideal1={algebra_mul(e1,w) for w in range(64)}
    assert len(ideal0)==4 and len(ideal1)==16 and {matrix_image(w) for w in ideal0}=={0}
    assert {matrix_image(w) for w in ideal1}==set(range(16))
    assert all(matrix_image(algebra_mul(a,b))==matrix_mul(matrix_image(a),matrix_image(b)) for a,b in product(ideal1,repeat=2))
    end=[m for m in range(16) if all(matrix_act(m,act(g,v))==act(g,matrix_act(m,v)) for g,v in product(range(6),range(4)))]
    assert end==[0,9]
    return {'schema':'rank-jump.standard-s3-class-module.v1','status':'PASS','group':[list(g) for g in G],
        'standard_action_columns':[[act(g,1),act(g,2)] for g in range(6)],'subgroups':records,
        'C3_trace_idempotent_mask':e0,'standard_idempotent_mask':e1,'trace_block_dimension':2,'standard_block_dimension':4,
        'standard_block_is_M2_F2':True,'standard_endomorphism_ring_size':len(end),
        'bindings':{str(Path(__file__).relative_to(r.ROOT)):r.digest(Path(__file__).read_bytes())},
        'boundary':'Finite module and normalized cochain arithmetic. The class-field identification and specialization applications require the mathematical argument in the note; no number-field class group is calculated.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS all six subgroup H1/H2 groups vanish; standard block M2(F2)')
