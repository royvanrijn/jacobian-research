#!/usr/bin/env sage-python
"""Complete height-six triangle gate with fixed old section P6.

Exact rational LDL enumeration, at most five million nodes and120 seconds.
Each triangle either matches a certified MW14/MW15 pencil or has two
additional independent rational vertical curves, forcing MW rank at most13.
The conclusion is restricted to triangles O+P6+Q with height(Q)=6.
"""
import argparse
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import signal
from sage.all import QQ,ZZ,matrix,vector,block_diagonal_matrix

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results'
SOURCE=ART/'elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
MW14=ART/'elkies-k3-curve302-anchor6-mw14-protocol-v1.json'
MW15=ART/'elkies-k3-curve302-six-mw15-triangles-protocol-v1.json'
OUT=ART/'elkies-k3-curve302-anchor6-triangle-gate-v1.json'

def exact_shell(H,norm):
    n=H.nrows(); L=matrix.identity(QQ,n); diagonal=[]
    for i in range(n):
        diagonal.append(QQ(H[i,i])-sum(L[i,k]**2*diagonal[k] for k in range(i)))
        assert diagonal[-1]>0
        for j in range(i+1,n):
            L[j,i]=(H[j,i]-sum(L[j,k]*L[i,k]*diagonal[k] for k in range(i)))/diagonal[i]
    assert L*matrix.diagonal(QQ,diagonal)*L.transpose()==H
    values=[0]*n; answer=[]; nodes=0
    def visit(i,remaining):
        nonlocal nodes
        nodes+=1; assert nodes<=5000000
        if i<0:
            if remaining==0: answer.append(tuple(values))
            return
        center=-sum((L[j,i]*values[j] for j in range(i+1,n)),QQ(0))
        radius=isqrt(int((remaining/diagonal[i]).floor()))+1
        for v in range(int(center.floor())-radius,int(center.ceil())+radius+1):
            cost=diagonal[i]*(v-center)**2
            if cost<=remaining:
                values[i]=v; visit(i-1,remaining-cost)
    visit(n-1,QQ(norm))
    return sorted(set(answer)),nodes

def build():
    G=matrix(ZZ,json.loads(SOURCE.read_text())['sections']['height_gram'])
    N=block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
    P=vector(ZZ,[int(i==5) for i in range(17)])
    def sec(w): return vector(ZZ,[1,(w*G*w)//2]+list(w))
    def key(w): return min(tuple(w),tuple(P-w))
    O=sec(vector(ZZ,17)); six,nodes6=exact_shell(G,6); four,nodes4=exact_shell(G,4)
    assert len(six)==53290 and len(four)==2626
    candidates=[vector(ZZ,w) for w in six if vector(ZZ,w)*G*P==3]
    triangles=sorted({key(w) for w in candidates})
    assert len(candidates)==1760 and len(triangles)==880
    recognized={}
    for path,rank in [(MW14,14),(MW15,15)]:
        for m in json.loads(path.read_text())['models']:
            k=key(vector(ZZ,m['Q_word'])); assert k not in recognized
            assert m['P_word']==list(P) and m['generic_rank']==rank
            recognized[k]=dict(protocol=str(path.relative_to(ROOT)),index=m['index'],generic_rank=rank)
    assert len(recognized)==94
    four=[vector(ZZ,w) for w in four]; records=[]
    for word in triangles:
        Q=vector(ZZ,word); D=O+sec(P)+sec(Q)
        assert D*N*D==0 and all(D*N*v==0 for v in [O,sec(P),sec(Q)])
        linear=(P+Q)*G
        vertical=[r for r in four if linear*r==6]
        if word in recognized:
            expected=15-recognized[word]['generic_rank']
            assert len(vertical)==expected
            record=dict(Q_word=list(map(int,Q)),kind='HIGH_RANK_PENCIL',**recognized[word])
        else:
            assert len(vertical)>=2
            chosen=vertical[:2]
            assert all(r*G*r==4 and D*N*sec(r)==0 for r in chosen)
            C=matrix(ZZ,[O,sec(P),sec(Q)]+list(map(sec,chosen)))
            pivots=list(C.pivots()); assert len(pivots)==5
            determinant=C.matrix_from_columns(pivots).det(); assert determinant
            record=dict(Q_word=list(map(int,Q)),kind='MW_AT_MOST_13_IF_JACOBIAN',
                        vertical_words=[list(map(int,r)) for r in chosen],
                        independence_columns=pivots,independence_determinant=int(determinant))
        records.append(record)
    assert sum(r['kind']=='MW_AT_MOST_13_IF_JACOBIAN' for r in records)==786
    return dict(schema='curve302.anchor6-triangle-gate.v1',status='COMPLETE_FIXED_P6_HEIGHT_SIX_GATE',
                input_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SOURCE,MW14,MW15]},
                shell_counts=dict(norm4=len(four),norm6=len(six)),
                shell_enumeration_nodes=dict(norm4=nodes4,norm6=nodes6),
                norm6_words_sha256=sha256(json.dumps(six,separators=(',',':')).encode()).hexdigest(),
                triangle_count=880,low_rank_count=786,MW14_count=88,MW15_count=6,records=records,
                boundary='All height-six Q with pairing(P6,Q)=3 on the pinned old fibration, modulo Q -> P6-Q. This does not classify other anchors, arbitrary degree-three divisors, or all elliptic fibrations. Inverse exclusions for the94 high-rank pencils are separate certificates.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    signal.alarm(120); result=build()
    if args.check: assert result==json.loads(OUT.read_text())
    else:
        assert not OUT.exists(); OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],result['triangle_count'],result['low_rank_count'],flush=True)
