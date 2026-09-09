#!/usr/bin/env sage-python
"""One exact affine ellipsoid: all original sections of alternate degree<=2.

Only generic Gram/pencil word and a previously saved unimodular reduction.
No target parameter, point coordinate, new subgroup-state or orbit census.
Each action <=25s, recursion<=100000 nodes, at most4096 source words.
"""
import argparse,hashlib,json,signal,time
from fractions import Fraction as Rat
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_low_degree_source_obstruction_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FRAME=ART/'det1092_genus1_picard_image_v1/frame.json'
REDUCTION=ART/'det1092_two_fibration_action_v1/degree-minima-v2-plus.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert read(path)==data,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    save('geometry-protocol.json',dict(classification='complete bounded generic low-degree source geometry',
        bound=16,alternate_degree_cap=2,node_cap=100000,source_cap=4096,seconds=25,
        definition='All v in M17 with degree D.S_v=v^tGv-w^tGv<=2; w is the fixed norm8 pencil word. Equivalently ||2v-w||^2<=16.',
        excluded_inputs=['first seed','carrier label','control address','catalogue rank','search output'],
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,FRAME,REDUCTION,Path(__file__)]}))
def enumerate_sources():
    protocol=read(OUT/'geometry-protocol.json')
    for name,digest in protocol['inputs'].items():assert sha(ROOT/name)==digest
    G=matrix(ZZ,read(PARENT)['generic_height_gram']);frame=read(FRAME)
    D=vector(QQ,frame['D']);assert list(D[:2])==[2,4]
    w=vector(ZZ,D[2:]);assert w*G*w==8
    U=matrix(ZZ,read(REDUCTION)['LLL_columns']);assert abs(U.det())==1
    inv=U.inverse();residue=[int(v%2) for v in inv*w]
    reduced=U.transpose()*G*U;n=17
    A=[[Rat(int(x)) for x in row] for row in reduced.rows()]
    L=[[Rat(int(i==j)) for j in range(n)] for i in range(n)];diag=[]
    for j in range(n):
        pivot=A[j][j]-sum(L[j][k]**2*diag[k] for k in range(j));assert pivot>0;diag.append(pivot)
        for i in range(j+1,n):L[i][j]=(A[i][j]-sum(L[i][k]*L[j][k]*diag[k] for k in range(j)))/pivot
    assert all(A[i][j]==sum(L[i][k]*diag[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
    z=[0]*n;nodes=0;found=[]
    def visit(j,used):
        nonlocal nodes
        nodes+=1
        if nodes>protocol['node_cap']:raise RuntimeError('UNKNOWN_NODE_CAP')
        if j<0:
            lifted=U*vector(ZZ,z);v=vector(ZZ,(lifted+w)/2)
            norm=lifted*G*lifted;degree=v*G*v-w*G*v
            assert norm==QQ(used.numerator)/used.denominator and norm==4*degree+8 and 0<=degree<=2
            found.append(dict(word=list(map(int,v)),degree=int(degree),centered_norm=int(norm)))
            if len(found)>protocol['source_cap']:raise RuntimeError('UNKNOWN_SOURCE_CAP')
            return
        shift=sum((L[k][j]*z[k] for k in range(j+1,n)),Rat(0))
        allowance=(Rat(16)-used)/diag[j]
        if allowance<0:return
        b=shift.denominator;a=shift.numerator
        radius=isqrt((allowance*b*b).numerator//(allowance*b*b).denominator)
        lower=-((radius+a)//b);upper=(radius-a)//b
        lower+=(residue[j]-lower)%2
        for value in range(lower,upper+1,2):
            z[j]=value;visit(j-1,used+diag[j]*(Rat(value)+shift)**2)
    try:visit(16,Rat(0))
    except BaseException as exc:
        save('geometry-failure.json',dict(status='UNKNOWN_INCOMPLETE',reason=str(exc),nodes=nodes,completed_words=found));raise
    found=sorted(found,key=lambda r:(r['degree'],r['word']))
    assert len({tuple(r['word']) for r in found})==len(found)
    counts={str(d):sum(r['degree']==d for r in found) for d in range(3)}
    save('sources.json',dict(status='PASS_COMPLETE_LOW_DEGREE_SOURCE_ELLIPSOID',
        counts=counts,nodes=nodes,word=list(map(int,w)),LLL_columns=[list(map(int,r)) for r in U.rows()],
        residue=residue,ldl_diagonal=list(map(str,diag)),sources=found,
        geometry_protocol_sha256=sha(OUT/'geometry-protocol.json'),
        boundary='This proves source completeness only. No first-seed carrier incidence has been evaluated.'))
    print('PASS_COMPLETE_LOW_DEGREE_SOURCE_ELLIPSOID',counts,'nodes',nodes,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','enumerate']);a=p.parse_args()
    signal.alarm(25);begun=time.monotonic()
    if a.action=='freeze':freeze()
    else:enumerate_sources()
    print('seconds',round(time.monotonic()-begun,3),flush=True)
