#!/usr/bin/env sage-python
"""Fixed generic parity sample outside a maximum-class span; no point search."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from sage.all import ZZ,matrix,pari
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
ROOT=Path(__file__).resolve().parents[2]
def run(family_name,output):
    source=ROOT/'artifacts/generated-results/elliptic-curves/r17_exact_maximum_parity_classes_v1.json'
    family=next(x for x in json.loads(source.read_text())['families'] if x['family']==family_name);piv={}
    def reduce(v):
        for k in sorted(piv,reverse=True):
            if v>>k&1:v^=piv[k]
        return v
    for row in family['classes']:
        v=reduce(row['mask'])
        if v:piv[v.bit_length()-1]=v
    if len(piv)>=17:raise ArithmeticError('no outside-span classes')
    masks=[];i=0
    while len(masks)<128:
        m=int.from_bytes(hashlib.sha256((family_name+'/outside-maximum-span/'+str(i)).encode()).digest(),'big')&((1<<17)-1);i+=1
        if reduce(m) and m not in masks:masks.append(m)
    g=2*matrix(ZZ,family['gram']);u=matrix(ZZ,pari(g).qflllgram()).transpose();assert abs(u.det())==1;inv=u.inverse();a=IntegerExactParity((u*g*u.transpose()).rows());b=ExactParity((u*g*u.transpose()).rows());rows=[]
    for mask in masks:
        word=matrix(ZZ,1,17,[(mask>>j)&1 for j in range(17)]);res=tuple(int(x)%2 for x in (word*inv).row(0));starts,_=a.babai(np.asarray([res],dtype=np.int64));start=tuple(map(int,starts[0]));proof=a.solve(res,start,2000000);assert proof==b.solve(res,start,2000000)
        rows.append({'mask':mask,'norm':proof['norm'],'proof':proof,'reduced_seed':start})
    largest=max(r['norm'] for r in rows);chosen=sorted((r for r in rows if r['norm']==largest),key=lambda r:r['mask'])[:16]
    result={'schema':'outside-maximum-seed-masks.v1','families':[{'family':family_name,'classes':[{'mask':r['mask']} for r in chosen]}],'generic_doubled_norm':largest,'maximum_span_rank':len(piv),'maximum_echelon_basis':list(piv.values()),'sample':rows,'gram':[list(map(int,row)) for row in g.rows()],'LLL':[list(map(int,row)) for row in u.rows()],'bindings':{str(source.relative_to(ROOT)):hashlib.sha256(source.read_bytes()).hexdigest()},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'claim_boundary':'Fixed generic parity coverage only; no point-existence or rank prediction.'}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print('PREPARED_OUTSIDE_MAXIMUM',family_name,len(chosen),largest,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',required=True);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.family,a.output.resolve())
