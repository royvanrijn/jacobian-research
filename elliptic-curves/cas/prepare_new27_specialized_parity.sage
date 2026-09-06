#!/usr/bin/env sage-python
"""Finite specialized parity sample using only the own certified27-point seed."""
import sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import matrix,ZZ,pari
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import new27_specialized_parity_six as control
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
mapper=SourceFileLoader('specialized_mapper',str(CAS/'prepare_fresh_r17_pari_batch.sage')).load_module()

def main():
    p=control.protocol();out=control.D/'maps.json'
    if out.exists():raise FileExistsError('preserve specialized geometry')
    seed=cert.read(control.SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    gram,asym=mapper.geometry.canonical_height_gram(model,points);rounded=mapper.geometry.rounded_gram(gram,1000000);g=matrix(ZZ,rounded);u=matrix(ZZ,pari(g).qflllgram()).transpose();inverse=u.inverse()
    if abs(u.det())!=1:raise ArithmeticError('unimodular metric change required')
    reduced=u*g*u.transpose();oracle=mapper.geometry.CosetOracle(reduced.rows());sample=[]
    data={'status':'RUNNING_SAMPLE','protocol_hash':digest(p),'metric_gram':[[str(v) for v in r] for r in gram],'maximum_gram_asymmetry':str(asym),'rounded_gram':[list(map(int,r)) for r in g.rows()],'change_of_basis':[list(map(int,r)) for r in u.rows()],'reduced_gram':[list(map(int,r)) for r in reduced.rows()],'sample':sample,'rows':[]};checkpoint(out,data)
    for i,mask in enumerate(control.masks(p)):
        residue=matrix(ZZ,1,27,[(mask>>j)&1 for j in range(27)]);target=[int(v)%2 for v in (residue*inverse).row(0)];norm,rep,error=oracle.solve(target);word=list(map(int,(matrix(ZZ,1,27,rep)*u).row(0)))
        if any((word[j]-(mask>>j))%2 for j in range(27)) or sum(word[j]*g[j,k]*word[k] for j in range(27) for k in range(27))!=norm:raise ArithmeticError('exact specialized parity/norm transport differs')
        sample.append({'parity':mask,'representative':word,'metric_norm':norm,'cvp_error':error,'reduced_representative':list(rep)})
        if (i+1)%127==0:checkpoint(out,data);print('SPECIALIZED SAMPLE',i+1,'of2048',flush=True)
    data['centres']=sorted(sample,key=lambda c:(-c['metric_norm'],c['parity']))[:49];data['status']='RUNNING_MAPS';checkpoint(out,data);mapper.pari.allocatemem(256000000,silent=True)
    for c in data['centres']:data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('FROZEN SPECIALIZED49 MAPS FROM2048 MASKS',flush=True)
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);a=p.parse_args();control.configure(a.index);main()
