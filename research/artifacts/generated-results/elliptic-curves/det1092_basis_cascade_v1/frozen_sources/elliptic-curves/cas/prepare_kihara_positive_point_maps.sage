#!/usr/bin/env sage-python
"""Exact coherent seed intake and prospective geometry for twelve fixed fibres."""
import argparse
import sys
from pathlib import Path
from dataclasses import asdict
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,EllipticCurve,matrix,pari
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import kihara_positive_point_pilot as control
import certify_compact_r17_candidates as cert
from mestre_parent_adapter import specialize
from research_runtime.store import checkpoint,digest
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from memory_rank_certificate import checked_rank

def maps(index):
    control.configure(index);p=control.protocol();out=control.D/'maps.json';rank=control.ROW['initial_rank']
    if out.exists():raise FileExistsError('preserve prospective maps')
    geometry=SourceFileLoader('mestre_geometry',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
    mapper=SourceFileLoader('mestre_mapper',str(CAS/'factor_free_pari_mapping.sage')).load_module()
    seed=cert.read(control.SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    gram,asym=geometry.canonical_height_gram(model,points);g=matrix(ZZ,geometry.rounded_gram(gram,1000000))
    u=matrix(ZZ,pari(g).qflllgram()).transpose();inverse=u.inverse()
    if abs(u.det())!=1:raise ArithmeticError('unimodular metric change required')
    reduced=u*g*u.transpose();oracle=geometry.CosetOracle(reduced.rows());sample=[]
    data={'status':'RUNNING_SAMPLE','protocol_hash':digest(p),'metric_gram':[[str(v) for v in r] for r in gram],
      'maximum_gram_asymmetry':str(asym),'rounded_gram':[list(map(int,r)) for r in g.rows()],
      'change_of_basis':[list(map(int,r)) for r in u.rows()],'reduced_gram':[list(map(int,r)) for r in reduced.rows()],
      'sample':sample,'rows':[]};checkpoint(out,data)
    for i,mask in enumerate(control.masks(p)):
        residue=matrix(ZZ,1,rank,[(mask>>j)&1 for j in range(rank)]);target=[int(v)%2 for v in (residue*inverse).row(0)]
        norm,rep,error=oracle.solve(target);word=list(map(int,(matrix(ZZ,1,rank,rep)*u).row(0)))
        if any((word[j]-(mask>>j))%2 for j in range(rank)) or sum(word[j]*g[j,k]*word[k] for j in range(rank) for k in range(rank))!=norm:raise ArithmeticError('exact parity/norm transport differs')
        sample.append({'parity':mask,'representative':word,'metric_norm':norm,'cvp_error':error,'reduced_representative':list(rep)})
        if (i+1)%511==0:checkpoint(out,data);print('PARITIES',i+1,'in fixed sample',flush=True)
    data['centres']=control.centres(sample,[list(map(int,row)) for row in g.rows()]);data['status']='RUNNING_MAPS';checkpoint(out,data)
    pari.allocatemem(256000000,silent=True)
    for c in data['centres']:
        data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data)
    print(control.ROW['id'],'FROZEN49 MAPS',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,required=True);a=parser.parse_args();maps(a.index)
