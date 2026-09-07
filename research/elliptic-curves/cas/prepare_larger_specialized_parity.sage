#!/usr/bin/env sage-python
"""Compute the declared enlarged parity sample in the existing frozen metric."""
import argparse,sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import matrix,ZZ
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import audit_larger_specialized_parity as control
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
geometry=SourceFileLoader('larger_parity_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
def main(index):
    p=control.protocol();case=p['rows'][index];old=cert.read(ROOT/case['maps']);out=control.D/case['id']/'sample.json'
    if out.exists():raise FileExistsError('preserve one enlarged numerical sample')
    g=matrix(ZZ,old['rounded_gram']);u=matrix(ZZ,old['change_of_basis']);inverse=u.inverse();rank=case['rank']
    if abs(u.det())!=1 or u*g*u.transpose()!=matrix(ZZ,old['reduced_gram']):raise ArithmeticError('frozen unimodular metric differs')
    oracle=geometry.CosetOracle(old['reduced_gram']);sample=list(old['sample'])
    data={'status':'RUNNING','protocol_sha256':cert.hashed(control.D/'protocol.json'),'case':case,'sample':sample};checkpoint(out,data)
    for i,mask in enumerate(control.masks(p,case)[2048:],start=2048):
        residue=matrix(ZZ,1,rank,[(mask>>j)&1 for j in range(rank)]);target=[int(v)%2 for v in (residue*inverse).row(0)]
        norm,rep,error=oracle.solve(target);word=list(map(int,(matrix(ZZ,1,rank,rep)*u).row(0)))
        if any((word[j]-(mask>>j))%2 for j in range(rank)) or sum(word[j]*g[j,k]*word[k] for j in range(rank) for k in range(rank))!=norm:raise ArithmeticError('exact parity/norm differs')
        sample.append({'parity':mask,'representative':word,'metric_norm':norm,'cvp_error':error,'reduced_representative':list(rep)})
        if (i+1)%8192==0:checkpoint(out,data);print('LARGER PARITY',case['id'],i+1,flush=True)
    data['centres']=sorted(sample,key=lambda r:(-r['metric_norm'],r['parity']))[:49];data['status']='COMPLETE_FIXED_SAMPLE';checkpoint(out,data)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);main(p.parse_args().index)
