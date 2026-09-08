#!/usr/bin/env python3
"""Odd-prime lower bounds for complete torsion-compatible pilot clouds."""
import argparse
from dataclasses import asdict
import kihara_positive_point_pilot as batch
import certify_compact_r17_candidates as cert
import audit_compact_r17_ambiguous as finite
from audit_retained_cloud_modl import insert
from research_runtime.store import checkpoint
def main(index):
    batch.configure(index);p=batch.protocol();path=batch.D/'result.json';data=cert.read(path);assert data['status']=='COMPLETE_DECLARED_POINT_ATTEMPT' and len(data['charts'])==49
    out=batch.ART/('kihara_positive_'+batch.ROW['id'].replace('-','_')+'_cloud_v1.json');assert not out.exists()
    model=tuple(map(cert.F,data['curve']));points=[tuple(map(cert.F,P)) for P in data['retained_points']];assert all(cert.is_on_weierstrass_curve(model,P) for P in points)
    result={'schema':'kihara-positive-full-cloud.v1','status':'RUNNING','id':batch.ROW['id'],'curve':data['curve'],'points':data['retained_points'],'initial_rank':batch.ROW['initial_rank'],'audits':[],'source_path':str(path.relative_to(batch.ROOT)),'source_sha256':cert.hashed(path)};checkpoint(out,result)
    for ell in (3,5):
        tp=finite.ml.find_no_rational_l_torsion_prime(model,modulus=ell);basis={};sigs=[]
        for prime in finite.ml._primes_up_to(997):
            if prime in (2,ell):continue
            try:sig=finite.signature(model,points,prime,ell)
            except ValueError:continue
            before=len(basis)
            for row in sig.rows:insert(basis,row,ell)
            if len(basis)>before:sigs.append(asdict(sig))
        result['audits'].append({'modulus':ell,'rank_lower_bound':len(basis),'independent_indices':sorted(basis),'no_rational_ell_torsion_prime':tp,'signatures':sigs});checkpoint(out,result)
        print(batch.ROW['id'],'mod',ell,'rank',len(basis),'cloud',len(points),flush=True)
    result['rank_lower_bound']=max(a['rank_lower_bound'] for a in result['audits']);assert result['rank_lower_bound']>=result['initial_rank'];result['discovered_rank_gain_lower_bound']=result['rank_lower_bound']-result['initial_rank'];result['status']='PASS';checkpoint(out,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);a=p.parse_args();main(a.index)
