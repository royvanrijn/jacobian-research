#!/usr/bin/env python3
"""Exact mod3/mod5 lower bounds on an immutable complete retained point cloud."""
import argparse,json
from dataclasses import asdict
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_compact_r17_ambiguous as finite
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ml=finite.ml


def insert(basis,row,ell):
    row=[int(v)%ell for v in row]
    for j,old in sorted(basis.items()):
        if row[j]:
            a=row[j];row=[(v-a*w)%ell for v,w in zip(row,old)]
    j=next((j for j,v in enumerate(row) if v),None)
    if j is None:return False
    inv=pow(row[j],-1,ell);basis[j]=[v*inv%ell for v in row];return True


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(finite.__file__).resolve(),Path(ml.__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')}


def build(path,output):
    if output.exists():raise FileExistsError('preserve mod-ell audit')
    data=cert.read(path);model=tuple(map(cert.F,data['curve']));points=[tuple(map(cert.F,p)) for p in data['points']]
    if data['status']!='COMPLETE_DECLARED_FINITE_AUDIT' or any(model[:3]) or any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('complete short-model point cloud required')
    result={'schema':'elliptic-curves.retained-cloud-modl.v1','status':'RUNNING','sources':sources(),'input_path':str(path.resolve().relative_to(ROOT)),'input_sha256':cert.hashed(path),'family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'points':data['points'],'original_lower_bound':data['rank_lower_bound'],'prime_bound':997,'moduli':[3,5],'audits':[],'claim_boundary':'Exact full finite column rank with no rational ell-torsion supplies a lower bound. A mod2 escape or a failure to increase finite rank is not an upper rank bound. No point search or new curve selection.'};checkpoint(output,result)
    for ell in result['moduli']:
        tp=ml.find_no_rational_l_torsion_prime(model,modulus=ell);basis={};audit={'modulus':ell,'no_rational_ell_torsion_prime':tp,'status':'RUNNING','finite_column_rank':0,'independent_column_indices':[],'signatures':[],'processed_primes':[]};result['audits'].append(audit)
        for p in ml._primes_up_to(997):
            if p in (2,ell):continue
            try:sig=finite.signature(model,points,p,ell)
            except ValueError:continue
            before=len(basis)
            for row in sig.rows:insert(basis,row,ell)
            if len(basis)>before:
                audit['signatures'].append(asdict(sig));audit.update(finite_column_rank=len(basis),independent_column_indices=sorted(basis));print('RETAINED CLOUD MODL',data['family'],ell,p,len(basis),flush=True)
            audit['processed_primes'].append(p);checkpoint(output,result)
        audit['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)
    result['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)


def check(path):
    data=cert.read(path)
    if data['sources']!=sources():raise ArithmeticError('mod-ell audit sources changed')
    model=tuple(map(cert.F,data['curve']));points=[tuple(map(cert.F,p)) for p in data['points']]
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('point is off curve')
    for audit in data['audits']:
        ell=audit['modulus'];rows=[]
        if not ml.no_rational_l_torsion_reduction_certificate(model,audit['no_rational_ell_torsion_prime'],ell):raise ArithmeticError('torsion witness failed')
        for old in audit['signatures']:
            actual=finite.signature(model,points,old['prime'],ell)
            if json.dumps(asdict(actual),sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite signature differs')
            rows.extend(actual.rows)
        indices=finite.pivots(rows,ell)
        if indices!=audit['independent_column_indices'] or len(indices)!=audit['finite_column_rank']:raise ArithmeticError('independent column certificate differs')
        print('REPLAYED CLOUD MODL',data['family'],ell,'rank >=',len(indices),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path);p.add_argument('--output',type=Path);p.add_argument('--check',type=Path);a=p.parse_args();check(a.check) if a.check else build(a.input,a.output)
