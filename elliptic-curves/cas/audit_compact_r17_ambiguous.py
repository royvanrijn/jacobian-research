#!/usr/bin/env python3
"""Bounded mod-3/mod-5 audit of every retained adaptive ambiguity, with replay.

An escape from an unsaturated subgroup is not automatically another rational
rank direction. Report the full finite column rank, with a torsion witness.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')];sys.set_int_max_str_digits(0)
import mod_l_reduction_independence as ml
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint


def signature(model,points,prime,modulus):
    # On a good integral short model a nonintegral affine point specializes
    # to O. The legacy routine rejects such denominators, so compute its
    # exact signature on the integral sublist and insert identity columns.
    included=[i for i,p in enumerate(points) if all(c.denominator%prime for c in p)]
    sig=ml.mod_l_reduction_signature(model,[points[i] for i in included],prime,modulus)
    rows=[]
    for source in sig.rows:
        row=[0]*len(points)
        for i,v in zip(included,source):row[i]=v
        rows.append(tuple(row))
    return ml.ModLReductionSignature(modulus,prime,sig.group_order,sig.multiple_subgroup_order,
        sig.quotient_dimension,tuple(rows))


def pivots(rows,modulus):
    a=[list(r) for r in rows];out=[];j=0
    for col in range(len(a[0]) if a else 0):
        k=next((k for k in range(j,len(a)) if a[k][col]%modulus),None)
        if k is None:continue
        a[j],a[k]=a[k],a[j];u=pow(a[j][col],-1,modulus);a[j]=[v*u%modulus for v in a[j]]
        for k in range(j+1,len(a)):
            if a[k][col]:
                q=a[k][col];a[k]=[(v-q*w)%modulus for v,w in zip(a[k],a[j])]
        out.append(col);j+=1
        if j==len(a):break
    return out


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in
        (Path(__file__).resolve(),CAS/'mod_l_reduction_independence.py',CAS/'mod2_reduction_independence.py',
         CAS/'certify_compact_r17_candidates.py',CAS/'elliptic_candidate_record.py',cert.MODEL,cert.SECTIONS)}


def build(input_path,output):
    if output.exists():raise FileExistsError('use a new immutable output')
    d=cert.read(input_path);model=tuple(map(F,d['curve']))
    basis=[tuple(map(F,p)) for p in d['final_state']['state']['reductions']['points']]
    cert.checked_rank(model,basis);cert.family_check(d['parameter'],model,basis)
    seen={(x,abs(y)) for x,y in basis};points=list(basis)
    for row in d['final_state']['state']['observations']:
        if row['status']!='AMBIGUOUS_FINITE_REDUCTIONS':continue
        x,y=map(F,row['point']);key=(x,abs(y))
        if key not in seen:seen.add(key);points.append((x,y))
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('off-curve observation')
    result={'schema':'elliptic-curves.compact-r17-ambiguity-audit.v1','status':'RUNNING','parameter':d['parameter'],
        'curve':d['curve'],'points':[list(map(str,p)) for p in points],'known_rank_lower_bound':len(basis),
        'input_path':str(input_path.resolve().relative_to(ROOT)),'input_sha256':cert.hashed(input_path),
        'sources':sources(),'protocol':{'moduli':[3,5],'prime_bound':997,'new_point_search':False,
            'workers':1,'wall_limit_seconds':600,'rss_limit_bytes':1610612736},'audits':[],
        'claim_boundary':'Finite column rank and no rational ell-torsion prove that lower bound. Non-escape is inconclusive. Never add marginal escape dimension to an unverified baseline.'}
    checkpoint(output,result)
    for ell in result['protocol']['moduli']:
        tp=ml.find_no_rational_l_torsion_prime(model,modulus=ell)
        audit={'modulus':ell,'no_rational_ell_torsion_prime':tp,'signatures':[],
            'finite_column_rank':0,'independent_column_indices':[],'processed_primes':[],'status':'RUNNING'}
        result['audits'].append(audit);selected=[]
        for p in ml._primes_up_to(997):
            if p in (2,ell):continue
            try:sig=signature(model,points,p,ell)
            except ValueError:continue
            new_rows=selected+list(sig.rows);columns=pivots(new_rows,ell)
            if len(columns)>audit['finite_column_rank']:
                selected=new_rows;audit['signatures'].append(asdict(sig));audit['finite_column_rank']=len(columns)
                audit['independent_column_indices']=columns
                print('AMBIGUITY',d['parameter'],'ell',ell,'prime',p,'rank',len(columns),flush=True)
            audit['processed_primes'].append(p);checkpoint(output,result)
        audit['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)
    result['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)


def check(path):
    d=cert.read(path)
    for p,h in d['sources'].items():
        if cert.hashed(ROOT/p)!=h:raise ArithmeticError('audit source changed: '+p)
    model=tuple(map(F,d['curve']));points=[tuple(map(F,p)) for p in d['points']]
    cert.family_check(d['parameter'],model,points)
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('off-curve point')
    for audit in d['audits']:
        ell=audit['modulus'];rows=[]
        if not ml.no_rational_l_torsion_reduction_certificate(model,audit['no_rational_ell_torsion_prime'],ell):
            raise ArithmeticError('torsion witness failed')
        for old in audit['signatures']:
            sig=signature(model,points,old['prime'],ell)
            if json.dumps(asdict(sig),sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite quotient changed')
            rows.extend(sig.rows)
        columns=pivots(rows,ell)
        if columns!=audit['independent_column_indices'] or len(columns)!=audit['finite_column_rank']:
            raise ArithmeticError('independent columns changed')
        print('REPLAYED AMBIGUITY',d['parameter'],'ell',ell,'rank >=',len(columns),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--input',type=Path);g.add_argument('--check',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
    if a.check:check(a.check)
    else:build(a.input,a.output)


if __name__=='__main__':main()
