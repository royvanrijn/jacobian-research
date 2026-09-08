#!/usr/bin/env python3
"""Bounded mod-3/mod-5 rank audit of retained prospective MW16 point clouds."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import audit_compact_r17_ambiguous as finite
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2];ml=finite.ml


def sources():
    paths=(Path(__file__).resolve(),Path(finite.__file__).resolve(),Path(cert.__file__).resolve(),Path(spec.__file__).resolve(),
           spec.ATLAS,ROOT/'elliptic-curves/cas/mod_l_reduction_independence.py',
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def family_check(data,model,points):
    family=next(f for f in cert.read(spec.ATLAS)['families'] if f['fibration_id']==data['family'])
    original,generic=spec.specialize(family,data['parameter'])
    if not cert.isomorphic(original,model):raise ArithmeticError('family model mismatch')
    a,b=cert.weierstrass_invariants(original),cert.weierstrass_invariants(model)
    u=cert.square_root(a['c6']*b['c4']/(a['c4']*b['c6']))
    if u is None:raise ArithmeticError('missing rational scale')
    for sign in (1,-1):
        v=sign*u
        if tuple((x/v**2,y/v**3) for x,y in generic)==tuple(points[:16]):return str(v)
    raise ArithmeticError('generic MW16 transport failed')


def build(path,output):
    if output.exists():raise FileExistsError('use a new immutable output')
    data=cert.read(path)
    if data['status'] not in ('COMPLETE_DECLARED_PILOT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY'):
        raise ArithmeticError('use a terminal immutable point-search result')
    model=tuple(map(cert.F,data['curve']));basis=[tuple(map(cert.F,p)) for p in data['final_state']['state']['reductions']['points']]
    cert.checked_rank(model,basis);family_check(data,model,basis)
    seen={(x,abs(y)) for x,y in basis};points=list(basis)
    for row in data['final_state']['state']['observations']:
        if row['status']!='AMBIGUOUS_FINITE_REDUCTIONS':continue
        x,y=map(cert.F,row['point']);key=(x,abs(y))
        if key not in seen:seen.add(key);points.append((x,y))
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('off-curve observation')
    result={'schema':'elliptic-curves.prospective-mw16-ambiguity-audit.v1','family':data['family'],'parameter':data['parameter'],
        'status':'RUNNING','curve':data['curve'],'points':[list(map(str,p)) for p in points],
        'known_rank_lower_bound':len(basis),'input_path':str(path.resolve().relative_to(ROOT)),'input_sha256':cert.hashed(path),
        'sources':sources(),'protocol':{'moduli':[3,5],'prime_bound':997,'new_point_search':False,'wall_limit_seconds':600,'rss_limit_bytes':1610612736},
        'audits':[],'claim_boundary':'Full finite column rank with no rational ell-torsion proves a lower bound. Failure to find another independent column is inconclusive; a quotient escape alone is not a rank gain.'}
    checkpoint(output,result)
    for ell in (3,5):
        tp=ml.find_no_rational_l_torsion_prime(model,modulus=ell)
        audit={'modulus':ell,'no_rational_ell_torsion_prime':tp,'signatures':[],'finite_column_rank':0,
            'independent_column_indices':[],'processed_primes':[],'status':'RUNNING'}
        result['audits'].append(audit);selected=[]
        for p in ml._primes_up_to(997):
            if p in (2,ell):continue
            try:sig=finite.signature(model,points,p,ell)
            except ValueError:continue
            rows=selected+list(sig.rows);columns=finite.pivots(rows,ell)
            if len(columns)>audit['finite_column_rank']:
                selected=rows;audit.update(finite_column_rank=len(columns),independent_column_indices=columns)
                audit['signatures'].append(asdict(sig));print('MW16 AMBIGUITY',data['family'],data['parameter'],'ell',ell,'rank',len(columns),flush=True)
            audit['processed_primes'].append(p);checkpoint(output,result)
        audit['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)
    result['status']='COMPLETE_BOUNDED_QUOTIENT_AUDIT';checkpoint(output,result)


def check(path):
    data=cert.read(path)
    if data['sources']!=sources():raise ArithmeticError('audit source binding changed')
    model=tuple(map(cert.F,data['curve']));points=[tuple(map(cert.F,p)) for p in data['points']]
    family_check(data,model,points)
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('off-curve point')
    for audit in data['audits']:
        ell=audit['modulus'];rows=[]
        if not ml.no_rational_l_torsion_reduction_certificate(model,audit['no_rational_ell_torsion_prime'],ell):raise ArithmeticError('torsion witness failed')
        for old in audit['signatures']:
            sig=finite.signature(model,points,old['prime'],ell)
            if json.dumps(asdict(sig),sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite signature changed')
            rows.extend(sig.rows)
        columns=finite.pivots(rows,ell)
        if columns!=audit['independent_column_indices'] or len(columns)!=audit['finite_column_rank']:raise ArithmeticError('finite column rank changed')
        print('REPLAYED MW16 AMBIGUITY',data['family'],data['parameter'],'ell',ell,'rank >=',len(columns),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--input',type=Path);g.add_argument('--check',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.input,a.output)
