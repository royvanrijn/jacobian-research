#!/usr/bin/env sage-python
"""Bounded exact polynomial-factor splitting of the two new rank25 discriminants."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import sys
from sage.all import ZZ,QQ,PolynomialRing,EllipticCurve,prime_range,prod,gcd
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as r17
import compact_mw16_specialization as mw16
from research_runtime.store import checkpoint

def run(kind,directory):
    spec=r17 if kind=='r17' else mw16
    certificate=ROOT/'artifacts/generated-results/elliptic-curves'/('compact_atlas_new_curves_v1.json' if kind=='r17' else 'prospective_mw16_results_v1.json')
    data=cert.read(certificate)
    row=next(r for r in data['curves'] if r.get('rank_lower_bound',r['rank_certificate']['rank_lower_bound'])==25)
    family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family' if kind=='r17' else 'fibration_id']==row['family'])
    model=tuple(QQ(q) for q in row['curve']);E=EllipticCurve(model);discriminant=E.discriminant()
    if discriminant.denominator()!=1:raise ArithmeticError('expected integral discriminant')
    n=abs(ZZ(discriminant));t=QQ(row['parameter']);R=PolynomialRing(QQ,'t')
    A=R([QQ(q) for q in family['A_coefficients_low_to_high']]);B=R([QQ(q) for q in family['B_coefficients_low_to_high']]);D=-16*(4*A**3+27*B**2)
    u=QQ(row['family_to_curve_scale_u'])
    if D(t)*t.denominator()**24/u**12!=discriminant:raise ArithmeticError('specialized discriminant transport failed')
    result={'kind':kind,'family':row['family'],'parameter':row['parameter'],'curve':row['curve'],'discriminant':str(discriminant),
        'certificate_sha256':cert.hashed(certificate),'atlas_sha256':cert.hashed(spec.ATLAS),'source_sha256':cert.hashed(Path(__file__).resolve()),
        'status':'POLYNOMIAL_FACTORIZATION_RUNNING','polynomial_degree':int(D.degree()),'polynomial_coefficients':list(map(str,D.list()))}
    checkpoint(directory/'result.json',result)
    factors=D.factor();assert factors.prod()==D
    result.update(polynomial_unit=str(factors.unit()),polynomial_factors=[{'coefficients':list(map(str,f.list())),'exponent':int(e),'degree':int(f.degree())} for f,e in factors],status='EXACT_POLYNOMIAL_FACTORIZATION')
    print('DISCRIMINANT FACTOR DEGREES',kind,[(int(f.degree()),int(e)) for f,e in factors],flush=True);checkpoint(directory/'result.json',result)
    remaining=n;small=[]
    for p in prime_range(2,10001):
        e=remaining.valuation(p)
        if e:small.append((ZZ(p),int(e)));remaining//=p**e
    pieces=[remaining] if remaining>1 else []
    for f,e in factors:
        v=abs(ZZ(f(t).numerator()));new=[]
        for q in pieces:
            h=gcd(q,v)
            if 1<h<q:new.extend([h,q//h])
            else:new.append(q)
        pieces=new
    assert prod(p**e for p,e in small)*prod(pieces)==n
    result.update(small_prime_factors=[[str(p),e] for p,e in small],remaining_pieces=list(map(str,pieces)),remaining_piece_bits=[int(q.nbits()) for q in pieces]);checkpoint(directory/'result.json',result)
    print('RESIDUAL PIECE BITS',kind,result['remaining_piece_bits'],flush=True)
    all_factors=list(small);unresolved=[]
    for q in pieces:
        if q.is_prime(proof=True):all_factors.append((q,1))
        elif q.nbits()<=192:all_factors.extend((p,int(e)) for p,e in q.factor(proof=True))
        else:unresolved.append(q)
    if unresolved:
        result.update(status='COMPOSITE_PIECE_EXCEEDS_192_BIT_GATE',unresolved_pieces=list(map(str,unresolved)));checkpoint(directory/'result.json',result);return
    collected={}
    for p,e in all_factors:collected[p]=collected.get(p,0)+e
    assert prod(p**e for p,e in collected.items())==n
    conductor=ZZ(1);local=[]
    for p,e in sorted(collected.items()):
        d=E.local_data(p,proof=True);v=int(d.conductor_valuation());conductor*=p**v
        local.append({'prime':str(p),'discriminant_valuation':e,'conductor_valuation':v,'kodaira_symbol':str(d.kodaira_symbol())})
    result.update(status='EXACT_CONDUCTOR',factorization=[[str(p),e] for p,e in sorted(collected.items())],local_data=local,conductor=str(conductor));checkpoint(directory/'result.json',result)
    print('EXACT CONDUCTOR',kind,conductor,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('kind',choices=['r17','mw16']);p.add_argument('--directory',type=Path,required=True);a=p.parse_args();run(a.kind,a.directory)
