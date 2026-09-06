#!/usr/bin/env python3
"""Independent finite-field characters and provenance for paired moment bounds."""
import argparse
from fractions import Fraction as F
from math import prod
from pathlib import Path
from sage.all import GF,PolynomialRing
import retrospective as r

HERE=Path(__file__).resolve().parent
SOURCE=r.OUT/'rank_jump_paired_character_moments_v1.json'
CASES=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
TRACES=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
GLOBAL=r.OUT/'rank_jump_global_pair_solubility_v1.json'
AD=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
CT=r.OUT/'rank_jump_labelled_carrier_ct_v2.json'
FIBRES=r.OUT/'rank_jump_solubility_first_v1.json'
OUTPUT=r.OUT/'rank_jump_paired_character_moments_verification_v1.json'


def verify():
    source=r.read(SOURCE);cases=r.read(CASES);base=r.read(TRACES)
    for data in (source,cases,base):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    roster={c['label']:c['form'] for case in cases['cases'] for c in case['covers']}
    fp=GF(131);Rp=PolynomialRing(fp,'z');z=Rp.gen();k=GF(131**2,'alpha',modulus=z*z-2);alpha=k.gen()
    chars={label:[] for label in roster}
    for n,a,b,mult,trace,unused in base['fibre_trace_ledger']:
        field=fp if n==1 else k;t=field(a) if n==1 else k(a)+b*alpha
        for label,(c0,c1,c2) in roster.items():
            value=field(c2) if a==-1 else c0+c1*t+c2*t*t
            chars[label].append(0 if not value else 1 if value.is_square() else -1)
    R=PolynomialRing(fp,'t');A,B,_=[R(c) for c in base['geometry']['modular_coefficients']]
    delta=-16*(4*A**3+27*B**2);fs=[f for f,e in delta.factor()];assert prod(fs)==delta.monic()
    local={label:[] for label in roster};original=[]
    for f in fs:
        field=fp if f.degree()==1 else GF(131**f.degree(),'beta',modulus=f)
        root=-f[0] if f.degree()==1 else field.gen();ring=PolynomialRing(field,'s')
        base_c6=ring(864*B)(root);original.append(-1 if base_c6.is_square() else 1)
        for label,form in roster.items():
            value=ring(form)(root);assert value
            local[label].append(1 if value.is_square() else -1)
    bounds=[]
    for row in source['rows']:
        labels=row['labels'];q=prod(R(roster[label]) for label in labels)
        assert q.degree()==2*len(labels) and q.is_squarefree() and q.gcd(delta)==1
        tr=[-sum(mult*trace*prod(chars[label][j] for label in labels) for j,(n,a,b,mult,trace,unused) in enumerate(base['fibre_trace_ledger']) if n==degree) for degree in(1,2)]
        assert row['status']=='PASS' and tr==row['Frobenius_traces']
        N=20+2*int(q.degree());s1=F(tr[0],131);s2=F(tr[1],131**2);c=F(row['quadratic_center'])
        bnd=((N+s2)/2-2*c*s1+N*c*c)/(1-c)**2
        assert str(bnd)==row['exact_moment_bound'];m=bnd.numerator//bnd.denominator
        W=prod(original)*prod(prod(local[label]) for label in labels)
        assert W==row['global_root_number']
        while (-1)**m!=W:m-=1
        assert m==row['arithmetic_generic_rank_upper_bound'];bounds.append(m)
    # Global solubility entries are inherited from their actual class certificates.
    old_global=r.read(GLOBAL);old_ad=r.read(AD);ct=r.read(CT)
    actual={row['id']:row['descent'] for row in old_global['rows']}
    for new,key in [('positive_FG','observed_positive'),('cross_group_FD','cross_group')]:
        row=next(x for x in source['paired_carrier_comparisons'] if x['id']==new)
        old=actual[key];assert old['global_carrier_solubility_proved'] and old['rank_lower_bound']==old['rank_upper_bound']==row['carrier_Jacobian_exact_rank']
        assert row['carrier_global_solubility']=='YES'
    row=next(x for x in source['paired_carrier_comparisons'] if x['id']=='obstructed_AD')
    old=old_ad['rows'][0]['descent'];assert old['rank_lower_bound']==old['rank_upper_bound']==row['carrier_Jacobian_exact_rank']==2
    assert ct['carrier_global_solubility']=='NO' and ct['carrier_nonzero_Sha_class'] and row['carrier_global_solubility']=='NO'
    fibres=r.read(FIBRES);outside=[]
    for row in source['rows']:
        if len(row['labels'])!=1 or row['system']=='obstructed_ABCD':continue
        fibre=next(x for x in fibres['rows'] if x['source_id']==row['system'])
        assert row['labels'][0] in [x['label'] for x in fibre['nonzero_square_hits']]
        outside.append({'source_id':row['system'],'label':row['labels'][0],
                        'generic_twist_upper_bound':row['arithmetic_generic_rank_upper_bound'],
                        'witness_directions_outside_entire_specialized_pullback_span_at_least':max(0,fibre['observed_quotient_rank']-row['arithmetic_generic_rank_upper_bound'])})
    return {'schema':'rank-jump.paired-character-moments-verification.v1','status':'PASS','characters_verified':len(bounds),
            'individual_cover_specialization_boundaries':outside,
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (SOURCE,CASES,TRACES,GLOBAL,AD,CT,FIBRES,Path(__file__),HERE/'retrospective.py')},
            'boundary':'Reuses the independently verified base fibre counts. Reweights them independently in actual finite fields; the global carrier solubility comparison is checked against the preceding exact class certificates.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=verify()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS',result['characters_verified'],'character bounds and three inherited carrier-solubility outcomes')
