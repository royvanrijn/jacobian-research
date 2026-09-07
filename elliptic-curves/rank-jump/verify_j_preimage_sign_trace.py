#!/usr/bin/env python3
"""CAS-free irreducibility, nonsquare and complete retained-panel replay."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import retrospective as r
import j_preimage_sign_trace as source
import second_generic_presentation as first
import verify_second_generic_presentation as rational
import verify_surface_discriminant_modular as finite

OUTPUT=r.OUT/'rank_jump_j_preimage_sign_trace_verification_v1.json'
COMPARISON=r.OUT/'rank_jump_j_preimage_sign_trace_comparison_v1.json'


def compute():
    data=r.read(source.INPUT);out=r.read(source.OUTPUT);old=r.read(first.OUTPUT);checks=0;cache={}
    for obj in (data,out,r.read(rational.OUTPUT)):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha,p
    targets={x['token']:x for x in data['cases']};families={x['family']:x for x in data['families']}
    inputs={(x['token'],x['family']):x for x in data['rows']};earlier={(x['token'],x['family']):x for x in old['rows']}
    assert len(out['rows'])==len(inputs)==112;rows=[]
    for row in out['rows']:
        assert row['status']=='PASS';key=row['token'],row['family'];inp=inputs[key];prior=earlier[key]
        q=list(map(Q,inp['polynomial_ascending']));n=len(q)-1;assert n==row['degree'] and n in (23,24) and q[-1]==1
        # Verify the exported residual polynomial against the original exact equation.
        product=q
        for x in prior['preimages']:
            assert x['parameter']!='infinity'
            product=rational.mul(product,rational.power([-Q(x['parameter']),Q(1)],x['multiplicity']))
        original=list(map(Q,prior['j_equation_ascending']))
        assert original==rational.scale(product,original[-1])
        possible=set(range(1,n));certificates={}
        for c in row['irreducibility_certificates']:
            p=c['prime'];assert rational.prime(p) and 3<=p<=251;certificates[p]=c
            qp=rational.mod(q,p);assert len(qp)==n+1 and qp[-1]==1
            parts=c['factors_ascending'];assert len({tuple(h) for h in parts})==len(parts)
            product=[1];degrees=[]
            for h in parts:
                assert h[-1]==1 and len(h)>1 and all(0<=x<p for x in h)
                k=p,tuple(h)
                if k not in cache:cache[k]=finite.irreducible(h,p)
                assert cache[k];checks+=1;product=finite.mul(product,h,p);degrees.append(len(h)-1)
            assert product==qp and degrees==c['factor_degrees']
            sums={0}
            for d in degrees:sums|={x+d for x in list(sums)}
            possible&=sums-{0,n};assert sorted(possible)==c['remaining_possible_factor_degrees']
        assert not possible and row['remaining_possible_factor_degrees']==[]
        w=row['nonsquare_isomorphism_class_witness'];p=w['prime'];c=certificates[p]
        h=c['factors_ascending'][w['factor_index']];assert h==w['residue_factor_ascending']
        t=targets[row['token']];f=families[row['family']];a,b=Q(t['a']),Q(t['b'])
        assert a*b and 4*a**3+27*b**2
        ap,bp=rational.mod([a,b],p);assert ap*bp%p
        A=rational.mod(list(map(Q,f['A'])),p);B=rational.mod(list(map(Q,f['B'])),p)
        num=finite.mod([bp*x%p for x in A],h,p);den=finite.mod([ap*x%p for x in B],h,p)
        assert num!=[0] and den!=[0]
        # Check the residue by multiplication, independently of worker inversion.
        residue=w['scaling_class_residue_ascending']
        assert len(residue)<len(h) and all(0<=x<p for x in residue)
        assert finite.mod(finite.mul(residue,den,p),h,p)==num
        d=len(h)-1;assert d==w['finite_field_degree'] and w['quadratic_character']==-1
        assert finite.power(residue,(p**d-1)//2,h,p)==[p-1]
        assert row['rational_trace_span_dimension']==0
        rows.append({'token':row['token'],'family':row['family'],'preimage_degree':n,
            'isomorphism_field_degree':2*n,'rational_irreducibility':'PASS',
            'nonsquare_prime':p,'nonsquare_residue_degree':d,'rational_trace_span_dimension':0})
    return {'schema':'rank-jump.j-preimage-sign-trace-verification.v1','status':'PASS','rows':rows,
        'finite_irreducible_factor_checks':checks,'distinct_cached_factors':len(cache),
        'bindings':first.bindings([Path(__file__),source.INPUT,source.OUTPUT,first.OUTPUT,rational.OUTPUT,
            Path(rational.__file__),Path(finite.__file__)]),
        'boundary':'Independent rational polynomial products, finite Rabin tests, subset-degree exclusions, and Euler quadratic characters. The trace-vanishing theorem is proved in the note. No field class group, elliptic points or rank labels enter this replay.'}


def compare(result):
    label_path=r.OUT/'rank_jump_fresh_governing_panel_manifest_v1.json';labels=r.read(label_path)
    capacity_path=r.OUT/'rank_jump_root_curve_capacity_comparison_v1.json';capacity=r.read(capacity_path)
    own={x['token']:x for x in r.read(first.INPUT)['cases']}
    rational_rows=r.read(rational.OUTPUT)['rows'];rows=[]
    for x in labels['rows']:
        token=x['token'];found=[y for y in result['rows'] if y['token']==token]
        assert len(found)==7
        c=next(y for y in capacity['rows'] if y['token']==token)
        rows.append({'token':token,'id':x['id'],'family':x['family'],'parameter':x['parameter'],
            'retained_rank_lower_bound':x['retained_rank_lower_bound'],'retained_gain_lower_bound':x['observed_quotient_rank'],
            'rational_presentations':sum(y['rational_preimages'] for y in rational_rows if y['token']==token),
            'rational_base_equivalence_classes':1,'nonrational_preimage_fields':7,
            'nonrational_isomorphism_field_degrees':[y['isomorphism_field_degree'] for y in found],
            'nonrational_source_section_rational_span_dimension':0,
            'constructor_dimensions_outside_global_pool':0,
            'retained_rank_forces_dimensions_outside_global_pool':c['rational_fibre_dimensions_outside_pool_lower_bound']})
    return {'schema':'rank-jump.j-preimage-sign-trace-comparison.v1','rows':rows,
        'bindings':first.bindings([Path(__file__),OUTPUT,label_path,capacity_path,rational.OUTPUT]),
        'boundary':'Retained labels joined after the equation-only exclusion. Every tested source section construction is confined to the inherited global pool. Censored low labels do not give exact low rank.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result);r.write_new(COMPARISON,compare(result))
    else:assert result==r.read(OUTPUT);assert compare(result)==r.read(COMPARISON)
    print('PASS',len(result['rows']),'preimage fields;',result['finite_irreducible_factor_checks'],'finite factor checks; all trace spans zero')
