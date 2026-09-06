#!/usr/bin/env python3
"""Explicit constant-cubic carrier and frozen matched-pair exclusion."""
import argparse
from pathlib import Path
from math import isqrt
import retrospective as r
import fixed_field_transfer_geometry as geometry
import fresh_governing_panel as panel

PROTOCOL=Path(__file__).with_name('CUBIC_ISOMORPHISM_CARRIER_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_cubic_isomorphism_carrier_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_cubic_isomorphism_carrier_v1.json'


def export():
    data=r.read(panel.INPUT);manifest=r.read(panel.MANIFEST)
    used={p[k] for p in manifest['pairs'] for k in ('high','low')}
    models=[{'token':x['token'],'model':x['model']} for x in data['cases'] if x['token'] in used]
    pairs=[{k:p[k] for k in ('high','low')} for p in manifest['pairs']]
    files=(Path(__file__),PROTOCOL,panel.INPUT,panel.MANIFEST,geometry.INPUT)
    r.write_new(INPUT,{'schema':'rank-jump.cubic-isomorphism-carrier-inputs.v1','models':models,'pairs':pairs,
        'families':r.read(geometry.INPUT)['families'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}})


def compute():
    from sage.all import QQ,PolynomialRing,matrix
    P=PolynomialRing(QQ,names=('a','b','v','w'));a,b,v,w=P.gens();R=PolynomialRing(P,'z');z=R.gen()
    f=z**3+a*z+b;q=v*z+w*(z**2+2*a/3)
    multiplication=matrix(P,3,3,lambda i,j:((q*z**j)%f)[i])
    AA=a*v**2+3*b*v*w-a**2*w**2/3
    BB=b*v**3-2*a**2*v**2*w/3-a*b*v*w**2-(2*a**3/27+b**2)*w**3
    char=multiplication.charpoly();assert list(char)==[BB,AA,0,1]
    basis=matrix(P,3,3,lambda i,j:((q**j)%f)[i]);index=v**3+a*v*w**2+b*w**3
    assert basis.det()==index and -4*AA**3-27*BB**2==(-4*a**3-27*b**2)*index**2
    inp=r.read(INPUT);families=[];G=PolynomialRing(QQ,'t')
    for family in inp['families']:
        A=G(family['A']);B=G(family['B']);t=QQ(family['irreducibility_witness_parameter']);a0,b0=A(t),B(t)
        assert -4*a0**3-27*b0**2
        families.append({'family':family['family'],'fixed_field_parameter':str(t),
            'fixed_cubic_ascending':list(map(str,[b0,a0,0,1])),
            'A_equation_v2_vw_w2':list(map(str,[a0,3*b0,-a0*a0/3])),
            'B_equation_v3_v2w_vw2_w3':list(map(str,[b0,-2*a0*a0/3,-a0*b0,-2*a0**3/27-b0*b0])),
            'index_v3_v2w_vw2_w3':list(map(str,[1,0,a0,b0])),
            'known_identity_lift':{'t':str(t),'v':'1','w':'0'}})
    models={x['token']:x['model'] for x in inp['models']};pairs=[]
    for pair in inp['pairs']:
        discriminants=[]
        for token in (pair['high'],pair['low']):
            short,_=r.short(models[token],[]);aa,bb=map(r.F,short[3:]);discriminants.append(-4*aa**3-27*bb**2)
        ratio=discriminants[0]/discriminants[1];n,d=ratio.numerator,ratio.denominator
        sn,sd=isqrt(abs(n)),isqrt(d);square=n>0 and sn*sn==n and sd*sd==d
        pairs.append({**pair,'cubic_discriminants':list(map(str,discriminants)),'discriminant_ratio':str(ratio),
            'absolute_numerator_floor_sqrt':str(sn),'denominator_floor_sqrt':str(sd),
            'same_quadratic_resolvent':square,'cubic_fields_isomorphic':False if not square else 'UNKNOWN'})
    files=(Path(__file__),PROTOCOL,INPUT,geometry.OUTPUT)
    return {'schema':'rank-jump.cubic-isomorphism-carrier.v1','status':'PASS','families':families,'pairs':pairs,
        'universal_trace_zero_generator':'v*eta+w*(eta^2+2*a/3)',
        'universal_A':'a*v^2+3*b*v*w-a^2*w^2/3',
        'universal_B':'b*v^3-2*a^2*v^2*w/3-a*b*v*w^2-(2*a^3/27+b^2)*w^3',
        'universal_basis_index':'v^3+a*v*w^2+b*w^3',
        'symbolic_characteristic_polynomial_identity':True,'symbolic_discriminant_identity':True,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Exact isomorphism-carrier equations and paired nonisomorphism obstructions. No rational-point enumeration on the carriers and no rank predictor.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);args=p.parse_args()
    if args.mode=='export':export()
    else:
        result=compute()
        if args.mode=='build':r.write_new(OUTPUT,result)
        else:assert result==r.read(OUTPUT)
        print('PASS',len(result['families']),'carriers',sum(not x['same_quadratic_resolvent'] for x in result['pairs']),'nonisomorphic pairs',flush=True)
