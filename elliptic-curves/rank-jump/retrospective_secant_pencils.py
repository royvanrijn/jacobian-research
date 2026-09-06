#!/usr/bin/env python3
"""One generic-point secant pencil for each of six retained high/low controls."""
import argparse
from math import isqrt
from pathlib import Path
import subprocess
import retrospective as r
import bad_prime_support as bad

PROTOCOL=Path(__file__).with_name('RETROSPECTIVE_SECANT_PENCIL_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_retrospective_secant_pencil_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_retrospective_secant_pencils_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-retrospective-secant-pencils-v1'


def rational_square(x):
    x=r.F(x)
    return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator


def export():
    rows=[]
    for i in r.read(PROTOCOL)['cases']:
        source=bad.cases()[i];model,points=r.short(source['model'],source['generic_points'][:2])
        profile,_,sigs=r.characterize(source)
        assert r.rank(sigs[:2])==2
        rows.append({'case_index':i,'id':source['id'],'short_model':model,'generic_pair':points,
            'generic_pair_fingerprints':sigs[:2],'original_generic_rank':len(source['generic_points']),
            'original_known_independent_rank':profile['certified_independent_subgroup_rank_exact'],
            'original_observed_quotient_rank':profile['certified_independent_quotient_rank_exact']})
    r.write_new(INPUT,{'schema':'rank-jump.retrospective-secant-pencil-inputs.v1',
        'source_sha256':r.digest(r.INPUT.read_bytes()),'rows':rows})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT)}


def compute(index):
    from sage.all import QQ,PolynomialRing,EllipticCurve
    from sage.version import version
    row=next(x for x in r.read(INPUT)['rows'] if x['case_index']==index)
    A,B=map(QQ,row['short_model'][3:]);(a,p),(b,q)=[tuple(map(QQ,P)) for P in row['generic_pair']]
    assert p*p==a**3+A*a+B and q*q==b**3+A*b+B and a!=b
    s=(q-p)/(b-a);intercept=p-s*a
    assert s!=0
    x0=-intercept/s;C=x0**3+A*x0+B;assert C!=0
    third=s*s-a-b;third_y=s*third+intercept
    assert third not in [a,b] and third_y*third_y==third**3+A*third+B
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field();X=PolynomialRing(K,'x');x=X.gen()
    f=x**3+A*x+B;L=s*x+intercept;g=f-L*L
    assert g==(x-a)*(x-b)*(x-third)
    E=EllipticCurve(K,[0,(t-1)*s*s,0,A+2*(t-1)*s*intercept,B+(t-1)*intercept**2])
    delta=R(E.discriminant());c4=R(E.c4())
    passed=delta.degree()==3 and c4.degree()==2 and delta.gcd(delta.derivative())==1 and delta.gcd(c4)==1
    if not passed:return {'bindings':bindings(),'case_index':index,'status':'UNKNOWN','reason':'generic fibre configuration gate not met'}
    assert delta(0)!=0 and delta(1)==-16*(4*A**3+27*B**2)
    N=PolynomialRing(QQ,'n');n=N.gen();En=EllipticCurve(N.fraction_field(),[0,(n*n-1)*s*s,0,
        A+2*(n*n-1)*s*intercept,B+(n*n-1)*intercept**2])
    P=En(a,n*p);Q=En(b,n*q);T=En(third,n*third_y)
    assert P+Q+T==En(0)
    delta_n=N(En.discriminant());c4_n=N(En.c4())
    assert delta_n==delta(n*n) and delta_n.degree()==6 and c4_n.degree()==4
    assert delta_n.gcd(delta_n.derivative())==1 and delta_n.gcd(c4_n)==1
    # Fixed point at x0 lies over the constant quadratic field sqrt(C).
    assert (x0**3+(t-1)*s*s*x0*x0+(A+2*(t-1)*s*intercept)*x0+B+(t-1)*intercept**2)==C
    square=rational_square(str(C));inherited=int(square)
    return {'bindings':bindings(),'case_index':index,'id':row['id'],'status':'PASS','software':{'sage':version},
        'original_generic_rank':row['original_generic_rank'],'original_known_independent_rank':row['original_known_independent_rank'],
        'original_observed_quotient_rank':row['original_observed_quotient_rank'],
        'selected_original_generic_indices':[0,1],'secant_slope':str(s),'secant_intercept':str(intercept),
        'residual_roots':list(map(str,[a,b,third])),'residual_cubic_coefficients':list(map(str,g.list())),
        'fixed_section_x':str(x0),'fixed_section_y_squared':str(C),'fixed_section_is_rational':square,
        'parent_discriminant_coefficients':list(map(str,delta)),'parent_c4_coefficients':list(map(str,c4)),
        'parent_fibres':['I1']*3+['I3*'],'parent_geometric_generic_rank':1,'parent_arithmetic_generic_rank':inherited,
        'fixed_section_height_lower_bound':'1/4 (disjoint from O; D7 correction at most 7/4)',
        'base_fibres':['I1']*6+['I6'],'base_geometric_generic_rank':3,'base_arithmetic_generic_rank':inherited+2,
        'new_base_change_rational_directions':2,
        'original_generic_quotient_contribution_of_selected_pair':0,
        'anchor_t':1,'anchor_n':1,
        'boundary':'The new pencil has generic rank 0 or 1, not the original 16 or 17. The selected two directions were already original generic directions; this construction supplies no original exceptional quotient.'}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in r.read(PROTOCOL)['cases']:
        if check:
            assert compute(i)==next(x for x in r.read(OUTPUT)['rows'] if x['case_index']==i)
            print('PASS secant pencil replay',i,flush=True);continue
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--index',str(i),
                        '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:failure='30-second timeout'
                if failure and not path.exists():r.write_new(path,{'bindings':bindings(),'case_index':i,'status':'UNKNOWN','reason':failure})
        record=r.read(path);assert record['bindings']==bindings();rows.append(record)
        print('checkpoint',i,record['status'],record.get('base_arithmetic_generic_rank'),flush=True)
    if not check:r.write_new(OUTPUT,{'schema':'rank-jump.retrospective-secant-pencils.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check']);p.add_argument('--index',type=int)
    p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index))
    else:capture(args.mode=='check')
