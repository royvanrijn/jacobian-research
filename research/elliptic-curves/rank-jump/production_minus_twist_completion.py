#!/usr/bin/env python3
"""Two bounded local root trees completing the frozen production twist audit."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import production_minus_twist as first
import scalar_cup

PROTOCOL=Path(__file__).with_name('PRODUCTION_MINUS_TWIST_COMPLETION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_production_minus_twist_completion_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_production_minus_twist_completion_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-production-minus-twist-completion-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,first.INPUT,first.OUTPUT,scalar_cup.OUTPUT)}


def compute(index,p):
    from sage.all import QQ, ZZ, PolynomialRing, pari, lcm
    from sage.version import version
    sys.path.insert(0,str(first.rem.bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    row=next(x for x in r.read(first.INPUT)['cases'] if x['case_index']==index)
    local=next(x for x in row['local'] if x['place']==p)
    assert not local['complete']
    primes=[x['place'] for x in row['local'] if x['place']!='infinity']
    R=PolynomialRing(QQ,'z');z=R.gen(); f=R(list(map(QQ,row['integral_cubic_ascending'])))
    nf=pari.nfinit([pari(f),primes]); theta=pari.Mod('z',pari(f)); chars=LocalSquareclasses(nf,p)
    twist=z**3+f[1]*z-f[0]; E=pari.ellinit([0,0,0,f[1],-f[0]])
    red=pari.elllocalred(E,p);u,shift=map(QQ,list(red[2])[:2])
    g=R(twist(u*u*z+shift)/u**6);h=g*lcm(c.denominator() for c in g)
    h/=ZZ(p)**min(c.valuation(p) for c in h if c)
    assert all(c.denominator()==1 for c in h)
    # Bind the local coordinate convention to every retained first-pass witness.
    for w in local['witnesses']:
        assert r.pack(chars.signature(pari(QQ(w['short_x']))+theta))==w['signature']
    basis=list(local['twist_basis']);witnesses=[];nodes=[ZZ(0)];modulus=ZZ(1);trials=0;trace=[]
    spec=r.read(PROTOCOL)['limits'];seen=set();complete=False
    for depth in range(1,spec['max_depth']+1):
        nextnodes=[]
        for center in nodes:
            for digit in range(p):
                t=center+digit*modulus
                if h(t)%(modulus*p)==0:nextnodes.append(t)
                if t in seen:continue
                seen.add(t);trials+=1
                if trials>spec['max_candidate_evaluations']:break
                x=u*u*t+shift;value=twist(x);square=first.local_square(value,p)
                if square is None or not square['square']:continue
                signature=r.pack(chars.signature(pari(x)+theta))
                if r.rank(basis+[signature])>len(basis):
                    basis.append(signature)
                    witnesses.append({'minimal_x':str(t),'short_x':str(x),'cubic_value':str(value),
                        'square_witness':square,'signature':signature,'depth':depth})
                if len(basis)==local['point_dimension']:complete=True;break
            if complete or trials>spec['max_candidate_evaluations']:break
        trace.append({'depth':depth,'root_nodes':len(nextnodes),'total_trials':trials,'basis_dimension':len(basis)})
        if complete or not nextnodes or len(nextnodes)>spec['max_root_nodes'] or trials>spec['max_candidate_evaluations']:break
        nodes=nextnodes;modulus*=p
    return {'bindings':bindings(),'case_index':index,'place':p,'software':{'sage':version,'pari':str(pari.version())},
        'complete':complete,'point_dimension':local['point_dimension'],'twist_basis':basis,
        'new_witnesses':witnesses,'root_polynomial_ascending':list(map(str,h)),
        'minimal_u':str(u),'minimal_shift':str(shift),'trace':trace}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i,p in r.read(PROTOCOL)['targets']:
        path=WORK/f'case-{i}-p-{p}.json'
        if not path.exists():
            with (WORK/f'case-{i}-p-{p}.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker',
                        '--index',str(i),'--prime',str(p),'--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,
                        timeout=r.read(PROTOCOL)['limits']['seconds_per_target'])
                    failure=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:failure='30-second timeout'
                if failure and not path.exists():r.write_new(path,{'bindings':bindings(),'case_index':i,'place':p,'complete':False,'reason':failure})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print('checkpoint',i,p,row['complete'],flush=True)
    r.write_new(INPUT,{'schema':'rank-jump.production-minus-twist-completion-inputs.v1','bindings':bindings(),'rows':rows})


def merged():
    data=r.read(INPUT);assert data['bindings']==bindings()
    cases=r.read(first.INPUT)['cases']
    for row in data['rows']:
        assert row['bindings']==bindings()
        if not row['complete']:continue
        case=next(x for x in cases if x['case_index']==row['case_index'])
        local=next(x for x in case['local'] if x['place']==row['place'])
        local['twist_basis']=row['twist_basis'];local['witnesses']+=row['new_witnesses'];local['complete']=True
        local['intersection_dimension']=len(first.intersection(local['original_basis'],local['twist_basis']))
    for case in cases:
        offset=0;product=[]
        for local in case['local']:
            product.extend(v<<offset for v in local['twist_basis']);offset+=local['width']
        case['twist_product_basis']=product;case['complete']=all(x['complete'] for x in case['local'])
    return cases


def build(check=False):
    rows=[]
    cup=r.read(scalar_cup.OUTPUT)
    for case in merged():
        if not case['complete']:
            rows.append({'case_index':case['case_index'],'status':'UNKNOWN'});continue
        boundary=case['joint_global_signatures']; original=case['original_product_basis']; twist=case['twist_product_basis']
        ell=case['relaxed_boundary_dimension'];n=case['witness_dimension'];k=case['known_strict_dimension']
        b=len(first.intersection(boundary,twist))
        c=len(first.intersection(first.intersection(boundary,original),twist))
        rows.append({'case_index':case['case_index'],'id':case['id'],'status':'PASS',
            'original_Selmer_constant':n,'twist_Selmer_constant':k+b,'common_Selmer_constant':k+c,
            'known_strict_dimension':k,'original_boundary_dimension':ell-1,'twist_boundary_dimension':b,
            'common_boundary_dimension':c,'dimension_change':b-ell+1,
            'local_intersection_codimension':sum(x['point_dimension']-x['intersection_dimension'] for x in case['local']),
            'changed_local_places':[x['place'] for x in case['local'] if x['point_dimension']>x['intersection_dimension']],
            'boundary':'All three Selmer constants have the same additive epsilon = dim Cl(O_K,S)/2 - k for this case; epsilon remains UNKNOWN.'})
    result={'schema':'rank-jump.production-minus-twist-completion.v1','bindings':bindings(),
            'input_sha256':r.digest(INPUT.read_bytes()),'rows':rows}
    if check:assert r.read(OUTPUT)==result;print('PASS completed production twist intersections')
    else:r.write_new(OUTPUT,result)
    for row in rows:print(row)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['capture','worker','build','check'])
    parser.add_argument('--index',type=int);parser.add_argument('--prime',type=int);parser.add_argument('--destination',type=Path)
    args=parser.parse_args()
    if args.mode=='capture':capture()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index,args.prime))
    else:build(args.mode=='check')
