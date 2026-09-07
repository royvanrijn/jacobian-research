#!/usr/bin/env python3
"""An equation-defined four-division class, with bounded local exclusion."""
import argparse
from collections import deque
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as panel
import fresh_governing_completion as completion
import fresh_strict_boundary_coordinates as derivative
import fresh_retained_factors as supplement
import fresh_governing_octics as octics

PROTOCOL=Path(__file__).with_name('FOUR_DIVISION_CLASS_CREATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_four_division_class_creation_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-four-division-class-creation-v1'
NEW=('case-08','case-10','case-11')


def bindings():
    paths=[Path(__file__),PROTOCOL,panel.INPUT,Path(panel.__file__),completion.OUTPUT,
           derivative.OUTPUT,supplement.OUTPUT,octics.OUTPUT,panel.LOCAL,Path(r.__file__)]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def sign_action(v,cycle=False):
    z,a,b,c=(v>>i&1 for i in range(4))
    return r.pack((z,c,z^a,z^b) if cycle else (z,z^a,c,b))


def beta_map(v):
    z,a,b,c=(v>>i&1 for i in range(4))
    return r.pack((z^a^b,a^c,z^b^c))


def module_certificate():
    # Coordinates: signs on i, sqrt(e1-e2), sqrt(e1-e3), sqrt(e2-e3).
    # Pullback by (12) and (123); either convention generates the same S3 action.
    S=[sign_action(v) for v in range(16)]
    T=[sign_action(v,True) for v in range(16)]
    phi=[beta_map(v) for v in range(16)]
    assert all(S[S[v]]==v and T[T[T[v]]]==v and S[T[S[v]]]==T[T[v]] for v in range(16))
    fixed=[v for v in range(16) if T[v]==v]
    standard=sorted({v^T[v] for v in range(16)})
    assert fixed==[v for v in range(16) if phi[v]==0]
    assert len(fixed)==len(standard)==4 and set(fixed)&set(standard)=={0}
    assert {phi[v] for v in standard}=={0,3,5,6}
    for v in range(16):
        x,y,z=(phi[v]>>i&1 for i in range(3))
        assert phi[S[v]]==r.pack((y,x,z))
        assert phi[T[v]]==r.pack((y,z,x))
    return {'sign_action_transposition':S,'sign_action_three_cycle':T,'derivative_map':phi,
            'three_cycle_fixed_subspace':fixed,'standard_subspace':standard,
            'derivative_kernel_dimension':2,'standard_multiplicity':1,
            'universal_four_division_Kummer_capacity':1}


def full_group(n):
    I=(1,0,0,1);gens=[(0,1,n-1,0),(1,1,0,1),(n-1,0,0,1)]
    if n>=8:gens.append((3,0,0,1))
    def mul(a,b):
        return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%n for i in range(2) for j in range(2))
    expr={I:(0,0)};pending=deque([I]);constraints=[]
    while pending:
        a=pending.popleft()
        for j,b in enumerate(gens):
            c=mul(a,b)
            v=tuple(expr[a][i]^((1<<(2*j)) if a[2*i]%2 else 0)^((1<<(2*j+1)) if a[2*i+1]%2 else 0) for i in range(2))
            if c not in expr:expr[c]=v;pending.append(c)
            else:constraints.extend(x^y for x,y in zip(expr[c],v))
    order=6*(n//2)**4;assert len(expr)==order
    rank=r.rank(constraints);zdim=2*len(gens)-rank;assert zdim==3
    solutions=[s for s in range(1<<(2*len(gens))) if all((s&v).bit_count()%2==0 for v in constraints)]
    assert len(solutions)==8
    return {'modulus':n,'group_order':order,'generators':list(map(list,gens)),
        'constraint_rank':rank,'Z1_dimension':zdim,'B1_dimension':2,'H1_dimension':1,
        'cocycle_generator_assignments':solutions,'application_to_panel':'UNKNOWN: actual image not computed'}


def local_worker(token):
    from sage.all import pari
    sys.path.insert(0,str(panel.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    saved=next(x for x in r.read(completion.OUTPUT)['rows'] if x['token']==token)
    f,pts,_=panel.model_data(token)
    primes=saved['tested_bad_primes']
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod('z',pari(f))
    beta=-pari(f.discriminant())*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==f.discriminant()**4
    trials=[]
    for p in primes:
        chars=LocalSquareclasses(nf,p)
        old=next(x for x in saved['local'] if x['place']==p)
        sigs=[list(chars.signature(pari(x)-theta)) for x,y in pts]
        assert sigs==old['signatures']
        dim=chars.point_kummer_dimension;rank=r.rank(list(map(r.pack,sigs)))
        ds=list(chars.signature(beta));outside=rank==dim and r.rank(list(map(r.pack,sigs))+[r.pack(ds)])>dim
        trials.append({'place':p,'point_dimension':dim,'generic_signatures':sigs,
                       'derivative_signature':ds,'outside_full_point_image':outside})
        if outside:
            return {'status':'PASS','token':token,'witness':trials[-1],'trials':trials,
                    'local_nf_primes':primes,'bindings':bindings()}
    return {'status':'UNKNOWN','token':token,'trials':trials,'bindings':bindings()}


def build():
    prior={x['token']:x for x in r.read(derivative.OUTPUT)['rows']}
    extra={x['token']:x['boundary'] for x in r.read(supplement.OUTPUT)['rows']}
    galois={x['token']:x['galois'] for x in r.read(octics.OUTPUT)['rows']}
    rows=[]
    for row in r.read(panel.INPUT)['cases']:
        token=row['token'];f,_,_=panel.model_data(token);delta=f.discriminant()
        assert galois[token]['galois_group']=='S3'
        beta=(-delta*f.derivative())%f
        assert f.resultant(beta)==delta**4
        result={'token':token,'cubic_ascending':list(map(str,f.list())),
                'discriminant':str(delta),'derivative_ascending':list(map(str,beta.list())),
                'four_division_class_dimension_upper_bound':1}
        if delta>0:
            witness={'place':'infinity','derivative_signs':[1,0,1],
                     'full_point_image':[[0,0,0],[0,1,1]],'source':'exact positive discriminant and ordered-root signs'}
        elif token in NEW:
            work=r.read(WORK/f'{token}.json');assert work['bindings']==bindings()
            witness=work.get('witness') if work['status']=='PASS' else None
        else:
            old=prior[token] if prior[token]['status']=='PASS' else extra[token]
            assert old['status']=='PASS' and old['witness'] is not None
            p=old['witness']['place'];loc=next(x for x in old['local'] if x['place']==p)
            assert loc['outside_full_point_image'] and loc['generic_local_rank']==loc['point_dimension']
            witness={'place':p,'retained_derivative_signature':loc['derivative_signature'],
                     'point_dimension':loc['point_dimension'],'source':'retained exact derivative exclusion'}
        result.update(status='PASS' if witness else 'UNKNOWN',witness=witness,
            four_division_class_dimension=1 if witness else 'UNKNOWN',
            four_division_Selmer_dimension=0 if witness else 'UNKNOWN',
            four_division_strict_dimension=0 if witness else 'UNKNOWN')
        rows.append(result)
    return {'schema':'rank-jump.four-division-class-creation.v1',
        'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'PARTIAL',
        'bindings':bindings(),'module':module_certificate(),
        'abstract_full_group_controls':[full_group(n) for n in r.read(PROTOCOL)['limits']['finite_group_moduli']],
        'rows':rows,'boundary':'Incidence capacity and local admissibility only. No new exceptional classes, points, CT values, or actual mod8/mod16 images computed. No whole-tower claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='worker':
        from sage.all import pari
        assert args.token in NEW
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        r.write_new(WORK/f'{args.token}.json',local_worker(args.token))
    else:
        WORK.mkdir(parents=True,exist_ok=True)
        for token in NEW:
            path=WORK/f'{token}.json'
            if not path.exists():
                with (WORK/f'{token}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,__file__,'worker','--token',token],stdout=log,stderr=log,timeout=30)
                        reason=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason:r.write_new(path,{'token':token,'status':'UNKNOWN','reason':reason,'bindings':bindings()})
            print(token,r.read(path)['status'],flush=True)
        result=build();r.write_new(OUTPUT,result);print(result['status'],flush=True)
