#!/usr/bin/env python3
"""Exact quadratic genus characters with S splitting; not elliptic incidence."""
import argparse
from fractions import Fraction as Q
from math import prod
from pathlib import Path
import retrospective as r
import fresh_governing_panel as panel
import fresh_retained_factors as supplement
import bounded_gain_reference as reference

PROTOCOL=Path(__file__).with_name('RESOLVENT_GENUS_COMPONENT_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_resolvent_genus_component_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_resolvent_genus_component_v1.json'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    old={x['token']:x for x in r.read(panel.OUTPUT)['rows']}
    for x in r.read(supplement.OUTPUT)['rows']:
        if x['factor']['status']=='PASS':old[x['token']]=x
    ref=r.read(reference.OUTPUT);old[ref['token']]={'token':ref['token'],**ref['stages']}
    rows=[]
    for token,x in old.items():
        f,loc=x['factor'],x['local']
        if f['status']!='PASS' or not loc.get('S_finite'):
            rows.append({'token':token,'status':'UNKNOWN','reason':'Retained complete factorization/local set unavailable; no retry'});continue
        rows.append({'token':token,'status':'READY','cubic_ascending':f['integral_cubic_ascending'],
            'elliptic_discriminant':f['elliptic_discriminant'],'discriminant_factors':f['factors'],
            'S_finite':loc['S_finite'],'cubic_irreducibility_prime':loc['galois']['irreducibility_prime']})
    assert len(rows)==17 and sum(x['status']=='READY' for x in rows)==12
    r.write_new(INPUT,{'schema':'rank-jump.resolvent-genus-component-inputs.v1','rows':rows,
        'bindings':bindings([Path(__file__),PROTOCOL,panel.OUTPUT,supplement.OUTPUT,reference.OUTPUT]),
        'whitelist':'Equation coefficients, factorization, equation-defined S and modular cubic irreducibility witness only.'})


def signature(a,p):
    a=Q(a)
    if p=='infinity':return int(a<0)
    n,d=a.numerator,a.denominator;e=0
    while n%p==0:n//=p;e+=1
    while d%p==0:d//=p;e-=1
    if p==2:
        u=n*pow(d,-1,8)%8
        # Unit generators -1 and 5, plus valuation generator 2.
        units={1:0,7:1,5:2,3:3}
        return (e%2)|(units[u]<<1)
    u=n*pow(d,-1,p)%p;leg=pow(u,(p-1)//2,p);assert leg in (1,p-1)
    return e%2 | (int(leg==p-1)<<1)


def kernel(rows,n):
    # Explicit elimination on constraints, preserving all free coefficient slots.
    a=[x for x in rows if x];pivots=[];k=0
    for j in range(n):
        i=next((i for i in range(k,len(a)) if a[i]>>j&1),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k]
        for i in range(len(a)):
            if i!=k and a[i]>>j&1:a[i]^=a[k]
        pivots.append(j);k+=1
    basis=[]
    for j in range(n):
        if j in pivots:continue
        v=1<<j
        for i,p in enumerate(pivots):
            if a[i]>>j&1:v|=1<<p
        assert all((v&x).bit_count()%2==0 for x in rows);basis.append(v)
    return basis


def genus(D,primes,S):
    atoms=[p if p%4==1 else -p for p in primes if p!=2]
    odd=prod(atoms);two=D//odd;assert D%odd==0 and two in (1,-4,8,-8)
    if two!=1:atoms.insert(0,two)
    assert prod(atoms)==D;count=len(atoms);constraints=[];places=[]
    for p in [*S,'infinity']:
        width=1 if p=='infinity' else 3 if p==2 else 2
        ds=signature(D,p);columns=[signature(a,p) for a in atoms]
        functionals=[u for u in range(1,1<<width) if (u&ds).bit_count()%2==0]
        rs=[r.pack((u&v).bit_count()%2 for v in columns) for u in functionals]
        constraints.extend(rs);places.append({'place':p,'local_dimension':width,'D_signature':ds,
            'prime_discriminant_signatures':columns,'annihilator_functionals':functionals,'constraint_rows':rs})
    k=kernel(constraints,count);trivial=(1<<count)-1
    assert all((trivial&v).bit_count()%2==0 for v in constraints)
    chosen=[];span=[trivial]
    for v in k:
        if r.rank(span+[v])>len(span):chosen.append(v);span.append(v)
    assert len(chosen)==len(k)-1
    signs=[r.pack(a<0 for a in atoms)] if D>0 else []
    ordinary=count-r.rank(signs)-1
    return {'quadratic_discriminant':str(D),'prime_discriminants':list(map(str,atoms)),
        'finite_unramified_genus_dimension':count-1,'ordinary_genus_dimension':ordinary,
        'places':places,'full_kernel_basis_masks':k,'trivial_total_product_mask':trivial,
        'S_genus_basis_masks':chosen,'S_genus_dimension':len(chosen),
        'S_genus_generators':[str(prod(a for i,a in enumerate(atoms) if v>>i&1)) for v in chosen],
        'contribution_to_standard_S3_component':0}


def compute():
    data=r.read(INPUT);rows=[]
    for x in data['rows']:
        if x['status']!='READY':rows.append(x);continue
        c=list(map(int,x['cubic_ascending']));assert len(c)==4 and c[2]==0 and c[3]==1
        delta=-4*c[1]**3-27*c[0]**2;assert 16*delta==int(x['elliptic_discriminant'])
        factors=x['discriminant_factors'];assert prod(p**e for p,e in factors)==abs(16*delta)
        d=(-1 if delta<0 else 1)*prod(p for p,e in factors if e%2)
        assert d!=1;D=d if d%4==1 else 4*d
        primes=sorted(p for p,e in factors if D%p==0);assert set(primes)<=set(x['S_finite'])
        rows.append({'token':x['token'],'status':'PASS',**genus(D,primes,x['S_finite'])})
    controls=[{'token':'control-'+str(D),'status':'PASS',**genus(D,ps,ps)}
        for D,ps in [(-56,[2,7]),(-84,[2,3,7]),(12,[2,3]),(136,[2,17])]]
    return {'schema':'rank-jump.resolvent-genus-component.v1','rows':rows,'small_controls':controls,
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT]),
        'boundary':'Exact quadratic resolvent genus component with S splitting; its S3 standard projection is zero by the theorem in the note. Additional cubic strict classes and their solubility are not computed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);args=p.parse_args()
    if args.mode=='export':export()
    else:
        result=compute()
        if args.mode=='build':r.write_new(OUTPUT,result)
        else:assert result==r.read(OUTPUT)
        print([(x['token'],x.get('ordinary_genus_dimension'),x.get('S_genus_dimension')) for x in result['rows']])
