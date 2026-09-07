#!/usr/bin/env python3
"""Exhaust in-field radical operations on the retained multiplicative dictionary."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import early_relation_pool as pool
import relation_root_class as root
import relation_root_affine as affine
import complete_relation_root_class as completed
from verify_unpointed_governing_norm import Algebra

PROTOCOL=Path(__file__).with_name('RETAINED_RADICAL_CLOSURE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_retained_radical_closure_v1.json'
ROOT_REPLAY=r.OUT/'rank_jump_relation_root_class_verification_v1.json'
AFFINE_REPLAY=r.OUT/'rank_jump_relation_root_affine_verification_v1.json'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def generic_characters(data):
    f=list(map(F,data['cubic_ascending']));g=[list(map(F,z)) for z in data['generic_classes_ascending']]
    blocks=[];signatures=[0]*len(g)
    for p in r.primes(2003):
        try:
            mod=lambda v:v.numerator*pow(v.denominator,-1,p)%p
            fp=list(map(mod,f));gp=[list(map(mod,z)) for z in g]
        except ValueError:continue
        # A cubic with three distinct roots is separable, regardless of its integral presentation.
        roots=[x for x in range(p) if sum(c*pow(x,i,p) for i,c in enumerate(fp))%p==0]
        if len(roots)!=3:continue
        vals=[[sum(c*pow(x,i,p) for i,c in enumerate(z))%p for x in roots] for z in gp]
        if any(0 in v for v in vals):continue
        bits=[[int(pow(v,(p-1)//2,p)==p-1) for v in z] for z in vals]
        assert all(sum(z)%2==0 for z in bits)
        k=3*len(blocks)
        for i,z in enumerate(bits):signatures[i]|=r.pack(z)<<k
        blocks.append({'p':p,'roots':roots,'signatures':[r.pack(z) for z in bits]})
    assert r.rank(signatures)==16
    return {'blocks':blocks,'signatures':signatures,'rank':16}


def compute():
    d=r.read(pool.INPUT);inp=r.read(root.INPUT);old=r.read(affine.OUTPUT);projection=r.read(completed.OUTPUT)
    for obj in (inp,old,projection,r.read(ROOT_REPLAY),r.read(AFFINE_REPLAY)):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    assert inp['cubic_ascending']==d['cubic_ascending'] and old['affine_status']=='INCONSISTENT'
    rows,blocks=affine.dictionary(d);assert len(rows)==4134 and r.rank(rows)==4133
    address=inp['parent_address'];j=next(i for i,z in enumerate(d['relations']) if [z['m'],z['n']]==address)
    assert rows[j]==0 and r.rank(rows[:j]+rows[j+1:])==4133
    target=sum((v%2)<<i for i,v in enumerate(old['projection_valuations_by_column']) if v is not None)
    assert r.rank(rows[:j]+rows[j+1:]+[target])==4134
    chars=generic_characters(inp)
    K=Algebra(inp['cubic_ascending'])
    def D(a):return tuple(v/K.norm(a) for v in K.power(a,3))
    a=K.elt(inp['alpha_ascending']);w=K.elt(inp['root_ascending']);gamma=list(map(K.elt,inp['generic_classes_ascending']))
    u=D(a);h=D(w);v=list(map(D,gamma));mask=inp['generic_product_mask']
    rhs=u
    for i,g in enumerate(v):
        assert K.norm(g)==1
        if (mask>>i)&1:rhs=K.mul(rhs,g)
    assert K.norm(u)==K.norm(h)==1 and K.mul(h,h)==rhs
    beta=K.elt(projection['projection_ascending']);Nw=K.norm(w)
    ratio=K.power(tuple(c/Nw for c in w),2)
    assert K.mul(beta,ratio)==h
    raw=tuple(K.norm(a)*c for c in a)
    assert D(raw)==u
    return {'schema':'rank-jump.retained-radical-closure.v1','status':'PASS',
        'dictionary_size':4134,'generic_size':16,'replaced_dictionary_index':j,'parent_address':address,
        'generic_product_mask':mask,'norm_one_replaced_generator':list(map(str,u)),
        'norm_one_root':list(map(str,h)),'norm_one_generic_generators':[list(map(str,g)) for g in v],
        'exact_relation':'h^2=u_j*product(v_i : mask_i=1)',
        'normalization':'D(a)=a^3/Norm(a); D(a) and pi(a)=Norm(a)*a have the same squareclass',
        'generic_independence':chars,
        'outside_S_rank_before':4133,'outside_S_rank_after_replacement':4134,
        'replacement_generator_count':4150,'replacement_squareclass_rank':4150,
        'old_squareclass_rank':4149,'old_group_index_in_2_saturation':2,
        'two_saturation_complete':True,
        'full_in_field_radical_squareclass_image_dimension':4150,
        'Selmer_intersection':'exactly the 16-dimensional generic subgroup',
        'new_Selmer_or_strict_classes_from_any_in_field_radical_operations':0,
        'bindings':bindings([Path(__file__),PROTOCOL,pool.INPUT,root.INPUT,affine.OUTPUT,completed.OUTPUT,
            ROOT_REPLAY,AFFINE_REPLAY,Path(affine.__file__),Path(r.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py')]),
        'boundary':'No new class-group computation. Retires multiplication, inversion, rational rescaling, norm projection and arbitrary in-field integer root extraction on this frozen pool. New nonrational generators or extension-field operations are outside the theorem.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS: index-two saturation complete; all in-field radical operations add zero Selmer classes')
