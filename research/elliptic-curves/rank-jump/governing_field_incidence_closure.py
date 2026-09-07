#!/usr/bin/env python3
"""Generic governing fields contain obstruction data, but no added S4 classes."""
import argparse
from collections import deque
from itertools import product
from pathlib import Path
import retrospective as r
import governing_cochain_gate as group_source
import fresh_governing_panel as panel
import fresh_governing_completion as completion
import fresh_governing_octics as octics
import fresh_retained_factors as supplement
import verify_fresh_governing_panel as panel_verifier

PROTOCOL=Path(__file__).with_name('GOVERNING_FIELD_INCIDENCE_CLOSURE_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_governing_field_incidence_closure_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_governing_field_incidence_closure_v1.json'
SUPPLEMENT_REPLAY=r.OUT/'rank_jump_fresh_retained_factor_verification_v1.json'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    base={x['token']:x for x in r.read(panel.OUTPUT)['rows']}
    patch={x['token']:x for x in r.read(completion.OUTPUT)['rows']}
    extra={x['token']:x for x in r.read(supplement.OUTPUT)['rows']}
    rows=[]
    for o in r.read(octics.OUTPUT)['rows']:
        token=o['token'];loc=base[token]['local'];p=patch.get(token,{})
        if extra.get(token,{}).get('factor',{}).get('status')=='PASS' and 'generic_strict_masks' in extra[token]['local']:
            loc=extra[token]['local'];p={}
        if 'local' in loc:
            local=loc['local'];masks=loc['generic_strict_masks']
            k=p.get('strict_generic_dimension',loc.get('strict_generic_dimension'))
            upper=p.get('strict_generic_dimension_upper_bound',k)
        else:
            local=p['local'];masks=p['tested_kernel_masks'];k=p['strict_generic_dimension'];upper=p['strict_generic_dimension_upper_bound']
        witnesses=[]
        for word in (1,2,3):
            for block in local:
                sig=0
                for i in range(2):
                    if (word>>i)&1:sig^=r.pack(block['signatures'][i])
                if sig:
                    witnesses.append({'word':word,'place':block['place'],
                        'pair_signatures':block['signatures'][:2],'nonzero_signature':sig});break
            else:raise AssertionError((token,'No local witness',word))
        assert not any(r.reduce(w,r.basis(masks))==0 for w in (1,2,3))
        assert o['status']=='PASS' and r.rank(o['independence_signatures'])==2 and o['galois']['galois_group']=='S3'
        rows.append({'token':token,'generic_dimension':base[token]['local'].get('generic_dimension',p.get('generic_dimension',loc.get('generic_dimension'))),
            'joint_class_field_degree':o['joint_class_field_degree'],'governing_field_degree':o['governing_field_degree'],
            'generic_pair_indices':o['generic_pair_indices'],'pair_independence_signatures':o['independence_signatures'],
            'local_witnesses':witnesses,'all_generic_strict_dimension_interval':[k,k] if k is not None else [0,upper]})
    assert len(rows)==16
    r.write_new(INPUT,{'schema':'rank-jump.governing-field-incidence-closure-inputs.v1','rows':rows,
        'bindings':bindings([Path(__file__),PROTOCOL,panel.OUTPUT,completion.OUTPUT,octics.OUTPUT,supplement.OUTPUT,
                            panel_verifier.OUTPUT,SUPPLEMENT_REPLAY]),
        'boundary':'Only inherited generic pair independence and local signatures; no exceptional points or rank labels.'})


def finite_group():
    G=[g for g in product((0,1),repeat=4) if (g[0]*g[3]-g[1]*g[2])%2]
    I=(1,0,0,1);unit=(0,0,I,0)
    def mul(x,y):
        a,b,g,z=x;c,d,h,w=y
        return (a^group_source.act(g,c),b^group_source.act(g,d),group_source.mm(g,h),z^w^group_source.weil(a,group_source.act(g,d)))
    elements=list(product(range(4),range(4),G,range(2)))
    gens=[(1,0,I,0),(2,0,I,0),(0,1,I,0),(0,2,I,0),(0,0,(0,1,1,0),0),(0,0,(0,1,1,1),0)]
    def linact(g,v):return tuple((v[0] if g[2*i] else 0)^(v[1] if g[2*i+1] else 0) for i in range(2))
    expressions={unit:(0,0)};pending=deque([unit]);constraints=[]
    while pending:
        x=pending.popleft()
        for j,g in enumerate(gens):
            y=mul(x,g);tail=linact(x[2],(1<<(2*j),1<<(2*j+1)))
            expr=tuple(a^b for a,b in zip(expressions[x],tail))
            if y not in expressions:expressions[y]=expr;pending.append(y)
            else:constraints.extend(a^b for a,b in zip(expressions[y],expr))
    assert len(expressions)==len(elements)==192
    rank=r.rank(constraints);solutions=[s for s in range(4096) if all((s&v).bit_count()%2==0 for v in constraints)]
    assert rank==8 and len(solutions)==16
    assignment=lambda fun:sum(fun(g)<<(2*i) for i,g in enumerate(gens))
    alpha=assignment(lambda g:g[0]);beta=assignment(lambda g:g[1])
    cob=[assignment(lambda g:group_source.act(g[2],v)^v) for v in range(4)]
    assert len(set(cob))==4 and set(solutions)=={c^(alpha if a else 0)^(beta if b else 0) for c in cob for a in (0,1) for b in (0,1)}
    center=(0,0,I,1);assert expressions[center]==(0,0)
    P=[x for x in elements if x[2]==I];assert len(P)==32
    inverse={x:next(y for y in P if mul(x,y)==unit and mul(y,x)==unit) for x in P}
    comm={mul(mul(mul(x,y),inverse[x]),inverse[y]) for x in P for y in P}
    assert comm=={unit,center}
    return {'group_order':192,'kernel_over_cubic_splitting_field_order':32,
        'kernel_commutator_order':2,'kernel_abelianization_order':16,
        'generator_count':6,'cochain_variable_count':12,'constraint_rank':rank,
        'Z1_dimension':4,'B1_dimension':2,'H1_dimension':2,
        'alpha_assignment':alpha,'beta_assignment':beta,'coboundary_assignments':cob,'cocycle_assignments':solutions,
        'central_value_expressions':list(expressions[center]),
        'basis_description':'Exactly the two inherited coordinate cocycles; the central governing bit contributes no H1 class.'}


def compute():
    inp=r.read(INPUT);rows=[]
    for row in inp['rows']:
        assert row['joint_class_field_degree']==96 and row['governing_field_degree']==192
        assert row['generic_pair_indices']==[0,1] and r.rank(row['pair_independence_signatures'])==2
        for w in row['local_witnesses']:
            value=0
            for i in range(2):
                if (w['word']>>i)&1:value^=r.pack(w['pair_signatures'][i])
            assert value==w['nonzero_signature']!=0
        rows.append({'token':row['token'],'pair_governing_total_class_dimension':2,
            'pair_governing_new_class_dimension':0,'pair_governing_strict_class_dimension':0,
            'pair_nontrivial_S4_quotients':3,'pair_strict_S4_quotients':0,
            'all_generic_pair_fields_new_class_dimension':0,
            'all_generic_pair_fields_strict_dimension_interval':row['all_generic_strict_dimension_interval']})
    return {'schema':'rank-jump.governing-field-incidence-closure.v1','status':'PASS','finite_group':finite_group(),'rows':rows,
        'third_independent_class_joint_field_degree':6*4**3,'retained_pair_governing_degree':192,
        'central_extension_invariance':'If C is central and acts trivially on V with V^S3=0, then H1(Gamma/C,V) -> H1(Gamma,V) is an isomorphism.',
        'generic_compositum_conclusion':'Every S4 cubic Kummer class whose Galois closure lies in the compositum of all generic pair governing fields belongs to the generic subgroup.',
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT,Path(group_source.__file__),Path(r.__file__)]),
        'boundary':'S4 Kummer quotients with the specified cubic resolvent. Other subfields may exist. Strict zero concerns the pair governing fields, not the whole curve. No new class or CT value computed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);args=p.parse_args()
    if args.mode=='export':export()
    else:
        result=compute()
        if args.mode=='build':r.write_new(OUTPUT,result)
        else:assert result==r.read(OUTPUT)
        print('PASS H1 dimension2; zero extra classes and zero strict classes in every retained pair governing field')
