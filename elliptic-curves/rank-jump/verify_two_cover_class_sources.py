#!/usr/bin/env python3
"""Replay the oracle class labels independently with Sage cubic arithmetic."""
import argparse
from pathlib import Path
import retrospective as r
import two_cover_pencil_geometry as pencils

OUTPUT=r.OUT/'rank_jump_two_cover_class_source_verification_v1.json'


def compute():
    from sage.all import QQ,PolynomialRing
    rows=[]
    inp=r.read(pencils.INPUT)
    for name,digest in inp['source_hashes'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    for data in inp['rows'][:6]:
        i=data['case_index'];source=pencils.bad.cases()[i]
        old=r.read(pencils.bad.INPUT)['cases'][i];block=r.read(pencils.strict.OUTPUT)['rows'][i]
        profile,_,_=r.characterize(source)
        assert profile['independent_input_indices']==old['selected_input_indices']
        assert data['id']==source['id'] and data['cubic_ascending']==old['integral_cubic_ascending']
        assert data['all_bad_places_complete']==block['all_bad_places_complete']
        assert data['generic_rank']==profile['generic_subgroup_rank_exact']
        assert data['known_independent_rank']==profile['certified_independent_subgroup_rank_exact']
        expected=[('generic',x) for x in block['generic_strict_kernel_masks']]
        expected += [('relative',x) for x in block['relative_strict_lift_masks']]
        assert [(x['kind'],x['witness_mask']) for x in data['classes']]==expected
        _,points=r.short(source['model'],source['generic_points']+source['points'])
        scale=QQ(old['elliptic_scaling_d'])
        R=PolynomialRing(QQ,'z');f=R(list(map(QQ,data['cubic_ascending'])))
        K=R.quotient(f,'theta');theta=K.gen()
        gammas=[];ys=[]
        for index in old['selected_input_indices']:
            x,y=map(QQ,points[index]);x*=scale**2;y*=scale**3
            assert f(x)==y*y
            gammas.append(x-theta);ys.append(y)
        for cls in data['classes']:
            beta=K(1);root=QQ(1)
            for j,gamma in enumerate(gammas):
                if cls['witness_mask']>>j&1:beta*=gamma;root*=ys[j]
            assert [str(beta.lift()[j]) for j in range(3)]==cls['beta']
            assert str(root)==cls['norm_root']
        m=data['generic_rank'];masks=[x['witness_mask'] for x in data['classes']]
        qdim=r.rank([mask>>m for mask in masks])
        assert r.rank(masks)==len(masks)
        assert qdim==len(block['relative_strict_lift_masks'])
        rows.append({'case_index':i,'id':source['id'],'status':'PASS','class_count':len(masks),
            'independent_class_dimension':len(masks),'relative_quotient_dimension':qdim,
            'original_generic_rank':m,'original_known_independent_rank':data['known_independent_rank'],
            'original_observed_quotient_rank':data['known_independent_rank']-m,
            'full_strict_locality_certified':data['all_bad_places_complete']})
    small=r.read(pencils.small.OUTPUT)['records']
    assert inp['rows'][6]['classes']==[{'kind':'small_strict','beta_index':x['beta_index'],
        'beta':x['beta'],'norm_root':'25'} for x in small]
    assert sum(x['class_count'] for x in rows)==48
    return {'schema':'rank-jump.two-cover-class-source-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),pencils.INPUT)},
        'rows':rows,'small_class_source_verified':True,
        'boundary':'Retrospective label verification. Incomplete local kernels are not full strict Selmer spaces; original known subgroup ranks are not curve-rank upper bounds.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS all 48 production class sources and two small controls')
    else:r.write_new(OUTPUT,data)
