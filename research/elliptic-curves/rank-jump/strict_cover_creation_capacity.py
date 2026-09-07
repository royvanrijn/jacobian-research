#!/usr/bin/env python3
"""Intersect global-pool capacity with strict Selmer necessities; no new search."""
import argparse
from pathlib import Path
import retrospective as r

POOL=r.OUT/'rank_jump_root_curve_capacity_comparison_v1.json'
BOUNDARY=r.OUT/'rank_jump_fresh_strict_boundary_coordinate_comparison_v1.json'
SUPPLEMENT=r.OUT/'rank_jump_fresh_retained_factor_comparison_v1.json'
OUTPUT=r.OUT/'rank_jump_strict_cover_creation_capacity_v1.json'


def compute():
    objects=[r.read(p) for p in (POOL,BOUNDARY,SUPPLEMENT)]
    for obj in objects:
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha,p
    boundary={x['token']:x for x in objects[1]['rows']}
    for x in objects[2]['rows']:
        if x['status']=='PASS':
            assert boundary[x['token']]['status']=='UNKNOWN';boundary[x['token']]=x
    rows=[]
    for x in objects[0]['rows']:
        b=boundary[x['token']];m=x['generic_dimension'];R=x['retained_rank_lower_bound']
        assert m==b['generic_rank'] and R==b['retained_rank_lower_bound'] and x['id']==b['id']
        d=x['arithmetic_global_pool_dimension_upper_bound'];q=d-m;assert q==2
        result={'token':x['token'],'id':x['id'],'family':x['family'],'generic_dimension':m,
            'retained_rank_lower_bound':R,'global_pool_excess_upper_bound':q,
            'all_rational_classes_outside_global_pool_lower_bound':max(0,R-d),
            'strict_boundary_status':b['status'],'explicit_additional_cover_basis':'UNKNOWN',
            'additional_cover_CT':'UNKNOWN'}
        if b['status']=='PASS':
            k=b['generic_strict_dimension'];a=b['additional_boundary_capacity_upper_bound']
            h=b['boundary_upper_bound'];assert h==m-k+a
            s=max(k,R-h);extra=max(0,R-m-a-q)
            assert extra==max(0,s-(k+q))
            result.update({'generic_strict_dimension':k,'additional_boundary_capacity_upper_bound':a,
                'total_strict_rational_dimension_lower_bound':s,
                'global_pool_intersection_with_strict_upper_bound':k+q,
                'strict_rational_classes_outside_global_pool_lower_bound':extra,
                'point_independent_strict_dimension_sufficient_to_exceed_global_pool':k+q+1})
        else:
            result['strict_rational_classes_outside_global_pool_lower_bound']=None
            result['reason']=b['reason']
        rows.append(result)
    # Fixed regression of the intersection bound on the seven completed high rows.
    assert {x['token']:x['strict_rational_classes_outside_global_pool_lower_bound'] for x in rows
            if x['retained_rank_lower_bound']>17 and x['strict_boundary_status']=='PASS'}=={
                'case-02':7,'case-04':3,'case-06':4,'case-09':6,'case-13':9,'case-14':2,'case-15':3}
    return {'schema':'rank-jump.strict-cover-creation-capacity.v1','status':'PASS','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),POOL,BOUNDARY,SUPPLEMENT]},
        'boundary':'Linear-algebra consequences of prior equation-only capacities joined with retained rank lower bounds. Not newly measured class groups, explicit new covers, or prospective rank features. Unknown boundary rows stay unknown.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',[(x['id'],x['strict_rational_classes_outside_global_pool_lower_bound']) for x in result['rows']])
