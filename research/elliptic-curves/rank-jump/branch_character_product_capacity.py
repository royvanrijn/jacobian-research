#!/usr/bin/env python3
"""Deduce all subset-cover capacities from the verified branch kernels."""
import argparse
from pathlib import Path
from math import comb
import retrospective as r
import verify_branch_divisibility_capacity as verified
import branch_divisibility_capacity as source

OUTPUT=r.OUT/'rank_jump_branch_character_product_capacity_v1.json'


def compute():
    v=r.read(verified.OUTPUT);assert v['status']=='PASS'
    for name,sha in v['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    assert len(v['rows'])==37 and len({x['finite_kernel_mask'] for x in v['rows']})==37
    assert v['pairwise_disjoint_branch_supports']==comb(37,2) and v['global_pool_upper_bound']==19
    rows=[]
    for k in range(1,38):
        # Every singleton character has capacity3. Every larger character
        # has joint branch rank17, hence capacity2.
        larger=sum(comb(k,j) for j in range(2,k+1))
        new=3*k+2*larger;assert new==k+2**(k+1)-2
        rows.append({'native_supports':k,'cover_degree':2**k,'cover_genus':1+2**(k-1)*(k-2),
            'singleton_characters':k,'larger_characters':larger,
            'generic_rank_gain_upper_bound':new,'full_generic_rank_upper_bound':17+new})
    assert rows[0]['generic_rank_gain_upper_bound']==3 and rows[1]['generic_rank_gain_upper_bound']==8
    assert rows[2]['cover_genus']==5
    return {'schema':'rank-jump.branch-character-product-capacity.v1','status':'PASS','rows':rows,
        'two_cover_pairs_covered':666,'single_twist_upper_bound':3,'product_twist_upper_bound':2,
        'native_supports_necessary_for_a_generic_gain_at_least_nine':3,
        'genus_necessary_for_a_gain_at_least_nine_in_this_full_native_fibre_product_construction':5,
        'boundary':'All nonzero rational scalar twists on these 37 supports. Bounds concern complete generic Mordell-Weil groups and allow integral gluing. They do not bound specialization jumps, prove a rational base point, or exclude other curves/carriers or other native supports.',
        'bindings':source.bindings([Path(__file__),verified.OUTPUT,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS all subset capacities; one cover <=3, two covers <=8 new generic directions')
