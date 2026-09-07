#!/usr/bin/env python3
"""Shared-class and carrier bounds after the exact Frobenius improvement."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r
import branch_divisibility_capacity as branch
import verify_branch_divisibility_capacity as old
import verify_native_root_frobenius_slack as frob
import verify_native_branch_half_section as halves

OUTPUT=r.OUT/'rank_jump_native_common_class_capacity_v1.json'


def compute():
    v=r.read(old.OUTPUT);f=r.read(frob.OUTPUT);h=r.read(halves.OUTPUT)
    assert all(x['status']=='PASS' for x in (v,f,h))
    for obj in (v,f,h):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    assert f['global_pool_dimension_interval']==[17,18]
    masks=[x['kernel_mask'] for x in h['rows']]
    assert len(set(masks))==37 and all(a^b^c for a,b,c in combinations(masks,3))
    rows=[]
    for k in range(1,38):
        gain=2*k+(2**k-1-k)
        rows.append({'supports':k,'genus':1+2**(k-1)*(k-2),
            'generic_gain_upper_bound':gain,'full_generic_rank_upper_bound':17+gain})
    # A_i=K_i\G is empty or a two-point affine line. Distinct trace
    # directions forbid parallel duplicate edges; no dependent triple
    # forbids a triangle. Enumerate every simple graph on at most8
    # vertices with at most4 edges to verify the nonconcurrent bound.
    edges=list(combinations(range(8),2));count=0;max_nonconcurrent=0;max_all=0
    degree_profiles=set()
    for n in range(5):
        for graph in combinations(edges,n):
            adj=[set() for _ in range(8)]
            for a,b in graph:adj[a].add(b);adj[b].add(a)
            if any(adj[a]&adj[b] for a,b in graph):continue
            degrees=[len(x) for x in adj];mixed=sum(2**x-1-x for x in degrees)
            count+=1;max_all=max(max_all,mixed)
            if max(degrees,default=0)<4:
                max_nonconcurrent=max(max_nonconcurrent,mixed)
                degree_profiles.add(tuple(sorted((x for x in degrees if x),reverse=True)))
    assert max_nonconcurrent==5 and max_all==11
    return {'schema':'rank-jump.native-common-class-capacity.v1','status':'PASS','rows':rows,
        'global_pool_dimension_interval':[17,18],'extra_global_class_dimension_upper_bound':1,
        'single_twist_rank_upper_bound':2,'multiple_twist_rank_upper_bound':1,
        'independent_trace_triples_checked':7770,'dependent_trace_triples':0,
        'triangle_free_graphs_checked':count,'nonconcurrent_degree_profiles':sorted(degree_profiles),
        'four_support_nonconcurrent_mixed_character_capacity':5,
        'four_support_nonconcurrent_total_gain_upper_bound':13,
        'four_support_gain_at_least14_requires_one_common_extra_class':True,
        'four_support_gain14_requires_common_class_rational_on_at_least_twists':10,
        'three_support_gain10_requires_common_class_rational_on_all_twists':7,
        'character_lattice_quotient_two_rank_if_gain_J_exceeds_k':'J-1',
        'boundary':'Necessary structure for large generic blocks on these retained native supports. The extra global class, its branch trivializations and its rational realizations on the twists remain UNKNOWN. No rank predictor or actual successful jump mechanism is claimed.',
        'bindings':branch.bindings([Path(__file__),old.OUTPUT,frob.OUTPUT,halves.OUTPUT,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS: gains <=2,5,10,19 on 1,2,3,4 supports;',result['triangle_free_graphs_checked'],'graphs;',result['independent_trace_triples_checked'],'trace triples')
