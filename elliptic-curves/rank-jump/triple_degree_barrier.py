#!/usr/bin/env python3
"""Check the geometric hypotheses of the panel-wide triple degree barrier."""
import argparse
from itertools import combinations
from pathlib import Path
from sage.all import QQ,PolynomialRing,matrix
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
VERIFIED=r.OUT/'rank_jump_degree_one_relation_panel_verification_v1.json'
TRIPLES=r.OUT/'rank_jump_trace_zero_triple_panel_v1.json'
OUTPUT=r.OUT/'rank_jump_triple_degree_barrier_v1.json'


def compute():
    inp=r.read(INPUT);v=r.read(VERIFIED);old=r.read(TRIPLES)
    assert v['status']=='PASS' and v['triple_parities_checked']==2853 and v['trace_even_triple_count']==0
    for path,sha in v['bindings'].items():assert hash_file(r.ROOT/path)==sha
    R=PolynomialRing(QQ,'t');A=R(inp['A']);B=R(inp['B']);delta=4*A**3+27*B**2
    assert delta.degree()==24 and delta.is_squarefree()
    qs={label:R(c['q']) for label,c in inp['covers'].items()}
    for q in qs.values():assert q.degree()==2 and q.is_squarefree() and q.gcd(delta)==1
    for labels in inp['pairs']:assert qs[labels[0]].gcd(qs[labels[1]])==1
    G=matrix(QQ,inp['gram']);assert G.is_positive_definite() and G[0,0]==4
    rows=[];count=0
    for key,block in inp['blocks'].items():
        n=0
        for labels in combinations(block['labels'],3):
            n+=1
            w=[sum(inp['covers'][label]['trace'][i] for label in labels) for i in range(17)]
            assert any(x%2 for x in w)
        if n:rows.append({'block_key':key,'triple_count':n,'all_signed_trace_sums_nonzero_mod_2M':True,
                          'total_intersection_degree_lower_bound':6})
        count+=n
    assert count==2853
    return {'schema':'rank-jump.triple-degree-barrier.v1','status':'PASS','triple_count':count,'rows':rows,
            'generic_hypotheses':{'rootless_nodal_fibres':24,'native_trace_height':10,'native_quadratic_maps':len(qs),
                'all_native_branch_divisors_disjoint_from_singular_fibres':True,'co_split_pair_branch_gcds_checked':len(inp['pairs']),
                'generic_lattice_minimum':4,'minimum_proof':'Independent LDL verification inherited from the bound verifier input.'},
            'class_identity':'B_tau=(2,2,tau), D_z=(h(z)/2+2,4,2z), B_tau.D_z=h(z-tau)+2; h(tau)=10.',
            'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,VERIFIED,TRIPLES,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
            'boundary':'Total proper intersection length is at least6 for all integer generic translates of the tested triples. This does not exclude degree-one rational components of a higher-degree intersection algebra. No claim about triples outside the frozen panel.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS 2853 triples; all generic translates have total intersection degree at least6')
