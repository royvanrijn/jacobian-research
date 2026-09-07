#!/usr/bin/env python3
"""Independent finite check of the abelian obstruction to one-prime twists."""
import argparse
from itertools import product
from pathlib import Path
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_governing_cochain_gate_v1.json'
OUTPUT=r.OUT/'rank_jump_governing_prime_compatibility_v1.json'


def compute():
    old=r.read(INPUT);assert old['status']=='PASS'
    for path,sha in old['bindings'].items():assert hash_file(r.ROOT/path)==sha
    G=[x for x in product((0,1),repeat=4) if (x[0]*x[3]-x[1]*x[2])%2]
    def action(g,v):return ((g[0]*v[0]+g[1]*v[1])%2,(g[2]*v[0]+g[3]*v[1])%2)
    def plus(v,w):return tuple((a+b)%2 for a,b in zip(v,w))
    def e(v,w):return (v[0]*w[1]+v[1]*w[0])%2
    def mat(g,h):return tuple(sum(g[2*i+k]*h[2*k+j] for k in range(2))%2 for i in range(2) for j in range(2))
    V=list(product((0,1),repeat=2));H=list(product(V,V,G,(0,1)));identity=((0,0),(0,0),(1,0,0,1),0)
    def mul(x,y):
        a,b,g,z=x;c,d,h,w=y
        return plus(a,action(g,c)),plus(b,action(g,d)),mat(g,h),(z+w+e(a,action(g,d)))%2
    inv={x:next(y for y in H if mul(x,y)==identity and mul(y,x)==identity) for x in H}
    comm={mul(mul(mul(x,y),inv[x]),inv[y]) for x in H for y in H}
    derived={identity};todo=[identity]
    while todo:
        x=todo.pop()
        for y in comm:
            z=mul(x,y)
            if z not in derived:derived.add(z);todo.append(z)
    assert len(derived)==96 and len(H)//len(derived)==2
    counts={0:0,1:0}
    for x in H:
        a,b,g,z=x
        if any(action(g,v)==v for v in V if v!=(0,0)):continue
        assert x in derived
        pp=next(v for v in V if plus(action(g,v),v)==a)
        counts[(e(pp,b)+z)%2]+=1
    assert counts=={0:32,1:32}
    rows=[]
    for row in old['rows']:
        k=row['strict_rational_block_dimension'];slots=k*(k-1)//2
        rows.append({'id':row['id'],'block_dimension':k,'pairing_slots':slots,
            'governing_extension_degree_formula_value':6*2**(2*k+slots),
            'conditional_prime_density_of_zero_restricted_CT':{'numerator':1,'denominator':2**slots},
            'density_type':'Theorem-derived second-descent distribution for a fixed block in congruence-conditioned inert prime twists; not a rational-rank distribution or an empirical count.'})
    return {'schema':'rank-jump.governing-prime-compatibility.v1','status':'PASS',
        'two_class_group':{'order':len(H),'commutator_subgroup_order':len(derived),'abelianization_order':2,
            'all_fixed_point_free_lifts_in_commutator_subgroup':True,'governing_values':{str(k):v for k,v in counts.items()}},
        'rows':rows,'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
        'boundary':'Finite n=2 algebra independently replayed. General-n abelianization and Chebotarev application are proved in the accompanying note. No actual governing polynomial, selected twist prime, or positive Mordell-Weil implication.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: independent order192 group, commutator order96, both governing values compatible with trivial abelian projection')
