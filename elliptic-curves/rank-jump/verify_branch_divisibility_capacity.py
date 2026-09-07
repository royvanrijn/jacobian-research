#!/usr/bin/env python3
"""Independent finite arithmetic, Sturm topology and F2 branch-rank replay."""
import argparse
from itertools import combinations
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import retrospective as r
import branch_divisibility_capacity as first
import complete_branch_divisibility_capacity as run

OUTPUT=r.OUT/'rank_jump_branch_divisibility_capacity_verification_v1.json'


def compute():
    from sage.all import QQ,GF,PolynomialRing,matrix
    d=r.read(first.INPUT);out=r.read(run.OUTPUT)
    for obj in (d,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    R=PolynomialRing(QQ,'t');A=R(d['A']);B=R(d['B']);D=-4*A**3-27*B**2
    section=[(R(s['x']),R(s['y'])) for s in d['sections']]
    assert len(section)==17 and all(y*y==x**3+A*x+B for x,y in section)
    assert [A.degree(),B.degree(),D.degree()]==[8,12,24] and D.gcd(D.derivative())==1 and D.gcd(A)==1
    def sturm(f):
        seq=[f,f.derivative()]
        while seq[-1]:
            rem=-(seq[-2]%seq[-1])
            if not rem:break
            seq.append(rem)
        return seq
    def variation(seq,t=None,side=1):
        values=[int(f(t).sign()) if t is not None else int(f.leading_coefficient().sign())*(-1 if side<0 and f.degree()%2 else 1) for f in seq]
        values=[x for x in values if x];return sum(a!=b for a,b in zip(values,values[1:]))
    ds,bs=sturm(D),sturm(B);top=out['topology'];events=top['events']
    assert variation(ds,side=-1)-variation(ds,side=1)==len(events)==18
    for i,e in enumerate(events):
        a,b=map(QQ,e['interval']);assert a<b and D(a)*D(b)<0
        assert variation(ds,a)-variation(ds,b)==1 and variation(bs,a)==variation(bs,b)
        assert int(B(a).sign())==int(B(b).sign())==e['B_sign'] and e['type']=='I1'
        if i:assert QQ(events[i-1]['interval'][1])<a
    samples=list(map(QQ,top['samples']));arities=[3 if D(t)>0 else 1 for t in samples]
    assert arities==top['arities']
    for i,t in enumerate(samples):
        assert t>QQ(events[i]['interval'][1])
        if i+1<len(events):assert t<QQ(events[i+1]['interval'][0])
    adjacency={(i,j):set() for i,n in enumerate(arities) for j in range(n)}
    def edge(a,b):adjacency[a].add(b);adjacency[b].add(a)
    for i,e in enumerate(events):
        left=(i-1)%len(events);right=i;assert sorted((arities[left],arities[right]))==[1,3]
        many=left if arities[left]==3 else right;one=right if many==left else left
        pair=(1,2) if e['B_sign']>0 else (0,1);single=0 if e['B_sign']>0 else 2
        edge((many,pair[0]),(many,pair[1]));edge((many,single),(one,0))
    unseen=set(adjacency);components=0
    while unseen:
        todo=[unseen.pop()];components+=1
        while todo:
            for v in adjacency[todo.pop()]:
                if v in unseen:unseen.remove(v);todo.append(v)
    assert components==top['real_components']==10 and top['genus']==10 and top['global_pool_upper_bound']==19
    def value(poly,t,p):
        z=0
        for c in reversed(poly.list()):
            c=Q(str(c));z=(z*t+c.numerator*pow(c.denominator,-1,p))%p
        return z
    block_checks=0;character_checks=0
    def check_block(block,q=None,parameter=None):
        nonlocal block_checks,character_checks
        p=block['p'];t=block['base_root'];assert 3<=p<=2003 and all(p%i for i in range(2,isqrt(p)+1))
        if q is not None:assert value(q,t,p)==0 and value(q.derivative(),t,p)!=0
        else:
            c=Q(parameter);assert t==c.numerator*pow(c.denominator,-1,p)%p
        a,b=value(A,t,p),value(B,t,p);assert (4*a**3+27*b*b)%p
        roots=sorted(x for x in range(p) if (x*x*x+a*x+b)%p==0)
        assert roots==block['cubic_roots'] and len(roots)==3
        squares={x*x%p for x in range(1,p)};signatures=[]
        for x,y in section:
            xx,yy=value(x,t,p),value(y,t,p);assert yy and yy*yy%p==(xx**3+a*xx+b)%p
            bits=[int((xx-z)%p not in squares) for z in roots];assert sum(bits)%2==0
            signatures.append(sum(bit<<i for i,bit in enumerate(bits)));character_checks+=3
        assert signatures==block['signatures'];block_checks+=1
        return [[(v>>i)&1 for v in signatures] for i in range(3)]
    global_rows=[]
    for block in out['generic_independence_blocks']:global_rows.extend(check_block(block,parameter=out['generic_independence_parameter']))
    assert matrix(GF(2),global_rows).rank()==out['generic_mod_two_dimension']==17
    results=[];polys=[];kernels=[]
    for inp,row in zip(d['covers'],out['rows']):
        assert inp['label']==row['label'];q=R(inp['q']);polys.append(q)
        assert q.degree()==2 and not q.discriminant().is_square() and q.gcd(D)==1
        w=row['irreducibility_witness'];assert w
        p,t=w['p'],w['base_root'];assert 3<=p<=2003 and all(p%i for i in range(2,isqrt(p)+1))
        assert value(q,t,p)==0 and value(q.derivative(),t,p)!=0
        a,b=value(A,t,p),value(B,t,p);assert (4*a**3+27*b*b)%p
        assert w['cubic_ascending']==[b,a,0,1]
        assert all((x*x*x+a*x+b)%p for x in range(p))
        chars=[]
        for block in row['blocks']:chars.extend(check_block(block,q=q))
        M=matrix(GF(2),chars);rank=int(M.rank());K=M.right_kernel()
        assert rank==row['branch_character_rank']==16 and K.dimension()==1
        mask=sum(int(x)<<i for i,x in enumerate(K.basis()[0]));kernels.append(mask)
        assert row['branch_divisibility_kernel_upper_bound']==1 and row['all_scalar_twists_rank_upper_bound']==19-rank==3
        results.append({'label':row['label'],'rank':rank,'finite_kernel_mask':mask,'twist_rank_upper_bound':3,'irreducibility_prime':p})
    assert len(results)==37 and len(set(kernels))==37
    for q,qq in combinations(polys,2):assert q.gcd(qq)==1
    return {'schema':'rank-jump.branch-divisibility-capacity-verification.v1','status':'PASS',
        'rows':results,'finite_blocks_verified':block_checks,'quadratic_character_bits_verified':character_checks,
        'generic_mod_two_dimension':17,'real_components':components,'global_pool_upper_bound':19,
        'distinct_branch_kernel_lines':37,'pairwise_disjoint_branch_supports':666,
        'all_single_scalar_twists_rank_upper_bound':3,
        'all_scalar_twists_of_products_of_at_least_two_distinct_supports_rank_upper_bound':2,
        'method':'Direct exhaustive finite cubic root tests and residue-square sets; Sage F2 ranks/kernels instead of bit elimination; exact rational Sturm isolation and independent graph traversal. No new sections or specializations are searched.',
        'bindings':first.bindings([Path(__file__),first.INPUT,run.OUTPUT,Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',len(result['rows']),'covers;',result['finite_blocks_verified'],'blocks;',result['quadratic_character_bits_verified'],'character bits; distinct kernels',result['distinct_branch_kernel_lines'])
