#!/usr/bin/env python3
"""Independent ideal-valuation and GF(2) replay of the retained batch gate."""
import argparse
from pathlib import Path
from math import prod
from fractions import Fraction as Q
import retrospective as r
import retained_norm_batch_capacity as source
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_retained_norm_batch_capacity_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,pari,GF,matrix,PolynomialRing
    pari.allocatemem(64000000,268435456,silent=True)
    data=r.read(source.INPUT);result=r.read(source.OUTPUT)
    for obj in (result,r.read(source.PROVENANCE)):
        for path,sha in obj['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    assert result['status']=='PASS'
    R=PolynomialRing(QQ,'x');f=R(data['cubic_ascending']);nf=pari.nfinit([pari(f),data['S_finite']]);th=pari.Mod('x',pari(f))
    assert str(nf.disc())==data['field_discriminant']
    b,a=map(QQ,data['affine_map_old_root_from_masked'])
    assert f(a*R.gen()+b)==a**3*R(data['masked_cubic_ascending'])
    cols=data['columns'];blocks={};ideals={}
    for i,c in enumerate(cols):blocks.setdefault(c['p'],[]).append(i)
    for p,indices in blocks.items():
        assert ZZ(p).is_prime(proof=True)
        dec=list(pari.idealprimedec(nf,p));assert len(dec)==len(indices)
        for i,P in zip(indices,dec):
            assert str(pari.idealhnf(nf,P))==cols[i]['hnf']
            assert int(P[2])==cols[i]['e'] and int(P[3])==cols[i]['f'];ideals[i]=P
    algebra=Algebra(data['cubic_ascending']);S=set(data['S_finite']);rows=[];checks=0
    assert len(result['rows'])==len(data['relations'])==296
    for rel,saved in zip(data['relations'],result['rows']):
        assert (rel['m'],rel['n'])==(saved['m'],saved['n'])
        alpha=pari(R(rel['alpha_ascending']))(th);N=int(saved['norm_alpha'])
        # Rational multiplication determinant, independent of PARI's field norm.
        assert algebra.norm(algebra.elt(rel['alpha_ascending']))==Q(N)
        vals=dict(rel['ideal_factorization']);norm_product=1;projection=pari(N)*alpha;sparse=[]
        for p in sorted({cols[i]['p'] for i in vals}):
            v=int(ZZ(N).valuation(p));assert v==sum(cols[i]['f']*vals.get(i,0) for i in blocks[p])
            norm_product*=p**v
            for i in blocks[p]:
                P=ideals[i];e=int(pari.idealval(nf,alpha,P));assert e==vals.get(i,0)
                ee=int(pari.idealval(nf,projection,P));assert ee==e+int(P[2])*v;checks+=2
                if p not in S and ee%2:sparse.append(i)
        assert norm_product==abs(N) # All other valuations vanish: alpha is integral.
        basis_coordinates=pari.nfalgtobasis(nf,alpha)
        assert all(QQ(x).denominator()==1 for x in basis_coordinates)
        assert sparse==saved['outside_parity_columns']
        rows.append(sparse)
    # Independent Sage elimination on the compact used-coordinate matrix.
    active=sorted({i for row in rows for i in row});position={c:i for i,c in enumerate(active)}
    M=matrix(GF(2),len(rows),len(active),sparse=True)
    for i,row in enumerate(rows):
        for j in row:M[i,position[j]]=1
    rank=int(M.rank());assert rank==result['outside_S_parity_rank']==296
    assert M.left_kernel().dimension()==result['coefficient_kernel_dimension']==0
    assert result['additional_strict_class_dimension_in_projected_dictionary']==0
    files=(Path(__file__),source.INPUT,source.OUTPUT,source.PROVENANCE,
           Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return {'schema':'rank-jump.retained-norm-batch-capacity-verification.v1','status':'PASS',
        'rational_norm_checks':len(rows),'independent_ideal_valuation_checks':checks,
        'complete_rational_prime_blocks':len(blocks),'used_odd_valuation_coordinates':len(active),
        'outside_parity_rank':rank,'coefficient_kernel_dimension':0,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Verifies the whole 296-element dictionary exclusion. Does not assert completeness of the dictionary or replay its original exhaustive smoothness sieve.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result,flush=True)
