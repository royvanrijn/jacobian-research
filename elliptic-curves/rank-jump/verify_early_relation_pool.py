#!/usr/bin/env python3
"""Independent valuation, matrix and square-root replay of the early pool."""
import argparse
from pathlib import Path
from math import prod,gcd
from fractions import Fraction as Q
import retrospective as r
import early_relation_pool as source
import retained_norm_inherited_hit as hit
import bounded_gain_reference as ref
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_early_relation_pool_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,pari,PolynomialRing,GF,matrix
    pari.allocatemem(64000000,268435456,silent=True)
    d=r.read(source.INPUT);out=r.read(source.OUTPUT)
    for path,sha in out['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);nf=pari.nfinit([pari(f),d['S_finite']]);th=pari.Mod('z',pari(f))
    assert str(nf.disc())==d['field_discriminant']
    cols=d['columns'];blocks={};ideals={}
    for i,c in enumerate(cols):blocks.setdefault(c['p'],[]).append(i)
    for p,indices in blocks.items():
        assert ZZ(p).is_prime(proof=True);dec=list(pari.idealprimedec(nf,p));assert len(indices)==len(dec)
        for i,P in zip(indices,dec):
            assert str(pari.idealhnf(nf,P))==cols[i]['hnf']
            assert int(P[2])==cols[i]['e'] and int(P[3])==cols[i]['f'];ideals[i]=P
    alg=Algebra(d['cubic_ascending']);S=set(d['S_finite']);rows=[];checks=0;seen=set()
    assert len(d['relations'])==len(out['norms_alpha'])==out['element_count']==4134
    for row,normstr in zip(d['relations'],out['norms_alpha']):
        key=(row['m'],row['n']);assert gcd(*key)==1 and key[1]>0 and key not in seen;seen.add(key)
        alpha=pari(R(row['alpha_ascending']))(th);N=int(normstr)
        assert alg.norm(alg.elt(row['alpha_ascending']))==Q(N)
        assert all(QQ(x).denominator()==1 for x in pari.nfalgtobasis(nf,alpha))
        vals=dict(row['ideal_factorization']);total=1;parity=[]
        for p in {cols[i]['p'] for i in vals}:
            v=int(ZZ(N).valuation(p));assert v==sum(cols[i]['f']*vals.get(i,0) for i in blocks[p]);total*=p**v
            for i in blocks[p]:
                e=int(pari.idealval(nf,alpha,ideals[i]));assert e==vals.get(i,0);checks+=1
                if p not in S and (e+cols[i]['e']*v)%2:parity.append(i)
        assert total==abs(N);rows.append(parity)
    active=sorted({i for row in rows for i in row});position={j:i for i,j in enumerate(active)}
    M=matrix(GF(2),len(rows),len(active),sparse=True)
    for i,row in enumerate(rows):
        for j in row:M[i,position[j]]=1
    assert M.rank()==out['outside_S_parity_rank']==4133
    kernel=M.left_kernel();supports=[[i for i,b in enumerate(v) if b] for v in kernel.basis()]
    assert supports==out['kernel_supports']==[[1769]]
    for phase in out['phases']:
        n=phase['cumulative_unique_elements'];rank=int(M[:n,:].rank())
        assert rank==phase['outside_S_parity_rank'] and n-rank==phase['kernel_dimension']
    # Reconstruct the sole dependency and its generic product in rational algebra.
    saved=r.read(hit.OUTPUT);row=d['relations'][1769]
    assert [row['m'],row['n']]==saved['address']==[20941,464]
    assert row['alpha_ascending']==saved['alpha_ascending'] and out['norms_alpha'][1769]==saved['norm_alpha']
    beta=alg.mul(alg.elt([saved['norm_alpha'],0,0]),alg.elt(row['alpha_ascending']))
    inp=r.read(ref.INPUT)['cases'][0];model,pts=r.short(inp['model'],inp['generic_sections'])
    scale=Q(r.read(ref.OUTPUT)['stages']['factor']['scale']);b,a=map(Q,d['affine_map_old_root_from_masked'])
    assert QQ(a).is_square()
    assert f(QQ(a)*R.gen()+QQ(b))==QQ(a)**3*R(d['masked_cubic_ascending'])
    mask=saved['generic_product_mask'];value=beta
    for i,(x,y) in enumerate(pts):
        if mask>>i&1:value=alg.mul(value,alg.elt([a*Q(x)*scale**2+b,-1,0]))
    root=alg.elt(saved['square_root_ascending']);assert alg.mul(root,root)==value==alg.elt(saved['product_ascending'])
    assert mask==out['inherited_generic_mask']==22674 and out['additional_strict_class_capacity']==0
    files=(Path(__file__),source.INPUT,source.OUTPUT,source.PROVENANCE,hit.OUTPUT,ref.INPUT,ref.OUTPUT,
           Path(__file__).with_name('verify_unpointed_governing_norm.py'),Path(r.__file__))
    return {'schema':'rank-jump.early-relation-pool-verification.v1','status':'PASS',
        'rational_norm_checks':len(rows),'independent_ideal_valuation_checks':checks,
        'complete_prime_blocks':len(blocks),'outside_parity_rank':4133,'kernel_dimension':1,
        'verified_generic_square_identities':1,'additional_strict_capacity':0,
        'phase_rank_replays':len(out['phases']),
        'original_adaptive_search_replayed':False,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();value=compute()
    if args.mode=='build':r.write_new(OUTPUT,value)
    else:assert value==r.read(OUTPUT)
    print({k:v for k,v in value.items() if k!='bindings'},flush=True)
