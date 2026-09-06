#!/usr/bin/env python3
"""Certify that a retained unramified norm-projection hit is inherited."""
import argparse
from pathlib import Path
import retrospective as r
import bounded_gain_reference as ref
import retained_norm_batch_capacity as batch
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_retained_norm_inherited_hit_v1.json'


def compute():
    from sage.all import QQ,ZZ,pari,PolynomialRing
    pari.allocatemem(64000000,268435456,silent=True)
    data=r.read(batch.INPUT);R=PolynomialRing(QQ,'z');f=R(data['cubic_ascending'])
    nf=pari.nfinit([pari(f),data['S_finite']]);th=pari.Mod('z',pari(f))
    # This address and mask came from explicitly exploratory retained-wave diagnostics.
    # No claim of a frozen prospective selection, or of whole-wave verification.
    m,n=20941,464;mask=88210
    alpha=pari(QQ(data['norm_generator']['fixed_a'])*m)+n*pari(R(data['norm_generator']['w_power_basis']))(th)
    N=pari.nfeltnorm(nf,alpha);beta=N*alpha
    ref.configure();fm,pts,scale=ref.base.model_data(ref.TOKEN)
    b,a=map(QQ,data['affine_map_old_root_from_masked'])
    assert f(a*R.gen()+b)==a**3*R(list(map(str,fm.list()))) and a.is_square()
    values=[pari(a*x+b)-th for x,y in pts]+[beta];product=pari.Mod(1,pari(f))
    for i,v in enumerate(values):
        if mask>>i&1:product*=v
    assert mask>>16==1
    # Extract the square root, retaining an independently checkable exact identity.
    roots=pari.nfroots(nf,pari('x')**2-product);assert len(roots)==2
    square_root=roots[0];assert square_root**2==product
    enc=lambda v:[str(pari.lift(v).polcoef(i)) for i in range(3)]
    algebra=Algebra(data['cubic_ascending']);value=algebra.elt([1,0,0])
    for i,v in enumerate(values):
        if mask>>i&1:value=algebra.mul(value,algebra.elt(enc(v)))
    assert algebra.mul(algebra.elt(enc(square_root)),algebra.elt(enc(square_root)))==value
    assert algebra.norm(algebra.elt(enc(alpha)))==r.F(str(N))
    files=(Path(__file__),batch.INPUT,ref.INPUT,Path(ref.__file__),Path(r.__file__),
           Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return {'schema':'rank-jump.retained-norm-inherited-hit.v1','status':'PASS',
        'address':[m,n],'alpha_ascending':enc(alpha),'norm_alpha':str(N),'norm_projection_ascending':enc(beta),
        'generic_product_mask':mask-(1<<16),'full_relation_mask':mask,
        'square_root_ascending':enc(square_root),'product_ascending':enc(product),
        'generic_dimension':16,'additional_squareclass_dimension':0,
        'selection':'Exploratory retained-wave diagnostic; the exact identity alone is certified here. No exhaustive early-wave replay.',
        'boundary':'The norm-projection class equals a product of generic Kummer classes. It is not an additional Selmer or rational direction. No exceptional point is an input and no elliptic point is searched for.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS inherited hit',result['address'],result['generic_product_mask'],flush=True)
