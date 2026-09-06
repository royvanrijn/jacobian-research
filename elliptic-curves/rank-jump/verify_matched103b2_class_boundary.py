#!/usr/bin/env python3
"""Independent cubic identities and exact class/boundary rank accounting."""
import argparse
from itertools import permutations
from fractions import Fraction as Q
from pathlib import Path
import retrospective as r
import matched103b2_class_boundary as source
import fresh_governing_panel as base
from verify_unpointed_governing_norm import Algebra
import verify_unpointed_governing_norm as algebra_source

OUTPUT=r.OUT/'rank_jump_matched103b2_class_boundary_verification_v1.json'


def det3(A):
    return sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))*A[0][p[0]]*A[1][p[1]]*A[2][p[2]] for p in permutations(range(3)))


def compute():
    from sage.all import pari,QQ,AA,PolynomialRing,EllipticCurve,GF,matrix
    data=r.read(source.OUTPUT)
    for name,sha in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    rows=[];labels={x['token']:x for x in r.read(base.MANIFEST)['rows']}
    for row in data['rows']:
        token=row['token'];b=row['boundary'];red=row['reduction'];assert b['status']==red['status']=='PASS'
        f,pts,_=base.model_data(token);original=Algebra(list(map(str,f.list())))
        reduced=Algebra(red['reduced_cubic_ascending']);theta=reduced.elt(red['original_root_in_reduced']);eta=original.elt(red['reduced_root_in_original'])
        assert reduced.evaluate(original.f,theta)==reduced.elt([])
        assert original.evaluate(reduced.f,eta)==original.elt([])
        assert original.evaluate(theta,eta)==original.elt([0,1])
        assert reduced.evaluate(eta,theta)==reduced.elt([0,1])
        basis=[reduced.elt(v) for v in red['transported_maximal_order_basis']]
        trace=lambda v:sum(reduced.matrix(v)[i][i] for i in range(3))
        gram=[[trace(reduced.mul(a,c)) for a in basis] for c in basis]
        assert det3(gram)==Q(red['field_discriminant'])==Q(b['field_discriminant'])
        delta=Q(b['polynomial_discriminant']);beta=original.elt(b['derivative_coefficients'])
        assert beta==original.elt([-delta*original.f[1],0,-3*delta])
        assert original.norm(beta)==delta**4
        # Independent local-prime coverage audit and omitted-support test.
        old=next(x for x in r.read(base.OUTPUT)['rows'] if x['token']==token)
        primes=[p for p,e in old['factor']['factors']];nf=pari.nfinit([pari(f),primes]);th=pari.Mod('z',pari(f))
        bb=pari(f.parent()(list(map(QQ,b['derivative_coefficients']))))(th)
        E=EllipticCurve([0,0,0,f[1],f[0]])
        expected=[2]+[p for p in primes if p!=2 and int(pari.elllocalred(pari.ellinit(E.ainvs()),p)[0])>0]
        assert expected==b['S_finite']
        for p in set(primes)-set(expected):
            assert p!=2
            for P in pari.idealprimedec(nf,p):assert int(pari.idealval(nf,bb,P))%2==0
        roots=f.roots(AA,multiplicities=False)
        assert [int(-QQ(delta)*f.derivative()(a)<0) for a in roots]==[1,0,1]
        # E(R) norm-square point classes 000 and011, signs in ascending-root order.
        assert all([int(x<a) for a in roots] in ([0,0,0],[0,1,1]) for x,y in pts)
        ell=sum(p['point_dimension'] for p in b['local']);assert ell==b['local_point_product_dimension']
        sigs=[0]*len(pts);width=0
        for loc in b['local']:
            for i,v in enumerate(loc['generic_signatures']):sigs[i]|=v<<width
            width+=loc['width']
        rank=int(matrix(GF(2),[[x>>j&1 for j in range(width)] for x in sigs]).rank())
        assert rank==17==b['generic_local_dimension'] and sigs==b['generic_joint_signatures']
        upper=ell-1;defect=upper-rank;assert defect==b['additional_boundary_dimension_upper_bound']
        # The interval is sharp for abstract isotropic boundary spaces with the
        # displayed G and one independent dual constraint. No arithmetic realization claimed.
        abstract=[]
        def pairing(x,y):
            mask=(1<<ell)-1
            return ((x&mask)&(y>>ell)).bit_count()%2 ^ ((y&mask)&(x>>ell)).bit_count()%2
        for e in range(defect+1):
            G=[1<<i for i in range(rank)];L=G+[1<<(ell+rank)]
            L += [1<<j if j<rank+1+e else 1<<(ell+j) for j in range(rank+1,ell)]
            assert r.rank(L)==ell and all(pairing(x,y)==0 for x in L for y in L)
            intersection_dim=ell-r.rank([x>>ell for x in L]);assert intersection_dim==rank+e
            abstract.append({'e':e,'boundary_dimension':ell,'Selmer_local_intersection_dimension':intersection_dim})
        # Labels are joined only here, after the equation/generic computation.
        label=labels[token];retained=label['retained_rank_lower_bound'];strict_min=max(0,retained-upper)
        rows.append({'token':token,'id':label['id'],'family':label['family'],'parameter':label['parameter'],
          'generic_rank':17,'retained_rank_lower_bound':retained,'Selmer_local_image_interval':[17,upper],
          'quotient_boundary_capacity':defect,'class_group_dimension':'UNKNOWN',
          'retrospective_rational_strict_dimension_lower_bound':strict_min,
          'retrospective_S_class_two_rank_lower_bound':strict_min,
          'retrospective_unramified_extension_degree_over_K_lower_bound':2**strict_min,
          'bound_for_arbitrary_rank_R':f'c_S >= R - {upper}',
          'abstract_boundary_completions':abstract,
          'label_dependency':'Rank-derived necessities only; not independent class-group measurements or explicit class representatives.'})
        print('PASS',token,'Selmer=c_S+17+e; e<=',defect,'rank-derived strict>=',strict_min,flush=True)
    files=(Path(__file__),source.OUTPUT,base.INPUT,base.MANIFEST,base.OUTPUT,Path(base.__file__),Path(r.__file__),Path(algebra_source.__file__))
    return {'schema':'rank-jump.matched103b2-class-boundary-verification.v1','status':'PASS',
      'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},'rows':rows,
      'boundary':'No additional representatives or CT entries computed. The low rank remains censored; the high strict-block lower bound uses its existing rank certificate only after arithmetic.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
