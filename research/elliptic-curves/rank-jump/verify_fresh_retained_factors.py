#!/usr/bin/env python3
"""Verify retained factor products, generic blocks and historic boundary supplement."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import sys
import retrospective as r
import fresh_retained_factors as source
import fresh_governing_panel as base
from verify_unpointed_governing_norm import Algebra
from verify_fresh_governing_panel import jacobi

OUTPUT=r.OUT/'rank_jump_fresh_retained_factor_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,AA,GF,matrix,pari
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    pari.allocatemem(64000000,r.read(source.PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    data=r.read(source.OUTPUT);report=r.read(source.REPORT);hints=r.read(source.HINTS)
    for artifact in (data,report,hints):
        for path,sha in artifact['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    for row in hints['rows']:assert r.digest((r.ROOT/row['source']).read_bytes())==row['source_sha256']
    results=[]
    for row,out in zip(data['rows'],report['rows']):
        token=row['token'];assert out['token']==token
        fac=row['factor'];f,pts,scale=base.model_data(token);D=int(16*f.discriminant())
        assert str(D)==fac['elliptic_discriminant']
        for p,e in fac['factors']:assert e>0 and ZZ(p).is_prime(proof=True)
        rem=int(fac['unresolved_cofactor']);assert prod(p**e for p,e in fac['factors'])*rem==abs(D)
        if fac['status']!='PASS':
            assert rem>1 and not ZZ(rem).is_prime(proof=True)
            assert row['local']['status']==row['boundary']['status']==out['status']=='UNKNOWN'
            results.append({'token':token,'status':'VERIFIED_INCOMPLETE','composite_cofactor_bits':rem.bit_length()});continue
        assert rem==1
        local=row['local'];b=row['boundary'];assert local['status']==b['status']=='PASS'
        nf=pari.nfinit([pari(f),[p for p,e in fac['factors']]]);th=pari.Mod('z',pari(f))
        E=pari.ellinit([0,0,0,pari(f[1]),pari(f[0])])
        S=[2]+[p for p,e in fac['factors'] if p!=2 and int(pari.elllocalred(E,p)[0])>0]
        assert S==local['S_finite']==b['S_finite']
        model=['0','0','0',str(f[1]),str(f[0])];blocks=local['independence_blocks']
        sigs=[r.point_signature(model,list(map(str,P)),blocks) for P in pts]
        assert sigs==local['independence_signatures'] and r.rank(sigs)==len(pts)
        columns=[];ell=0
        for place in local['local']:
            if place['place']=='infinity':
                roots=f.roots(AA,multiplicities=False);sigs=[[int(x<a) for a in roots] for x,y in pts];dim=int(len(roots)==3)
            else:
                chars=LocalSquareclasses(nf,place['place']);sigs=[list(chars.signature(pari(x)-th)) for x,y in pts]
                dim=len(pari.idealprimedec(nf,place['place']))-1+int(place['place']==2)
            assert sigs==place['signatures'] and dim==place['point_kummer_dimension']
            columns.extend(list(zip(*sigs)));ell+=dim
        kernel=matrix(GF(2),columns).right_kernel();k=int(kernel.dimension())
        assert k==local['strict_generic_dimension']==b['generic_strict_dimension']
        assert ell==b['local_point_product_dimension'];m=len(pts);g=m-k
        masks=local['generic_strict_masks'];assert len(masks)==k
        power_checks=0;ideal_checks=0;jacobi_checks=0
        for mask,record in zip(masks,local['class_records']):
            assert matrix(GF(2),[[mask>>i&1 for i in range(m)]]).row(0) in kernel
            beta=pari.Mod(1,pari(f))
            for i,(x,y) in enumerate(pts):
                if mask>>i&1:
                    d=ZZ(x.denominator()).sqrt();beta*=pari(ZZ(x*d*d))-pari(d*d)*th
            assert [str(pari.lift(beta).polcoef(i)) for i in range(3)]==record['beta_ascending']
            for p in S:
                for P in pari.idealprimedec(nf,p):
                    assert pari.nfislocalpower(nf,P,beta,2)==1;power_checks+=1
            for a in f.roots(AA,multiplicities=False):assert f.parent()(list(map(QQ,record['beta_ascending'])))(a)>0
            I=pari(matrix(QQ,record['half_ideal_hnf']));assert pari.idealpow(nf,I,2)==pari.idealhnf(nf,beta);ideal_checks+=1
        A=local['Artin_matrix'];M=local['minus_twist_CT_matrix']
        assert M==[[A[i][j]^A[j][i] for j in range(k)] for i in range(k)]
        assert r.rank(list(map(r.pack,M)))==local['minus_twist_CT_rank']==0
        for j,c in enumerate(local['artin_columns']):
            reduced=pari(matrix(QQ,c['reduced_ideal_hnf']));alpha=pari(f.parent()(list(map(QQ,c['multiplier_ascending']))))(th)
            I=pari(matrix(QQ,local['class_records'][j]['half_ideal_hnf']))
            assert pari.idealmul(nf,reduced,alpha)==I
            good=reduced
            for p in S:
                for P in pari.idealprimedec(nf,p):
                    e=int(pari.idealval(nf,good,P))
                    if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
            H=pari(matrix(QQ,c['coprime_ideal_hnf']));assert pari.idealhnf(nf,good)==H
            N=int(c['norm']);assert int(pari.idealnorm(nf,H))==N and all(N%p for p in S)
            assert c['cyclic'] and H[0,0]==N and H[1,1]==H[2,2]==1
            for i,e in enumerate(c['evaluations']):
                beta=pari(f.parent()(list(map(QQ,local['class_records'][i]['beta_ascending']))))(th)
                v=pari.nfalgtobasis(nf,beta);residue=int(v[0]-H[0,1]*v[1]-H[0,2]*v[2])%N
                assert residue==int(e['residue'])
                assert int(jacobi(residue,N)==-1)==e['artin_bit']==A[i][j];jacobi_checks+=1
        K=Algebra(list(map(str,f.list())));delta=Q(b['polynomial_discriminant']);derivative=K.elt(b['derivative_coefficients'])
        assert delta==Q(str(f.discriminant())) and derivative==K.elt([-delta*K.f[1],0,-3*delta])
        assert K.norm(derivative)==delta**4==Q(b['derivative_norm'])
        beta=pari(f.parent()(list(map(QQ,b['derivative_coefficients']))))(th)
        expected_omitted=[]
        for p in sorted(set(p for p,e in fac['factors'])-set(S)):
            vals=[int(pari.idealval(nf,beta,P)) for P in pari.idealprimedec(nf,p)]
            assert all(v%2==0 for v in vals);expected_omitted.append({'prime':p,'valuations':vals})
        assert expected_omitted==b['omitted_good_prime_valuations']
        witness=b['witness'];assert witness is not None
        if witness['place']=='infinity':
            assert [int(-QQ(str(delta))*f.derivative()(a)<0) for a in f.roots(AA,multiplicities=False)]==[1,0,1]
        else:
            p=witness['place'];chars=LocalSquareclasses(nf,p)
            saved=next(x for x in local['local'] if x['place']==p);sigs=saved['signatures'];ds=list(chars.signature(beta))
            rank=int(matrix(GF(2),sigs).rank());assert rank==chars.point_kummer_dimension
            assert matrix(GF(2),sigs+[ds]).rank()==rank+1
        h=ell-1;a=h-g;R=out['retained_rank_lower_bound']
        assert b['Selmer_boundary_dimension_interval']==[g,h] and b['additional_boundary_capacity_upper_bound']==a
        assert out['necessary_additional_strict_rational_dimension']==max(0,R-h-k)
        results.append({'token':token,'status':'PASS','proved_factor_count':len(fac['factors']),
            'generic_independence_dimension':m,'generic_strict_dimension':k,'local_product_dimension':ell,
            'additional_boundary_capacity':a,'rank_label_dependent_additional_strict_minimum':max(0,R-h-k),
            'local_power_checks':power_checks,'ideal_square_checks':ideal_checks,'Jacobi_checks':jacobi_checks})
    files=(Path(__file__),source.OUTPUT,source.REPORT,source.HINTS,base.INPUT,base.LOCAL,
           Path(__file__).with_name('verify_unpointed_governing_norm.py'),Path(__file__).with_name('verify_fresh_governing_panel.py'))
    return {'schema':'rank-jump.fresh-retained-factor-verification.v1','status':'PASS','rows':results,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'scope':'Prime products, composite residuals, generic local coordinates and good-prime independence, separate local-power tests, full half-ideal and Jacobi identities, derivative norm/support and boundary arithmetic. No exceptional quotient certificate.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],result['rows'],flush=True)
