#!/usr/bin/env python3
"""Replay equation/support bounds and distinguish arithmetic from joined rank labels."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
import sys
import retrospective as r
import fresh_strict_boundary_coordinates as source
import fresh_governing_panel as base
import matched103b2_class_boundary as prior
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_fresh_strict_boundary_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,AA,pari,GF,matrix
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    data=r.read(source.OUTPUT);report=r.read(source.REPORT)
    for artifact in (data,report):
        for path,sha in artifact['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,r.read(source.PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    rows=[];labels={x['token']:x for x in r.read(base.MANIFEST)['rows']}
    direct={x['token']:x for x in r.read(r.OUT/'rank_jump_fresh_strict_boundary_v1.json')['rows']}
    for row,out in zip(data['rows'],report['rows']):
        token=row['token'];assert out['token']==token
        if row['status']!='PASS':
            assert out['status']=='UNKNOWN';rows.append({'token':token,'status':'UNKNOWN','reason':row['reason']});continue
        old,f,pts,primes,nf=prior.setup(token);K=Algebra(list(map(str,f.list())))
        delta=Q(row['polynomial_discriminant']);beta=K.elt(row['derivative_coefficients'])
        assert delta==Q(str(f.discriminant())) and beta==K.elt([-delta*K.f[1],0,-3*delta])
        assert K.norm(beta)==delta**4==Q(row['derivative_norm'])
        th=pari.Mod('z',pari(f));bb=pari(f.parent()(list(map(QQ,row['derivative_coefficients']))))(th)
        E=pari.ellinit([0,0,0,pari(f[1]),pari(f[0])])
        S=[2]+[p for p in primes if p!=2 and int(pari.elllocalred(E,p)[0])>0]
        assert S==row['S_finite']
        omitted=[]
        for p in sorted(set(primes)-set(S)):
            vals=[int(pari.idealval(nf,bb,P)) for P in pari.idealprimedec(nf,p)]
            assert all(v%2==0 for v in vals);omitted.append({'prime':p,'valuations':vals})
        assert omitted==row['omitted_good_prime_valuations']
        roots=f.roots(AA,multiplicities=False);real_dim=int(len(roots)==3)
        signs=[int(-QQ(str(delta))*f.derivative()(a)<0) for a in roots]
        assert signs==row['real_derivative_signs'];witness=False
        if row['witness'] and row['witness']['place']=='infinity':
            assert real_dim==1 and signs==[1,0,1] and row['complete_real_point_basis']==[[0,1,1]]
            witness=True
        blocks=[];ell=real_dim;local_replays=0
        for loc in row['local']:
            p=loc['place'];saved=next(x for x in old['local']['local'] if x['place']==p)
            P=list(pari.idealprimedec(nf,p));d=len(P)-1+int(p==2);ell+=d
            assert d==loc['point_dimension'];sigs=saved['signatures'];rank=int(matrix(GF(2),sigs).rank())
            assert rank==loc['generic_local_rank'] and rank<=d
            blocks.extend(list(zip(*sigs)))
            if loc['tested_derivative']:
                chars=LocalSquareclasses(nf,p);bs=list(chars.signature(bb));assert r.pack(bs)==loc['derivative_signature']
                # Same coordinate backend, independent Sage rank calculation.
                assert [list(chars.signature(pari(x)-th)) for x,y in pts]==sigs
                outside=int(matrix(GF(2),sigs+[bs]).rank())>rank
                assert rank==d and outside==loc['outside_full_point_image'];local_replays+=1
                if row['witness'] and row['witness']['place']==p:
                    assert outside;witness=True
        real_sigs=[[int(x<a) for a in roots] for x,y in pts];blocks.extend(list(zip(*real_sigs)))
        g=int(matrix(GF(2),blocks).rank());m=len(pts);k=m-g;b=int(witness)
        assert ell==row['local_point_product_dimension'] and g==row['generic_local_dimension']
        assert k==row['generic_strict_dimension'] and b==row['reciprocity_constraint_rank_lower_bound']
        h=ell-b;a=h-g;assert a>=0 and row['Selmer_boundary_dimension_interval']==[g,h]
        assert row['additional_boundary_capacity_upper_bound']==a
        y=labels[token];R=y['retained_rank_lower_bound'];assert m==y['generic_rank']
        cmin=max(k,R-h);extra=cmin-k
        assert out['necessary_total_rational_strict_dimension']==cmin
        assert out['necessary_additional_rational_strict_dimension']==extra==max(0,R-m-a)
        assert out['S_class_dimension_lower_bound_from_generic_classes']==k
        assert out['S_class_dimension_lower_bound_after_rank_label']==cmin
        # Direct pairing outcomes, where completed, are an independent witness route.
        drow=direct[token];agreement=None
        if drow['status']=='PASS':
            assert drow['Selmer_boundary_dimension_interval']==row['Selmer_boundary_dimension_interval']
            agreement=True
        rows.append({'token':token,'status':'PASS','generic_dimension':m,
            'local_point_product_dimension':ell,'generic_strict_dimension':k,
            'boundary_upper_bound':h,'additional_boundary_capacity':a,
            'rank_label_dependent_additional_strict_minimum':extra,
            'coordinate_place_replays':local_replays,'agrees_with_completed_direct_Hilbert_method':agreement})
    assert len(rows)==len(labels)==16
    files=(Path(__file__),source.OUTPUT,source.REPORT,base.MANIFEST,base.OUTPUT,
           Path(prior.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py'),
           r.OUT/'rank_jump_fresh_strict_boundary_v1.json',base.LOCAL)
    return {'schema':'rank-jump.fresh-strict-boundary-verification.v1','status':'PASS',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},'rows':rows,
        'scope':'Rational norm identities; independent PARI bad-prime coverage versus Sage capture; prime splitting dimensions; Sage matrix accounting. Coordinate witnesses reuse LocalSquareclasses, with direct Hilbert comparison where completed. Missing rows and additional CT remain UNKNOWN.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],result['rows'],flush=True)
