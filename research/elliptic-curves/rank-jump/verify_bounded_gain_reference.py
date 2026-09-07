#!/usr/bin/env python3
"""Replay masked local/ideal/octic data and audit the separate conditional join."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import sys
import retrospective as r
import bounded_gain_reference as source
import bounded_gain_reference_completion as complete
import verify_explicit_governing_octic as octic
from verify_fresh_governing_panel import jacobi
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_bounded_gain_reference_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,pari,GF,matrix,AA
    source.configure();pari.allocatemem(64000000,268435456,silent=True)
    sys.path.insert(0,str(source.base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    raw=r.read(source.OUTPUT);data=r.read(complete.OUTPUT);report=r.read(source.REPORT)
    for value in (raw,data,report,r.read(source.PROVENANCE)):
        for path,sha in value['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    inputs=r.read(source.INPUT)
    assert set(inputs)=={'schema','cases','equation_factor_hints'}
    assert set(inputs['cases'][0])=={'token','model','generic_sections'}
    d=data['stages'];loc=d['local'];b=d['boundary'];fac=d['factor'];o=d['octic']
    f,pts,scale=source.base.model_data(source.TOKEN);m=len(pts);assert m==16
    assert prod(p**e for p,e in fac['factors'])==abs(16*f.discriminant())
    assert all(e>0 and ZZ(p).is_prime(proof=True) for p,e in fac['factors'])
    assert fac['unresolved_cofactor']=='1' and fac['status']=='PASS'
    primes=[p for p,e in fac['factors']];nf=pari.nfinit([pari(f),primes]);th=pari.Mod('z',pari(f))
    assert str(nf.disc())==loc['field_discriminant']
    E=pari.ellinit([0,0,0,pari(f[1]),pari(f[0])])
    S=[2]+[p for p in primes if p!=2 and pari.elllocalred(E,p)[0]>0]
    assert S==loc['S_finite']==b['S_finite']
    model=['0','0','0',str(f[1]),str(f[0])]
    assert r.galois(model)==loc['galois'] and loc['galois']['galois_group']=='S3'
    sigs=[r.point_signature(model,list(map(str,P)),loc['independence_blocks']) for P in pts]
    assert sigs==loc['independence_signatures'] and r.rank(sigs)==m
    roots=f.roots(AA,multiplicities=False);columns=[];ell=0;character_checks=0
    for place in loc['local']:
        p=place['place']
        if p=='infinity':
            sigs=[[int(x<a) for a in roots] for x,y in pts];dim=int(len(roots)==3)
        else:
            chars=LocalSquareclasses(nf,p);sigs=[list(chars.signature(pari(x)-th)) for x,y in pts]
            dim=len(pari.idealprimedec(nf,p))-1+int(p==2);character_checks+=m
        assert sigs==place['signatures'] and dim==place['point_kummer_dimension']
        columns.extend(zip(*sigs));ell+=dim
    kernel=matrix(GF(2),list(columns)).right_kernel();k=int(kernel.dimension())
    masks=loc['generic_strict_masks'];assert r.rank(masks)==len(masks)==k==6
    beta=[];ideals=[];power_checks=0
    for mask,rec in zip(masks,loc['class_records']):
        assert matrix(GF(2),[[mask>>i&1 for i in range(m)]]).row(0) in kernel
        value=pari.Mod(1,pari(f))
        for i,(x,y) in enumerate(pts):
            if mask>>i&1:
                den=ZZ(x.denominator()).sqrt();assert den in ZZ
                value*=pari(ZZ(x*den**2))-pari(den**2)*th
        assert rec['generic_mask']==mask
        assert [str(pari.lift(value).polcoef(i)) for i in range(3)]==rec['beta_ascending']
        for p in S:
            for P in pari.idealprimedec(nf,p):
                assert pari.nfislocalpower(nf,P,value,2)==1;power_checks+=1
        for a in roots:assert f.parent()(list(map(QQ,rec['beta_ascending'])))(a)>0
        I=pari(matrix(QQ,rec['half_ideal_hnf']));assert pari.idealpow(nf,I,2)==pari.idealhnf(nf,value)
        beta.append(value);ideals.append(I)
    A=loc['Artin_matrix'];M=loc['minus_twist_CT_matrix'];repairs={(v['row'],v['column']):v for v in data['repair']['repairs']}
    jacobi_checks=0;repair_checks=0
    for j,c in enumerate(loc['artin_columns']):
        reduced=pari(matrix(QQ,c['reduced_ideal_hnf']));alpha=pari(f.parent()(list(map(QQ,c['multiplier_ascending']))))(th)
        assert pari.idealmul(nf,reduced,alpha)==ideals[j]
        good=reduced
        for p in S:
            for P in pari.idealprimedec(nf,p):
                e=int(pari.idealval(nf,good,P))
                if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
        H=pari(matrix(QQ,c['coprime_ideal_hnf']));N=int(c['norm'])
        assert pari.idealhnf(nf,good)==H and int(pari.idealnorm(nf,H))==N
        assert c['cyclic'] and H[0,0]==N and H[1,1]==H[2,2]==1 and all(N%p for p in S)
        for i,entry in enumerate(c['evaluations']):
            v=pari.nfalgtobasis(nf,beta[i]);value=int(v[0]-H[0,1]*v[1]-H[0,2]*v[2])
            if 'artin_bit' in entry:
                assert value%N==int(entry['residue'])
                bit=int(jacobi(value,N)==-1);assert bit==entry['artin_bit'];jacobi_checks+=1
            else:
                rec=repairs[i,j];cofactor=N;bit=0
                for pl in rec['local']:
                    p=pl['prime'];assert ZZ(p).is_prime(proof=True) and f.discriminant()%p
                    assert ZZ(N).valuation(p)==pl['norm_exponent'];cofactor//=p**pl['norm_exponent']
                    contributions=[]
                    for P in pari.idealprimedec(nf,p):
                        e=int(pari.idealval(nf,H,P))
                        if not e:continue
                        assert int(P[3])==1
                        valuation=int(pari.idealval(nf,beta[i],P));assert valuation%2==0
                        frob=int(pari.nfislocalpower(nf,P,beta[i],2)==0)
                        tc=pari.nfalgtobasis(nf,th);root=int((tc[0]-H[0,1]*tc[1]-H[0,2]*tc[2])%p)
                        assert pari.idealval(nf,th-root,P)>0
                        point_bit=0
                        for n,Pnt in enumerate(pts):
                            if masks[i]>>n&1:point_bit^=r.point_signature(model,list(map(str,Pnt)),[(p,[root])])
                        assert point_bit==frob;bit^=(e%2)*frob
                        contributions.append({'exponent':e,'valuation':valuation,'frobenius_bit':frob,'independent_generic_character_bit':point_bit,'root':root})
                    assert contributions==pl['contributions']
                assert cofactor==int(rec['cofactor']) and value%cofactor==int(rec['residue'])
                jj=jacobi(value,cofactor);assert jj==rec['jacobi'] and jj in (-1,1)
                bit^=int(jj==-1);assert bit==rec['artin_bit'];repair_checks+=1
            assert bit==A[i][j]
    assert M==[[A[i][j]^A[j][i] for j in range(k)] for i in range(k)]
    assert matrix(GF(2),M).rank()==loc['minus_twist_CT_rank']==6
    alg=Algebra(list(map(str,f.list())));delta=Q(str(f.discriminant()));deriv=alg.elt(b['derivative_coefficients'])
    assert deriv==alg.elt([-delta*alg.f[1],0,-3*delta]) and alg.norm(deriv)==delta**4==Q(b['derivative_norm'])
    assert b['witness']['place']=='infinity'
    assert [int(-f.discriminant()*f.derivative()(a)<0) for a in roots]==b['real_derivative_signs']==[1,0,1]
    assert b['complete_real_point_basis']==[[0,1,1]] and ell==b['local_point_product_dimension']==11
    assert b['omitted_good_prime_valuations']==[] and set(primes)==set(S)
    assert b['Selmer_boundary_dimension_interval']==[m-k,ell-1]==[10,10]
    assert b['additional_boundary_capacity_upper_bound']==0
    x,y=map(Q,map(str,pts[0]));u,v=map(Q,map(str,pts[1]));c=u-x
    expect=[c**6,0,-4*c**3*(v-y),0,6*c*c*(x+u),0,-4*(y+v),0,1]
    assert list(map(Q,o['rational_octic_ascending']))==expect
    scale=Q(o['integral_octic_root_scale']);H=list(map(int,o['integral_octic_ascending']))
    assert H==[expect[i]*scale**(8-i) for i in range(9)]
    assert octic.discriminant(H)==int(o['integral_octic_discriminant'])
    transform={'cubic_ascending':list(map(str,f.list())),'integral_octic_ascending':o['integral_octic_ascending'],'scaled_points':[list(map(str,P)) for P in pts[:2]]}
    for e in o['inert_prime_table']:
        octic.prime_replay(transform,{'prime':e['prime'],'octic_factor_degrees':e['factor_degrees'],'radical_norm_mod_p':e['radical_norm'],'psi':e['psi'],'independent_radical_psi':e['psi']})
    assert o['governing_field_degree']==192 and o['joint_class_field_degree']==96
    # This audits only the source's stated label and hashes, not its extensive proof.
    label=r.read(source.LABEL);out=report['result']
    assert label['status']=='PASS' and label['conditional_on_grh_exact_rank']==label['conditional_on_grh_selmer_dimension']==22
    assert label['field_discriminant']==loc['field_discriminant']
    assert out['assumption']==label['assumption'] and out['generic_dimension']==m
    assert out['localized_class_dimension_conditional_interval']==[12,12]
    assert out['additional_strict_rational_dimension_conditional_interval']==[6,6]
    assert out['label_Sha2_dimension_conditional']==0 and out['inherited_minus_twist_CT_rank']==6
    assert out['additional_quotient_CT_computed_point_blind']=='UNKNOWN'
    files=(Path(__file__),source.OUTPUT,complete.OUTPUT,source.INPUT,source.REPORT,source.LABEL,Path(octic.__file__),Path(r.__file__),source.base.LOCAL)
    return {'schema':'rank-jump.bounded-gain-reference-verification.v1','status':'PASS',
        'generic_mod2_independence':m,'generic_strict_dimension':k,'local_character_checks':character_checks,
        'local_power_checks':power_checks,'half_ideal_checks':k,'direct_Jacobi_checks':jacobi_checks,
        'nonunit_Artin_repairs':repair_checks,'CT_switch_rank':6,'octic_prime_replays':len(o['inert_prime_table']),
        'boundary_dimension_exact':10,'conditional_localized_class_dimension':12,
        'conditional_label_full_proof_replayed':False,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result,flush=True)
