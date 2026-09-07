#!/usr/bin/env python3
"""Exact geometry, character counts, isogenies and bounded descent replays."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,EllipticCurve,Jacobian,matrix,pari
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_native_genus_five_lift_inputs_v1.json'
RESULT=r.OUT/'rank_jump_native_genus_five_lift_v1.json'
PAIRS=r.OUT/'rank_jump_native_pair_factor_descent_v1.json'
OLD=r.OUT/'rank_jump_minimal_native_block_carrier_v1.json'
OUTPUT=r.OUT/'rank_jump_native_genus_five_lift_verification_v1.json'


def verify():
    inp=r.read(INPUT);res=r.read(RESULT);pairs=r.read(PAIRS);old=r.read(OLD)
    for data in (inp,res,pairs,old):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');h,f,g=qs=[R(c['form']) for c in inp['covers']]
    assert all(q.degree()==2 and q.is_squarefree() for q in qs)
    assert all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
    t0=QQ(inp['retained_lift']['t']);roots=list(map(QQ,inp['retained_lift']['roots']))
    assert all(q(t0)==u*u and u!=0 for q,u in zip(qs,roots,strict=True))
    geo=res['geometry'];assert [geo[k] for k in ('pair_carrier_degree','pair_carrier_genus','triple_carrier_degree','triple_carrier_genus','relative_cover_degree','geometric_ramification_points','Prym_dimension')]==[4,1,8,5,2,8,4]
    # Independent elimination: norm of z-t-u-v, first over u,v and then over t.
    Z=PolynomialRing(QQ,'z');z=Z.gen();T=PolynomialRing(Z,'t');t=T.gen()
    hh,ff,gg=[T(q.list()) for q in qs]
    norm_uv=((z-t)**2-ff-gg)**2-4*ff*gg
    norm_all=Z(hh.resultant(norm_uv));norm_all=norm_all.monic()
    cp=Z(geo['ramification_primitive_polynomial']);assert cp==norm_all and cp.degree()==8 and cp.is_squarefree()
    M=matrix(QQ,geo['ramification_multiplication_matrix']);powers=matrix(QQ,geo['ramification_primitive_power_matrix'])
    assert Z(M.charpoly().list())==cp and powers.rank()==8
    # Each norm is a nonsquare in Q, proving f,g,fg nonsquare in Q[t]/h.
    assert not h.discriminant().is_square()
    norms=[h.resultant(q)/h.leading_coefficient()**q.degree() for q in (f,g,f*g)]
    assert all(n and not n.is_square() for n in norms) and norms[2]==norms[0]*norms[1]
    # Thus this etale algebra is one degree-eight field, not a rational splitting.
    assert cp.is_irreducible()
    factors=[h*f,h*g,f*g,h*f*g]
    assert [int((q.degree()-2)//2) for q in factors]==[1,1,1,2]
    counts=[]
    for p in (131,137):
        F=GF(p);rp=PolynomialRing(F,'t');reduced=[rp(q) for q in qs]
        assert all(q.degree()==2 for q in reduced)
        assert rp(h*f*g).degree()==6 and rp(h*f*g).is_squarefree()
        def nr(x):return 1 if x==0 else 2 if x.is_square() else 0
        total=sum(nr(reduced[0](x))*nr(reduced[1](x))*nr(reduced[2](x)) for x in F)
        total+=nr(reduced[0][2])*nr(reduced[1][2])*nr(reduced[2][2])
        qcounts=[sum(nr(rp(q)(x)) for x in F)+nr(F(q.leading_coefficient())) for q in factors]
        assert total==sum(qcounts)-3*(p+1)
        counts.append({'prime':p,'triple_points':int(total),'character_quotient_points':list(map(int,qcounts))})
    # Rational origins parameterize the first conic on each pair carrier.
    pairrows=pairs['pairs']+[old];models=[];rankrows=[]
    pari.allocatemem(67108864,silent=True)
    RR=PolynomialRing(QQ,names=('t','y'));tt,yy=RR.gens()
    for i,indices in enumerate(((0,1),(0,2),(1,2))):
        row=pairrows[i];pg=row['geometry'];pd=row['descent'];product=res['elliptic_factors'][i]['descent']
        q1,q2=[qs[j] for j in indices]
        anchor=pg['anchor'];assert QQ(anchor['t'])==t0
        assert QQ(anchor['u'])**2==q1(t0) and QQ(anchor['v'])**2==q2(t0)
        N=R(pg['conic_parameter_numerator']);D=R(pg['conic_parameter_denominator']);U=R(pg['conic_root_numerator'])
        assert U*U==sum(q1[j]*N**j*D**(2-j) for j in range(3))
        quartic=R(pg['quartic_coefficients']);scale=QQ(pg['quartic_square_scale'])
        assert quartic==scale**2*sum(q2[j]*N**j*D**(2-j) for j in range(3))
        assert quartic.degree()==4 and quartic.is_squarefree()
        E=EllipticCurve(QQ,list(map(QQ,pg['minimal_Jacobian_model'])))
        assert Jacobian(yy**2-sum(quartic[j]*tt**j for j in range(5))).is_isomorphic(E)
        J=EllipticCurve(QQ,list(map(QQ,product['Jacobian_model'])))
        pol=q1*q2
        assert Jacobian(yy**2-sum(pol[j]*tt**j for j in range(5))).is_isomorphic(J)
        maps=[phi for phi in E.isogenies_prime_degree(2) if phi.codomain().is_isomorphic(J)]
        assert len(maps)==1
        for curve,desc in ((E,pd),(J,product)):
            for P in desc['points']:assert curve(list(map(QQ,P)))
            ans=pari.ellrank(pari.ellinit(curve.a_invariants()),0)
            assert list(map(int,ans[:3]))==[desc[k] for k in ('rank_lower_bound','rank_upper_bound','CT_Sha2_mod_2Sha4_dimension')]
            td={0:0,1:1,3:2}[len(curve.division_polynomial(2).roots(QQ,multiplicities=False))]
            assert td==desc['rational_2_torsion_dimension']
            assert desc['full_2_Selmer_dimension']==desc['rank_upper_bound']+td+desc['CT_Sha2_mod_2Sha4_dimension']
        rank=pd['rank_lower_bound'];assert rank==pd['rank_upper_bound']
        sha=product['full_2_Selmer_dimension']-product['rational_2_torsion_dimension']-rank
        rankrows.append({'labels':[inp['covers'][j]['label'] for j in indices],'exact_rank':rank,
            'pair_Jacobian_model':pg['minimal_Jacobian_model'],'product_Jacobian_model':product['Jacobian_model'],
            'isogeny_kernel_polynomial':list(map(str,maps[0].kernel_polynomial().list())),
            'pair_Sha2_dimension':pd['full_2_Selmer_dimension']-pd['rational_2_torsion_dimension']-rank,
            'product_Sha2_dimension':sha,'product_2Sha4_dimension':sha-product['CT_Sha2_mod_2Sha4_dimension']})
    assert [row['exact_rank'] for row in rankrows]==[3,2,3]
    paths=[INPUT,RESULT,PAIRS,OLD,Path(__file__),HERE/'retrospective.py']
    return {'schema':'rank-jump.native-genus-five-lift-verification.v1','status':'PASS',
        'ramification_algebra_field_degree':8,'h_discriminant':str(h.discriminant()),
        'nonsquare_branch_norms_f_g_fg':list(map(str,norms)),
        'independent_resultant_equals_primitive_polynomial':True,'finite_field_character_checks':counts,
        'elliptic_factors':rankrows,'Jacobian_rank_lower_bound':8,'genus_two_Jacobian_rank':'UNKNOWN',
        'ordinary_rank_below_genus_Chabauty_gate':'FAIL: rank at least 8, genus 5',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'boundary':'Isogenies and branch algebra checked exactly. Rank replays use the same PARI descent implementation. No genus-two descent, complete rational-point enumeration, or original-fibre rank prediction.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    print('PASS: branch field degree8; elliptic ranks3,2,3; genus5 Jacobian rank>=8')
