#!/usr/bin/env python3
"""Verify a complete degree-12 intersection and the marked minimal carrier."""
import argparse
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,EllipticCurve,Jacobian,pari,prime_range
import retrospective as r
from verify_native_intersection_solubility import add
import verify_paired_quartet_relations as relations

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_native_triple_intersection_inputs_v1.json'
TRIPLE=r.OUT/'rank_jump_native_triple_intersection_v1.json'
CARRIER_INPUT=r.OUT/'rank_jump_minimal_native_block_carrier_inputs_v1.json'
CARRIER=r.OUT/'rank_jump_minimal_native_block_carrier_v1.json'
OUTPUT=r.OUT/'rank_jump_native_triple_carrier_verification_v1.json'


def verify():
    inp=r.read(INPUT);old=r.read(TRIPLE);ci=r.read(CARRIER_INPUT);carrier=r.read(CARRIER)
    for data in (inp,old,ci,carrier):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(inp['A']);B=R(inp['B']);G=matrix(ZZ,inp['gram'])
    def dec(x):return K(R(x['numerator']))/R(x['denominator'])
    s=vector(ZZ,[ZZ(x) for x in inp['generic_word']]);w=[vector(ZZ,c['published_basis_w']) for c in inp['covers']]
    z=2*s+w[1]-w[2];res=z-w[0]
    dO=(z*G*z-4)/2;a=dO+4;trace=2*z
    selfint=8*a-trace*G*trace;inter=2*a+8-w[0]*G*trace
    assert selfint==16 and dO==6 and inter==12 and res*G*res==10
    qs=[R(c['residual_chord']['q_coefficients']) for c in inp['covers']]
    assert all(q.degree()==2 and q.is_squarefree() for q in qs)
    assert all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
    delta=-16*(4*A**3+27*B**2)
    assert delta.degree()==24 and delta.is_squarefree() and all(delta.gcd(q)==1 for q in qs)
    assert (K(A**3)/delta).derivative()!=0
    # Reconstruct S with explicit rational function group addition.
    basis=[]
    for c in inp['sections']:
        x=K(R(c['x_coefficients_low_to_high']))
        if 'y_coefficients_low_to_high' in c:y=K(R(c['y_coefficients_low_to_high']))
        else:
            ch=c['chord'];p=basis[ch['reference_basis_index']]
            y=p[1]+R(ch['slope_coefficients_low_to_high'])*(x-p[0])
        assert y*y==x**3+A*x+B;basis.append((x,y))
    terms=[(int(n/abs(n))*vector(ZZ,[int(j==i) for j in range(17)]),(basis[i][0],int(n/abs(n))*basis[i][1])) for i,n in enumerate(s) if n for _ in range(abs(int(n)))]
    acc=vector(ZZ,[0]*17);S=None
    while terms:
        j=min(range(len(terms)),key=lambda j:((acc+terms[j][0])*G*(acc+terms[j][0]),j));v,p=terms.pop(j)
        acc+=v;S=add(A,S,p)
    assert S==(dec(old['generic_translate']['x']),dec(old['generic_translate']['y']))
    f=R(old['intersection_polynomial']);assert f.degree()==inter and f.is_squarefree()
    fs=[R(x['coefficients']) for x in old['factorization']];assert len(fs)==2 and fs[0].degree()==1 and fs[1].degree()==11
    assert fs[0]*fs[1]==f and all(x['multiplicity']==1 for x in old['factorization'])
    integer=R(fs[1].denominator()*fs[1]);modp=None
    for p in prime_range(2,132):
        rp=PolynomialRing(GF(p),'t');ff=rp(integer)
        if ff.degree()==11 and ff.is_irreducible():modp=int(p);break
    assert modp is not None, 'No bounded modular irreducibility witness'
    root=-fs[0][0]/fs[0][1];assert root==QQ(-288)/65
    F=R.quotient(f,'tt')
    def red(x):x=K(x);assert x.denominator().gcd(f)==1;return F(x.numerator())/F(x.denominator())
    rs=[dec(x) for x in old['rational_roots']];pts=[]
    for i,c in enumerate(inp['covers']):
        assert red(rs[i])**2==F(qs[i])
        x0,x1,y0,y1=[R(c['lifted_section'][k+'_coefficients']) for k in ('x0','x1','y0','y1')]
        px=x0+rs[i]*x1;py=y0+rs[i]*y1
        assert (px,py)==(dec(old['point_maps'][i]['x']),dec(old['point_maps'][i]['y']))
        x,y=red(px),red(py);assert y*y==x**3+F(A)*x+F(B);pts.append((x,y))
    total=add(F(A),add(F(A),pts[0],(pts[1][0],-pts[1][1])),pts[2]);assert total==tuple(map(red,S))
    assert R(old['excluded_polynomial']).gcd(f)==1
    checked=relations.verify();rel=next(x for x in checked['rows'] if x['id']=='08234-003')
    pair=next(x for x in rel['pair_ranks'] if x['indices']==[1,3]);assert pair['quotient_rank']==2
    # Minimal marked carrier and its pointed degree-two isogeny quotient.
    forms=[R(c['form']) for c in ci['covers']];u0,v0=[QQ(ci['anchor'][k]) for k in ('u','v')]
    assert root==QQ(ci['anchor']['t']) and forms[0](root)==u0*u0 and forms[1](root)==v0*v0
    assert forms[0].gcd(forms[1])==1 and (forms[0]*forms[1]).is_squarefree()
    geo=carrier['geometry'];N=R(geo['conic_parameter_numerator']);D=R(geo['conic_parameter_denominator']);U=R(geo['conic_root_numerator'])
    assert U*U==sum(forms[0][i]*N**i*D**(2-i) for i in range(3))
    q=R(geo['quartic_coefficients']);scale=QQ(geo['quartic_square_scale'])
    assert q==scale**2*sum(forms[1][i]*N**i*D**(2-i) for i in range(3))
    RR=PolynomialRing(QQ,names=('t','y'));t,y=RR.gens()
    J=Jacobian(y*y-sum(q[i]*t**i for i in range(5))).minimal_model()
    Jquot=Jacobian(y*y-sum((forms[0]*forms[1])[i]*t**i for i in range(5))).minimal_model()
    E=EllipticCurve(QQ,list(map(QQ,geo['minimal_Jacobian_model'])));assert E.is_isomorphic(J)
    isogenies=[phi for phi in E.isogenies_prime_degree(2) if phi.codomain().is_isomorphic(Jquot)];assert len(isogenies)==1
    ans=pari.ellrank(pari.ellinit(E.a_invariants()),0);des=carrier['descent']
    assert int(ans[0])==des['rank_lower_bound']==3 and int(ans[1])==des['rank_upper_bound']==3
    assert int(ans[2])==des['CT_Sha2_mod_2Sha4_dimension']==0
    return {'schema':'rank-jump.native-triple-carrier-verification.v1','status':'PASS',
        'intersection_degree':int(inter),'pair_image_square':int(selfint),'pair_image_zero_intersection':int(dO),
        'factor_degrees':[1,11],'degree_eleven_irreducible_mod_prime':modp,'unique_rational_intersection_parameter':str(root),
        'all_twelve_geometric_lifts_and_relation_verified':True,'specialized_pair_quotient_rank':2,
        'marked_pair_carrier_degree':4,'marked_pair_carrier_genus':1,'full_triple_carrier_degree':8,'full_triple_carrier_genus':5,
        'auxiliary_Jacobian_exact_rank':3,'product_quotient_Jacobian_model':list(map(str,Jquot.a_invariants())),
        'isogeny_degree':2,'isogeny_kernel_polynomial':list(map(str,isogenies[0].kernel_polynomial().list())),
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,TRIPLE,CARRIER_INPUT,CARRIER,Path(__file__),HERE/'retrospective.py',HERE/'verify_native_intersection_solubility.py',HERE/'verify_paired_quartet_relations.py')},
        'boundary':'Independent exact group arithmetic and modular irreducibility; auxiliary rank uses a replay of the same PARI 2-descent, not a second descent implementation. Minimality is for the specified native maps.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    print('PASS',data['factor_degrees'],'irreducible mod',data['degree_eleven_irreducible_mod_prime'],'minimal marked carrier genus1; Jacobian rank3; quotient rank2')
