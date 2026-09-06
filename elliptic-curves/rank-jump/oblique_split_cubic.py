#!/usr/bin/env python3
"""A split residual cubic gives two new directions beyond a horizontal norm obstruction."""
import argparse
from pathlib import Path
import retrospective as r

PROTOCOL=Path(__file__).with_name('OBLIQUE_SPLIT_CUBIC_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_oblique_split_cubic_v1.json'


def compute():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,EllipticCurve
    from sage.version import version
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field();X=PolynomialRing(K,'x');x=X.gen()
    g=(x+1)*(x-2)*(x-QQ(1)/2);f=g+t*x*x
    E=EllipticCurve(K,[0,t-QQ(3)/2,0,-QQ(3)/2,1])
    delta=R(E.discriminant());c4=R(E.c4())
    assert f==x**3+(t-QQ(3)/2)*x*x-QQ(3)/2*x+1
    assert delta.degree()==3 and delta.gcd(delta.derivative())==1 and delta.gcd(c4)==1
    assert c4.degree()==2 and delta(0)!=0
    # Finite 3 I1 and infinity (v(c4),v(Delta))=(2,9), hence I3*.
    parent_rank=10-2-7
    N=PolynomialRing(QQ,'n');n=N.gen();Kn=N.fraction_field()
    En=EllipticCurve(Kn,[0,n*n-QQ(3)/2,0,-QQ(3)/2,1])
    delta_n=N(En.discriminant());c4_n=N(En.c4())
    assert delta_n==delta(n*n) and c4_n==c4(n*n)
    assert delta_n.degree()==6 and delta_n.gcd(delta_n.derivative())==1 and delta_n.gcd(c4_n)==1
    assert c4_n.degree()==4
    # Finite 6 I1 and infinity (v(c4),v(Delta))=(0,6), hence I6.
    base_rank=10-2-5
    R0=En(0,1);P=En(-1,-n);Q=En(2,2*n);T=En(QQ(1)/2,n/2)
    assert P+Q+T==En(0)
    assert parent_rank==1 and base_rank==3
    n0=QQ(r.read(PROTOCOL)['fixed_anchor_n']);a2=n0*n0-QQ(3)/2
    E0=EllipticCurve(QQ,[0,a2,0,-QQ(3)/2,1])
    # Integral model x_int=4*x, y_int=8*y avoids denominator-two reductions.
    Ei=EllipticCurve(QQ,[0,4*a2,0,-24,64])
    points=[Ei(0,8),Ei(-4,-8*n0),Ei(8,16*n0)]
    Z=PolynomialRing(QQ,'z');z=Z.gen();pol=z**3+4*a2*z*z-24*z+64
    assert Z.change_ring(GF(7))(pol).is_irreducible()
    realroots=pol.roots(AA,multiplicities=False);signatures=[r.pack(int(P[0]<a) for a in realroots) for P in points]
    real_signatures=list(signatures);offset=len(realroots);records=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['largest_fingerprint_prime']):
        if pol.discriminant()%p==0:continue
        fp=PolynomialRing(GF(p),'z')(pol);roots=sorted(int(a) for a in fp.roots(multiplicities=False))
        if len(roots)!=3:continue
        local=[]
        for i,P in enumerate(points):
            bits=[]
            for a in roots:
                value=(int(P[0])-a)%p
                if not value:value=int(fp.derivative()(a))
                bits.append(int(pow(value,(p-1)//2,p)==p-1))
            sig=r.pack(bits);local.append(sig);signatures[i]|=sig<<offset
        records.append({'prime':p,'roots':roots,'signatures':local});offset+=3
    assert r.rank(signatures)==3
    minus_A=E0.c4()/48
    assert minus_A==QQ(2227)/12 and minus_A.valuation(17)==1 and ZZ(-3).kronecker(17)==-1
    # The common-cover relation is a special abscissa identity, not duplication.
    assert 2*En(-1,n)!=En(2,2*n)
    return {'schema':'rank-jump.oblique-split-cubic.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL)},
        'software':{'sage':version},'residual_cubic_coefficients':list(map(str,g.list())),
        'parent_discriminant_coefficients':list(map(str,delta)),'parent_c4_coefficients':list(map(str,c4)),
        'parent_fibres':['I1','I1','I1','I3*'],'parent_geometric_generic_rank':parent_rank,'parent_arithmetic_generic_rank':1,
        'base_discriminant_coefficients':list(map(str,delta_n)),'base_c4_coefficients':list(map(str,c4_n)),
        'base_fibres':['I1']*6+['I6'],'base_geometric_generic_rank':base_rank,'base_arithmetic_generic_rank':3,
        'new_generic_rational_directions':2,'line_relation':'(-1,-n)+(2,2*n)+(1/2,n/2)=O',
        'deck_action_on_R_P_Q':[1,-1,-1],
        'anchor':{'n':int(n0),'original_ainvariants':list(map(str,E0.ainvs())),
            'integral_ainvariants':list(map(str,Ei.ainvs())),
            'independent_points':[list(map(str,P)) for P in points],
            'real_signatures':real_signatures,'finite_fingerprints':records,'joint_fingerprints':signatures,
            'retained_subgroup_rank':3,'whole_curve_rank':'UNKNOWN here',
            'minus_short_A':str(minus_A),'horizontal_obstruction_prime':17,'valuation_minus_A':1},
        'boundary':'Exact generic ranks and a fixed independent subgroup, not a production rank result or a guarantee for every specialization.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS oblique split-cubic mechanism')
    else:r.write_new(OUTPUT,data);print('PASS generic rank 1 -> 3 despite horizontal norm obstruction')
