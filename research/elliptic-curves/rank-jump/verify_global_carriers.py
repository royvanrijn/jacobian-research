#!/usr/bin/env python3
"""Replay complete local proofs, auxiliary descents and the unresolved class audits."""
import argparse
from fractions import Fraction as F
from math import prod
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
OUT=r.OUT/'rank_jump_global_carrier_verification_v1.json'
MAIN=r.OUT/'rank_jump_global_pair_solubility_v1.json'
EXTRA=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
CLASS=r.OUT/'rank_jump_carrier_sha_class_v1.json'
ISOG=r.OUT/'rank_jump_carrier_isogeny_obstruction_v1.json'


def local_square(n,p):
    if not n:return True
    valuation=0
    while n%p==0:n//=p;valuation+=1
    return valuation%2==0 and (n%8==1 if p==2 else pow(n%p,(p-1)//2,p)==1)


def verify():
    from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,Jacobian,pari
    from sage.version import version
    pari.allocatemem(67108864,silent=True)
    sources=[r.read(p) for p in (MAIN,EXTRA,CLASS,ISOG)]
    for data in sources:
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'v');v=R.gen();checks=[]
    for row in sources[0]['rows']+sources[1]['rows']:
        g,lc,des=row['geometry'],row['local'],row['descent'];f,h=[R(q) for q in g['forms']]
        N,D,U=[R(g[k]) for k in ('parameter_numerator','parameter_denominator','conic_root_numerator')]
        t,u=map(QQ,g['single_conic_anchor']);assert u*u==f(t)
        assert N/D-t==(f.derivative()(t)-2*u*v)/(v*v-f[2])
        assert U*U==f[0]*D*D+f[1]*N*D+f[2]*N*N
        quartic=R(g['quartic_coefficients']);den=QQ(g['quartic_square_scaling'])
        assert quartic==den*den*(h[0]*D*D+h[1]*N*D+h[2]*N*N)
        assert quartic.degree()==4 and quartic.discriminant()!=0 and N.gcd(D).degree()==0
        e,d,c,b,a=[quartic[i] for i in range(5)]
        I,J=12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
        assert [str(I),str(J)]==g['invariants_I_J']
        E=EllipticCurve(QQ,list(map(QQ,g['minimal_Jacobian_model'])))
        E0=EllipticCurve(QQ,[-27*I,-27*J]);assert E.is_isomorphic(E0)
        S=PolynomialRing(QQ,names=('v','w'));vv,ww=S.gens()
        assert E.is_isomorphic(Jacobian(ww*ww-sum(quartic[i]*vv**i for i in range(5))))
        bad=abs(ZZ(2*f.discriminant()*h.discriminant()*f.resultant(h)))
        factors=[(ZZ(p),int(k)) for p,k in g['bad_support_factorization']]
        assert bad==ZZ(g['bad_reduction_integer'])==prod(p**k for p,k in factors)
        assert len({p for p,k in factors})==len(factors) and all(p.is_prime(proof=True) for p,k in factors)
        assert {str(p) for p,k in factors}=={x['prime'] for x in lc['finite_bad_places']}
        for x in lc['finite_bad_places']:
            T,Z=x['base_point'];assert ZZ(T).gcd(ZZ(Z))==1
            vals=[int(q[0]*Z*Z+q[1]*T*Z+q[2]*T*T) for q in (f,h)]
            assert list(map(str,vals))==x['exact_form_values']
            assert all(local_square(n,int(x['prime'])) for n in vals)
        T,Z=lc['real_witness'];assert all(q[0]*Z*Z+q[1]*T*Z+q[2]*T*T>0 for q in (f,h))
        assert lc['everywhere_locally_soluble']
        points=[[QQ(x) for x in P] for P in des['points']]
        assert all(E(P) for P in points)
        tors=[P for P in E.torsion_points() if P];assert len(tors)==1 and 2*tors[0]==E(0)
        short,pts=r.short(g['minimal_Jacobian_model'],des['points']+[list(map(str,tors[0][:2]))])
        blocks=[(p,r.roots_at(short[3],short[4],p)) for p in r.primes(503)
                if r.roots_at(short[3],short[4],p) is not None and len(r.roots_at(short[3],short[4],p))==3]
        signatures=[r.point_signature(short,P,blocks) for P in pts]
        finite_lower=r.rank(signatures)-1
        ans=pari.ellrank(pari.ellinit(E.a_invariants()),0,points)
        lo,hi,s=map(int,ans[:3]);assert [lo,hi,s]==[des['rank_lower_bound'],des['rank_upper_bound'],des['CT_Sha2_mod_2Sha4_dimension']]
        assert des['rational_2_torsion_dimension']==1 and des['full_2_Selmer_dimension']==hi+1+s
        assert lo==hi and des['Sha_2_dimension']==s
        assert des['global_carrier_solubility_proved']==(s==0)
        if row['id']!='complementary_AD':assert finite_lower==hi
        checks.append({'id':row['id'],'rank_interval_replayed':[lo,hi],'Sha_2_dimension':s,
                       'finite_Kummer_free_rank_lower_bound':finite_lower,
                       'complete_local_support':list(map(str,(p for p,k in factors))),
                       'global_solubility_proved_by_Sha_vanishing':s==0})
    # Audit the degree-four map and its cubic-algebra square identity independently.
    cl=sources[2];m=cl['mapping'];f,N,G=[R(m[k]) for k in ('quartic','N','G')]
    E=EllipticCurve(QQ,list(map(QQ,m['Jacobian_model'])))
    assert G*G==N**3+E.a4()*N*f*f+E.a6()*f**3
    assert f.gcd(N).degree()==0 and N.degree()==f.degree()==4
    Q=PolynomialRing(QQ,'theta');theta=Q.gen();A=Q.quotient(theta**3+E.a4()*theta+E.a6(),'theta');th=A.gen()
    V=PolynomialRing(A,'v');beta=A(Q(m['beta']))
    square=V([A(Q(c)) for c in m['square_polynomial_coefficients']])
    assert V(N)-th*V(f)==beta*square*square
    model=m['Jacobian_model'];blocks=[(p,r.roots_at(model[3],model[4],p)) for p in cl['proof_primes']]
    assert all(len(roots)==3 for p,roots in blocks)
    points=cl['transported_rational_points_including_torsion'];r.short(model,points)
    sigs=[r.point_signature(model,P,blocks) for P in points];bs=0;col=0
    for p,roots in blocks:
        for root in roots:
            bval=sum(r.mod(c,p)*pow(root,i,p) for i,c in enumerate(m['beta']))%p
            assert bval;bs|=int(pow(bval,(p-1)//2,p)==p-1)<<col;col+=1
    assert sigs==cl['rational_Kummer_fingerprints'] and bs==cl['beta_fingerprint']
    assert r.rank(sigs)==cl['rational_Kummer_image_dimension']==2
    assert r.rank(sigs+[bs])==cl['dimension_after_beta']==3
    assert not cl['rational_Kummer_space_complete'] and cl['carrier_global_solubility']=='UNKNOWN'
    # The scalar alternative also failed its predeclared completeness gate.
    iso=sources[3];a,b,bp=map(QQ,(iso['a'],iso['b'],iso['b_prime']));e=QQ(iso['Jacobian_root'])
    assert a==3*e and b==3*e*e+E.a4() and e**3+E.a4()*e+E.a6()==0 and bp==a*a-4*b
    xp=v+a+b/v;mult=1-b/v**2
    assert (v**3+a*v*v+b*v)*mult**2==xp**3-2*a*xp*xp+bp*xp
    Ep=EllipticCurve(QQ,list(map(QQ,iso['isogenous_model'])))
    assert list(Ep.a_invariants())==[0,-2*a,0,bp,0]
    assert all(Ep([QQ(x) for x in P]) for P in iso['isogenous_points_including_2torsion'])
    for item in (iso['original_scalar_classes_with_beta_last'],iso['isogenous_scalar_classes']):
        recomputed=[]
        for x in item['values']:
            val=F(x);signature=0
            for i,p in enumerate(item['primes']):
                assert val.numerator%p and val.denominator%p
                signature|=int(pow(r.mod(val,p),(p-1)//2,p)==p-1)<<i
            recomputed.append(signature)
        assert recomputed==item['fingerprints'] and r.rank(recomputed)==item['rank']
    assert iso['isogenous_scalar_classes']['rank']==1 and not iso['original_scalar_image_complete']
    assert iso['carrier_global_solubility']=='UNKNOWN'
    return {'schema':'rank-jump.global-carrier-verification.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (MAIN,EXTRA,CLASS,ISOG,Path(__file__),HERE/'retrospective.py')},
            'software':{'sage':version,'pari':str(pari('version()'))},'carrier_checks':checks,
            'specific_AD_class':'UNKNOWN: both attempted completeness gates remain unsatisfied',
            'boundary':'This replay proves the two global nonemptiness results and preserves the unresolved A,D torsor class. It does not certify rank jumps on new original fibres.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args()
    result=verify()
    if a.mode=='build':r.write_new(OUT,result)
    else:assert r.read(OUT)==result
    print(result['status']);print(result['carrier_checks']);print(result['specific_AD_class'])
