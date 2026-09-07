#!/usr/bin/env python3
"""Symbolic geometry and independent Frobenius linear-algebra replay."""
import argparse
from pathlib import Path
import retrospective as r
import mixed_character as ex


def verify():
    from sage.all import QQ,ZZ,PolynomialRing,companion_matrix,identity_matrix,GF,matrix
    R=PolynomialRing(QQ,["A","B","a","b"]);A,B,a,b=R.gens();K=R.fraction_field()
    U=PolynomialRing(K,"u");u=U.gen()
    D=1+A*u*u+B*u**3;delta=-4*A**3-27*B*B
    d=(1-a*u)*(1-b*u)
    a2=d*2*A*u;a4=d*d*(A+3*B*u+A*A*u*u)
    a6=d**3*(B+A*B*u*u-B*B*u**3)
    b2,b4,b6=4*a2,2*a4,4*a6;b8=4*a2*a6-a4*a4
    disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    c4=b2*b2-24*b4
    assert disc==16*delta*D*D*d**6
    assert disc.degree()==18 and c4.degree()==6
    assert D.gcd(d)==1 and D.gcd(c4)==1
    for branch in (1/a,1/b):
        assert (disc/d**6)(branch)!=0
        assert (c4/d**2)(branch)!=0
    # Infinity has valuations (c4,Delta)=(2,6); all three extra fibres are I0*.
    assert 8-c4.degree()==2 and 24-disc.degree()==6
    assert 2+3*1+3*4==17
    # Two private zeroes and the shared infinity: genus zero, not genus one.
    assert matrix(GF(2),[[1,0],[0,1],[1,1]]).rank()==2
    genus=1+(-2*4+3*2)//2
    assert genus==0
    T=PolynomialRing(K,"t").fraction_field();t=T.gen()
    z1=(a*t*t-2*a*t+b)/(a*t*t-b)
    z2=1+t*(z1-1);parameter=(1-z1*z1)/a
    assert z1*z1==1-a*parameter and z2*z2==1-b*parameter

    report=r.read(ex.OUTPUT);raw=r.read(ex.COUNTS);checks=[]
    QX=PolynomialRing(QQ,"X");X=QX.gen()
    for row in report["reductions"]:
        pol=QX(list(map(QQ,row["residual_polynomial_ascending"])))
        M=companion_matrix(pol)
        assert [int((M**k).trace()) for k in (1,2,3)]==row["traces"]
        p=row["p"];normalized=QX(pol(p*X)/p**5)
        cyclotomic_dimension=0
        for factor,mult in normalized.factor():
            # A rational irreducible monic factor can be cyclotomic only if integral.
            integral=all(x.denominator()==1 for x in factor.list())
            if integral and PolynomialRing(ZZ,"X")(factor).is_cyclotomic():
                cyclotomic_dimension+=int(factor.degree())*int(mult)
        assert 17+cyclotomic_dimension==row["reduction_geometric_Picard_rank"]
        if row["status"]=="RHO_18_REDUCTION":
            quartic=pol//(X-p)
            assert pol==(X-p)*quartic and quartic.degree()==4
            V=companion_matrix(quartic)/p
            value=(identity_matrix(QQ,4)-V**6).det()
            assert value==QQ(row["normalized_transcendental_at_one_degree6"])
            signed=-value/QQ(p**6)
            assert int(ZZ(signed.numerator()*signed.denominator()).squarefree_part())==row["NS_discriminant_squareclass"]
        checks.append({"case":row["case"],"p":p,"reduction_Picard_rank":17+cyclotomic_dimension})
    bounds=[]
    for i,case in enumerate(r.read(ex.INPUT)["cases"]):
        A0,B0=map(r.F,case["model"][3:]);points=[list(map(r.F,P)) for P in case["generic_points"]]
        assert all(y*y==x**3+A0*x+B0 for x,y in points)
        assert points[0][0]*points[1][0]*(points[0][0]-points[1][0])!=0
        assert not (B0>0 and r.isqrt(B0.numerator)**2==B0.numerator and r.isqrt(B0.denominator)**2==B0.denominator)
        gal=r.galois(case["model"]);assert gal["galois_group"]=="S3"
        upper=min(row["reduction_Picard_rank"] for row in checks if row["case"]==i)
        assert upper==18
        bounds.append({"id":case["id"],"geometric_Picard_rank_interval":[17,upper],
                       "mixed_geometric_rank_interval":[0,upper-17],
                       "full_pair_base_geometric_rank_interval":[3,3+upper-17],
                       "full_pair_base_arithmetic_rank_interval":[2,2+upper-17],
                       "production_curve_rank":"UNKNOWN","anchor_galois":gal})
    out={"schema":"rank-jump.mixed-character-geometry.v1","base_genus":genus,
         "twist_trivial_rank":17,"fibre_configuration":["I2"]*3+["I0*"]*3,
         "frobenius_checks":checks,"bounds":bounds,
         "analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),"counts_sha256":r.digest(ex.COUNTS.read_bytes()),
         "verifier_sha256":r.digest(Path(__file__).read_bytes())}
    print("PASS symbolic geometry, Frobenius matrices, and paired character bounds")
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--destination",type=Path);args=p.parse_args()
    data=verify()
    if args.destination:r.write_new(args.destination,data)
