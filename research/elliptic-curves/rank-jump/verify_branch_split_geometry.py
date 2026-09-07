#!/usr/bin/env python3
"""Symbolic branch-cover equations and all eight twist fibre configurations."""
import itertools
import argparse
from pathlib import Path
import retrospective as r
import branch_split_picard as experiment
from cubic_bridge import Cubic


def verify():
    from sage.all import QQ,GF,PolynomialRing,matrix,prod,EllipticCurve,gcd
    R=PolynomialRing(QQ,["r","s"]);F=R.fraction_field()
    r0,s=map(F,R.gens());roots=[r0,s,-r0-s]
    A=sum(roots[i]*roots[j] for i in range(3) for j in range(i))
    B=-prod(roots)
    U=PolynomialRing(F,"u");u=U.gen()
    gamma=[1-t*u for t in roots]
    a2=2*A*u;a4=A+3*B*u+A*A*u*u;a6=B+A*B*u*u-B*B*u**3
    delta=-4*A**3-27*B*B;D=prod(gamma)
    rows=[]
    for mask in range(8):
        d=prod(gamma[i] for i in range(3) if mask>>i&1)
        c2,c4,c6=d*a2,d*d*a4,d**3*a6
        b2,b4,b6=4*c2,2*c4,4*c6;b8=4*c2*c6-c4*c4
        disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
        inv4=b2*b2-24*b4
        assert disc==16*delta*D*D*d**6
        finite=[]
        for i,t in enumerate(roots):
            root=1/t;power=8 if mask>>i&1 else 2
            f=disc
            for _ in range(power):
                assert f(root)==0
                f=f//(u-root)
            assert f(root)!=0
            assert (inv4(root)==0)==bool(mask>>i&1)
            finite.append("I2*" if mask>>i&1 else "I2")
        chi=1 if mask.bit_count()<=1 else 2
        infinity=12*chi-disc.degree()
        assert infinity==(6 if mask.bit_count()%2==0 else 0)
        root_rank=sum(6 if x=="I2*" else 1 for x in finite)+(4 if infinity else 0)
        rows.append({"mask":mask,"finite":finite,"infinity":"I0*" if infinity else "I0",
                     "chi":chi,"trivial_rank":2+root_rank})

    C=PolynomialRing(QQ,["A","B","a","b","c","h"])
    AA,BB,a,b,c,h=C.gens();T=PolynomialRing(C,"theta");theta=T.gen()
    z=a+b*theta+c*theta*theta
    rem=z*z%(theta**3+AA*theta+BB)
    assert rem[0]==a*a-2*BB*b*c
    assert rem[1]==2*a*b-2*AA*b*c-BB*c*c
    assert rem[2]==b*b+2*a*c-AA*c*c
    assert matrix(QQ,[[2,0,0,-2],[0,0,2,0]]).rank()==2
    # Three private zeroes and one common pole: four inertia vectors, rank three.
    inertia=matrix(GF(2),[[1,0,0],[0,1,0],[0,0,1],[1,1,1]])
    assert inertia.rank()==3
    genus=1+(-2*8+4*4)//2
    assert genus==1
    cert=r.read(experiment.OUTPUT)
    assert cert["generic_geometric_Picard_rank"]==19
    ranks=[(10 if row["chi"]==1 else (19 if row["mask"].bit_count()==2 else 20))
           -row["trivial_rank"] for row in rows]
    assert ranks==[1,0,0,0,0,0,0,0]
    anchor=r.read(experiment.lc.INPUT)["anchor"]
    aa,bb=map(QQ,anchor["short_model_ainvariants"][3:])
    rational_a,rational_b=map(lambda x:r.F(str(x)),(aa,bb))
    K=Cubic(rational_a,rational_b);parameter=4*bb/(aa*aa)
    z=K.add(K.one,K.scale(K.square(K.theta),2/rational_a))
    gamma=K.sub(K.one,K.scale(K.theta,r.F(str(parameter))))
    assert K.square(z)==gamma
    assert K.norm(z)==1+8*rational_b*rational_b/(rational_a**3)
    norm_z=QQ(str(K.norm(z)))
    assert norm_z!=0 and 1+aa*parameter**2+bb*parameter**3==norm_z**2
    auxiliary=EllipticCurve([0,aa,0,0,bb*bb]);Q=auxiliary([0,bb])
    assert list((3*Q)[:2])==[4*bb*bb/(aa*aa),-bb-8*bb**3/(aa**3)]
    assert (3*Q)[1]!=0
    orders=[]
    for p in (23,59):
        Ep=auxiliary.change_ring(GF(p));order=int(Ep([0,GF(p)(bb)]).order())
        orders.append({"p":p,"point_order":order})
    assert gcd([x["point_order"] for x in orders])==6
    print("PASS all eight character fibre configurations, genus-one splitting cover, rank sum",sum(ranks))
    return {"schema":"rank-jump.branch-split-geometry.v1","character_geometry":rows,
            "character_ranks":ranks,"cover_genus":genus,
            "split_control":{"parameter":str(parameter),"square_root_coordinates":list(map(str,z)),
                             "norm_square_root":str(K.norm(z)),"D":str(K.norm(z)**2),
                             "specialized_curve_rank":"UNKNOWN"},
            "auxiliary_point_reduction_orders":orders,"auxiliary_point_is_nontorsion":True,
            "third_multiple_coordinates":list(map(str,(3*Q)[:2])),
            "input_sha256":r.digest(experiment.lc.INPUT.read_bytes()),
            "picard_certificate_sha256":r.digest(experiment.OUTPUT.read_bytes()),
            "verifier_sha256":r.digest(Path(__file__).read_bytes())}


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--destination",type=Path)
    args=p.parse_args();data=verify()
    if args.destination:r.write_new(args.destination,data)
