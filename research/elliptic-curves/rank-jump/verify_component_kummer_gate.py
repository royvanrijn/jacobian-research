#!/usr/bin/env python3
"""Independent universal identities and additive-class witness checks."""
import argparse
from pathlib import Path
import retrospective as r
import component_kummer_gate as ex


def verify():
    from sage.all import QQ,GF,PolynomialRing,matrix,EllipticCurve
    R=PolynomialRing(QQ,["A","B","a","u"]);A,B,a,u=R.gens()
    T=PolynomialRing(R,"theta");theta=T.gen();f=theta**3+A*theta+B
    gamma=1-u*theta;kappa=1+u*theta+u*u*(A+theta*theta);D=1+A*u*u+B*u**3
    assert (gamma*kappa-D)%f==0
    assert f.resultant(kappa)==D*D
    fa=a**3+A*a+B
    specialized_kappa=a*a+a*theta+A+theta*theta
    assert ((a-theta)*specialized_kappa-fa)%f==0
    assert (theta*theta*(A+theta*theta)+B*theta)%f==0
    alpha=theta+u*theta*theta
    columns=[T(1),alpha,alpha*alpha%f]
    change=matrix(R,[[p[i] for p in columns] for i in range(3)])
    assert change.det()==D
    infinity=matrix(R,[[p[i] for p in [T(1),theta*theta,theta**4%f]] for i in range(3)])
    assert infinity.det()==B
    # The known singleton section realizes the canonical nonconstant class.
    X=a+(A*a+B)*u*u;g=1-a*u
    assert (X-g*alpha-(a-theta)*kappa)%f==0
    # Smooth reduced points of the twisted fibre at a rational I0* branch
    # have X=s^2, Y=s^3, so all three Kummer residues are squares.
    S=PolynomialRing(QQ,"s");s=S.gen()
    assert (s**3)**2==(s*s)**3

    data=r.read(ex.OUTPUT);inputs=r.read(ex.INPUT)["cases"];verified=[]
    for row,input_row in zip(data["cases"],inputs):
        A,B=map(QQ,input_row["model"][3:]);E=EllipticCurve([0,0,0,A,B])
        points=[E(list(map(QQ,P))) for P in input_row["generic_points"]]
        RQ=PolynomialRing(QQ,"t");t=RQ.gen();pol=t**3+A*t+B
        for sample in row["samples"]:
            p=sample["p"];F=GF(p);root=F(sample["root"])
            assert root**3+F(A)*root+F(B)==0
            residues=[F(P[0])-root for P in points]+[-F(B)*root]
            assert all(residues)
            assert [int(not x.is_square()) for x in residues]==sample["bits"]
        for j,control in zip((1,2),row["same_class_control"]["derived_points"]):
            target=points[0]+2*points[j]
            assert list(map(str,target[:2]))==control["point"]
            sigma=RQ(list(map(QQ,control["square_root_of_class_ratio"])))
            assert ((points[0][0]-t)*sigma*sigma-(target[0]-t))%pol==0
            line1=RQ(list(map(QQ,control["first_chord"])))
            line2=RQ(list(map(QQ,control["second_chord"])))
            assert (line1*sigma-line2)%pol==0
        verified.append({"id":row["id"],"prime_root_characters":len(row["samples"]),
                         "same_class_square_roots_checked":2})
    assert matrix(QQ,[[1,0,0],[1,2,0],[1,0,2]]).det()==4
    out={"schema":"rank-jump.component-kummer-gate-verification.v1","status":"PASS",
         "universal_identities":["gamma*kappa=D","Norm(kappa)=D^2","a^2*kappa(1/a)*(a-theta)=f(a)",
                                 "theta^2*kappa_leading=-B*theta","residue algebra determinants D and B",
                                 "singleton Kummer class=(a-theta)*kappa"],
         "cases":verified,"analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),
         "input_sha256":r.digest(ex.INPUT.read_bytes()),"verifier_sha256":r.digest(Path(__file__).read_bytes())}
    print("PASS universal identities, finite squareclass witnesses, and independent elliptic additions")
    print(verified)
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--destination",type=Path);args=p.parse_args()
    data=verify()
    if args.destination:r.write_new(args.destination,data)
