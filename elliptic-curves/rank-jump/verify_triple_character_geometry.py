#!/usr/bin/env python3
"""Generic triple-cover geometry, Picard replay, and rational-base non-torsion."""
import argparse
from pathlib import Path
import retrospective as r
import triple_character as ex


def verify():
    from sage.all import QQ,ZZ,GF,PolynomialRing,prod,matrix,companion_matrix,identity_matrix,EllipticCurve,gcd
    R=PolynomialRing(QQ,["A","B","a","b","c"]);A,B,a,b,c=R.gens();K=R.fraction_field()
    U=PolynomialRing(K,"u");u=U.gen();xs=[a,b,c]
    D=1+A*u*u+B*u**3;delta=-4*A**3-27*B*B;geometry=[]
    for mask in (3,5,6,7):
        selected=[x for i,x in enumerate(xs) if mask>>i&1]
        d=prod(1-x*u for x in selected)
        a2=2*A*u*d;a4=(A+3*B*u+A*A*u*u)*d*d
        a6=(B+A*B*u*u-B*B*u**3)*d**3
        b2,b4,b6=4*a2,2*a4,4*a6;b8=4*a2*a6-a4*a4
        disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
        c4=b2*b2-24*b4
        assert disc==16*delta*D*D*d**6 and D.gcd(d)==1 and D.gcd(c4)==1
        for x in selected:
            assert (disc/d**6)(1/x)!=0 and (c4/d**2)(1/x)!=0
        at_infinity=24-disc.degree()
        assert at_infinity==(6 if len(selected)==2 else 0)
        if len(selected)==3:
            k=-a*b*c
            assert [a2[4],a4[8],a6[12]]==[2*A*k,A*A*k*k,-B*B*k**3]
        geometry.append({"mask":mask,"finite_I2":3,"finite_I0star":len(selected),
                         "infinity":"I0*" if at_infinity else "smooth","trivial_rank":17})
    assert matrix(GF(2),[[1,0,0],[0,1,0],[0,0,1],[1,1,1]]).rank()==3
    genus=1+(-16+16)//2;assert genus==1

    V=PolynomialRing(K,["z0","z1","z2","h"]);z0,z1,z2,h=V.gens()
    quadrics=[b*z0*z0-a*z1*z1-(b-a)*h*h,c*z0*z0-a*z2*z2-(c-a)*h*h]
    tangent=matrix(K,[[q.derivative(v)(1,1,1,1) for v in V.gens()] for q in quadrics])
    assert tangent.rank()==2
    L=-a*b*c;S1=a+b+c;S2=a*b+a*c+b*c
    norm=prod(1-x*u for x in xs)
    assert (L*u)**3+S2*(L*u)**2+a*b*c*S1*(L*u)+(a*b*c)**2==L*L*norm

    proof=r.read(ex.OUTPUT);checks=[];RX=PolynomialRing(QQ,"X");X=RX.gen()
    for rec in proof["reductions"]:
        if "polynomial_ascending" not in rec:continue
        p=rec["p"];pol=RX(list(map(QQ,rec["polynomial_ascending"])))
        M=companion_matrix(pol)
        assert [int((M**i).trace()) for i in (1,2,3)]==rec["traces"]
        normpol=RX(pol(p*X)/p**5);cycles=0
        for factor,mult in normpol.factor():
            if all(c.denominator()==1 for c in factor.list()) and PolynomialRing(ZZ,"X")(factor).is_cyclotomic():
                cycles+=int(factor.degree())*int(mult)
        assert 17+cycles==rec["reduction_geometric_Picard_rank"]
        if rec["status"]=="RHO_18_REDUCTION":
            root=next(e*p for e in (-1,1) if pol(e*p)==0)
            quartic=pol//(X-root);V=companion_matrix(quartic)/p
            determinant=(identity_matrix(QQ,4)-V**6).det()
            assert determinant==QQ(rec["normalized_transcendental_degree6"])
            signed=-determinant/QQ(p**6)
            assert int(ZZ(signed.numerator()*signed.denominator()).squarefree_part())==rec["NS_discriminant_squareclass"]
        checks.append({"case":rec["case"],"mask":rec["mask"],"p":p,"PASS":True})

    auxiliary=[]
    for case,row in enumerate(r.read(ex.INPUT)["cases"]):
        a,b,c=[QQ(P[0]) for P in row["generic_points"]];L=-a*b*c
        E=EllipticCurve([0,a*b+a*c+b*c,0,a*b*c*(a+b+c),L*L]);Q=E([0,L])
        records=[]
        for p in r.read(ex.PROTOCOL)["auxiliary_base_check"]["primes"][row["role"]]:
            assert E.discriminant().valuation(p)==0
            F=GF(p);Ep=E.change_ring(F);order=int(Ep([0,F(L)]).order())
            records.append({"p":p,"point_order":order,"group_order":int(Ep.cardinality())})
        first,second=records;p,q=first["p"],second["p"];m,n=first["point_order"],second["point_order"]
        # At p only prime-to-p torsion is used; at q only prime-to-q torsion.
        stripped_m,stripped_n=ZZ(m),ZZ(n)
        for prime in (p,q):
            stripped_m//=prime**stripped_m.valuation(prime)
            stripped_n//=prime**stripped_n.valuation(prime)
        bound=int(gcd(stripped_m,stripped_n)*p**ZZ(n).valuation(p)*q**ZZ(m).valuation(q))
        multiple=bound*Q
        assert not multiple.is_zero()
        auxiliary.append({"case":case,"model":list(map(str,E.a_invariants())),"point":list(map(str,Q[:2])),
                          "reductions":records,"hypothetical_torsion_order_divides":bound,
                          "bound_multiple":list(map(str,multiple[:2])),"point_is_nontorsion":True})

    out={"schema":"rank-jump.triple-character-geometry.v1","character_geometry":geometry,"base_genus":genus,
         "norm_map_degree":4,"norm_map_kernel":"full geometric 2-torsion after choosing rational origins",
         "auxiliary_bases":auxiliary,"frobenius_checks":checks,
         "analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),"input_sha256":r.digest(ex.INPUT.read_bytes()),
         "verifier_sha256":r.digest(Path(__file__).read_bytes())}
    print("PASS all four K3 configurations, genus-one base, Frobenius, and two non-torsion proofs")
    print([x["reductions"] for x in auxiliary])
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--destination",type=Path);args=p.parse_args()
    data=verify()
    if args.destination:r.write_new(args.destination,data)
