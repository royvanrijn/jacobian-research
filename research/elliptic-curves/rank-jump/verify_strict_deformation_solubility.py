#!/usr/bin/env python3
"""Independent local-power replay of strict Sha blocks, with exact support checks."""
import argparse
import subprocess
from pathlib import Path

import retrospective as r
import local_collision as lc
import affine_selmer as af
import strict_deformation_solubility as experiment


def verify(u):
    from sage.all import AA,QQ,ZZ,PolynomialRing,GF,matrix,pari
    raw=r.read(af.INPUT)
    record=next(x for x in raw["cases"] if x["u"]==u)
    row=next(x for x in r.read(experiment.OUTPUT)["single_deformations"] if x["u"]==u)
    source=r.read(lc.INPUT);anchor=source["anchor"]
    r.short(anchor["short_model_ainvariants"],anchor["known_points_on_short_model"])
    R=PolynomialRing(QQ,"z");z=R.gen()
    f=R(list(map(QQ,anchor["base_polynomial_ascending"])))
    assert f.is_irreducible()
    product=ZZ(1)
    for entry in anchor["base_discriminant_factorization"]:
        p,e=entry["prime"],entry["exponent"]
        assert ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert product==abs(f.discriminant())
    A,B=map(QQ,anchor["short_model_ainvariants"][3:])
    D=1+A*u*u+B*u**3
    fu=z**3+2*A*u*z*z+(A+3*B*u+A*A*u*u)*z+B+A*B*u*u-B*B*u**3
    assert fu.discriminant()==D*D*f.discriminant()
    primes=[entry["place"] for entry in record["local"] if entry["place"]!="infinity"]
    assert 2 in primes
    for p in primes:assert ZZ(p).is_prime(proof=True)
    for value in (abs(ZZ(f.discriminant())),abs(ZZ(D))):
        for p in primes:
            while value%p==0:value//=p
        assert value==1
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod("z",pari(f))
    betas=[sum(pari(QQ(c))*theta**i for i,c in enumerate(beta)) for beta in record["anchor_beta_coordinates"]]
    K=row["strict_anchor_basis"]
    # Independent Sage binary kernel; no reuse of the repository elimination routine.
    equations=[r.pack(x) for local in record["local"] for x in zip(*local["class_signature_rows"][:20])]
    M=matrix(GF(2),[[(x>>j)&1 for j in range(20)] for x in equations])
    kernel=M.right_kernel()
    assert kernel.dimension()==len(K)
    assert all(vector in kernel for vector in matrix(GF(2),[[(mask>>j)&1 for j in range(20)] for mask in K]).rows())
    real_roots=f.roots(AA,multiplicities=False)
    checked=0
    for mask in K:
        beta=pari.Mod(1,pari(f))
        for i,b in enumerate(betas):
            if mask>>i&1:beta*=b
        denominator=pari.denominator(pari.nfalgtobasis(nf,beta))
        beta*=denominator**2
        for p in primes:
            for P in pari.idealprimedec(nf,p):
                assert pari.nfislocalpower(nf,P,beta,2)
                checked+=1
        bpoly=R(pari.lift(beta))
        assert all(bpoly(root)>0 for root in real_roots)
    original=next(x for x in source["rows"] if int(x["parameter_u"])==u)
    W=original["W_u_basis"]
    CT=next(x["matrix"] for x in source["ct"] if x["u"]==u)
    T=matrix(GF(2),CT)
    assert T==T.transpose() and all(T[i,i]==0 for i in range(T.nrows()))
    C=matrix(GF(2),[[(word>>j)&1 for j in range(len(W))]
                    for word in row["CT_cross_report"]["strict_basis_in_inherited_coordinates"]])
    assert (C*T).rank()==len(K)==row["CT_cross_report"]["cross_pairing_rank"]
    for i,mask in enumerate(K):
        assert lc.lift(r.pack(map(int,C.row(i))),W)==mask
    print("PASS u",u,"strict dimension",len(K),";",checked,
          "independent finite local-square checks; real positivity; complete support; injective CT cross-map")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--u",type=int)
    p.add_argument("--all",action="store_true")
    args=p.parse_args()
    if args.all:
        for u in r.read(experiment.PROTOCOL)["parameters"]:
            subprocess.run(["sage","-python",str(Path(__file__).resolve()),"--u",str(u)],
                           check=True,timeout=30,cwd=r.ROOT)
    else:
        assert args.u in r.read(experiment.PROTOCOL)["parameters"]
        verify(args.u)
