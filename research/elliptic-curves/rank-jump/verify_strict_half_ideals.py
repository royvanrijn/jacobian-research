#!/usr/bin/env python3
"""Replay half-ideal identities by exact lattice products, without ideal arithmetic."""
import argparse
import copy

import retrospective as r
import strict_half_ideals as h
import strict_class_blocks as strict
import remaining_bad_primes as rem


def check_lattice(K, order_basis, H, gamma, expected_norm):
    from sage.all import QQ, ZZ, matrix, vector
    B=matrix(QQ,[list(x) for x in order_basis]).transpose()
    Binv=B.inverse()
    assert H.det()==expected_norm and expected_norm>0
    assert all(x in ZZ for x in H.list())
    lattice_basis=[K(list(B*H.column(i))) for i in range(3)]
    Hinv=H.inverse()
    # Multiplication by every integral-basis element preserves the lattice.
    for x in lattice_basis:
        for y in order_basis:
            assert all(c in ZZ for c in Hinv*Binv*vector(QQ,list(x*y)))
    # Six products generate J^2. Containment and equal index imply equality.
    for i in range(3):
        for j in range(i,3):
            assert all(c in ZZ for c in Binv*vector(QQ,list(lattice_basis[i]*lattice_basis[j]/gamma)))
    assert gamma.norm()==expected_norm**2


def verify(index, negative_controls=False):
    from sage.all import AA, QQ, ZZ, NumberField, PolynomialRing, matrix, pari
    row=next(x for x in r.read(h.INPUT)["cases"] if x["case_index"]==index)
    assert row["bindings"]==h.bindings()
    R=PolynomialRing(QQ,"z");f=R(list(map(QQ,row["integral_cubic_ascending"])))
    assert f.is_irreducible()
    K=NumberField(f,"t");theta=K.gen()
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    product=ZZ(1)
    for p,e in factor["factors"]:
        assert ZZ(p).is_prime(proof=True)
        product*=ZZ(p)**e
    assert factor["factorization_complete"] and product==abs(16*f.discriminant())
    nf=pari.nfinit([pari(f),[p for p,e in factor["factors"]]])
    assert list(map(str,nf.nf_get_zk()))==row["maximal_order_basis"]
    order_basis=[K(R(x)) for x in row["maximal_order_basis"]]
    assert matrix(QQ,3,3,lambda i,j:(order_basis[i]*order_basis[j]).trace()).det()==ZZ(row["field_discriminant"])
    assert f.discriminant()==ZZ(row["field_discriminant"])*ZZ(row["polynomial_order_index"])**2
    real_roots=f.roots(AA,multiplicities=False)
    assert [len(real_roots),(3-len(real_roots))//2]==row["signature"]
    old=r.read(rem.bad.INPUT)["cases"][index]
    source=rem.bad.cases()[index]
    _,points=r.short(source["model"],source["generic_points"]+source["points"])
    scale=QQ(old["elliptic_scaling_d"])
    gammas=[];signs=[];numerators=[]
    for position,p in enumerate(row["points"]):
        assert p["position"]==position
        i=old["selected_input_indices"][position]
        assert i==p["input_index"]
        a,b,d=(ZZ(p[k]) for k in ("a","b","d"))
        x,y=map(QQ,points[i])
        assert (a/d**2,b/d**3)==(x*scale**2,y*scale**3)
        assert b*b==a**3+QQ(f[1])*a*d**4+QQ(f[0])*d**6
        gamma=a-d*d*theta
        assert list(gamma)==list(map(QQ,p["gamma_coordinates"]))
        gammas.append(gamma);numerators.append(abs(b))
        signs.append([int(AA(a-d*d*t).sign()) for t in real_roots])
    block=r.read(strict.OUTPUT)["rows"][index]
    masks=block["generic_strict_kernel_masks"]+block["relative_strict_lift_masks"]
    assert r.rank(masks)==len(masks)
    assert [x["point_mask"] for x in row["half_ideals"]]==masks
    assert row["generic_strict_count"]==len(block["generic_strict_kernel_masks"])
    for position,certificate in enumerate(row["half_ideals"]):
        mask=certificate["point_mask"];gamma=K(1);norm=ZZ(1)
        positive=[1]*len(real_roots)
        for i in range(len(gammas)):
            if mask>>i&1:
                gamma*=gammas[i];norm*=numerators[i]
                positive=[s*t for s,t in zip(positive,signs[i])]
        assert positive==[1]*len(real_roots)
        H=matrix(QQ,certificate["half_ideal_hnf"])
        assert norm==ZZ(certificate["half_ideal_norm"])
        check_lattice(K,order_basis,H,gamma,norm)
        if negative_controls and position==0:
            for changed in (2*H,copy.copy(H)):
                if changed==H:
                    changed[0,2]+=1  # Preserve determinant, disturb containment.
                try:
                    check_lattice(K,order_basis,changed,gamma,norm)
                except AssertionError:
                    pass
                else:
                    raise AssertionError("altered half-ideal passed")
    print("PASS",row["id"],len(masks),"ideal squares; positivity; maximal-order and point transports",
          "; two altered ideals rejected" if negative_controls else "")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--negative-controls",action="store_true")
    args=p.parse_args()
    verify(args.index,args.negative_controls)
