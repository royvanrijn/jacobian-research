#!/usr/bin/env python3
"""Independent exact lattice/ring and integer-Jacobi replay of half-ideal Artin data."""
import argparse
import retrospective as r
import half_ideal_artin as art
import half_ideal_artin_completion as complete
import strict_half_ideals as h
import remaining_bad_primes as rem


def jacobi(a,n):
    assert n>0 and n%2
    a%=n;sign=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5):sign=-sign
        a,n=n,a
        if a%4==n%4==3:sign=-sign
        a%=n
    return sign if n==1 else 0


def verify(index):
    from sage.all import QQ,ZZ,NumberField,PolynomialRing,matrix,vector,pari
    assert r.read(art.INPUT)["bindings"]==art.bindings()
    assert r.read(complete.INPUT)["bindings"]==complete.bindings()
    row=next(x for x in r.read(art.INPUT)["cases"] if x["case_index"]==index)
    source=next(x for x in r.read(h.INPUT)["cases"] if x["case_index"]==index)
    completed=next(x for x in r.read(complete.INPUT)["cases"] if x["case_index"]==index)
    assert completed==complete.compute(index)  # Each local repair also checks a separate point-signature implementation.
    repairs={(x["column"],x["character"]):x for x in completed["repairs"]}
    R=PolynomialRing(QQ,"z");f=R(list(map(QQ,source["integral_cubic_ascending"])))
    K=NumberField(f,"t")
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    product=ZZ(1)
    for p,e in factor["factors"]:
        assert ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert factor["factorization_complete"] and product==abs(16*f.discriminant())
    nf=pari.nfinit([pari(f),[p for p,e in factor["factors"]]])
    assert list(map(str,nf.nf_get_zk()))==source["maximal_order_basis"]
    basis=[K(R(x)) for x in source["maximal_order_basis"]]
    B=matrix(QQ,[list(x) for x in basis]).transpose();Binv=B.inverse()
    gamma_points=[K(list(map(QQ,P["gamma_coordinates"]))) for P in source["points"]]
    gammas=[]
    for mask in row["character_point_masks"]:
        value=K(1)
        for i,g in enumerate(gamma_points):
            if mask>>i&1:value*=g
        gammas.append(value)
    prime_ideals=[(p,j,P) for p,e in factor["factors"] for j,P in enumerate(pari.idealprimedec(nf,p))]
    checked=0
    for j,column in enumerate(row["columns"]):
        original=matrix(QQ,source["half_ideals"][j]["half_ideal_hnf"])
        reduced=matrix(QQ,column["reduced_ideal_hnf"])
        alpha=K(list(map(QQ,column["principal_multiplier_coordinates"])))
        T=matrix(QQ,[list(Binv*vector(QQ,list(alpha*b))) for b in basis]).transpose()
        change=original.inverse()*T*reduced
        assert all(c in ZZ for c in change.list()) and abs(change.det())==1
        H=matrix(QQ,column["coprime_ideal_hnf"]);N=ZZ(column["norm"])
        product=pari(H)
        assert len(column["removed_S_prime_ideals"])==len(prime_ideals)
        for stored,(p,k,P) in zip(column["removed_S_prime_ideals"],prime_ideals):
            assert (stored["prime"],stored["prime_index"])==(p,k)
            e=stored["exponent"]
            assert e==int(pari.idealval(nf,pari(reduced),P))
            if e:product=pari.idealmul(nf,product,pari.idealpow(nf,P,e))
        assert product==pari(reduced)
        assert N==H.det()>0 and all(N%p for p,e in factor["factors"])
        assert H[0,0]==N and H[1,1]==H[2,2]==1
        assert all(H[i,k]==0 for i in range(3) for k in range(i))
        residues=vector(QQ,[1,-H[0,1],-H[0,2]])
        # This proves that the cyclic lattice quotient is the ring Z/N.
        for a in range(3):
            for b in range(3):
                coords=Binv*vector(QQ,list(basis[a]*basis[b]))
                assert all(c in ZZ for c in coords)
                assert ZZ(residues.dot_product(coords)-residues[a]*residues[b])%N==0
        for i,gamma in enumerate(gammas):
            coordinates=Binv*vector(QQ,list(gamma))
            assert all(c in ZZ for c in coordinates)
            value=ZZ(residues.dot_product(coordinates))
            entry=column["evaluations"][i]
            if "artin_bit" in entry:
                assert value%N==ZZ(entry["residue"])
                symbol=jacobi(int(value),int(N))
                assert symbol in (-1,1)
                assert symbol==entry["jacobi_symbol"] and (symbol==-1)==entry["artin_bit"]
            else:
                repair=repairs[j,i];p=repair["prime"];e=repair["prime_exponent"]
                assert value.gcd(N)==ZZ(entry["gcd"])==p
                assert N.valuation(p)==e
                M=N//ZZ(p)**e
                assert M==ZZ(repair["cofactor"])
                assert value%M==ZZ(repair["cofactor_residue"])
                symbol=jacobi(int(value),int(M))
                assert symbol in (-1,1) and symbol==repair["cofactor_jacobi_symbol"]
                assert repair["artin_bit"]==int(symbol==-1)^((e%2)*repair["local_frobenius_bit"])
            checked+=1
    print("PASS",row["id"],len(row["columns"]),"exact ideal transports and cyclic rings;",
          checked,"Artin entries;",len(repairs),"independently checked local repairs")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--index",type=int,required=True)
    args=parser.parse_args()
    verify(args.index)
