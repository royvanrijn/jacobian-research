#!/usr/bin/env python3
"""Independent lattice, Jacobi and CT checks for the strict Sha elementary factor."""
from sage.all import QQ,ZZ,PolynomialRing,NumberField,matrix,vector,pari,GF
import retrospective as r
import local_collision as lc
import strict_Sha_Artin as experiment
import strict_deformation_solubility as strict
from verify_strict_half_ideals import check_lattice
from verify_half_ideal_artin import jacobi


def verify():
    data=r.read(experiment.INPUT);assert data["bindings"]==experiment.bindings()
    R=PolynomialRing(QQ,"z");f=R(list(map(QQ,data["cubic_ascending"])))
    K=NumberField(f,"t");theta=K.gen()
    nf=pari.nfinit([pari(f),data["S_finite"]])
    assert list(map(str,nf.nf_get_zk()))==data["maximal_order_basis"]
    basis=[K(R(x)) for x in data["maximal_order_basis"]]
    B=matrix(QQ,[list(x) for x in basis]).transpose();Binv=B.inverse()
    anchor=r.read(lc.INPUT)["anchor"]
    gammas=[]
    for p,old in zip(data["points"],anchor["known_points_on_short_model"]):
        a,b,d=(ZZ(p[k]) for k in ("a","b","d"))
        assert (a/d**2,b/d**3)==tuple(map(QQ,old))
        assert b*b==a**3+f[1]*a*d**4+f[0]*d**6
        gamma=a-d*d*theta
        assert list(gamma)==list(map(QQ,p["gamma_coordinates"]))
        gammas.append(gamma)
    classes=[];norms=[]
    for mask in data["point_masks"]:
        beta=K(1);norm=ZZ(1)
        for i,gamma in enumerate(gammas):
            if mask>>i&1:
                beta*=gamma;norm*=abs(ZZ(data["points"][i]["b"]))
        classes.append(beta);norms.append(norm)
    Ps=[(p,j,P) for p in data["S_finite"] for j,P in enumerate(pari.idealprimedec(nf,p))]
    for j,column in enumerate(data["ideals"]):
        J=matrix(QQ,column["half_ideal_hnf"])
        assert ZZ(column["half_ideal_norm"])==norms[j]
        check_lattice(K,basis,J,classes[j],norms[j])
        reduced=matrix(QQ,column["reduced_ideal_hnf"])
        a=K(list(map(QQ,column["principal_multiplier_coordinates"])))
        T=matrix(QQ,[list(Binv*vector(QQ,list(a*b))) for b in basis]).transpose()
        change=J.inverse()*T*reduced
        assert all(c in ZZ for c in change.list()) and abs(change.det())==1
        H=matrix(QQ,column["coprime_ideal_hnf"]);N=ZZ(column["norm"])
        product=pari(H)
        for stored,(p,k,P) in zip(column["removed_S_prime_ideals"],Ps):
            assert (stored["prime"],stored["prime_index"])==(p,k)
            e=stored["exponent"]
            assert e==int(pari.idealval(nf,pari(reduced),P))
            if e:product=pari.idealmul(nf,product,pari.idealpow(nf,P,e))
        assert product==pari(reduced)
        assert N==H.det()>0 and all(N%p for p in data["S_finite"])
        assert H[0,0]==N and H[1,1]==H[2,2]==1
        assert all(H[i,k]==0 for i in range(3) for k in range(i))
        residues=vector(QQ,[1,-H[0,1],-H[0,2]])
        for i in range(3):
            for k in range(3):
                coords=Binv*vector(QQ,list(basis[i]*basis[k]))
                assert all(c in ZZ for c in coords)
                assert ZZ(residues.dot_product(coords)-residues[i]*residues[k])%N==0
        for i,beta in enumerate(classes):
            coords=Binv*vector(QQ,list(beta))
            value=ZZ(residues.dot_product(coords))
            entry=column["evaluations"][i]
            # All 25 recorded entries took the unit branch; no exception was used.
            assert value.gcd(N)==1 and value%N==ZZ(entry["residue"])
            symbol=jacobi(int(value),int(N))
            assert symbol==entry["jacobi_symbol"] and (symbol==-1)==entry["bit"]
    result=r.read(experiment.OUTPUT)["result"]
    M=matrix(GF(2),result["Artin_matrix_rows"])
    assert M.rank()==3
    sub=M.matrix_from_rows_and_columns(result["selected_character_indices"],result["selected_half_ideal_indices"])
    assert sub.det()==1
    record=next(x for x in r.read(strict.OUTPUT)["single_deformations"] if x["u"]==-1)
    cross=record["CT_cross_report"]["cross_pairing_rows"]
    selected=[cross[i] for i in result["selected_character_indices"]]
    assert matrix(GF(2),[[(v>>j)&1 for j in range(18)] for v in selected]).rank()==3
    assert all(lc.lift(mask,selected)!=0 for mask in range(1,8))
    print("PASS five half-ideal squares and transports; 25 independent Jacobi entries;")
    print("PASS an elementary S-class factor of dimension three whose seven nonzero dual characters are CT-obstructed on E_-1")


if __name__=="__main__":
    verify()
