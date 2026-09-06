#!/usr/bin/env sage -python
"""Independent prime-ideal and Frobenius checks for the retained Artin matrix."""
import argparse
from sage.all import QQ,ZZ,PolynomialRing,pari
import retrospective as r
import strict_artin as art
import remaining_bad_primes as rem


def verify(index):
    data = r.read(art.OUTPUT)
    assert data["bindings"] == art.bindings()
    row = next(x for x in data["rows"] if x["case_index"]==index)
    assert row == art.calculate(index)
    old = r.read(rem.bad.INPUT)["cases"][index]
    source = rem.bad.cases()[index]
    _,points = r.short(source["model"],source["generic_points"]+source["points"])
    selected = [points[i] for i in old["selected_input_indices"]]
    d = QQ(old["elliptic_scaling_d"])
    R = PolynomialRing(QQ,"z")
    pol = pari(R(list(map(QQ,old["integral_cubic_ascending"]))))
    primes = sorted(set(x["prime"] for x in row["ideal_dictionary"]))
    nf = pari.nfinit([pol,primes])
    theta = pari.Mod("z",pol)
    betas = [pari(QQ(P[0])*d*d)-theta for P in selected]
    classes = []
    for mask in row["character_masks"]:
        beta = pari.Mod(1,pol)
        for i,b in enumerate(betas):
            if mask>>i&1:
                beta *= b
        denominator = pari.denominator(pari.nfalgtobasis(nf,beta))
        classes.append(beta*denominator**2)
    checked = 0
    for j,ideal in enumerate(row["ideal_dictionary"]):
        p,a = ideal["prime"],ideal["theta_residue"]
        assert ZZ(p).is_prime(proof=True) and p <= r.read(art.PROTOCOL)["limits"]["rational_prime_bound"]
        assert int(pari.poldisc(pol))%p != 0
        matching = [P for P in pari.idealprimedec(nf,p) if int(P[3])==1 and pari.idealval(nf,theta-a,P)>0]
        assert len(matching)==1
        P = matching[0]
        assert pari.idealnorm(nf,P)==p==ideal["ideal_norm"]
        for i,beta in enumerate(classes):
            assert int(pari.idealval(nf,beta,P))%2==0
            bit = int(pari.nfislocalpower(nf,P,beta,2)==0)
            assert bit==row["artin_matrix_rows"][i][j]
            checked += 1
    print("PASS",row["id"],len(row["ideal_dictionary"]),"labelled prime ideals;",
          checked,"independent Frobenius checks; matrix rank",row["artin_rank"],flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index",type=int,required=True)
    args = parser.parse_args()
    verify(args.index)
