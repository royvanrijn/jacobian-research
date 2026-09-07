#!/usr/bin/env sage -python
"""Replay local signatures and independently check strict-kernel local powers."""
import argparse
import sys
from sage.all import QQ,ZZ,AA,PolynomialRing,pari
import retrospective as r
import remaining_bad_primes as rem
import strict_class_blocks as blocks

sys.path.insert(0,str(rem.bad.LOCAL_SOURCE.parents[1]))
from research_runtime.local_kummer import LocalSquareclasses


def verify(index):
    case = r.read(rem.INPUT)["cases"][index]
    factor,record = case["factor"],case["local"]
    assert factor["bindings"] == record["bindings"] == rem.bindings()
    initial = rem.initial_factors(index)
    for key in initial:
        assert factor[key] == initial[key]
    assert all(ZZ(p).is_prime(proof=True) and e>0 for p,e in factor["factors"])
    product = ZZ(1)
    for p,e in factor["factors"]:
        product *= ZZ(p)**e
    disc = abs(ZZ(factor["model_discriminant"]))
    assert disc%product == 0
    if factor["factorization_complete"]:
        assert product == disc
    else:
        assert factor["factors"] == factor["trial_factors"]
    old = r.read(rem.bad.INPUT)["cases"][index]
    source = rem.bad.cases()[index]
    model,allpoints = r.short(source["model"],source["generic_points"]+source["points"])
    assert old["short_model"] == model
    points = [allpoints[i] for i in old["selected_input_indices"]]
    R = PolynomialRing(QQ,"z")
    pol = pari(R(list(map(QQ,old["integral_cubic_ascending"]))))
    primes = [p for p,e in factor["factors"] if p>13]
    nf = pari.nfinit([pol,primes]) if primes else None
    th = pari.Mod("z",pol)
    d = QQ(old["elliptic_scaling_d"])
    betas = [pari(QQ(P[0])*d*d)-th for P in points]
    assert record["local_order_basis"] == (list(map(str,nf.nf_get_zk())) if nf else [])
    assert [c["prime"] for c in record["local"]] == primes
    B,A,_,_ = map(QQ,old["integral_cubic_ascending"])
    assert A == QQ(model[3])*d**4 and B == QQ(model[4])*d**6
    E = pari.ellinit([0,0,0,A,B])
    for local in record["local"]:
        p = local["prime"]
        chars = LocalSquareclasses(nf,p)
        assert local["point_signature_rows"] == [list(map(int,chars.signature(b))) for b in betas]
        assert local["point_kummer_dimension"] == chars.point_kummer_dimension
        assert local["prime_decomposition"] == [{"ramification_index":int(P[2]),"residue_degree":int(P[3])} for P in chars.primes]
        red = pari.elllocalred(E,p)
        minimal = pari.ellchangecurve(E,red[2])
        expected = {"conductor_exponent":int(red[0]),"kodaira_code":int(red[1]),
                    "minimal_change":list(map(str,red[2])),"tamagawa_number":int(red[3]),
                    "minimal_discriminant_valuation":int(pari.valuation(minimal.disc(),p))}
        assert local["local_reduction"] == expected
    # Check that the complete discriminant support really includes every bad
    # place, including the small primes whose signatures were reused.
    all_local = old["local"]+r.read(rem.dy.INPUT)["cases"][index]["local"]+record["local"]
    for p,e in factor["factors"]:
        local = next(c for c in all_local if c["prime"] == p)
        red = pari.elllocalred(E,p)
        assert int(red[0]) == local["local_reduction"]["conductor_exponent"]
    summary = blocks.calculate(index)
    assert summary == r.read(blocks.OUTPUT)["rows"][index]
    # Cross-check every displayed strict-kernel generator by PARI's separate
    # local-power routine, not the signature implementation.
    places = [p for p in summary["tested_places"] if p!="infinity"]
    nf = pari.nfinit([pol,places])
    prime_ideals = [P for p in places for P in pari.idealprimedec(nf,p)]
    roots = R(list(map(QQ,old["integral_cubic_ascending"]))).roots(AA,multiplicities=False)
    checks = 0
    for mask in summary["witness_strict_kernel_masks"]:
        beta = pari.Mod(1,pol)
        for i,b in enumerate(betas):
            if mask>>i&1:
                beta *= b
        denominator = pari.denominator(pari.nfalgtobasis(nf,beta))
        beta *= denominator**2
        for P in prime_ideals:
            assert pari.nfislocalpower(nf,P,beta,2) == 1
            checks += 1
        for root in roots:
            assert sum(int(QQ(points[i][0])*d*d<root) for i in range(len(points)) if mask>>i&1)%2 == 0
    print("PASS",old["id"],len(primes),"new local primes;",checks,
          "independent local-power checks; complete",factor["factorization_complete"],flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index",type=int,required=True)
    args = parser.parse_args()
    verify(args.index)
