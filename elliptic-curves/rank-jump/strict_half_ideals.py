#!/usr/bin/env python3
"""Factorization-free half-ideal certificates for strict cubic Kummer classes."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import strict_class_blocks as strict
import remaining_bad_primes as rem

PROTOCOL=Path(__file__).with_name("STRICT_HALF_IDEAL_PROTOCOL.json")
INPUT=r.OUT/"rank_jump_strict_half_ideal_inputs_v2.json"
OUTPUT=r.OUT/"rank_jump_strict_half_ideals_v2.json"
WORK=r.ROOT/"artifacts/local/rank-jump-strict-half-ideals-v2"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,strict.OUTPUT,rem.INPUT,rem.bad.INPUT,r.INPUT)}


def matrix_record(M):
    return [[str(M[i,j]) for j in range(3)] for i in range(3)]


def compute(index):
    from sage.all import QQ,ZZ,PolynomialRing,pari
    from sage.version import version
    old=r.read(rem.bad.INPUT)["cases"][index]
    block=r.read(strict.OUTPUT)["rows"][index]
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    assert factor["factorization_complete"] and block["all_bad_places_complete"]
    primes=[p for p,e in factor["factors"]]
    assert set(primes)<=set(block["tested_places"])
    R=PolynomialRing(QQ,"z")
    f=R(list(map(QQ,old["integral_cubic_ascending"])))
    disc=ZZ(f.discriminant())
    product=ZZ(1)
    for p,e in factor["factors"]:
        assert ZZ(p).is_prime(proof=True)
        product*=ZZ(p)**e
    assert product==abs(16*disc)
    nf=pari.nfinit([pari(f),primes])
    theta=pari.Mod("z",pari(f))
    source=rem.bad.cases()[index]
    _,allpoints=r.short(source["model"],source["generic_points"]+source["points"])
    scale=QQ(old["elliptic_scaling_d"])
    gammas=[];ideals=[];point_rows=[];numerators=[]
    for position,i in enumerate(old["selected_input_indices"]):
        x,y=map(QQ,allpoints[i])
        x*=scale**2;y*=scale**3
        d=ZZ(x.denominator()).sqrt()
        assert d in ZZ and x*d*d in ZZ and y*d**3 in ZZ
        d=ZZ(d);a=ZZ(x*d*d);b=ZZ(y*d**3)
        assert a.gcd(d)==1 and b.gcd(d)==1 and b!=0
        gamma=pari(a)-pari(d*d)*theta
        assert pari.nfeltnorm(nf,gamma)==b*b
        I=pari.idealadd(nf,pari(b),gamma)
        gammas.append(gamma);ideals.append(I);numerators.append(abs(b))
        point_rows.append({"position":position,"input_index":i,"a":str(a),"b":str(b),"d":str(d),
                           "gamma_coordinates":[str(pari.lift(gamma).polcoef(j)) for j in range(3)],
                           "gcd_ideal_hnf":matrix_record(I)})
    prime_ideals=[(p,j,P) for p in primes for j,P in enumerate(pari.idealprimedec(nf,p))]
    rows=[]
    masks=block["generic_strict_kernel_masks"]+block["relative_strict_lift_masks"]
    for mask in masks:
        gamma=pari.Mod(1,pari(f));I=pari.idealhnf(nf,1);norm=ZZ(1)
        for i in range(len(gammas)):
            if mask>>i&1:
                gamma*=gammas[i]
                I=pari.idealmul(nf,I,ideals[i])
                norm*=numerators[i]
        corrections=[];J=I
        for p,j,P in prime_ideals:
            vg=int(pari.idealval(nf,gamma,P))
            vi=int(pari.idealval(nf,I,P))
            assert vg%2==0
            e=vg//2-vi
            corrections.append({"prime":p,"prime_index":j,"gamma_valuation":vg,
                                "initial_ideal_valuation":vi,"correction_exponent":e})
            if e: J=pari.idealmul(nf,J,pari.idealpow(nf,P,e))
        J=pari.idealhnf(nf,J)
        principal=pari.idealhnf(nf,gamma)
        assert pari.idealpow(nf,J,2)==principal
        assert pari.idealnorm(nf,J)==norm
        rows.append({"point_mask":mask,"corrections":corrections,
                     "half_ideal_hnf":matrix_record(J),"half_ideal_norm":str(norm),
                     "principal_product_hnf":matrix_record(principal)})
    r1,r2=map(int,nf.nf_get_sign())
    assert r1+2*r2==3
    field_disc=ZZ(nf.disc())
    order_index=QQ(disc/field_disc).sqrt()
    assert order_index in ZZ
    return {"bindings":bindings(),"case_index":index,"id":old["id"],
            "software":{"sage":version,"pari":str(pari.version())},
            "integral_cubic_ascending":old["integral_cubic_ascending"],
            "maximal_order_basis":list(map(str,nf.nf_get_zk())),
            "field_discriminant":str(field_disc),"polynomial_order_index":str(order_index),
            "signature":[r1,r2],"generic_strict_count":len(block["generic_strict_kernel_masks"]),
            "points":point_rows,"half_ideals":rows}


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    rows=[]
    for i in r.read(PROTOCOL)["cases"]:
        path=WORK/f"case-{i}.json"
        if not path.exists():
            with (WORK/f"case-{i}.log").open("x") as log:
                try:
                    p=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                      "--index",str(i),"--destination",str(path)],
                                     cwd=r.ROOT,stdout=log,stderr=log,timeout=40)
                    failure=None if p.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:
                    failure="40-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"case_index":i,"status":"UNKNOWN","reason":failure})
        row=r.read(path);assert row["bindings"]==bindings()
        rows.append(row)
        print("checkpoint",i,row.get("status",len(row.get("half_ideals",[]))),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.strict-half-ideal-inputs.v2","bindings":bindings(),"cases":rows})


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    rows=[]
    for row in data["cases"]:
        assert row["bindings"]==bindings()
        if row.get("status")=="UNKNOWN":
            rows.append(row);continue
        k=len(row["half_ideals"]);g=row["generic_strict_count"]
        unit_bound=sum(row["signature"])-1
        assert row["signature"] in ([3,0],[1,1])
        rows.append({"case_index":row["case_index"],"id":row["id"],
                     "known_strict_dimension":k,"generic_strict_dimension":g,
                     "totally_positive_unit_squareclass_dimension_upper_bound":unit_bound,
                     "ordinary_half_ideal_span_dimension_lower_bound":max(0,k-unit_bound),
                     "generic_ordinary_half_ideal_span_dimension_lower_bound":max(0,g-unit_bound),
                     "relative_ordinary_half_ideal_span_dimension_lower_bound":max(0,k-g-unit_bound),
                     "field_discriminant":row["field_discriminant"],
                     "boundary":"Only the displayed ordinary ideal-class 2-torsion span, modulo its generic strict image. No specified independent subset, full class rank, localized S-class upper bound or new rational class."})
    result={"schema":"rank-jump.strict-half-ideals.v2","bindings":bindings(),
            "input_sha256":r.digest(INPUT.read_bytes()),"rows":rows}
    if check:
        assert r.read(OUTPUT)==result
        print("PASS unit-kernel and ideal-span bounds")
    else:r.write_new(OUTPUT,result)
    for row in rows:
        print(row.get("id",row["case_index"]),row.get("ordinary_half_ideal_span_dimension_lower_bound",row.get("status")),
              row.get("relative_ordinary_half_ideal_span_dimension_lower_bound"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","build","check","verify"))
    parser.add_argument("--index",type=int);parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="capture":capture()
    elif args.mode=="worker":r.write_new(args.destination,compute(args.index))
    elif args.mode=="verify":
        row=next(x for x in r.read(INPUT)["cases"] if x["case_index"]==args.index)
        assert row==compute(args.index)
        print("PASS exact ideal-square reconstruction",row["id"])
    else:build(args.mode=="check")
