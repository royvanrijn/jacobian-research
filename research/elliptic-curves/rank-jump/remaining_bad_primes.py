#!/usr/bin/env python3
"""Bounded completion of bad-place quotient support on the frozen paired panel."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import bad_prime_support as bad
import dyadic_real_support as dy

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"REMAINING_BAD_PRIMES_PROTOCOL.json"
INPUT = r.OUT/"rank_jump_remaining_bad_prime_inputs_v1.json"
OUTPUT = r.OUT/"rank_jump_remaining_bad_primes_v1.json"
WORK = r.ROOT/"artifacts/local/rank-jump-remaining-bad-primes-v1"


def bindings():
    paths = [Path(__file__),PROTOCOL,bad.INPUT,dy.INPUT,r.INPUT,bad.LOCAL_SOURCE,
             Path(bad.__file__),Path(dy.__file__)]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def initial_factors(index):
    old = r.read(bad.INPUT)["cases"][index]
    B,A,_,_ = map(int,old["integral_cubic_ascending"])
    disc = -16*(4*A**3+27*B*B)
    remainder = abs(disc)
    factors = []
    primes = [2]+r.primes(r.read(PROTOCOL)["limits"]["trial_prime_bound"])
    for p in primes:
        e = 0
        while remainder%p == 0:
            e += 1
            remainder //= p
        if e:
            factors.append([p,e])
    return {"bindings":bindings(),"case_index":index,"id":old["id"],
            "model_discriminant":str(disc),"trial_factors":factors,
            "remaining_cofactor":str(remainder)}


def factor_worker(index,path):
    from sage.all import ZZ
    row = r.read(WORK/f"trial-{index}.json")
    assert row["bindings"] == bindings()
    factors = row["trial_factors"] + [[int(p),int(e)] for p,e in ZZ(row["remaining_cofactor"]).factor(proof=True)]
    assert all(ZZ(p).is_prime(proof=True) for p,e in factors)
    r.write_new(path,{**row,"factorization_complete":True,"factors":sorted(factors)})


def local_worker(index,path):
    from sage.all import QQ,PolynomialRing,pari
    from sage.version import version
    sys.path.insert(0,str(bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old = r.read(bad.INPUT)["cases"][index]
    factor = r.read(WORK/f"factor-{index}.json")
    primes = [p for p,e in factor["factors"] if p>13]
    row = bad.cases()[index]
    model,allpoints = r.short(row["model"],row["generic_points"]+row["points"])
    points = [allpoints[i] for i in old["selected_input_indices"]]
    R = PolynomialRing(QQ,"z")
    pol = pari(R(list(map(QQ,old["integral_cubic_ascending"]))))
    nf = pari.nfinit([pol,primes]) if primes else None
    th = pari.Mod("z",pol)
    d = QQ(old["elliptic_scaling_d"])
    betas = [pari(QQ(P[0])*d*d)-th for P in points]
    B,A,_,_ = map(QQ,old["integral_cubic_ascending"])
    E = pari.ellinit([0,0,0,A,B])
    rows = []
    for p in primes:
        print("LOCAL",index,p,flush=True)
        chars = LocalSquareclasses(nf,p)
        sigs = [list(map(int,chars.signature(b))) for b in betas]
        red = pari.elllocalred(E,p)
        minimal = pari.ellchangecurve(E,red[2])
        rows.append({"prime":p,"point_signature_rows":sigs,
                     "point_kummer_dimension":chars.point_kummer_dimension,
                     "prime_decomposition":[{"ramification_index":int(P[2]),"residue_degree":int(P[3])} for P in chars.primes],
                     "local_reduction":{"conductor_exponent":int(red[0]),"kodaira_code":int(red[1]),
                                        "minimal_change":list(map(str,red[2])),"tamagawa_number":int(red[3]),
                                        "minimal_discriminant_valuation":int(pari.valuation(minimal.disc(),p))}})
    r.write_new(path,{"bindings":bindings(),"case_index":index,"id":old["id"],
                     "factorization_complete":factor["factorization_complete"],"local":rows,
                     "software":{"sage":version,"pari":str(pari.version())},
                     "local_order_basis":list(map(str,nf.nf_get_zk())) if nf else []})


def run(kind,index,seconds):
    path = WORK/f"{kind}-{index}.json"
    if not path.exists():
        with (WORK/f"{kind}-{index}.log").open("x") as log:
            try:
                proc = subprocess.run(["sage","-python",str(Path(__file__).resolve()),kind,
                                       "--index",str(index),"--destination",str(path)],
                                      cwd=r.ROOT,stdout=log,stderr=log,timeout=seconds)
                failure = None if proc.returncode == 0 else "worker failure"
            except subprocess.TimeoutExpired:
                failure = str(seconds)+"-second timeout"
            if failure:
                log.write("\nUNKNOWN: "+failure+"\n")
                if kind == "factor":
                    trial = r.read(WORK/f"trial-{index}.json")
                    r.write_new(path,{**trial,"factorization_complete":False,
                                      "factors":trial["trial_factors"],"reason":failure})
                else:
                    r.write_new(path,{"bindings":bindings(),"case_index":index,"status":"UNKNOWN","reason":failure})
    row = r.read(path)
    assert row["bindings"] == bindings()
    return row


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    rows = []
    limits = r.read(PROTOCOL)["limits"]
    for i in range(6):
        trial = WORK/f"trial-{i}.json"
        if not trial.exists():
            r.write_new(trial,initial_factors(i))
        factor = run("factor",i,limits["factor_seconds_per_curve"])
        print("FACTORS",i,factor["factorization_complete"],len(factor["factors"]),flush=True)
        local = run("local",i,limits["local_seconds_per_curve"])
        rows.append({"factor":factor,"local":local})
        print("CHECKPOINT",i,local.get("status","local arithmetic complete"),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.remaining-bad-prime-inputs.v1","bindings":bindings(),"cases":rows})


def build(check=False):
    data = r.read(INPUT)
    assert data["bindings"] == bindings()
    rows = []
    old = r.read(bad.INPUT)["cases"]
    two = r.read(dy.INPUT)["cases"]
    for i,case in enumerate(data["cases"]):
        fac,local = case["factor"],case["local"]
        assert fac["bindings"] == bindings() and local["bindings"] == bindings()
        product = 1
        for p,e in fac["factors"]:
            product *= p**e
        disc = abs(int(fac["model_discriminant"]))
        assert disc%product == 0
        if fac["factorization_complete"]:
            assert product == disc
        if local.get("status") == "UNKNOWN":
            rows.append({"id":old[i]["id"],"status":"UNKNOWN","reason":local["reason"]})
            continue
        assert [x["prime"] for x in local["local"]] == [p for p,e in fac["factors"] if p>13]
        merged = {**old[i],"local":old[i]["local"]+two[i]["local"]+local["local"]}
        profile = bad.characterize(merged)
        bad_places = [p for p in merged["local"] if p["prime"] in (2,"infinity") or p["local_reduction"].get("conductor_exponent",0)>0]
        bad_profile = bad.characterize({**old[i],"local":bad_places})
        rows.append({"id":old[i]["id"],"factorization_complete":fac["factorization_complete"],
                     "tested_all_bad_places":fac["factorization_complete"],"extended_dictionary":profile,
                     "bad_places_and_two_real":bad_profile})
    report = {"schema":"rank-jump.remaining-bad-primes.v1","bindings":bindings(),
              "input_sha256":r.digest(INPUT.read_bytes()),"rows":rows,
              "boundary":"Only the finite bad-place dictionary plus 2 and infinity. No full Selmer, global divisibility, exact full rank or selector claim."}
    if check:
        assert r.read(OUTPUT) == report
    else:
        r.write_new(OUTPUT,report)
    for row in rows:
        if row.get("status") == "UNKNOWN":
            print(row["id"],"UNKNOWN")
            continue
        p = row["bad_places_and_two_real"]
        print(row["id"],"complete",row["tested_all_bad_places"],"places",len(p["local"]),
              "joint",p["joint_generic_image_dimension"],p["joint_witness_image_dimension"],
              "full product",p["full_product_point_image_dimension"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","factor","local","build","check"))
    parser.add_argument("--index",type=int)
    parser.add_argument("--destination",type=Path)
    args = parser.parse_args()
    if args.mode == "capture":
        capture()
    elif args.mode == "factor":
        factor_worker(args.index,args.destination)
    elif args.mode == "local":
        local_worker(args.index,args.destination)
    else:
        build(args.mode == "check")
