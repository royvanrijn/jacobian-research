#!/usr/bin/env python3
"""Artin evaluation on strict half ideals via cyclic residue rings and Jacobi symbols."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import strict_half_ideals as h
import remaining_bad_primes as rem

PROTOCOL=Path(__file__).with_name("HALF_IDEAL_ARTIN_PROTOCOL.json")
INPUT=r.OUT/"rank_jump_half_ideal_artin_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_half_ideal_artin_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-half-ideal-artin-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,h.INPUT,rem.INPUT)}


def compute(index):
    from sage.all import QQ,ZZ,PolynomialRing,matrix,pari
    from sage.version import version
    source=next(x for x in r.read(h.INPUT)["cases"] if x["case_index"]==index)
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    assert factor["factorization_complete"]
    R=PolynomialRing(QQ,"z");f=pari(R(list(map(QQ,source["integral_cubic_ascending"]))))
    primes=[p for p,e in factor["factors"]]
    nf=pari.nfinit([f,primes])
    assert list(map(str,nf.nf_get_zk()))==source["maximal_order_basis"]
    assert str(nf.nf_get_zk()[0])=="1"
    prime_ideals=[(p,j,P) for p in primes for j,P in enumerate(pari.idealprimedec(nf,p))]
    gammas=[pari.Mod(R(list(map(QQ,p["gamma_coordinates"]))),f) for p in source["points"]]
    classes=[]
    for c in source["half_ideals"]:
        gamma=pari.Mod(1,f)
        for i,g in enumerate(gammas):
            if c["point_mask"]>>i&1:gamma*=g
        coordinates=list(map(ZZ,pari.nfalgtobasis(nf,gamma)))
        classes.append(coordinates)
    columns=[]
    for position,c in enumerate(source["half_ideals"]):
        J=pari(matrix(QQ,c["half_ideal_hnf"]))
        reduced,alpha=pari.idealred(nf,[J,1])
        assert pari.idealmul(nf,reduced,alpha)==J
        multiplier=pari.nfbasistoalg(nf,alpha)
        good=reduced;removed=[]
        for p,j,P in prime_ideals:
            e=int(pari.idealval(nf,reduced,P))
            assert e>=0
            removed.append({"prime":p,"prime_index":j,"exponent":e})
            if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
        good=pari.idealhnf(nf,good)
        N=ZZ(pari.idealnorm(nf,good))
        assert N>0 and N%2==1
        assert all(N%p for p in primes)
        cyclic=good[1,1]==good[2,2]==1 and good[0,0]==N
        evaluations=[]
        if cyclic:
            for coordinates in classes:
                residue=ZZ(sum(coordinates[j]*([1,-good[0,1],-good[0,2]][j]) for j in range(3)))%N
                gcd=residue.gcd(N)
                if gcd!=1:
                    evaluations.append({"status":"UNKNOWN","reason":"nonunit residue","gcd":str(gcd)})
                else:
                    symbol=int(pari.kronecker(residue,N))
                    assert symbol in (-1,1)
                    evaluations.append({"residue":str(residue),"jacobi_symbol":symbol,"artin_bit":int(symbol==-1)})
        columns.append({"position":position,"point_mask":c["point_mask"],
                        "reduced_ideal_hnf":h.matrix_record(reduced),
                        "principal_multiplier_coordinates":[str(pari.lift(multiplier).polcoef(j)) for j in range(3)],
                        "removed_S_prime_ideals":removed,"coprime_ideal_hnf":h.matrix_record(good),
                        "norm":str(N),"cyclic_residue_ring":bool(cyclic),"evaluations":evaluations,
                        "status":"complete" if cyclic and all("artin_bit" in e for e in evaluations) else "UNKNOWN"})
    return {"bindings":bindings(),"case_index":index,"id":source["id"],
            "software":{"sage":version,"pari":str(pari.version())},
            "generic_strict_count":source["generic_strict_count"],
            "character_point_masks":[c["point_mask"] for c in source["half_ideals"]],
            "columns":columns}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);cases=[]
    for index in r.read(PROTOCOL)["cases"]:
        path=WORK/f"case-{index}.json"
        if not path.exists():
            with (WORK/f"case-{index}.log").open("x") as log:
                try:
                    process=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                            "--index",str(index),"--destination",str(path)],
                                           cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if process.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:
                    failure="30-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"case_index":index,"status":"UNKNOWN","reason":failure})
        row=r.read(path);assert row["bindings"]==bindings()
        cases.append(row)
        print("checkpoint",index,row.get("status",[x["status"] for x in row["columns"]]),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.half-ideal-artin-inputs.v1","bindings":bindings(),"cases":cases})


def summarize(row):
    if row.get("status")=="UNKNOWN":return row
    full=[c for c in row["columns"] if c["status"]=="complete"]
    packed=[r.pack(e["artin_bit"] for e in c["evaluations"]) for c in full]
    generic=[c for c in full if c["position"]<row["generic_strict_count"]]
    gcols=[r.pack(e["artin_bit"] for e in c["evaluations"]) for c in generic]
    generic_complete=len(generic)==row["generic_strict_count"]
    generic_upper=r.rank(gcols) if generic_complete else row["generic_strict_count"]
    return {"case_index":row["case_index"],"id":row["id"],
            "known_strict_dimension":len(row["character_point_masks"]),
            "generic_strict_dimension":row["generic_strict_count"],
            "complete_columns":len(full),"generic_columns_complete":generic_complete,
            "artin_rank":r.rank(packed),"generic_artin_rank":r.rank(gcols),
            "S_class_mod_two_half_ideal_image_dimension_lower_bound":r.rank(packed),
            "relative_S_class_mod_two_half_ideal_image_dimension_lower_bound":max(0,r.rank(packed)-generic_upper),
            "matrix_rows":[[c["evaluations"][i]["artin_bit"] if c["status"]=="complete" else None
                            for c in row["columns"]] for i in range(len(row["character_point_masks"]))],
            "boundary":"Retrospective pairing of known strict characters with known half ideals; no full class group or CT computation. Zero on known characters is not a triviality or divisibility certificate."}


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    rows=[summarize(row) for row in data["cases"]]
    report={"schema":"rank-jump.half-ideal-artin.v1","bindings":bindings(),
            "input_sha256":r.digest(INPUT.read_bytes()),"rows":rows}
    if check:
        assert r.read(OUTPUT)==report
        print("PASS Artin matrix accounting")
    else:r.write_new(OUTPUT,report)
    for row in rows:
        print(row.get("id",row["case_index"]),row.get("complete_columns"),
              row.get("artin_rank"),row.get("relative_S_class_mod_two_half_ideal_image_dimension_lower_bound"))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("capture","worker","build","check"))
    p.add_argument("--index",type=int);p.add_argument("--destination",type=Path)
    args=p.parse_args()
    if args.mode=="capture":capture()
    elif args.mode=="worker":r.write_new(args.destination,compute(args.index))
    else:build(args.mode=="check")
