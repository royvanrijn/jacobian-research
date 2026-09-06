#!/usr/bin/env python3
"""Resolve only frozen nonunit Jacobi entries by four small local prime ideals."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import local_collision as lc
import half_ideal_artin as art
import strict_half_ideals as h
import remaining_bad_primes as rem

PROTOCOL=Path(__file__).with_name("HALF_IDEAL_ARTIN_COMPLETION_PROTOCOL.json")
INPUT=r.OUT/"rank_jump_half_ideal_artin_completion_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_half_ideal_artin_completion_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-half-ideal-artin-completion-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,art.INPUT,h.INPUT,rem.INPUT)}


def compute(index):
    from sage.all import QQ,ZZ,PolynomialRing,matrix,pari
    row=next(c for c in r.read(art.INPUT)["cases"] if c["case_index"]==index)
    source=next(c for c in r.read(h.INPUT)["cases"] if c["case_index"]==index)
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    R=PolynomialRing(QQ,"z");f=pari(R(list(map(QQ,source["integral_cubic_ascending"]))))
    nf=pari.nfinit([f,[p for p,e in factor["factors"]]])
    theta=pari.Mod("z",f)
    gamma_points=[pari.Mod(R(list(map(QQ,P["gamma_coordinates"]))),f) for P in source["points"]]
    point_model=["0","0","0",source["integral_cubic_ascending"][1],source["integral_cubic_ascending"][0]]
    rational_points=[[str(QQ(P["a"])/QQ(P["d"])**2),str(QQ(P["b"])/QQ(P["d"])**3)] for P in source["points"]]
    repairs=[]
    for exception in r.read(PROTOCOL)["frozen_exceptions"]:
        if exception["case"]!=index:continue
        j=exception["column"];p=exception["prime"];column=row["columns"][j]
        assert ZZ(p).is_prime(proof=True) and int(pari.poldisc(f))%p
        H=pari(matrix(QQ,column["coprime_ideal_hnf"]));N=ZZ(column["norm"])
        e=int(N.valuation(p));M=N//ZZ(p)**e
        assert e>0 and M%p
        matching=[P for P in pari.idealprimedec(nf,p) if pari.idealval(nf,H,P)>0]
        assert len(matching)==1
        P=matching[0]
        assert int(P[2])==int(P[3])==1 and int(pari.idealval(nf,H,P))==e
        tc=pari.nfalgtobasis(nf,theta)
        theta_residue=int((tc[0]-H[0,1]*tc[1]-H[0,2]*tc[2])%p)
        assert pari.idealval(nf,theta-theta_residue,P)>0
        point_bits=[r.point_signature(point_model,Q,[(p,[theta_residue])]) for Q in rational_points]
        for i in exception["characters"]:
            old=column["evaluations"][i]
            assert old["status"]=="UNKNOWN" and int(old["gcd"])==p
            mask=row["character_point_masks"][i];beta=pari.Mod(1,f)
            for k,gamma in enumerate(gamma_points):
                if mask>>k&1:beta*=gamma
            coords=pari.nfalgtobasis(nf,beta)
            residue=ZZ(coords[0]-H[0,1]*coords[1]-H[0,2]*coords[2])%M
            assert residue.gcd(M)==1
            symbol=int(pari.kronecker(residue,M))
            valuation=int(pari.idealval(nf,beta,P));assert valuation%2==0
            local_bit=int(pari.nfislocalpower(nf,P,beta,2)==0)
            point_bit=(r.pack(point_bits)&mask).bit_count()%2
            assert point_bit==local_bit
            bit=int(symbol==-1)^((e%2)*local_bit)
            repairs.append({"column":j,"character":i,"prime":p,"prime_exponent":e,
                            "theta_residue":theta_residue,"cofactor":str(M),"cofactor_residue":str(residue),
                            "cofactor_jacobi_symbol":symbol,"beta_local_valuation":valuation,
                            "local_frobenius_bit":local_bit,"independent_point_bit":point_bit,"artin_bit":bit})
    return {"bindings":bindings(),"case_index":index,"id":row["id"],"repairs":repairs}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for index in r.read(PROTOCOL)["cases"]:
        path=WORK/f"case-{index}.json"
        if not path.exists():
            with (WORK/f"case-{index}.log").open("x") as log:
                try:
                    p=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                      "--index",str(index),"--destination",str(path)],
                                     cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if p.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:failure="30-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"case_index":index,"status":"UNKNOWN","reason":failure})
        row=r.read(path);assert row["bindings"]==bindings();rows.append(row)
        print("checkpoint",index,row.get("status",len(row.get("repairs",[]))),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.half-ideal-artin-completion-inputs.v1","bindings":bindings(),"cases":rows})


def calculate(row):
    original=next(c for c in r.read(art.INPUT)["cases"] if c["case_index"]==row["case_index"])
    if row.get("status")=="UNKNOWN":return row
    repairs={(x["column"],x["character"]):x for x in row["repairs"]}
    n=len(original["character_point_masks"]);g=original["generic_strict_count"]
    columns=[]
    for j,c in enumerate(original["columns"]):
        bits=[]
        for i,entry in enumerate(c["evaluations"]):
            if "artin_bit" in entry:bits.append(entry["artin_bit"])
            else:bits.append(repairs[j,i]["artin_bit"])
        assert len(bits)==n
        columns.append(r.pack(bits))
    rank=r.rank(columns);grank=r.rank(columns[:g])
    dual_words=[lc.coordinates(1<<i,columns) for i in range(n)] if rank==n else []
    return {"case_index":row["case_index"],"id":row["id"],"character_dimension":n,
            "generic_strict_dimension":g,"artin_rank":rank,"generic_artin_rank":grank,
            "relative_S_class_mod_two_image_dimension_lower_bound":rank-grank,
            "columns":columns,"matrix_rows":[[(c>>i)&1 for c in columns] for i in range(n)],
            "full_rank":rank==n,"dual_half_ideal_words":dual_words,
            "direct_factor_dimension":n if rank==n else None,
            "boundary":"All ideals are certified 2-torsion. Full Artin rank certifies their direct factor in the S-class group; deficient rank is only a lower bound. The remaining class factor and unknown rational solubility are unresolved."}


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    report={"schema":"rank-jump.half-ideal-artin-completion.v1","bindings":bindings(),
            "input_sha256":r.digest(INPUT.read_bytes()),"rows":[calculate(row) for row in data["cases"]]}
    if check:
        assert r.read(OUTPUT)==report;print("PASS completed matrices and dual half-ideal words")
    else:r.write_new(OUTPUT,report)
    for row in report["rows"]:
        print(row.get("id",row["case_index"]),"rank",row.get("artin_rank"),
              "relative",row.get("relative_S_class_mod_two_image_dimension_lower_bound"),
              "direct",row.get("direct_factor_dimension"))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("capture","worker","build","check","verify"))
    p.add_argument("--index",type=int);p.add_argument("--destination",type=Path)
    args=p.parse_args()
    if args.mode=="capture":capture()
    elif args.mode=="worker":r.write_new(args.destination,compute(args.index))
    elif args.mode=="verify":
        row=next(x for x in r.read(INPUT)["cases"] if x["case_index"]==args.index)
        assert row==compute(args.index);print("PASS local completion",args.index)
    else:build(args.mode=="check")
