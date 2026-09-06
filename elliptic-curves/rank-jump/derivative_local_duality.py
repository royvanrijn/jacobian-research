#!/usr/bin/env python3
"""Local Tate-duality certificate for the canonical derivative boundary class."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import derivative_reciprocity as first
import strict_class_blocks as strict
import remaining_bad_primes as rem

PROTOCOL=Path(__file__).with_name("DERIVATIVE_LOCAL_DUALITY_PROTOCOL.json")
INPUT=r.OUT/"rank_jump_derivative_local_duality_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_derivative_local_duality_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-derivative-local-duality-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,first.INPUT,first.OUTPUT,strict.OUTPUT,rem.INPUT,rem.OUTPUT,
             rem.bad.INPUT,rem.dy.INPUT,r.INPUT,rem.bad.LOCAL_SOURCE)}


def compute(index):
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari
    from sage.version import version
    sys.path.insert(0,str(rem.bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old=r.read(rem.bad.INPUT)["cases"][index]
    source=rem.bad.cases()[index]
    block=r.read(strict.OUTPUT)["rows"][index]
    factor=r.read(rem.INPUT)["cases"][index]["factor"]
    assert block["all_bad_places_complete"] and factor["factorization_complete"]
    primes=[p for p in block["tested_places"] if p!="infinity"]
    assert {p for p,e in factor["factors"]} <= set(primes)
    _,allpoints=r.short(source["model"],source["generic_points"]+source["points"])
    points=[allpoints[i] for i in old["selected_input_indices"]]
    m=old["generic_dimension"]
    R=PolynomialRing(QQ,"z")
    f=R(list(map(QQ,old["integral_cubic_ascending"])))
    delta=f.discriminant()
    assert delta!=0 and 16*delta==ZZ(factor["model_discriminant"])
    assert all(ZZ(p).is_prime(proof=True) for p,e in factor["factors"])
    product=ZZ(1)
    for p,e in factor["factors"]: product*=ZZ(p)**e
    assert product==abs(16*delta)
    nf=pari.nfinit([pari(f),primes])
    theta=pari.Mod("z",pari(f))
    beta=-pari(delta)*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==delta**4
    d=QQ(old["elliptic_scaling_d"])
    gammas=[pari(QQ(P[0])*d*d)-theta for P in points]
    roots=f.roots(AA,multiplicities=False)
    beta_signs=[int(-delta*f.derivative()(a)<0) for a in roots]
    assert beta_signs==([1,0,1] if delta>0 else [0])
    rows=[]
    for local in strict.local_rows(index):
        p=local["prime"]
        if p=="infinity":
            actual=[r.pack(int(QQ(P[0])*d*d<a) for a in roots) for P in points]
            signature=r.pack(beta_signs)
            width=len(roots)
            chars=None
        else:
            chars=LocalSquareclasses(nf,p)
            actual=[r.pack(chars.signature(g)) for g in gammas]
            raw=chars.signature(beta)
            signature=r.pack(raw)
            width=len(raw)
        selected=[]
        basis=[]
        for i in range(m):
            if r.rank(basis+[actual[i]])>len(basis):
                selected.append(i);basis.append(actual[i])
        assert len(basis)==local["point_kummer_dimension"]<=r.read(PROTOCOL)["limits"]["max_point_image_dimension"]
        assert all(r.reduce(v,r.basis(basis))==0 for v in actual)
        residual=r.reduce(signature,r.basis(basis))
        member=residual==0
        separator=None
        word=None
        if member:
            word=lc.coordinates(signature,basis)
        else:
            separator=next(a for a in lc.orthogonal(basis,width) if (a&signature).bit_count()%2)
            assert all((a&separator).bit_count()%2==0 for a in basis)
        # Independent local-square comparison of every possible correction.
        local_square_results=[]
        for mask in range(1<<len(basis)):
            expected=(signature==lc.lift(mask,basis))
            if p=="infinity":
                value=expected
            else:
                b=beta
                for j,i in enumerate(selected):
                    if mask>>j&1: b*=gammas[i]
                den=pari.denominator(pari.nfalgtobasis(nf,b))
                b*=den**2
                value=all(pari.nfislocalpower(nf,P,b,2)==1 for P in chars.primes)
                assert value==expected
            local_square_results.append(bool(value))
        assert member==any(local_square_results)
        rows.append({"place":p,"width":width,"point_dimension":len(basis),"basis_generic_indices":selected,
                     "basis_signatures":basis,"derivative_signature":signature,"is_in_point_image":member,
                     "membership_word":word,"separating_character":separator,
                     "local_square_comparison_results":local_square_results})
    return {"bindings":bindings(),"case_index":index,"id":old["id"],
            "software":{"sage":version,"pari":str(pari.version())},
            "delta":str(delta),"norm_beta":str(delta**4),
            "beta_coordinates":[str(pari.lift(beta).polcoef(i)) for i in range(3)],
            "beta_real_sign_bits":beta_signs,"local":rows,
            "outside_places":[x["place"] for x in rows if not x["is_in_point_image"]]}


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    cases=[]
    for i in r.read(PROTOCOL)["cases"]:
        path=WORK/f"case-{i}.json"
        if not path.exists():
            with (WORK/f"case-{i}.log").open("x") as log:
                try:
                    p=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                      "--index",str(i),"--destination",str(path)],
                                     cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if p.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:
                    failure="30-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"case_index":i,"status":"UNKNOWN","reason":failure})
        record=r.read(path);assert record["bindings"]==bindings()
        cases.append(record)
        print("checkpoint",i,record.get("outside_places",record.get("status")),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.derivative-local-duality-inputs.v1","bindings":bindings(),"cases":cases})


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    rows=[]
    for row in data["cases"]:
        assert row["bindings"]==bindings()
        if row.get("status")=="UNKNOWN":
            rows.append(row);continue
        outside=[]
        for local in row["local"]:
            b=local["basis_signatures"];sig=local["derivative_signature"]
            member=r.reduce(sig,r.basis(b))==0
            assert member==local["is_in_point_image"]
            if member:
                assert lc.lift(local["membership_word"],b)==sig
            else:
                a=local["separating_character"]
                assert (a&sig).bit_count()%2==1
                assert all((a&v).bit_count()%2==0 for v in b)
                outside.append(local["place"])
        assert outside==row["outside_places"]
        old=r.read(strict.OUTPUT)["rows"][row["case_index"]]
        profile=r.read(rem.OUTPUT)["rows"][row["case_index"]]["bad_places_and_two_real"]
        ell=profile["full_product_point_image_dimension"]
        a=profile["joint_witness_image_dimension"]
        assert ell-a==1
        resolved=bool(outside)
        n=old["witness_dimension"];m=old["generic_dimension"]
        rows.append({"case_index":row["case_index"],"id":row["id"],"outside_places":outside,
                     "boundary_bit":0 if resolved else None,
                     "full_Selmer_local_image_dimension":a if resolved else None,
                     "local_product_dimension":ell,
                     "full_Selmer_dimension_formula":f"{n} + epsilon"+("" if resolved else " + b"),
                     "full_relative_Selmer_dimension_formula":f"{n-m} + epsilon"+("" if resolved else " + b"),
                     "epsilon_definition":f"dim_F2 Cl(O_K,S)/2 - {len(old['witness_strict_kernel_masks'])} >= 0",
                     "boundary":"No upper bound on epsilon and no new rational point or global solubility claim."})
    result={"schema":"rank-jump.derivative-local-duality.v1","bindings":bindings(),
            "input_sha256":r.digest(INPUT.read_bytes()),"rows":rows}
    if check:
        assert r.read(OUTPUT)==result
        print("PASS exact local nonmembership and reciprocity consequence")
    else: r.write_new(OUTPUT,result)
    for row in rows: print(row.get("id",row["case_index"]),row.get("full_Selmer_dimension_formula",row.get("status")))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","build","check","verify"))
    parser.add_argument("--index",type=int);parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="capture": capture()
    elif args.mode=="worker": r.write_new(args.destination,compute(args.index))
    elif args.mode=="verify":
        record=next(x for x in r.read(INPUT)["cases"] if x["case_index"]==args.index)
        assert record==compute(args.index)
        print("PASS exact local square comparisons",record["id"])
    else: build(args.mode=="check")
