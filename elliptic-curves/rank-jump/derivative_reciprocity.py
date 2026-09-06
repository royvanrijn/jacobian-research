#!/usr/bin/env python3
"""Canonical cubic derivative class and its Selmer reciprocity hyperplane."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import strict_class_blocks as strict
import remaining_bad_primes as rem

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"DERIVATIVE_RECIPROCITY_PROTOCOL.json"
INPUT = r.OUT/"rank_jump_derivative_reciprocity_inputs_v1.json"
OUTPUT = r.OUT/"rank_jump_derivative_reciprocity_v1.json"
WORK = r.ROOT/"artifacts/local/rank-jump-derivative-reciprocity-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,strict.OUTPUT,rem.INPUT,rem.OUTPUT,rem.bad.INPUT,rem.dy.INPUT,r.INPUT,rem.bad.LOCAL_SOURCE)}


def compute(index):
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari
    from sage.version import version
    sys.path.insert(0,str(rem.bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old = r.read(rem.bad.INPUT)["cases"][index]
    row = rem.bad.cases()[index]
    block = r.read(strict.OUTPUT)["rows"][index]
    factor = r.read(rem.INPUT)["cases"][index]["factor"]
    assert block["all_bad_places_complete"] and factor["factorization_complete"]
    places = block["tested_places"]
    primes = [p for p in places if p!="infinity"]
    assert {p for p,e in factor["factors"]} <= set(primes)
    _,allpoints = r.short(row["model"],row["generic_points"]+row["points"])
    points = [allpoints[i] for i in old["selected_input_indices"]]
    n,m = old["witness_dimension"],old["generic_dimension"]
    R = PolynomialRing(QQ,"z")
    f = R(list(map(QQ,old["integral_cubic_ascending"])))
    delta = f.discriminant()
    assert delta>0
    assert 16*delta == ZZ(factor["model_discriminant"])
    assert all(ZZ(p).is_prime(proof=True) for p,e in factor["factors"])
    product = ZZ(1)
    for p,e in factor["factors"]:
        product *= ZZ(p)**e
    assert product == 16*delta
    nf = pari.nfinit([pari(f),primes])
    theta = pari.Mod("z",pari(f))
    beta = -pari(delta)*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==delta**4
    d = QQ(old["elliptic_scaling_d"])
    gammas = [pari(QQ(P[0])*d*d)-theta for P in points]
    roots = f.roots(AA,multiplicities=False)
    assert len(roots)==3
    derivative_signs = [int(-delta*f.derivative()(a)<0) for a in roots]
    assert derivative_signs==[1,0,1]
    records = []
    all_coordinates = [0]*n
    functional = 0
    offset = 0
    for local in strict.local_rows(index):
        place = local["prime"]
        signatures = list(map(r.pack,local["point_signature_rows"]))
        basis_indices = []
        basis_signatures = []
        for i in range(m):
            if r.rank(basis_signatures+[signatures[i]])>len(basis_signatures):
                basis_indices.append(i)
                basis_signatures.append(signatures[i])
        dim = local["point_kummer_dimension"]
        assert len(basis_indices)==dim
        coordinates = [lc.coordinates(s,basis_signatures) for s in signatures]
        hilbert = []
        if place=="infinity":
            for i in basis_indices:
                signs = [int(QQ(points[i][0])*d*d<a) for a in roots]
                terms = [-1 if a and b else 1 for a,b in zip(derivative_signs,signs)]
                hilbert.append(terms)
            local_beta_signature = derivative_signs
        else:
            chars = LocalSquareclasses(nf,place)
            actual = list(map(r.pack,[chars.signature(g) for g in gammas]))
            actual_basis = [actual[i] for i in basis_indices]
            assert r.rank(actual_basis)==dim
            assert all(lc.lift(c,actual_basis)==v for c,v in zip(coordinates,actual))
            local_beta_signature = list(chars.signature(beta))
            for i in basis_indices:
                # Always pass a prime ideal: no global Hilbert routine or
                # global factorization is called.
                hilbert.append([int(pari.nfhilbert(nf,beta,gammas[i],P)) for P in chars.primes])
        bits = [sum(int(s==-1) for s in terms)%2 for terms in hilbert]
        for i,c in enumerate(coordinates):
            all_coordinates[i] |= c<<offset
        functional |= r.pack(bits)<<offset
        records.append({"place":place,"dimension":dim,"basis_generic_indices":basis_indices,
                        "point_coordinates":coordinates,"hilbert_symbols_by_basis":hilbert,
                        "functional_bits":bits,"beta_local_signature":local_beta_signature})
        offset += dim
    assert functional != 0
    assert all((functional&v).bit_count()%2==0 for v in all_coordinates)
    assert r.rank(all_coordinates)==offset-1
    assert next(x for x in records if x["place"]=="infinity")["functional_bits"]==[1]
    return {"bindings":bindings(),"case_index":index,"id":old["id"],
            "software":{"sage":version,"pari":str(pari.version())},
            "cubic_discriminant":str(delta),"derivative_beta_coordinates":[str(pari.lift(beta).polcoef(i)) for i in range(3)],
            "norm_beta":str(delta**4),"beta_real_sign_bits":derivative_signs,
            "local_order_basis":list(map(str,nf.nf_get_zk())),"local":records,
            "local_product_dimension":offset,"functional_mask":functional,
            "known_point_local_coordinates":all_coordinates,
            "known_image_dimension":r.rank(all_coordinates),
            "generic_image_dimension":r.rank(all_coordinates[:m]),
            "known_image_equals_reciprocity_hyperplane":True}


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    records = []
    for i in r.read(PROTOCOL)["cases"]:
        path = WORK/f"case-{i}.json"
        if not path.exists():
            with (WORK/f"case-{i}.log").open("x") as log:
                try:
                    proc = subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                           "--index",str(i),"--destination",str(path)],
                                          cwd=r.ROOT,stdout=log,stderr=log,timeout=40)
                    failure = None if proc.returncode==0 else "worker failure"
                except subprocess.TimeoutExpired:
                    failure = "40-second timeout"
                if failure:
                    log.write("\nUNKNOWN: "+failure+"\n")
                    r.write_new(path,{"bindings":bindings(),"case_index":i,"status":"UNKNOWN","reason":failure})
        record = r.read(path)
        assert record["bindings"]==bindings()
        records.append(record)
        print("checkpoint",i,record.get("status","hyperplane verified"),flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.derivative-reciprocity-inputs.v1","bindings":bindings(),"cases":records})


def build(check=False):
    data = r.read(INPUT)
    assert data["bindings"]==bindings()
    results = []
    for row in data["cases"]:
        assert row["bindings"]==bindings()
        if row.get("status")=="UNKNOWN":
            results.append(row)
            continue
        ell = row["local_product_dimension"]
        coords = row["known_point_local_coordinates"]
        functional = row["functional_mask"]
        assert 0<functional<1<<ell
        assert all((v&functional).bit_count()%2==0 for v in coords)
        assert r.rank(coords)==ell-1==row["known_image_dimension"]
        old = r.read(strict.OUTPUT)["rows"][row["case_index"]]
        n,m = old["witness_dimension"],old["generic_dimension"]
        k = len(old["witness_strict_kernel_masks"])
        assert n==k+ell-1
        results.append({"case_index":row["case_index"],"id":row["id"],
                        "local_product_dimension":ell,"full_Selmer_local_image_dimension":ell-1,
                        "full_Selmer_local_image_equals_known_point_image":True,
                        "resolved_boundary_bit":0,"known_strict_dimension":k,
                        "generic_strict_dimension":len(old["generic_strict_kernel_masks"]),
                        "full_Selmer_dimension_formula":f"{n} + epsilon",
                        "full_relative_Selmer_dimension_formula":f"{n-m} + epsilon",
                        "epsilon_definition":f"dim_F2 Cl(O_K,S)/2 - {k} >= 0",
                        "boundary":"Full S-class-group dimension and rational solubility of any unseen strict classes remain UNKNOWN."})
    result = {"schema":"rank-jump.derivative-reciprocity.v1","bindings":bindings(),
              "input_sha256":r.digest(INPUT.read_bytes()),"rows":results}
    if check:
        assert r.read(OUTPUT)==result
        print("PASS reciprocity hyperplane and full Selmer image accounting")
    else:
        r.write_new(OUTPUT,result)
    for row in results:
        print(row.get("id",row["case_index"]),row.get("full_Selmer_dimension_formula",row.get("status")))


def verify(index):
    row = next(x for x in r.read(INPUT)["cases"] if x["case_index"]==index)
    assert row==compute(index)
    print("PASS exact derivative class, local Hilbert symbols and image",row["id"])


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","build","check","verify"))
    parser.add_argument("--index",type=int)
    parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="capture":
        capture()
    elif args.mode=="worker":
        r.write_new(args.destination,compute(args.index))
    elif args.mode=="verify":
        verify(args.index)
    else:
        build(args.mode=="check")
