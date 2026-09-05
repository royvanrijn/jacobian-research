#!/usr/bin/env python3
"""Six fixed affine Selmer intersections. Reuse local witnesses; no searches."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
from ct_variation import solve

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"AFFINE_SELMER_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_affine_selmer_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_affine_selmer_v1.json"
CHECKPOINT=r.ROOT/"artifacts/local/rank-jump-affine-selmer-v1"
LOCAL=r.ROOT/"elliptic-curves/cas/research_runtime/local_kummer.py"


def source():
    path=r.ROOT/r.read(PROTOCOL)["source"]
    raw=path.read_bytes()
    old=r.read(lc.INPUT)
    assert r.digest(raw)==old["bindings"][str(path.relative_to(r.ROOT))]
    return json.loads(raw)


def bindings():
    paths=[PROTOCOL,Path(__file__),lc.INPUT,LOCAL,HERE/"local_collision.py",
           HERE/"ct_variation.py",HERE/"retrospective.py"]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def local_square(value,p):
    value=r.F(value)
    if not value:
        return True
    n,d=value.numerator,value.denominator
    val=0
    while n%p==0:
        n//=p;val+=1
    while d%p==0:
        d//=p;val-=1
    if val%2:
        return False
    modulus=8 if p==2 else p
    unit=(n*pow(d,-1,modulus))%modulus
    return unit==1 if p==2 else pow(unit,(p-1)//2,p)==1


def compute(index):
    from sage.all import QQ,AA,PolynomialRing,pari
    from sage.version import version
    sys.path.insert(0,str(r.ROOT/"elliptic-curves/cas"))
    from research_runtime.local_kummer import LocalSquareclasses
    original=source()
    u=r.read(PROTOCOL)["parameters"][index]
    row=next(x for x in original["runs"] if int(x["parameter_u"])==u)
    anchor=original["anchor"]
    A,B=map(QQ,anchor["short_model_ainvariants"][3:])
    polyring=PolynomialRing(QQ,"z")
    z=polyring.gen();f=z**3+A*z+B
    primes=[x["prime"] for x in row["finite_local_conditions"]]
    # Listed primes only: no global class group/maximal-order campaign.
    nf=pari.nfinit([pari(f),primes])
    th=pari.Mod(pari(z),pari(f))
    D=1+A*u*u+B*u**3
    alpha=th+u*th**2
    eta=pari(D)*(1-u*th)
    betas=[pari(QQ(P[0]))-th for P in anchor["known_points_on_short_model"]]
    local=[]
    for pinfo in row["finite_local_conditions"]:
        p=pinfo["prime"];chars=LocalSquareclasses(nf,p)
        xs=[QQ(w["raw_x"]) for w in pinfo["basis_witnesses"]]
        for x in xs:
            val=x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3
            assert local_square(str(val),p)
        gens=[pari(x)-alpha for x in xs]
        signatures=[list(chars.signature(value)) for value in gens+betas+[eta]]
        m=len(gens)
        assert r.rank(map(r.pack,signatures[:m]))==m==chars.point_kummer_dimension
        valuations=[[int(pari.idealval(nf,value,P)) for P in chars.primes]
                    for value in betas+[eta]]
        local.append({"place":p,"point_generator_x":list(map(str,xs)),
                      "point_dimension":m,"point_signature_rows":signatures[:m],
                      "class_signature_rows":signatures[m:],
                      "class_ideal_valuations":valuations,
                      "prime_decomposition":[[int(P[2]),int(P[3])] for P in chars.primes],
                      "new_relative_to_anchor":p in row["newly_bad_primes_relative_to_anchor"],
                      "old_constraint_rows":pinfo["known_span_quotient_rows"]})
    roots=sorted(f.roots(AA,multiplicities=False))
    avals=[a+u*a*a for a in roots]
    if len(roots)==3:
        order=sorted(range(3),key=lambda i:avals[i])
        signs=[1]*3;signs[order[0]]=0
        pg=[signs]
    else:
        order=[0];pg=[]
    realrows=[[int(AA(QQ(P[0]))-a<0) for a in roots]
              for P in anchor["known_points_on_short_model"]]
    realrows.append([int(AA(D)*(1-u*a)<0) for a in roots])
    local.append({"place":"infinity","point_dimension":len(pg),
                  "point_signature_rows":pg,"class_signature_rows":realrows,
                  "old_constraint_rows":row["real_local_condition"]["known_span_quotient_rows"],
                  "alpha_root_order":order})
    return {"schema":"rank-jump.affine-selmer-case.v1","bindings":bindings(),"u":u,
            "software":{"sage":version,"pari":str(pari("version()"))},
            "anchor_model":anchor["short_model_ainvariants"],
            "anchor_beta_coordinates":anchor["known_kummer_basis_beta_power_coordinates"],
            "D":str(D),"eta_coordinates":[str(D),str(-D*u),"0"],
            "old_inherited_basis":[x["mask"] for x in row["W_u_basis"]],
            "local":local}


def constraints(entry):
    points=list(map(r.pack,entry["point_signature_rows"]))
    sigs=list(map(r.pack,entry["class_signature_rows"]))
    n=len(entry["class_signature_rows"][0])
    ann=lc.orthogonal(points,n)
    return lc.canonical([r.pack([(v&a).bit_count()%2 for v in sigs]) for a in ann])


def affine_solve(equations,n):
    """Retain row-provenance for a checkable 0=1 certificate."""
    pivots={}
    for i,row in enumerate(equations):
        a=row&((1<<n)-1);b=row>>n;word=1<<i
        while a:
            j=a.bit_length()-1
            if j not in pivots:
                pivots[j]=(a,b,word)
                break
            x,y,w=pivots[j];a^=x;b^=y;word^=w
        else:
            if b:
                assert lc.lift(word,equations)==1<<n
                return {"consistent":False,"inconsistent_row_combination":word}
    particular=solve([row&((1<<n)-1) for row in equations],
                     r.pack([row>>n for row in equations]),n)
    return {"consistent":True,"particular_anchor_mask":particular}


def analyze(record):
    equations=[];places=[];independence=[]
    for entry in record["local"]:
        rows=constraints(entry)
        old=lc.canonical(r.transpose(entry["old_constraint_rows"]))
        new=lc.canonical([v&((1<<20)-1) for v in rows])
        assert old==new
        start=len(equations);equations.extend(rows)
        places.append({"place":entry["place"],"constraint_rows":rows,
                       "global_row_indices":list(range(start,len(equations))),
                       "affine_local_solution":affine_solve(rows,20)})
        if entry.get("new_relative_to_anchor"):
            valuations=entry["class_ideal_valuations"]
            for j in range(len(valuations[0])):
                if valuations[-1][j]%2 and all(v[j]%2==0 for v in valuations[:-1]):
                    independence.append({"prime":entry["place"],"prime_ideal_index":j,
                                         "eta_valuation":valuations[-1][j],
                                         "all_anchor_valuations_even":True})
    assert independence
    kernel=lc.orthogonal([v&((1<<20)-1) for v in equations],20)
    assert kernel==lc.canonical(record["old_inherited_basis"])
    solution=affine_solve(equations,20)
    combined=lc.orthogonal(equations,21)
    expected=len(kernel)+int(solution["consistent"])
    assert len(combined)==expected
    out={"u":record["u"],"places":places,"all_constraint_rows":equations,
         "eta_outside_inherited_space_witnesses":independence,
         "inherited_selmer_dimension":len(kernel),"inherited_selmer_basis":kernel,
         "affine_solution":solution,"extended_selmer_basis":combined,
         "extended_selmer_dimension":len(combined)}
    if solution["consistent"]:
        mask=solution["particular_anchor_mask"]
        corrections=[]
        for entry in record["local"]:
            signature=lc.lift(mask,list(map(r.pack,entry["class_signature_rows"][:20]))) ^ r.pack(entry["class_signature_rows"][20])
            coefficients=lc.coordinates(signature,list(map(r.pack,entry["point_signature_rows"])))
            corrections.append({"place":entry["place"],"point_generator_mask":coefficients})
        out["admissible_representative_local_corrections"]=corrections
    return out


def capture():
    CHECKPOINT.mkdir(parents=True,exist_ok=True)
    records=[]
    for i,u in enumerate(r.read(PROTOCOL)["parameters"]):
        dest=CHECKPOINT/f"case-{i}.json"
        if not dest.exists():
            result=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                   "--case",str(i),"--destination",str(dest)],
                                  capture_output=True,text=True,timeout=30,cwd=r.ROOT)
            with (CHECKPOINT/f"case-{i}.log").open("x") as f:
                f.write(result.stdout+result.stderr)
            if result.returncode or not dest.exists():
                raise RuntimeError(f"worker {i} incomplete; inspect checkpoint log")
        record=r.read(dest)
        assert record["bindings"]==bindings()
        records.append(record);print("checkpoint",u,flush=True)
    r.write_new(INPUT,{"schema":"rank-jump.affine-selmer-inputs.v1",
                       "bindings":bindings(),"cases":records})


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    rows=[analyze(record) for record in data["cases"]]
    out={"schema":"rank-jump.affine-selmer.v1","bindings":bindings(),
         "input_sha256":r.digest(INPUT.read_bytes()),"cases":rows,
         "boundary":"Complete local intersection only inside W+<eta>; neither full Selmer computation nor rational point/rank certification."}
    if check:
        assert r.read(OUTPUT)==out
        print("PASS affine Selmer replay")
    else:
        r.write_new(OUTPUT,out)
    for row in rows:
        print(row["u"],row["inherited_selmer_dimension"],row["extended_selmer_dimension"],row["affine_solution"])


def verify():
    from sage.all import QQ,PolynomialRing,pari,AA
    data=r.read(INPUT);out=r.read(OUTPUT)
    assert data["bindings"]==bindings()
    for i,(record,row) in enumerate(zip(data["cases"],out["cases"])):
        assert compute(i)==record
        if not row["affine_solution"]["consistent"]:
            assert lc.lift(row["affine_solution"]["inconsistent_row_combination"],
                           row["all_constraint_rows"])==1<<20
            print("PASS recomputation and exact inconsistency",row["u"])
            continue
        A,B=map(QQ,record["anchor_model"][3:]);u=row["u"];D=QQ(record["D"])
        R=PolynomialRing(QQ,"z");z=R.gen();pol=z**3+A*z+B
        primes=[x["place"] for x in record["local"] if x["place"]!="infinity"]
        nf=pari.nfinit([pari(pol),primes]);th=pari.Mod(pari(z),pari(pol))
        alpha=th+u*th*th;value=pari(D)*(1-u*th)
        mask=row["affine_solution"]["particular_anchor_mask"]
        for j,beta in enumerate(record["anchor_beta_coordinates"]):
            if mask>>j&1:
                value*=sum(pari(QQ(b))*th**k for k,b in enumerate(beta))
        for local,correction in zip(record["local"],row["admissible_representative_local_corrections"]):
            if local["place"]=="infinity":
                continue
            val=value
            for j,x in enumerate(local["point_generator_x"]):
                if correction["point_generator_mask"]>>j&1:
                    val*=pari(QQ(x))-alpha
            den=pari.denominator(pari.nfalgtobasis(nf,val))
            val*=den**2
            assert all(pari.nfislocalpower(nf,P,val,2) for P in pari.idealprimedec(nf,local["place"]))
        print("PASS recomputation and independent local-power check",u)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("capture","worker","build","check","verify"))
    parser.add_argument("--case",type=int);parser.add_argument("--destination",type=Path)
    args=parser.parse_args()
    if args.mode=="worker":
        r.write_new(args.destination,compute(args.case))
    elif args.mode=="capture":
        capture()
    elif args.mode=="verify":
        verify()
    else:
        build(args.mode=="check")
