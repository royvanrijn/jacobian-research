#!/usr/bin/env python3
"""A fixed five-character Artin experiment on a certified strict Sha block."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import local_collision as lc
import affine_selmer as af
import strict_deformation_solubility as strict
import strict_half_ideals as half

PROTOCOL=Path(__file__).with_name("STRICT_SHA_ARTIN_PROTOCOL.json")
INPUT=r.OUT/"rank_jump_strict_Sha_Artin_inputs_v1.json"
OUTPUT=r.OUT/"rank_jump_strict_Sha_Artin_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-strict-Sha-Artin-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,strict.OUTPUT,af.INPUT,lc.INPUT)}


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,matrix,pari
    from sage.version import version
    protocol=r.read(PROTOCOL);u=protocol["u"];masks=protocol["point_masks"]
    sr=next(x for x in r.read(strict.OUTPUT)["single_deformations"] if x["u"]==u)
    assert masks==sr["strict_anchor_basis"]
    assert sr["CT_cross_report"]["cross_pairing_rank"]==len(masks)==5
    local=next(x for x in r.read(af.INPUT)["cases"] if x["u"]==u)
    anchor=r.read(lc.INPUT)["anchor"]
    r.short(anchor["short_model_ainvariants"],anchor["known_points_on_short_model"])
    R=PolynomialRing(QQ,"z");f=R(list(map(QQ,anchor["base_polynomial_ascending"])))
    assert f.is_irreducible()
    primes=[x["place"] for x in local["local"] if x["place"]!="infinity"]
    product=ZZ(1)
    for e in anchor["base_discriminant_factorization"]:
        p=e["prime"]
        assert p in primes and ZZ(p).is_prime(proof=True)
        product*=ZZ(p)**e["exponent"]
    assert product==abs(f.discriminant())
    nf=pari.nfinit([pari(f),primes]);theta=pari.Mod("z",pari(f))
    assert str(nf.nf_get_zk()[0])=="1"
    Ps=[(p,j,P) for p in primes for j,P in enumerate(pari.idealprimedec(nf,p))]
    points=[];gammas=[];gcd_ideals=[]
    for raw in anchor["known_points_on_short_model"]:
        x,y=map(QQ,raw);d=ZZ(x.denominator()).sqrt();assert d in ZZ
        d=ZZ(d);a=ZZ(x*d*d);b=ZZ(y*d**3)
        gamma=pari(a)-pari(d*d)*theta
        assert pari.nfeltnorm(nf,gamma)==b*b
        I=pari.idealadd(nf,pari(b),gamma)
        points.append({"a":str(a),"b":str(b),"d":str(d),
                       "gamma_coordinates":[str(pari.lift(gamma).polcoef(j)) for j in range(3)],
                       "gcd_ideal_hnf":half.matrix_record(I)})
        gammas.append(gamma);gcd_ideals.append(I)
    classes=[];ideals=[]
    for mask in masks:
        beta=pari.Mod(1,pari(f));I=pari.idealhnf(nf,1);norm=ZZ(1)
        for j,gamma in enumerate(gammas):
            if mask>>j&1:
                beta*=gamma;I=pari.idealmul(nf,I,gcd_ideals[j]);norm*=abs(ZZ(points[j]["b"]))
        J=I;corrections=[]
        for p,j,P in Ps:
            vb=int(pari.idealval(nf,beta,P));vi=int(pari.idealval(nf,I,P))
            assert vb%2==0
            e=vb//2-vi
            corrections.append({"prime":p,"prime_index":j,"beta_valuation":vb,
                                "gcd_product_valuation":vi,"exponent":e})
            if e:J=pari.idealmul(nf,J,pari.idealpow(nf,P,e))
        assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,beta)
        assert pari.idealnorm(nf,J)==norm
        reduced,alpha=pari.idealred(nf,[J,1])
        assert pari.idealmul(nf,reduced,alpha)==J
        good=reduced;removed=[]
        for p,j,P in Ps:
            e=int(pari.idealval(nf,reduced,P))
            assert e>=0
            removed.append({"prime":p,"prime_index":j,"exponent":e})
            if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
        good=pari.idealhnf(nf,good);N=ZZ(pari.idealnorm(nf,good))
        assert N>0 and all(N%p for p in primes)
        a=pari.nfbasistoalg(nf,alpha)
        ideals.append({"point_mask":mask,"half_ideal_hnf":half.matrix_record(J),"half_ideal_norm":str(norm),
                       "corrections":corrections,"reduced_ideal_hnf":half.matrix_record(reduced),
                       "principal_multiplier_coordinates":[str(pari.lift(a).polcoef(j)) for j in range(3)],
                       "removed_S_prime_ideals":removed,"coprime_ideal_hnf":half.matrix_record(good),
                       "norm":str(N),"cyclic":bool(good[0,0]==N and good[1,1]==good[2,2]==1)})
        classes.append(beta)
    for column in ideals:
        H=pari(matrix(QQ,column["coprime_ideal_hnf"]));N=ZZ(column["norm"])
        evaluations=[]
        if column["cyclic"]:
            for beta in classes:
                coords=pari.nfalgtobasis(nf,beta)
                value=ZZ(coords[0]-H[0,1]*coords[1]-H[0,2]*coords[2])
                gcd=value.gcd(N)
                if gcd==1:
                    symbol=int(pari.kronecker(value,N));assert symbol in (-1,1)
                    evaluations.append({"residue":str(value%N),"jacobi_symbol":symbol,"bit":int(symbol==-1)})
                elif gcd<=1000 and gcd.is_prime(proof=True):
                    p=int(gcd);e=int(N.valuation(p));M=N//ZZ(p)**e
                    if value.gcd(M)!=1:
                        evaluations.append({"status":"UNKNOWN","reason":"nonunit cofactor","gcd":str(gcd)});continue
                    matching=[P for P in pari.idealprimedec(nf,p) if pari.idealval(nf,H,P)>0]
                    assert len(matching)==1
                    P=matching[0]
                    assert int(P[3])==1 and int(pari.idealval(nf,H,P))==e
                    local_bit=int(pari.nfislocalpower(nf,P,beta,2)==0)
                    symbol=int(pari.kronecker(value,M));assert symbol in (-1,1)
                    evaluations.append({"gcd_prime":p,"prime_exponent":e,"cofactor":str(M),
                                        "cofactor_residue":str(value%M),"cofactor_jacobi_symbol":symbol,
                                        "local_bit":local_bit,"bit":int(symbol==-1)^((e%2)*local_bit)})
                else:evaluations.append({"status":"UNKNOWN","reason":"gcd outside frozen rule","gcd":str(gcd)})
        column["evaluations"]=evaluations
        column["complete"]=column["cyclic"] and len(evaluations)==5 and all("bit" in x for x in evaluations)
    return {"bindings":bindings(),"u":u,"software":{"sage":version,"pari":str(pari.version())},
            "cubic_ascending":anchor["base_polynomial_ascending"],"S_finite":primes,
            "maximal_order_basis":list(map(str,nf.nf_get_zk())),"point_masks":masks,
            "points":points,"ideals":ideals}


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    if not INPUT.exists():
        with (WORK/"worker.log").open("x") as log:
            try:
                p=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker"],
                                 cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                reason=None if p.returncode==0 else "worker failure"
            except subprocess.TimeoutExpired:reason="30-second timeout"
            if reason:
                log.write("\nUNKNOWN: "+reason+"\n")
                r.write_new(INPUT,{"bindings":bindings(),"status":"UNKNOWN","reason":reason})
    data=r.read(INPUT);assert data["bindings"]==bindings()
    print(data.get("status",[(x["complete"],x["norm"]) for x in data.get("ideals",[])]))


def build(check=False):
    data=r.read(INPUT);assert data["bindings"]==bindings()
    if data.get("status")=="UNKNOWN":
        result=data
    else:
        complete_columns=[(i,r.pack(x["bit"] for x in column["evaluations"]))
                          for i,column in enumerate(data["ideals"]) if column["complete"]]
        selected=[];columns=[]
        for i,col in complete_columns:
            if r.rank(columns+[col])>len(columns):selected.append(i);columns.append(col)
        rows=[r.pack((col>>i)&1 for col in columns) for i in range(5)]
        selected_rows=[];row_basis=[]
        for i,row in enumerate(rows):
            if r.rank(row_basis+[row])>len(row_basis):selected_rows.append(i);row_basis.append(row)
        d=len(columns)
        subcolumns=[r.pack((col>>i)&1 for i in selected_rows) for col in columns]
        dual=[lc.lift(lc.coordinates(1<<i,subcolumns),[1<<j for j in selected]) for i in range(d)]
        result={"u":data["u"],"strict_Sha_dimension":5,"complete_Artin_columns":len(complete_columns),
                "Artin_rank":d,"elementary_S_class_direct_factor_dimension":d,
                "selected_half_ideal_indices":selected,"selected_character_indices":selected_rows,
                "dual_half_ideal_words":dual,
                "Artin_matrix_rows":[[column["evaluations"][i]["bit"] if column["complete"] else None
                                      for column in data["ideals"]] for i in range(5)],
                "boundary":"Selected strict characters are rational on E0 but inject into Sha(E_-1)[2]. Their elementary S-class direct factor and Artin data are identical on both curves. No full class-group or rank upper bound."}
    report={"schema":"rank-jump.strict-Sha-Artin.v1","bindings":bindings(),
            "input_sha256":r.digest(INPUT.read_bytes()),"result":result}
    if check:assert r.read(OUTPUT)==report;print("PASS strict Sha Artin accounting")
    else:r.write_new(OUTPUT,report)
    print(result)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("capture","worker","build","check"))
    args=p.parse_args()
    if args.mode=="capture":capture()
    elif args.mode=="worker":r.write_new(INPUT,compute())
    else:build(args.mode=="check")
