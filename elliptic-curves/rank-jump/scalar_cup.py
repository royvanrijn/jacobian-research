#!/usr/bin/env python3
"""Scalar cup-product experiment; independent norm control over one small cubic."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import local_collision as lc
import half_ideal_artin_completion as completed
import strict_Sha_Artin as art

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"SCALAR_CUP_PROTOCOL.json"
CONTROL=r.OUT/"rank_jump_scalar_cup_control_v1.json"
OUTPUT=r.OUT/"rank_jump_scalar_cup_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-scalar-cup-v1"


def symmetrize(M):
    return [[M[i][j]^M[j][i] for j in range(len(M))] for i in range(len(M))]


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,completed.OUTPUT,art.OUTPUT)}


def control():
    from sage.all import QQ,PolynomialRing,NumberField,RealIntervalField,pari
    from sage.version import version
    spec=r.read(PROTOCOL)["small_control"]
    R=PolynomialRing(QQ,"y");f=R(spec["polynomial_ascending"])
    bnf=pari.bnfinit(pari(f),1);nf=pari.nfinit(pari(f))
    assert pari.bnfcertify(bnf)==1
    assert list(map(int,bnf.bnf_get_cyc()))==spec["expected_class_group"]
    K=NumberField(f,"theta");embeddings=K.embeddings(RealIntervalField(128))
    dec=lambda a:K([QQ(pari.lift(a).polcoef(i)) for i in range(3)])
    enc=lambda a:[str(pari.lift(a).polcoef(i)) for i in range(3)]
    small_basis=[pari.Mod(-1,pari(f)),*bnf.bnf_get_fu()]
    for J in bnf.bnf_get_gen():
        JJ=pari.idealpow(nf,J,2)
        relation=pari.bnfisprincipal(bnf,JJ,1)
        assert all(x==0 for x in relation[0])
        beta=pari.nfbasistoalg(nf,relation[1])
        assert pari.idealhnf(nf,beta)==JJ
        small_basis.append(beta)
    assert len(small_basis)==5
    places=[P for p in spec["S_finite"] for P in pari.idealprimedec(nf,p)]
    masks=[];betas=[];checks=[]
    for mask in range(1,32):
        beta=pari.Mod(1,pari(f))
        for j,a in enumerate(small_basis):
            if mask>>j&1:beta*=a
        signs=[bool(e(dec(beta))>0) for e in embeddings]
        local=[int(pari.nfislocalpower(nf,P,beta,2)) for P in places]
        good=all(signs) and all(local)
        checks.append({"mask":mask,"positive":signs,"S_local_square":local,"strict":good})
        if good and r.rank(masks+[mask])>len(masks):
            masks.append(mask);betas.append(beta)
    assert len(masks)==2
    print("CONTROL strict basis",masks,flush=True)
    ideals=[]
    for beta in betas:
        fact=pari.idealfactor(nf,beta);J=pari.idealhnf(nf,1)
        for P,e in zip(fact[0],fact[1]):
            assert int(e)%2==0
            J=pari.idealmul(nf,J,pari.idealpow(nf,P,int(e)//2))
        assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,beta)
        ideals.append(J)

    def artin(beta,I):
        fact=pari.idealfactor(nf,I);bit=0;records=[]
        for P,e in zip(fact[0],fact[1]):
            p=int(P[0])
            if p in spec["S_finite"]:continue
            local=int(pari.nfislocalpower(nf,P,beta,2)==0)
            bit^=(int(e)%2)*local
            records.append({"prime":p,"prime_hnf":str(pari.idealhnf(nf,P)),
                            "exponent":int(e),"unramified_frobenius_bit":local})
        return bit,records

    A=[[artin(beta,J)[0] for J in ideals] for beta in betas]
    target=symmetrize(A)
    relative=pari.rnfinit(nf,pari("x^2+1"),1)
    absolute=pari.nfinit(relative)
    norm_setup=pari.rnfisnorminit(nf,pari("x^2+1"),1)
    witnesses=[];observed=[]
    for beta in betas:
        z,remainder=pari.rnfisnorm(norm_setup,beta)
        assert remainder==1
        a,b=[z.lift().polcoef(i) for i in range(2)]
        assert a*a+b*b==beta
        za=pari.rnfeltreltoabs(relative,z)
        fact=pari.idealfactor(absolute,za)
        primes=sorted({int(P[0]) for P in fact[0]}-set(spec["S_finite"]))
        I=pari.idealhnf(nf,1);valuations=[]
        for p in primes:
            for P in pari.idealprimedec(nf,p):
                over=pari.rnfidealprimedec(relative,P)
                exps=[int(pari.idealval(absolute,za,Q)) for Q in over]
                assert len({e%2 for e in exps})==1
                parity=exps[0]%2
                if parity:I=pari.idealmul(nf,I,P)
                valuations.append({"prime":p,"base_prime_hnf":str(pari.idealhnf(nf,P)),
                                   "norm_witness_valuations_above":exps,"descended_parity":parity})
        values=[artin(psi,I) for psi in betas]
        row=[v[0] for v in values];observed.append(row)
        witnesses.append({"beta":enc(beta),"norm_a":enc(a),"norm_b":enc(b),
                          "parity_ideal_hnf":str(I),"valuation_parities":valuations,
                          "Artin_evaluations":[v[1] for v in values]})
        print("CONTROL norm and cup row",row,flush=True)
    result={"bindings":bindings(),"software":{"sage":version,"pari":str(pari.version())},
            "polynomial_ascending":spec["polynomial_ascending"],"S_finite":spec["S_finite"],
            "certified_class_group":list(map(int,bnf.bnf_get_cyc())),
            "squareclass_generators":[enc(a) for a in small_basis],"local_checks":checks,
            "strict_masks":masks,"half_ideal_hnfs":list(map(str,ideals)),
            "Artin_matrix":A,"predicted_scalar_cup_matrix":target,
            "independent_norm_cup_matrix":observed,"norm_witnesses":witnesses,
            "status":"PASS" if observed==target else "REFUTED"}
    r.write_new(CONTROL,result)


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    if not CONTROL.exists():
        with (WORK/"control.log").open("x") as log:
            try:
                proc=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker"],
                                    cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                reason=None if proc.returncode==0 else "worker failed"
            except subprocess.TimeoutExpired:reason="30-second timeout"
            if reason:
                r.write_new(CONTROL,{"bindings":bindings(),"status":"UNKNOWN","reason":reason,
                                    "transcript":(WORK/"control.log").read_text()})
    print(r.read(CONTROL)["status"])


def build(check=False):
    control=r.read(CONTROL)
    assert control["bindings"]==bindings()
    rows=[]
    for old in r.read(completed.OUTPUT)["rows"]:
        if old["case_index"] not in r.read(PROTOCOL)["production_cases"]:continue
        A=old["matrix_rows"];M=symmetrize(A);n=len(M);rank=r.rank(list(map(r.pack,M)))
        rows.append({"case_index":old["case_index"],"id":old["id"],
                     "strict_dimension":n,"Artin_matrix":A,"scalar_cup_matrix":M,
                     "detected_scalar_cup_rank":rank,
                     "retained_space_necessary_twist_solubility_dimension_upper_bound":n-rank,
                     "annihilator_masks":lc.orthogonal(list(map(r.pack,M)),n),
                     "claim_status":"THEOREM_AND_INDEPENDENT_CONTROL" if control["status"]=="PASS" else "WITHHELD"})
    A=r.read(art.OUTPUT)["result"]["Artin_matrix_rows"]
    M=symmetrize(A)
    out={"schema":"rank-jump.scalar-cup.v1","bindings":bindings(),
         "control_sha256":r.digest(CONTROL.read_bytes()),"control_status":control["status"],
         "production_cases":rows,"fixed_cubic_strict_case":{
             "Artin_matrix":A,"scalar_cup_matrix":M,"detected_scalar_cup_rank":r.rank(list(map(r.pack,M))),
             "boundary":"This is scalar -1 cup, not gamma=1+theta for the pencil parameter u=-1."},
         "boundary":r.read(PROTOCOL)["boundary"]}
    if check:assert r.read(OUTPUT)==out;print("PASS retained scalar cup replay")
    else:r.write_new(OUTPUT,out)
    for row in rows:print(row["id"],row["strict_dimension"],row["detected_scalar_cup_rank"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("worker","capture","build","check"))
    mode=p.parse_args().mode
    if mode=="worker":control()
    elif mode=="capture":capture()
    else:build(mode=="check")
