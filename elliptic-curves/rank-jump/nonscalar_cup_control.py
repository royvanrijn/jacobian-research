#!/usr/bin/env python3
"""Two fixed nonscalar norm/ideal controls in the retained small cubic."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import scalar_cup as old

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"NONSCALAR_CUP_CONTROL_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_nonscalar_cup_control_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-nonscalar-cup-control-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,old.CONTROL)}


def worker(index):
    from sage.all import QQ,PolynomialRing,pari
    from sage.version import version
    spec=r.read(PROTOCOL);source=r.read(old.CONTROL)
    assert source["status"]=="PASS" and source["polynomial_ascending"]==spec["field_polynomial_ascending"]
    R=PolynomialRing(QQ,"y");f=pari(R(spec["field_polynomial_ascending"]));nf=pari.nfinit(f)
    gamma=pari.Mod(pari(R(spec["multipliers_ascending"][index])),f)
    assert abs(pari.nfeltnorm(nf,gamma))==1
    relative_polynomial=pari("x")**2-gamma
    relative=pari.rnfinit(nf,relative_polynomial,1);absolute=pari.nfinit(relative)
    norm_setup=pari.rnfisnorminit(nf,relative_polynomial,1)
    enc=lambda a:[str(pari.lift(a).polcoef(i)) for i in range(3)]
    betas=[pari.Mod(pari(R(w["beta"])),f) for w in source["norm_witnesses"]]
    witnesses=[]
    for j,beta in enumerate(betas):
        z,remainder=pari.rnfisnorm(norm_setup,beta)
        assert remainder==1
        a,b=[z.lift().polcoef(i) for i in range(2)]
        assert a*a-gamma*b*b==beta
        za=pari.rnfeltreltoabs(relative,z);fact=pari.idealfactor(absolute,za)
        primes=sorted({int(P[0]) for P in fact[0]}-set(spec["S_finite"]))
        parity_records=[];values=[0]*len(betas)
        for p in primes:
            for P in pari.idealprimedec(nf,p):
                over=pari.rnfidealprimedec(relative,P)
                exponents=[int(pari.idealval(absolute,za,Q)) for Q in over]
                assert len({e%2 for e in exponents})==1
                parity=exponents[0]%2
                bits=[int(pari.nfislocalpower(nf,P,psi,2)==0) for psi in betas]
                values=[v^(parity*bit) for v,bit in zip(values,bits)]
                parity_records.append({"p":p,"prime_hnf":str(pari.idealhnf(nf,P)),"residue_degree":int(P[3]),
                                       "valuations_above":exponents,"parity":parity,"character_bits":bits})
        witness={"beta":enc(beta),"a":enc(a),"b":enc(b),"parity_records":parity_records,"cup_row":values}
        witnesses.append(witness)
        r.write_new(WORK/f"multiplier-{index}-class-{j}.json",{"bindings":bindings(),"gamma":enc(gamma),**witness})
        print("checkpoint",index,j,values,flush=True)
    r.write_new(WORK/f"multiplier-{index}.json",{"bindings":bindings(),"status":"COUNTED","index":index,
                "gamma":enc(gamma),"witnesses":witnesses,"matrix":[w["cup_row"] for w in witnesses],
                "software":{"sage":version,"pari":str(pari.version())}})


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for index in range(2):
        path=WORK/f"multiplier-{index}.json"
        if not path.exists():
            with (WORK/f"multiplier-{index}.log").open("x") as log:
                try:
                    result=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker","--index",str(index)],
                                          cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    error=None if result.returncode==0 else "worker failed"
                except subprocess.TimeoutExpired:error="30-second timeout"
            if error:r.write_new(path,{"bindings":bindings(),"status":"UNKNOWN","index":index,"reason":error,
                                      "transcript":(WORK/f"multiplier-{index}.log").read_text()})
        row=r.read(path);assert row["bindings"]==bindings();rows.append(row)
        print("captured",index,row["status"],flush=True)
    expected=r.read(old.CONTROL)["independent_norm_cup_matrix"]
    observed=None
    if all(v["status"]=="COUNTED" for v in rows):
        observed=[[rows[0]["matrix"][i][j]^rows[1]["matrix"][i][j] for j in range(2)] for i in range(2)]
    r.write_new(OUTPUT,{"schema":"rank-jump.nonscalar-cup-control.v1","bindings":bindings(),"records":rows,
                       "expected_sum":expected,"observed_sum":observed,
                       "bilinearity_status":"UNKNOWN" if observed is None else ("PASS" if observed==expected else "REFUTED")})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["worker","capture"])
    parser.add_argument("--index",type=int);args=parser.parse_args()
    if args.mode=="worker":worker(args.index)
    else:capture()
