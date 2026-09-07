#!/usr/bin/env python3
"""Independent rational arithmetic and Hensel replay; no number-field CAS."""
import argparse
from fractions import Fraction as F
from itertools import permutations
from math import lcm
from pathlib import Path
import retrospective as r
import nonscalar_cup_control as ex
import verify_scalar_cup_control as local

OUTPUT=r.OUT/"rank_jump_nonscalar_cup_control_verification_v1.json"
POLY=[-1,-14,-11,1]


def mul(a,b):
    out=[F(0)]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    for degree in (4,3):
        for j in range(3):out[degree-3+j]-=out[degree]*POLY[j]
    return out[:3]


def norm(a):
    columns=[mul(a,[F(i==j) for i in range(3)]) for j in range(3)]
    value=F(0)
    for perm in permutations(range(3)):
        term=F((-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))
        for j in range(3):term*=columns[j][perm[j]]
        value+=term
    return value


def prime_factors(n):
    n=abs(n);p=2;out=[]
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out


def verify(check=False):
    data=r.read(ex.OUTPUT);source=r.read(ex.old.CONTROL);spec=r.read(ex.PROTOCOL)
    assert data["bindings"]==ex.bindings() and data["bilinearity_status"]=="PASS"
    assert source["certified_class_group"]==[2,2]
    betas=[w["beta"] for w in source["norm_witnesses"]]
    precision=spec["limits"]["independent_valuation_precision"]
    records=[];matrices=[]
    for record in data["records"]:
        gamma=list(map(F,record["gamma"]))
        assert abs(norm(gamma))==1 and all(x.denominator==1 for x in gamma)
        rows=[];witnesses=[]
        for witness in record["witnesses"]:
            beta,a,b=[list(map(F,witness[key])) for key in ("beta","a","b")]
            aa=mul(a,a);gbb=mul(gamma,mul(b,b))
            assert [x-y for x,y in zip(aa,gbb)]==beta and norm(beta)==625
            denominator=lcm(*(x.denominator for x in a+b))
            primes=sorted(set(prime_factors(denominator)+[5])-set(spec["S_finite"]))
            # Integral coefficients and norm(beta)=5^4 exclude support elsewhere.
            values=[0,0];details=[]
            for p in primes:
                roots=[local.hensel(POLY,t,p,precision) for t in range(p)
                       if sum(c*t**i for i,c in enumerate(POLY))%p==0]
                assert len(roots)==3 and p not in (2,163)
                for theta in roots:
                    v_beta,_=local.local_value(beta,theta,p,precision)
                    assert v_beta%2==0
                    g=sum(int(c)*pow(theta,i,p**precision) for i,c in enumerate(gamma))%(p**precision)
                    assert g%p
                    square_roots=[s for s in range(p) if s*s%p==g%p]
                    if square_roots:
                        sqrt_g=[local.hensel([-g,0,1],s,p,precision) for s in square_roots]
                        valuations=[local.norm_value(a,b,theta,s,p,precision) for s in sqrt_g]
                        assert len(valuations)==2 and sum(valuations)==v_beta
                    else:
                        # Unramified quadratic extension: norm valuation is twice valuation.
                        valuations=[v_beta//2]
                    assert len({v%2 for v in valuations})==1
                    parity=valuations[0]%2;bits=[]
                    for psi in betas:
                        v,unit=local.local_value(psi,theta,p,precision)
                        assert v%2==0
                        bits.append(int(pow(unit,(p-1)//2,p)==p-1))
                    values=[v^(parity*b) for v,b in zip(values,bits)]
                    details.append({"p":p,"theta_mod_p":theta%p,"extension":"split" if square_roots else "inert",
                                    "valuations":sorted(valuations),"parity":parity,"character_bits":bits})
            canonical=lambda p,v,e,b:(p,tuple(sorted(v)),e,tuple(b))
            independently=sorted(canonical(v["p"],v["valuations"],v["parity"],v["character_bits"]) for v in details)
            retained=sorted(canonical(v["p"],v["valuations_above"],v["parity"],v["character_bits"]) for v in witness["parity_records"])
            assert independently==retained and values==witness["cup_row"]
            rows.append(values);witnesses.append({"support_primes":primes,"places":details,"cup_row":values})
        assert rows==record["matrix"]
        matrices.append(rows);records.append({"index":record["index"],"gamma_norm":str(norm(gamma)),"witnesses":witnesses,"matrix":rows})
    summed=[[matrices[0][i][j]^matrices[1][i][j] for j in range(2)] for i in range(2)]
    assert summed==source["independent_norm_cup_matrix"]==data["expected_sum"]
    out={"schema":"rank-jump.nonscalar-cup-control-verification.v1","status":"PASS","records":records,
         "matrix_sum":summed,"analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),
         "verifier_sha256":r.digest(Path(__file__).read_bytes()),
         "Hensel_source_sha256":r.digest(Path(local.__file__).read_bytes()),
         "source_control_sha256":r.digest(ex.old.CONTROL.read_bytes()),"precision":precision,
         "boundary":"Complete parity support verified by integrality and norm, not a bounded list of guessed primes."}
    if check:assert out==r.read(OUTPUT)
    else:r.write_new(OUTPUT,out)
    print("PASS all four rational norm identities, complete parity supports, and nonscalar cup bilinearity")
    print(matrices)
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    verify(parser.parse_args().check)
