#!/usr/bin/env python3
"""Bounded mixed-character K3 counts on one matched production pair."""
import argparse
import base64
import struct
import subprocess
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"MIXED_CHARACTER_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_mixed_character_pair_inputs_v1.json"
COUNTS=r.OUT/"rank_jump_mixed_character_counts_v2.json"
OUTPUT=r.OUT/"rank_jump_mixed_character_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-mixed-character-v2"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT)}


def field_modulus(p,degree):
    if degree==1:return [0,1]
    if degree==2:
        n=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1)
        return [(-n)%p,0,1]
    for c0 in range(1,p):
        for c1 in range(p):
            if all((x**3+c1*x+c0)%p for x in range(p)):
                return [c0,c1,0,1]
    raise AssertionError("no cubic found")


def count(case,p):
    from sage.all import GF,PolynomialRing,EllipticCurve
    from sage.version import version
    data=r.read(INPUT)["cases"][case];A,B=map(r.F,data["model"][3:])
    a,b=[r.F(P[0]) for P in data["generic_points"]]
    A,B,a,b=[r.mod(v,p) for v in (A,B,a,b)]
    assert B*(-4*A**3-27*B*B)*a*b*(a-b)*(a**3+A*a+B)*(b**3+A*b+B)%p
    rows=[]
    for degree in (1,2,3):
        modulus=field_modulus(p,degree);Fp=GF(p);R=PolynomialRing(Fp,"z")
        F=Fp if degree==1 else GF(p**degree,name="w",modulus=R(modulus))
        w=F(0) if degree==1 else F.gen();q=p**degree
        def element(i):
            v=F(0);power=F(1)
            for _ in range(degree):
                v+=F(i%p)*power;i//=p;power*=w
            return v
        def encode(v):
            if degree==1:return int(v)
            return sum(int(c)*p**i for i,c in enumerate(v.polynomial().list()))
        roots=PolynomialRing(F,"x")([B,A,0,1]).roots(multiplicities=False)
        traces=[None]*q;calls=0;singular=0
        for i in range(q):
            if traces[i] is not None:continue
            u=element(i);D=1+F(A)*u*u+F(B)*u**3
            if D==0:
                theta=next(t for t in roots if u*t==1)
                single=2*theta;double=-F(A)*u-theta
                assert single!=double
                value=-1 if (double-single).is_square() else 1;singular+=1
            else:
                E=EllipticCurve(F,[0,2*F(A)*u,0,F(A)+3*F(B)*u+F(A*A)*u*u,
                                  F(B)+F(A*B)*u*u-F(B*B)*u**3])
                value=int(E.cardinality(algorithm="pari")-q-1);calls+=1
            orbit=u
            while traces[encode(orbit)] is None:
                traces[encode(orbit)]=value;orbit=orbit**p
            if calls and calls%512==0:print("checkpoint",case,p,degree,calls,flush=True)
        total=0
        for i,value in enumerate(traces):
            u=element(i);d=(1-F(a)*u)*(1-F(b)*u)
            chi=0 if d==0 else (1 if d.is_square() else -1)
            total+=chi*value
        raw=struct.pack("<"+"h"*q,*traces)
        row={"p":p,"degree":degree,"q":q,"modulus_ascending":modulus,"roots_of_f":len(roots),
             "residual_H2_trace":total,"surface_point_count":q*q+1+q*(5+4*len(roots))+total,
             "trace_values_i16_le_base64":base64.b64encode(raw).decode(),
             "smooth_cardinality_calls":calls,"singular_checks":singular}
        r.write_new(WORK/f"case{case}-p{p}-degree{degree}.json",{"bindings":bindings(),**row})
        rows.append(row);print("completed",case,p,degree,total,flush=True)
    return {"case":case,"p":p,"status":"COUNTED","bindings":bindings(),"sage":version,"fields":rows}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for case,data in enumerate(r.read(INPUT)["cases"]):
        for p in r.read(PROTOCOL)["primes"][data["role"]]:
            path=WORK/f"case{case}-p{p}.json"
            if not path.exists():
                logpath=WORK/f"case{case}-p{p}.log"
                with logpath.open("x") as log:
                    try:
                        proc=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                             "--case",str(case),"--prime",str(p)],cwd=r.ROOT,
                                            stdout=log,stderr=log,timeout=40)
                        reason=None if proc.returncode==0 else "worker failure"
                    except subprocess.TimeoutExpired:reason="40-second cap"
                if reason:r.write_new(path,{"case":case,"p":p,"status":"UNKNOWN","reason":reason,
                                           "bindings":bindings(),"transcript":logpath.read_text()})
            row=r.read(path);assert row["bindings"]==bindings();rows.append(row)
            print("checkpoint",case,p,row["status"],flush=True)
    r.write_new(COUNTS,{"schema":"rank-jump.mixed-character-counts.v1","bindings":bindings(),"reductions":rows})


def analyze(check=False):
    from sage.all import QQ,PolynomialRing,cyclotomic_polynomial,ZZ
    R=PolynomialRing(QQ,"X");X=R.gen();rows=[]
    data=r.read(COUNTS)
    producer=r.read(r.OUT/"rank_jump_mixed_character_count_source_v1.json")
    assert r.digest(producer["source"].encode())==producer["source_sha256"]
    expected=bindings();expected[str(Path(__file__).relative_to(r.ROOT))]=producer["source_sha256"]
    assert data["bindings"]==expected
    for raw in data["reductions"]:
        row={"case":raw["case"],"p":raw["p"],"status":"UNKNOWN"}
        if raw["status"]=="COUNTED":
            p=raw["p"];t1,t2,t3=[x["residual_H2_trace"] for x in raw["fields"]];options=[]
            for sign in (-1,1):
                s1=t1-sign*p;s2=t2-p*p;c2=QQ(s1*s1-s2)/2
                if t3-sign*p**3==s1**3-3*s1*c2+3*s1*p*p:
                    quartic=X**4-s1*X**3+c2*X*X-s1*p*p*X+p**4
                    options.append((sign,quartic))
            row["traces"]=[t1,t2,t3]
            if len(options)==1:
                sign,quartic=options[0];pol=(X-sign*p)*quartic
                residual=pol;cyclotomic=[]
                for n in (1,2,3,4,5,6,8,10,12):
                    c=R(cyclotomic_polynomial(n));c=R(c(X/p)*p**c.degree())
                    while residual%c==0:
                        residual=residual//c;cyclotomic.append(n)
                row.update({"residual_polynomial_ascending":list(map(str,pol.list())),
                            "cyclotomic_orders":cyclotomic,"noncyclotomic_factor_ascending":list(map(str,residual.list())),
                            "reduction_geometric_Picard_rank":int(17+5-residual.degree())})
                if residual.degree()==4 and len(cyclotomic)==1:
                    # Over F_(p^6), all trivial classes and the +/-p line are fixed.
                    # The product (1-alpha^6/p^6) is the normalized resultant.
                    determinant=quartic.resultant(X**6-p**6)/QQ(p**24)
                    signed_disc=-determinant/QQ(p**6)
                    sf=int(ZZ(signed_disc.numerator()*signed_disc.denominator()).squarefree_part())
                    row.update({"status":"RHO_18_REDUCTION","NS_discriminant_squareclass":sf,
                                "normalized_transcendental_at_one_degree6":str(determinant)})
        rows.append(row)
    cases=[]
    for i,case in enumerate(r.read(INPUT)["cases"]):
        usable=[x for x in rows if x["case"]==i and x["status"]=="RHO_18_REDUCTION"]
        witness=next(([a,b] for j,a in enumerate(usable) for b in usable[j+1:]
                      if a["NS_discriminant_squareclass"]!=b["NS_discriminant_squareclass"]),None)
        cases.append({"id":case["id"],"role":case["role"],"two_place_witness":witness,
                      "geometric_Picard_rank":17 if witness else "UNKNOWN",
                      "mixed_character_geometric_MW_rank":0 if witness else "UNKNOWN",
                      "production_curve_rank":"UNKNOWN"})
    out={"schema":"rank-jump.mixed-character.v1","bindings":bindings(),
         "counts_sha256":r.digest(COUNTS.read_bytes()),"reductions":rows,"cases":cases,
         "boundary":"Function-field character ranks only; no production fibre rank upper bound."}
    if check:assert r.read(OUTPUT)==out;print("PASS mixed-character analysis replay")
    else:r.write_new(OUTPUT,out)
    print(out)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["worker","capture","build","check"])
    p.add_argument("--case",type=int);p.add_argument("--prime",type=int);args=p.parse_args()
    if args.mode=="worker":
        row=count(args.case,args.prime);r.write_new(WORK/f"case{args.case}-p{args.prime}.json",row)
    elif args.mode=="capture":capture()
    else:analyze(args.mode=="check")
