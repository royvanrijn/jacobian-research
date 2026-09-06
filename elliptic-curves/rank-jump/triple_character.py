#!/usr/bin/env python3
"""All four mixed characters of the fixed three-generic-point genus-one cover."""
import argparse
import base64
from pathlib import Path
import struct
import subprocess
import retrospective as r
import mixed_character as previous

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"TRIPLE_CHARACTER_PROTOCOL.json"
INPUT=r.OUT/"rank_jump_triple_character_inputs_v1.json"
RAW=r.OUT/"rank_jump_triple_character_counts_v1.json"
OUTPUT=r.OUT/"rank_jump_triple_character_v1.json"
WORK=r.ROOT/"artifacts/local/rank-jump-triple-character-v1"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,INPUT,previous.COUNTS,Path(previous.__file__))}


def new_count(case,p):
    from sage.all import GF,PolynomialRing,EllipticCurve
    from sage.version import version
    row=r.read(INPUT)["cases"][case];A,B=[r.mod(x,p) for x in row["model"][3:]]
    assert p>3 and B*(-4*A**3-27*B*B)%p
    fields=[]
    for degree in (1,2,3):
        Fp=GF(p);R=PolynomialRing(Fp,"w");modulus=previous.field_modulus(p,degree)
        F=Fp if degree==1 else GF(p**degree,name="w",modulus=R(modulus))
        w=F(0) if degree==1 else F.gen();q=p**degree
        def elt(i):return sum(F((i//p**k)%p)*w**k for k in range(degree))
        def encode(v):
            return int(v) if degree==1 else sum(int(c)*p**i for i,c in enumerate(v.polynomial().list()))
        roots=PolynomialRing(F,"x")([B,A,0,1]).roots(multiplicities=False)
        traces=[None]*q;calls=0;nodes=0
        for i in range(q):
            if traces[i] is not None:continue
            u=elt(i);D=1+F(A)*u*u+F(B)*u**3
            if D==0:
                t=next(x for x in roots if u*x==1)
                double=-F(A)*u-t;single=2*t;assert double!=single
                value=-1 if (double-single).is_square() else 1;nodes+=1
            else:
                E=EllipticCurve(F,[0,2*F(A)*u,0,F(A)+3*F(B)*u+F(A*A)*u*u,
                                  F(B)+F(A*B)*u*u-F(B*B)*u**3])
                value=int(E.cardinality(algorithm="pari")-q-1);calls+=1
            orbit=u
            while traces[encode(orbit)] is None:
                traces[encode(orbit)]=value;orbit=orbit**p
        field={"p":p,"degree":degree,"q":q,"modulus_ascending":modulus,"roots_of_f":len(roots),
               "trace_values_i16_le_base64":base64.b64encode(struct.pack("<"+"h"*q,*traces)).decode(),
               "smooth_cardinality_calls":calls,"singular_checks":nodes}
        fields.append(field)
        r.write_new(WORK/f"case{case}-p{p}-degree{degree}.json",{"bindings":bindings(),**field})
        print("checkpoint",case,p,degree,flush=True)
    return {"case":case,"p":p,"status":"COUNTED","fields":fields,"sage":version,"bindings":bindings()}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);records=[]
    old=r.read(previous.COUNTS)
    assert r.read(INPUT)["previous_counts_sha256"]==r.digest(previous.COUNTS.read_bytes())
    for case,row in enumerate(r.read(INPUT)["cases"]):
        for p in r.read(PROTOCOL)["primes"][row["role"]]:
            retained=next((x for x in old["reductions"] if x["case"]==case and x["p"]==p),None)
            if retained is not None:
                assert retained["status"]=="COUNTED"
                records.append({"case":case,"p":p,"status":"COUNTED","fields":retained["fields"],
                                "origin":"reused previous exact untwisted arrays",
                                "previous_counts_sha256":r.digest(previous.COUNTS.read_bytes())})
                continue
            assert p in r.read(PROTOCOL)["new_count_primes"]
            path=WORK/f"case{case}-p{p}.json"
            if not path.exists():
                with (WORK/f"case{case}-p{p}.log").open("x") as log:
                    try:
                        proc=subprocess.run(["sage","-python",str(Path(__file__).resolve()),"worker",
                                             "--case",str(case),"--prime",str(p)],cwd=r.ROOT,
                                            stdout=log,stderr=log,timeout=40)
                        error=None if proc.returncode==0 else "worker failure"
                    except subprocess.TimeoutExpired:error="40-second cap"
                if error:r.write_new(path,{"case":case,"p":p,"status":"UNKNOWN","reason":error,"bindings":bindings(),
                                          "transcript":(WORK/f"case{case}-p{p}.log").read_text()})
            value=r.read(path);assert value["bindings"]==bindings();records.append(value)
            print("captured",case,p,value["status"],flush=True)
    r.write_new(RAW,{"schema":"rank-jump.triple-character-counts.v1","bindings":bindings(),"records":records})


def good(row,p,mask):
    try:
        A,B=[r.mod(x,p) for x in row["model"][3:]]
        xs=[r.mod(P[0],p) for i,P in enumerate(row["generic_points"]) if mask>>i&1]
    except ValueError:return False
    return bool(B*(-4*A**3-27*B*B)%p and all(x*(x**3+A*x+B)%p for x in xs) and len(set(xs))==len(xs))


def recover(p,traces):
    from sage.all import QQ,ZZ,PolynomialRing,cyclotomic_polynomial
    R=PolynomialRing(QQ,"X");X=R.gen();t1,t2,t3=traces;options=[]
    for e in (-1,1):
        s=t1-e*p;c=QQ(s*s-(t2-p*p))/2
        if t3-e*p**3==s**3-3*s*c+3*s*p*p:
            options.append((e,X**4-s*X**3+c*X*X-s*p*p*X+p**4))
    out={"traces":traces,"status":"UNKNOWN"}
    if len(options)!=1:return out
    e,quartic=options[0];poly=(X-e*p)*quartic;rest=poly;orders=[]
    for n in (1,2,3,4,5,6,8,10,12):
        cyc=R(cyclotomic_polynomial(n));cyc=R(cyc(X/p)*p**cyc.degree())
        while rest%cyc==0:rest=rest//cyc;orders.append(n)
    out.update({"status":"RECOVERED","polynomial_ascending":list(map(str,poly.list())),
                "noncyclotomic_factor_ascending":list(map(str,rest.list())),
                "cyclotomic_orders":orders,"reduction_geometric_Picard_rank":int(22-rest.degree())})
    if rest.degree()==4:
        assert len(orders)==1
        value=quartic.resultant(X**6-p**6)/QQ(p**24)
        signed=-value/QQ(p**6)
        out.update({"status":"RHO_18_REDUCTION","NS_discriminant_squareclass":
                    int(ZZ(signed.numerator()*signed.denominator()).squarefree_part()),
                    "normalized_transcendental_degree6":str(value)})
    return out


def build(check=False):
    from sage.all import GF,PolynomialRing,EllipticCurve,prod
    data=r.read(RAW);assert data["bindings"]==bindings();rows=[];cases=r.read(INPUT)["cases"]
    for raw in data["records"]:
        case,p=raw["case"],raw["p"];row=cases[case]
        for mask in r.read(PROTOCOL)["masks"]:
            if raw["status"]!="COUNTED" or not good(row,p,mask):
                rows.append({"case":case,"p":p,"mask":mask,"status":"INELIGIBLE" if not good(row,p,mask) else "UNKNOWN"})
                continue
            A,B=[r.mod(x,p) for x in row["model"][3:]]
            xs=[r.mod(P[0],p) for i,P in enumerate(row["generic_points"]) if mask>>i&1];fields=[]
            for retained in raw["fields"]:
                degree=retained["degree"];q=retained["q"];Fp=GF(p);R=PolynomialRing(Fp,"w")
                F=Fp if degree==1 else GF(q,name="w",modulus=R(retained["modulus_ascending"]))
                w=F(0) if degree==1 else F.gen()
                trace=struct.unpack("<"+"h"*q,base64.b64decode(retained["trace_values_i16_le_base64"]))
                total=0
                for i,value in enumerate(trace):
                    u=sum(F((i//p**j)%p)*w**j for j in range(degree))
                    d=prod(1-F(x)*u for x in xs)
                    chi=0 if d==0 else (1 if d.is_square() else -1)
                    total+=chi*value
                infinity=0
                if len(xs)==3:
                    k=-prod(F(x) for x in xs)
                    E=EllipticCurve(F,[0,2*F(A)*k,0,F(A*A)*k*k,-F(B*B)*k**3])
                    infinity=int(E.cardinality(algorithm="pari")-q-1)
                fields.append({"degree":degree,"finite_sum":total,"infinity_sum":infinity,
                               "residual_H2_trace":total+infinity,
                               "surface_point_count":q*q+1+q*(5+4*retained["roots_of_f"])+total+infinity})
            rec=recover(p,[x["residual_H2_trace"] for x in fields])
            rows.append({"case":case,"p":p,"mask":mask,"fields":fields,**rec})
    bounds=[]
    for case,row in enumerate(cases):
        chars=[]
        for mask in r.read(PROTOCOL)["masks"]:
            relevant=[x for x in rows if x["case"]==case and x["mask"]==mask]
            upper=min([3]+[x["reduction_geometric_Picard_rank"]-17 for x in relevant if "reduction_geometric_Picard_rank" in x])
            rank18=[x for x in relevant if x["status"]=="RHO_18_REDUCTION"]
            witness=next(([a["p"],b["p"]] for j,a in enumerate(rank18) for b in rank18[j+1:]
                          if a["NS_discriminant_squareclass"]!=b["NS_discriminant_squareclass"]),None)
            if witness:upper=0
            chars.append({"mask":mask,"geometric_rank_interval":[0,int(upper)],"zero_rank_witness_primes":witness})
        extra=sum(x["geometric_rank_interval"][1] for x in chars)
        bounds.append({"id":row["id"],"role":row["role"],"characters":chars,
                       "full_base_geometric_rank_interval":[4,4+extra],
                       "full_base_arithmetic_rank_interval":[3,3+extra],
                       "production_curve_rank":"UNKNOWN"})
    result={"schema":"rank-jump.triple-character.v1","bindings":bindings(),"counts_sha256":r.digest(RAW.read_bytes()),
            "reductions":rows,"bounds":bounds,
            "boundary":"Function-field bounds for the new genus-one base; not rank bounds on the production specializations."}
    if check:assert result==r.read(OUTPUT);print("PASS triple-character replay")
    else:r.write_new(OUTPUT,result)
    for row in bounds:print(row)
    for row in rows:print(row["case"],row["mask"],row["p"],row["status"],row.get("reduction_geometric_Picard_rank"),row.get("NS_discriminant_squareclass"))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["worker","capture","build","check"])
    p.add_argument("--case",type=int);p.add_argument("--prime",type=int);args=p.parse_args()
    if args.mode=="worker":r.write_new(WORK/f"case{args.case}-p{args.prime}.json",new_count(args.case,args.prime))
    elif args.mode=="capture":capture()
    else:build(args.mode=="check")
