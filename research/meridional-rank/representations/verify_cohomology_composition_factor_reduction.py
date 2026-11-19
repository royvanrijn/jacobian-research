#!/usr/bin/env python3
"""Exact finite controls for composition-factor reduction of the GGLW bound."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys

if not __debug__:
    raise SystemExit("Run without -O: assertion controls are required.")
HERE=Path(__file__).resolve().parent

def digest(p):return sha256(p.read_bytes()).hexdigest()
def eye(n):return [[int(i==j)for j in range(n)]for i in range(n)]
def add(A,B,p,sign=1):return [[(a+sign*b)%p for a,b in zip(r,s)]for r,s in zip(A,B)]
def mul(A,B,p):return [[sum(a*b for a,b in zip(row,col))%p for col in zip(*B)]for row in A]
def flat(A):return tuple(x for row in A for x in row)

def rref(A,p):
    A=[[v%p for v in row]for row in A]
    if not A:return A,[]
    piv=[];r=0
    for c in range(len(A[0])):
        j=next((j for j in range(r,len(A))if A[j][c]),None)
        if j is None:continue
        A[r],A[j]=A[j],A[r]
        inv=pow(A[r][c],-1,p);A[r]=[(v*inv)%p for v in A[r]]
        for j in range(len(A)):
            if j!=r and A[j][c]:
                z=A[j][c];A[j]=[(x-z*y)%p for x,y in zip(A[j],A[r])]
        piv.append(c);r+=1
        if r==len(A):break
    return A,piv

def rank(A,p):return len(rref(A,p)[1])
def kernel(A,p,n=None):
    if n is None:n=len(A[0])
    R,piv=rref(A,p);out=[]
    for c in range(n):
        if c in piv:continue
        v=[0]*n;v[c]=1
        for j,k in enumerate(piv):v[k]=(-R[j][c])%p
        assert all(sum(a*b for a,b in zip(row,v))%p==0 for row in A)
        out.append(v)
    return out

def certificate(delta,q):
    return 1 if not q else max(1,1+(-((-delta)//q)))

def module(A,B,p):
    d=len(A);I=eye(d)
    assert all(len(row)==d for row in A+B)
    assert rank(A,p)==rank(B,p)==d
    assert mul(mul(A,B,p),A,p)==mul(mul(B,A,p),B,p)
    # Derivatives of aba=bab: no braid propagation implementation is used.
    F=add(add(I,mul(A,B,p),p),B,p,-1)
    G=add(add(A,I,p,-1),mul(B,A,p),p,-1)
    equations=[r+s for r,s in zip(F,G)]
    basis=kernel(equations,p)
    zdim=len(basis);delta=zdim-d;q=rank(add(A,I,p,-1),p)
    return {"dimension":d,"zdim":zdim,"delta":delta,"q":q,
            "certificate":certificate(delta,q),"kernel":basis}

def extension(A,B,C,D,T,U):
    a=len(A);c=len(C)
    E=[row+tail for row,tail in zip(A,T)]+[[0]*a+row for row in C]
    F=[row+tail for row,tail in zip(B,U)]+[[0]*a+row for row in D]
    return E,F

def projection_check(v,w,e,p):
    dv,dw=v["dimension"],w["dimension"];de=dv+dw
    projected=[z[dv:de]+z[de+dv:2*de]for z in e["kernel"]]
    image_rank=rank(projected,p)
    assert e["zdim"]-image_rank==v["zdim"]
    assert image_rank<=w["zdim"]
    assert e["delta"]<=v["delta"]+w["delta"]
    assert e["q"]>=v["q"]+w["q"]
    assert e["certificate"]<=max(v["certificate"],w["certificate"])
    return image_rank

def irreducible_two(A,B,p):
    # Every possible proper nonzero submodule in dimension2 is a line.
    lines=[(1,t)for t in range(p)]+[(0,1)]
    for v in lines:
        images=[[sum(a*b for a,b in zip(row,v))%p for row in X]for X in [A,B]]
        if all(rank([list(v),u],p)==1 for u in images):return False
    return True

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo-root",type=Path,default=Path.cwd())
    parser.add_argument("--output",type=Path,default=HERE/"cohomology_composition_factor_verification.json")
    args=parser.parse_args()
    inputs_path=HERE/"cohomology_composition_factor_inputs.json"
    ledger_path=HERE/"cohomology_composition_factor_source_ledger.json"
    inputs=json.loads(inputs_path.read_text(encoding="utf-8-sig"))
    ledger=json.loads(ledger_path.read_text(encoding="utf-8-sig"))
    pins=[]
    for source in ledger["sources"]:
        path=args.repo_root/source["path"]
        assert digest(path)==source["sha256"]
        pins.append({"path":source["path"],"sha256":digest(path),"role":source["role"]})
    counts={"character_extensions":0,"natural_trivial_extensions":0,
            "natural_natural_extensions":0,"strict_rank_superadditivity":0,
            "strict_cocycle_projection":0,"strict_certificate_loss":0,"strict_normalized_defect_loss":0,
            "zero_defect_all_trivial_nonzero_rank":0,"root_affine_preserved_rank":0}
    irreducible_records=[];jordan=None
    all_modules=[];factor_modules={};rounding_example=None
    def check(A,B,C,D,E,F,p,kind):
        nonlocal rounding_example
        v,w,e=module(A,B,p),module(C,D,p),module(E,F,p)
        image=projection_check(v,w,e,p)
        all_modules.extend([v,w,e])
        factor_modules[(p,flat(A),flat(B))]=v
        factor_modules[(p,flat(C),flat(D))]=w
        if e["delta"]==1 and e["q"]==2 and rounding_example is None:
            rounding_example={"p":p,"A":E,"B":F,"delta":1,"q":2,"certificate":2}
        counts[kind]+=1
        counts["strict_rank_superadditivity"]+=int(e["q"]>v["q"]+w["q"])
        counts["strict_cocycle_projection"]+=int(image<w["zdim"])
        counts["strict_certificate_loss"]+=int(e["certificate"]<max(v["certificate"],w["certificate"]))
        ratio=lambda x:max(Fraction(0),Fraction(x["delta"],x["q"])) if x["q"] else Fraction(0)
        counts["strict_normalized_defect_loss"]+=int(ratio(e)<max(ratio(v),ratio(w)))
        if v["q"]==w["q"]==0 and e["q"]>0 and e["delta"]==0:
            counts["zero_defect_all_trivial_nonzero_rank"]+=1
        return v,w,e
    for p in inputs["primes"]:
        for lam,mu,s,t in product(range(1,p),range(1,p),range(p),range(p)):
            A=B=[[lam]];C=D=[[mu]]
            E,F=extension(A,B,C,D,[[s]],[[t]])
            if mul(mul(E,F,p),E,p)!=mul(mul(F,E,p),F,p):continue
            v,w,e=check(A,B,C,D,E,F,p,"character_extensions")
            if p==2 and lam==mu==s==t==1:
                jordan={"p":p,"A":E,"B":F,**{k:e[k]for k in ["dimension","zdim","delta","q","certificate"]}}
        A=[[1,1],[0,1]];B=[[1,0],[p-1,1]];C=D=[[1]]
        assert irreducible_two(A,B,p)
        v=module(A,B,p)
        irreducible_records.append({"p":p,"A":A,"B":B,**{k:v[k]for k in ["dimension","zdim","delta","q","certificate"]}})
        for coeff in product(range(p),repeat=v["zdim"]):
            z=[sum(c*b[i]for c,b in zip(coeff,v["kernel"]))%p for i in range(4)]
            T=[[z[0]],[z[1]]];U=[[z[2]],[z[3]]]
            E,F=extension(A,B,C,D,T,U)
            vv,w,e=check(A,B,C,D,E,F,p,"natural_trivial_extensions")
            assert e["delta"]<=vv["delta"] and e["certificate"]<=vv["certificate"]
            if e["q"]==vv["q"]:counts["root_affine_preserved_rank"]+=1
    p=2;A=[[1,1],[0,1]];B=[[1,0],[1,1]]
    for vals in product(range(2),repeat=8):
        T=[list(vals[:2]),list(vals[2:4])];U=[list(vals[4:6]),list(vals[6:8])]
        E,F=extension(A,B,A,B,T,U)
        if mul(mul(E,F,p),E,p)!=mul(mul(F,E,p),F,p):continue
        check(A,B,A,B,E,F,p,"natural_natural_extensions")
    assert jordan and jordan["delta"]==0 and jordan["q"]==1 and jordan["zdim"]==2
    assert all(counts[k]>0 for k in counts if k!="strict_certificate_loss")
    # Explicit false strengthening: delta>=0 does not force a positive-rank factor.
    assert jordan["delta"]>=0 and jordan["q"]>0
    jordan["composition_factor_ranks"]=[0,0]
    # A numerical division control, not an assertion of a realized knot module.
    assert Fraction(-2,3)>Fraction(-1,1)
    assert certificate(-2,3)==certificate(-1,1)==1
    caps=[]
    for cap_text in inputs["uniform_cap_values"]:
        cap=Fraction(cap_text);integer_cap=cap.numerator//cap.denominator
        all_bound=all(v["certificate"]<=cap for v in all_modules)
        factor_bound=all(v["certificate"]<=cap for v in factor_modules.values())
        numeric=(cap>=1 and all(v["delta"]<=(integer_cap-1)*v["q"]
                     for v in factor_modules.values() if v["q"]>0))
        assert all_bound==factor_bound==numeric
        caps.append({"B":str(cap),"floor_B":integer_cap,
                     "all_fixtures_bound":all_bound,
                     "all_irreducible_factors_bound":factor_bound,
                     "exact_numeric_criterion":numeric})
    trivial_caps=[]
    trivial=next(v for v in factor_modules.values() if v["q"]==0)
    assert trivial["dimension"]==1 and trivial["delta"]==0 and trivial["certificate"]==1
    for cap_text in inputs["uniform_cap_values"]:
        cap=Fraction(cap_text)
        all_bound=(jordan["certificate"]<=cap and trivial["certificate"]<=cap)
        irred_bound=trivial["certificate"]<=cap
        assert all_bound==irred_bound==(cap>=1)
        trivial_caps.append({"B":str(cap),"all_trivial_factor_modules_bound":all_bound,
                             "trivial_irreducible_bound":irred_bound,
                             "baseline_guard":cap>=1})
    assert not(trivial["certificate"]<=Fraction(0))  # empty positive-factor condition is insufficient
    assert rounding_example is not None
    cap=Fraction(3,2)
    assert Fraction(rounding_example["delta"],rounding_example["q"])<=cap-1
    assert rounding_example["certificate"]>cap
    assert not(rounding_example["delta"]<=(int(cap)-1)*rounding_example["q"])
    # Corrupt matrix fixture must fail the representation relation.
    try:module([[1,1],[0,1]],[[1,1],[1,0]],2)
    except AssertionError:corrupt_rejected=True
    else:raise AssertionError("Nonrepresentation corruption accepted")
    packet_names=["COHOMOLOGY_COMPOSITION_FACTOR_REDUCTION.md",Path(__file__).name,
                  inputs_path.name,ledger_path.name]
    receipt={"status":"PASS","theorem":"CohomologyCompositionFactorReduction",
       "scope":"Every field; finite-dimensional modules of finitely generated G normally generated by m with G_ab=Z. Only the clipped GGLW cohomology certificate is compared.",
       "finite_control_group":"<a,b | aba=bab>; abelianization Z and a normally generates",
       "counts":counts,"irreducible_fixtures":irreducible_records,
       "zero_defect_counterexample":jordan,
       "uniform_cap_controls":caps,"trivial_factor_uniform_cap_controls":trivial_caps,
       "noninteger_cap_counterexample":rounding_example,
       "negative_controls":{"nonrepresentation_rejected":corrupt_rejected,
         "rank_additivity_refuted":counts["strict_rank_superadditivity"]>0,
         "universal_cocycle_surjectivity_refuted":counts["strict_cocycle_projection"]>0,
         "zero_defect_positive_factor_strengthening_refuted":True,
         "noninteger_cap_without_floor_refuted":True,
         "vacuous_positive_factor_condition_below_one_refuted":True,
         "negative_ratio_division_guard":"-2/3 > -1; both clipped certificates1; numerical only"},
       "source_pins":pins,"packet_sha256":{name:digest(HERE/name)for name in packet_names},
       "runtime":{"python":sys.version.split()[0],"standard_library_only":True},
       "assurance":"Exact bounded controls supplement a self-contained unrestricted proof. No formal verification or external human review."}
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","counts":counts,"output":str(args.output)},sort_keys=True))

if __name__=="__main__":main()


