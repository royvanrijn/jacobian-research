#!/usr/bin/env python3
"""Direct rational verification of two cubic index-four fibre obstructions."""
if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertion checks are required.")
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import permutations
import argparse, copy, hashlib, json

HERE=Path(__file__).resolve().parent

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def validate(case):
    assert set(case)=={"label","dimension","coordinate_grades","H_terms"}
    n=case["dimension"]
    assert type(n) is int and n in (7,9)
    g=case["coordinate_grades"]
    assert len(g)==n and all(type(x) is int and x>=0 for x in g)
    assert g[0]==0 and g.count(0)==1 and g==sorted(g) and len(set(g))==n
    seen=set()
    for t in case["H_terms"]:
        assert set(t)=={"output","inputs","coefficient"}
        o=t["output"]; ids=t["inputs"]
        assert type(o) is int and 0<=o<n
        assert isinstance(ids,list) and len(ids)==3
        assert all(type(i) is int and 0<=i<o for i in ids)
        assert ids==sorted(ids) and (o,tuple(ids)) not in seen
        seen.add((o,tuple(ids)))
        assert isinstance(t["coefficient"],str) and F(t["coefficient"])!=0
        assert g[o]==1+sum(g[i] for i in ids)
    return n,[2*x+1 for x in g]

def jacobian(terms):
    out={}
    for t in terms:
        ids=t["inputs"]
        for i,m in Counter(ids).items():
            mon=list(ids);mon.remove(i)
            p=out.setdefault((t["output"],i),{})
            key=tuple(mon);p[key]=p.get(key,F(0))+F(t["coefficient"])*m
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}

def multiply(a,b):
    byrow={}
    for (r,c),p in b.items():byrow.setdefault(r,[]).append((c,p))
    out={}
    for (i,k),p in a.items():
        for j,q in byrow.get(k,[]):
            r=out.setdefault((i,j),{})
            for u,c in p.items():
                for v,d in q.items():
                    m=tuple(sorted(u+v));r[m]=r.get(m,F(0))+c*d
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}

def mat_eval(a,point,n):
    out=[[F(0) for _ in range(n)] for _ in range(n)]
    for (i,j),p in a.items():
        for mon,c in p.items():
            v=c
            for k in mon:v*=point[k]
            out[i][j]+=v
    return out

def act(a,v):
    return [sum((c*x for c,x in zip(row,v)),F(0)) for row in a]

def mu(terms,a,b,c,n):
    out=[F(0)]*n
    for t in terms:
        # Each ordered slot assignment has coefficient cubic_coefficient/6.
        val=sum((a[i]*b[j]*c[k] for i,j,k in permutations(t["inputs"])),F(0))/6
        out[t["output"]]+=F(t["coefficient"])*val
    return out

def run(case):
    n,w=validate(case);terms=case["H_terms"]
    j=jacobian(terms);j2=multiply(j,j);j3=multiply(j2,j);j4=multiply(j3,j)
    assert j3 and not j4
    e0=[F(int(i==0)) for i in range(n)];e3=[F(int(i==3)) for i in range(n)]
    actual=act(mat_eval(j3,e0,n),e0)
    assert actual==[F(27*int(i==3)) for i in range(n)]
    l=e0
    for _ in range(3):l=mu(terms,e0,e0,l,n)
    assert l==e3
    product=mu(terms,e3,e3,e0 if n==7 else e3,n)
    assert product==[F((1 if n==7 else -1)*int(i==n-1)) for i in range(n)]
    v=[F(1)]+[F(0)]*(n-1)
    for o in range(1,n):
        for t in terms:
            if t["output"]==o:
                a,b,c=t["inputs"];v[o]+=F(t["coefficient"])*v[a]*v[b]*v[c]
    expected=([1,1,3,9,27,-81,0] if n==7 else [1,1,3,9,27,27,-162,243,729])
    assert v==list(map(F,expected))
    if n==7:
        assert w[-1]==15 and w.count(1)==1 and w.count(15)==1 and not v[-1]
        assert w[-2]==13 and v[-2]==-81 and all(x%2 for x in w)
    else:assert w[-1]==21 and v[-1]==729
    perturbed=copy.deepcopy(case)
    # Change one actual coefficient in each literal map.
    k=2 if n==7 else 0
    perturbed["H_terms"][k]["coefficient"]=str(F(perturbed["H_terms"][k]["coefficient"])+1)
    p=jacobian(perturbed["H_terms"]);p2=multiply(p,p);p4=multiply(multiply(p2,p),p)
    assert p4
    return {"label":case["label"],"dimension":n,"cubic_terms":len(terms),
            "jacobian_cube_nonzero_entries":len(j3),"jacobian_fourth_power_entries":len(j4),
            "L_e0_cubed_e0":[str(x) for x in l],"obstruction_product":[str(x) for x in product],
            "inverse_line_coefficients":[str(x) for x in v],"inverse_degree":13 if n==7 else 21,
            "perturbed_fourth_power_nonzero_entries":len(p4)}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=HERE/"engel_four_fibre_obstructions.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"engel_four_fibre_obstructions_source_ledger.json")
    ap.add_argument("--output",type=Path,default=HERE/"engel_four_fibre_obstructions_verification.json")
    args=ap.parse_args()
    ledger=json.loads(args.ledger.read_text())
    assert ledger["schema_version"]==1 and len(ledger["sources"])==3
    expected={"proof","literal_input","checker"};assert {s["role"] for s in ledger["sources"]}==expected
    for s in ledger["sources"]:
        p=(args.ledger.parent/s["path"]).resolve()
        assert p.is_file() and digest(p)==s["sha256"],"Source pin mismatch: "+str(p)
        if s["role"]=="literal_input":assert p==args.input.resolve()
        if s["role"]=="checker":assert p==Path(__file__).resolve()
    data=json.loads(args.input.read_text())
    assert set(data)=={"schema_version","cases"} and data["schema_version"]==1
    assert len(data["cases"])==2 and [x["dimension"] for x in data["cases"]]==[7,9]
    results=[run(c) for c in data["cases"]]
    rejected=0
    for mutate in ("duplicate","zero","index","grade"):
        c=copy.deepcopy(data["cases"][0])
        if mutate=="duplicate":c["H_terms"].append(copy.deepcopy(c["H_terms"][0]))
        elif mutate=="zero":c["H_terms"][0]["coefficient"]="0"
        elif mutate=="index":c["H_terms"][0]["inputs"][0]=7
        else:c["coordinate_grades"][3]=4
        try:validate(c)
        except (AssertionError,ValueError):rejected+=1
    assert rejected==4
    receipt={"schema_version":1,"status":"PASS","sources_sha256":digest(args.ledger),
             "cases":results,"malformed_inputs_rejected":rejected}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()
