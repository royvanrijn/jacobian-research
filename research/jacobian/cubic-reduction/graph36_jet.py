#!/usr/bin/env python3
"""Quartic own-variable jets can obstruct target corrections of every degree."""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json
import sympy as s
from graph36 import solve_or_separate
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def projected(f,X,i):
    return {e:c for e,c in s.Poly(s.expand(f),*X).terms() if sum(e)==4 and e[i]>0}
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text());X=s.symbols(d["variables"])
    K=[s.expand(sum(s.Rational(t["coefficient"])*s.prod(v**n for v,n in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    folder=HERE/"jet36";folder.mkdir(exist_ok=True);cases=[]
    for j in range(3,12):
        graph=s.expand(X[j]-K[j]);assert not graph.has(X[j])
        kept={i:s.expand(f.subs(X[j],graph)) for i,f in enumerate(K) if i!=j}
        assert all(f.subs(dict.fromkeys(X,0))==0 for f in kept.values())
        for i,f in kept.items():
            assert [s.diff(f,v).subs(dict.fromkeys(X,0)) for v in X]==[int(k==i) for k in range(12)]
        tests=[]
        for i,f in kept.items():
            target=projected(f,X,i)
            if not target:continue
            allowed=[a for a in kept if a!=i]
            tuples=[ix for degree in range(1,4) for ix in combinations_with_replacement(allowed,degree)]
            columns=[projected(s.prod(kept[a] for a in ix),X,i) for ix in tuples]
            result,nrows=solve_or_separate(columns,target)
            tests.append(dict(component=i+1,allowed_outputs=[a+1 for a in allowed],coefficient_rows=nrows,columns=len(columns),**result))
            print(json.dumps(dict(pivot=j+1,component=i+1,status=result["status"],rank=result["rank"],augmented_rank=result["augmented_rank"],functional_terms=len(result.get("functional",[])))),flush=True)
            if result["status"]=="separated":break
        result=dict(status="arbitrary-degree-stage-obstructed" if any(t["status"]=="separated" for t in tests) else "jet4-inconclusive",
            pivot=j+1,tests=tests,source_sha256=sha(Path(__file__)),linear_solver_sha256=sha(HERE/"graph36.py"),input_sha256=sha(source))
        (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
        cases.append(dict(pivot=j+1,status=result["status"],file="pivot-%02d.json"%(j+1)))
        receipt=dict(status="running",source_sha256=sha(Path(__file__)),linear_solver_sha256=sha(HERE/"graph36.py"),
            input_sha256=sha(source),cases=cases,
            scope="A separated quartic coefficient functional involving the missing raw coordinate annihilates all target polynomials: degree>=4 has no such quartic terms. This obstructs one raw-output target stage, not general sequential changes.")
        (folder/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    receipt["status"]="complete";(folder/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

