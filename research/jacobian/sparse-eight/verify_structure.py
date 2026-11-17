#!/usr/bin/env python3
"""Classify the complete eight-term cover by three necessary coordinate gates."""
from collections import Counter
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
def support_id(rows):
    content=json.dumps(rows,sort_keys=True,separators=(",",":"))
    return "s-"+hashlib.sha256(content.encode()).hexdigest()[:16]

def classify(rows):
    univariate=[i for i,row in enumerate(rows) if row and all(
        all(e[j]==0 for j in range(3) if j!=i) for e in row)]
    if univariate:
        return "univariate-coordinate",univariate[0]
    shear=[i for i,row in enumerate(rows) if all(e[i]==0 for e in row)]
    if shear:
        return "elementary-shear-coordinate",shear[0]
    hyperplane=[i for i,row in enumerate(rows) if row and
                all(e[i]>=1 for e in row) and any(e[i]==1 for e in row)]
    if hyperplane:
        return "hyperplane-normal-factor",hyperplane[0]
    return "unresolved",None

def main():
    path=HERE/"results/supports.json"
    data=json.loads(path.read_text())
    assert data["status"]=="complete" and (data["degree"],data["bound"])==(7,8)
    counts=Counter()
    small=[]
    raw_sizes=Counter()
    for rows in data["supports"]:
        size=sum(map(len,rows))
        tag,coordinate=classify(rows)
        counts[size,tag]+=1
        if tag!="univariate-coordinate":
            entry=dict(support_id=support_id(rows),size=size,gate=tag,
                       coordinate=coordinate,support=rows)
            if tag=="elementary-shear-coordinate":
                delta=max([1]+[sum(e) for e in rows[coordinate]])
                degree=max([1]+[sum(e) for row in rows for e in row])
                entry["fiber_plane_degree_bound"]=degree*delta
                assert degree*delta<=49
            small.append(entry)
        raw_sizes[size]+=1
    expected={(6,"univariate-coordinate"):450,
              (7,"univariate-coordinate"):1875,
              (7,"elementary-shear-coordinate"):1,
              (8,"univariate-coordinate"):17435,
              (8,"elementary-shear-coordinate"):46,
              (8,"hyperplane-normal-factor"):1}
    assert dict(counts)==expected,dict(counts)
    assert len(small)==48 and all(row["gate"]!="unresolved" for row in small)
    residual=[entry for entry in small if entry["gate"]=="hyperplane-normal-factor"]
    assert len(residual)==1
    assert residual[0]["support"]==[
        [[0,0,2],[1,0,1],[3,0,0]],
        [[1,1,0]],
        [[0,0,2],[1,0,1],[2,0,0],[3,0,0]],
    ]
    receipt=dict(
        status="pass",schema="sparse-eight-coordinate-gates.v1",
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supports_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        counts=[dict(size=size,gate=gate,count=n) for (size,gate),n in sorted(counts.items())],
        rare_supports=small,
        rational_algebra_dependency=False,
        external_dependency="Moh plane Keller injectivity through maximum component degree 100",
        scope="Complete retained raw-degree-at-most-seven support cover through eight nonlinear occurrences.",
    )
    (HERE/"structural-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!="rare_supports"}),flush=True)
if __name__=="__main__":
    main()

