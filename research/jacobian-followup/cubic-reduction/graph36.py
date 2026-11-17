#!/usr/bin/env python3
"""Exact graph-deletion screen in the new36-term cubic frame, with dual certificates."""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json,time
import sympy as s
from flint import fmpq_mat
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def high(f,X):return {e:c for e,c in s.Poly(s.expand(f),*X).terms() if sum(e)>3}
def pivots(rr,rank):
    return [next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(rank)]
def solve_or_separate(columns,target):
    monomials=sorted(set(target).union(*(col.keys() for col in columns)))
    rows=[[str(col.get(e,0)) for col in columns]+[str(target.get(e,0))] for e in monomials]
    M=fmpq_mat(rows);rr,rank_aug=M.rref();pp=pivots(rr,rank_aug);count=len(columns)
    if count not in pp:
        solution=[s.Rational(0)]*count
        for i,j in enumerate(pp):solution[j]=s.Rational(str(rr[i,count]))
        return dict(status="span-member",rank=rank_aug,augmented_rank=rank_aug,
          solution=[str(c) for c in solution]),len(monomials)
    # Independent rows of the augmented matrix support a small left-null witness.
    tr,rank_t=M.transpose().rref();selected=pivots(tr,rank_t);assert rank_t==rank_aug
    sub=fmpq_mat([rows[j] for j in selected]);coeff=fmpq_mat([row[:-1] for row in [rows[j] for j in selected]])
    ct,r=coeff.transpose().rref();cp=pivots(ct,r)
    assert r+1==rank_aug
    free=next(j for j in range(len(selected)) if j not in cp)
    lam=[s.Rational(0)]*len(selected);lam[free]=1
    for i,j in enumerate(cp):lam[j]=-s.Rational(str(ct[i,free]))
    value=sum(q*s.Rational(rows[j][-1]) for q,j in zip(lam,selected));assert value
    lam=[q/value for q in lam]
    functional=fmpq_mat([[str(c) for c in lam]])
    product=functional*sub
    assert all(product[0,j]==0 for j in range(count)) and product[0,count]==1
    certificate=[dict(exponent=list(monomials[j]),weight=str(q)) for j,q in zip(selected,lam) if q]
    return dict(status="separated",rank=r,augmented_rank=rank_aug,functional=certificate),len(monomials)
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text())
    X=s.symbols(d["variables"])
    K=[s.expand(sum(s.Rational(t["coefficient"])*s.prod(v**n for v,n in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    folder=HERE/"graph36";folder.mkdir(exist_ok=True);cases=[];start=time.monotonic()
    for j in range(3,12):
        t0=time.monotonic();graph=s.expand(X[j]-K[j]);assert not graph.has(X[j])
        kept={i:s.expand(f.subs(X[j],graph)) for i,f in enumerate(K) if i!=j}
        bad=[i for i,f in kept.items() if high(f,X)];assert bad
        result=dict(pivot=j+1,graph=str(graph),bad_components=[i+1 for i in bad],tests=[])
        for i in bad:
            allowed=[a for a in kept if a!=i]
            indices=[ix for degree in range(1,4) for ix in combinations_with_replacement(allowed,degree)]
            columns=[high(s.prod(kept[a] for a in ix),X) for ix in indices]
            answer,nrows=solve_or_separate(columns,high(kept[i],X))
            row=dict(component=i+1,target_degree=3,allowed_outputs=[a+1 for a in allowed],
                columns=len(columns),coefficient_rows=nrows,**answer)
            if answer["status"]=="span-member":
                row["solution"]=[dict(indices=[a+1 for a in ix],coefficient=c) for ix,c in zip(indices,answer["solution"]) if s.Rational(c)]
            result["tests"].append(row)
            print(json.dumps(dict(pivot=j+1,component=i+1,status=answer["status"],rank=answer["rank"],augmented_rank=answer["augmented_rank"],functional_terms=len(answer.get("functional",[])))),flush=True)
            if answer["status"]=="separated":break
        result["status"]="raw-component-obstructed" if any(row["status"]=="separated" for row in result["tests"]) else "individual-completions-found"
        result["source_sha256"]=sha(Path(__file__));result["input_sha256"]=sha(source)
        result["seconds"]=round(time.monotonic()-t0,3)
        (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
        cases.append(dict(pivot=j+1,status=result["status"],file="pivot-%02d.json"%(j+1)))
        checkpoint=dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,seconds=round(time.monotonic()-start,3),
            scope="Nine fixed graph slices in the36-term frame; target corrections of degree<=3 in other raw retained outputs. No sequential-stage or global dimension exclusion.")
        (folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
    checkpoint["status"]="complete"
    (folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

