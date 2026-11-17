#!/usr/bin/env python3
"""Bounded cubic discovery with reuse of existing output coordinates as factors."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse,hashlib,json,time
from verify import add,mul,scale,variable,evaluate
from search import initial
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def lt(p):return max(p,key=lambda e:(sum(e),e))
def padded(p,d):return {e+(0,)*d:c for e,c in p.items()}
def metric(s):
    rows=s["rows"];high={e for p in rows for e in p if sum(e)>3}
    return (max(sum(e) for p in rows for e in p),
      sum((sum(e)-3)**2 for e in high),len(high),sum(len(p) for p in rows))
def key(s):return tuple(tuple(sorted(p.items())) for p in s["rows"])
def linear_reduce(s):
    rows=[dict(p) for p in s["rows"]];ops=[]
    while True:
        groups={}
        for i,p in enumerate(rows):groups.setdefault(lt(p),[]).append(i)
        duplicate=next(((e,g) for e,g in sorted(groups.items(),key=lambda t:(sum(t[0]),t[0]),reverse=True) if len(g)>1),None)
        if duplicate is None:break
        e,g=duplicate;base=min(g,key=lambda j:(len(rows[j]),j))
        for i in g:
            if i==base:continue
            c=rows[i][e]/rows[base][e]
            rows[i]=add(rows[i],rows[base],-c);ops.append(dict(kind="linear",target=i,factor=base,coefficient=str(c)))
    return dict(rows=rows,points=s["points"],moves=s["moves"]+ops)
def product_reduce(s,a,b,ia,ib):
    # ia/ib are output rows, or None to introduce w-a or w-b.
    rows=[dict(p) for p in s["rows"]];points=[list(p) for p in s["points"]];moves=list(s["moves"])
    aa,bb=a,b
    def extend(e):
        nonlocal rows,points,moves
        n=len(rows);rows=[padded(p,1) for p in rows]
        e=e+(0,)*(n-len(e));ep=e+(0,)
        rows.append(add(variable(n,n+1),{ep:Q(1)},Q(-1)))
        points=[p+[evaluate({e:Q(1)},p)] for p in points]
        moves.append(dict(kind="carrier",exponent=list(e)))
        return n
    if ia is None:ia=extend(a)
    if ib is None:
        if a==b:ib=ia
        else:ib=extend(b)
    n=len(rows);a=a+(0,)*(n-len(a));b=b+(0,)*(n-len(b))
    m=tuple(x+y for x,y in zip(a,b))
    assert lt(rows[ia])==a and lt(rows[ib])==b
    prod=mul(rows[ia],rows[ib]);den=rows[ia][a]*rows[ib][b];targets=[]
    for i,p in enumerate(rows):
        if i in (ia,ib):continue
        c=p.get(m,Q(0))/den
        if c:rows[i]=add(p,prod,-c);targets.append([i,str(c)])
    assert targets
    moves.append(dict(kind="product",factors=[ia,ib],targets=targets))
    # Clear the two existing leading monomials by elementary target row operations.
    for j,e in ((ia,a),(ib,b)):
        for i,p in enumerate(rows):
            if i in (ia,ib):continue
            c=p.get(e,Q(0))/rows[j][e]
            if c:rows[i]=add(p,rows[j],-c);moves.append(dict(kind="linear",target=i,factor=j,coefficient=str(c)))
    return linear_reduce(dict(rows=rows,points=points,moves=moves))
def zero_closure(s):
    s=linear_reduce(s)
    for step in range(200):
        if metric(s)[0]<=3:return s
        m=max((e for p in s["rows"] for e in p),key=lambda e:(sum(e),e))
        leads=[lt(p) for p in s["rows"]];pairs=[]
        for i,a in enumerate(leads):
            if sum(a)<2:continue
            for j in range(i,len(leads)):
                b=leads[j]
                if sum(b)>=2 and tuple(x+y for x,y in zip(a,b))==m:
                    if any(k not in (i,j) and m in p for k,p in enumerate(s["rows"])):
                        pairs.append((len(s["rows"][i])*len(s["rows"][j]),i,j,a,b))
        if not pairs:return s
        _,i,j,a,b=min(pairs)
        s=product_reduce(s,a,b,i,j)
        if sum(len(p) for p in s["rows"])>3000:return None
    return None
def serialize(s):
    return dict(dimension=len(s["rows"]),degree=metric(s)[0],
       coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in sorted(p.items())] for p in s["rows"]],
       points=[[str(q) for q in p] for p in s["points"]],moves=s["moves"])
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--beam",type=int,default=60);ap.add_argument("--max-dimension",type=int,default=17);args=ap.parse_args()
    t=time.monotonic();q={3:[zero_closure(initial())]};logs=[];processed=0
    for n in range(3,args.max_dimension+1):
        states=q.pop(n,[]);unique={}
        for s in states:
            if s is not None:unique.setdefault(key(s),s)
        states=sorted(unique.values(),key=metric)[:args.beam]
        if not states:continue
        row=dict(dimension=n,candidates=len(unique),retained=len(states),best_metric=metric(states[0]),seconds=round(time.monotonic()-t,3));logs.append(row);print(json.dumps(row),flush=True)
        receipt=dict(status="running",source_sha256=sha(Path(__file__)),
            helper_sha256={"verify.py":sha(HERE/"verify.py"),"search.py":sha(HERE/"search.py")},
            parameters=vars(args),log=logs,processed=processed,
            boundary="Bounded beam with greedy zero-cost closure; no optimality claim.")
        complete=[s for s in states if metric(s)[0]<=3]
        if complete:
            receipt.update(status="complete",witness=serialize(min(complete,key=lambda s:sum(len(p) for p in s["rows"]))))
            (HERE/"reuse-witness.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
            print("CUBIC",n,flush=True);return
        for s in states:
            processed+=1;rows=s["rows"];m=max((e for p in rows for e in p),key=lambda e:(sum(e),e))
            leads={}
            for i,p in enumerate(rows):leads.setdefault(lt(p),[]).append(i)
            for a in product(*(range(d+1) for d in m)):
                b=tuple(d-e for d,e in zip(m,a))
                if min(sum(a),sum(b))<2 or a>b:continue
                for ia in leads.get(a,[])+[None]:
                    for ib in leads.get(b,[])+[None]:
                        extra=int(ia is None)+int(ib is None and (a!=b or ia is not None))
                        if not extra or n+extra>args.max_dimension:continue
                        new=zero_closure(product_reduce(s,a,b,ia,ib))
                        if new is not None:q.setdefault(len(new["rows"]),[]).append(new)
        (HERE/"reuse-checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    receipt["status"]="bounded-search-exhausted"
    (HERE/"reuse-checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

