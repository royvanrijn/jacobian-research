#!/usr/bin/env python3
"""Bounded discovery search for a small explicit cubic stable reduction.

Every move is an exact tame source/target equivalence. Beam pruning concerns
discovery only and supplies no lower bound on dimension or support.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse,hashlib,json,time
HERE=Path(__file__).resolve().parent
def clean(p):return {e:c for e,c in p.items() if c}
def add(a,b,scale=Q(1)):
    r=dict(a)
    for e,c in b.items():r[e]=r.get(e,Q(0))+scale*c
    return clean(r)
def mon(n,positions):
    e=[0]*n
    for j,d in positions.items():e[j]=d
    return tuple(e)
def initial():
    # Tangent-normalized foundational map: (F3/2,F2,F1).
    n=3
    rows=[
      {(1,0,0):Q(1),(2,1,0):Q(-3,2),(3,0,1):Q(-1,2)},
      {(0,1,0):Q(1),(1,0,1):Q(3),(2,1,1):Q(6),(3,2,1):Q(3),(1,2,0):Q(12),(2,3,0):Q(9)},
      {(0,0,1):Q(1),(1,1,1):Q(3),(2,2,1):Q(3),(3,3,1):Q(1),(0,2,0):Q(4),(1,3,0):Q(7),(2,4,0):Q(3)}]
    return dict(rows=rows,moves=[],points=[[Q(0),Q(-2),Q(-16)],[Q(1),Q(-2),Q(8)]])
def evaluate_monomial(e,p):return Q(1)*prod(q**d for q,d in zip(p,e))
def prod(values):
    r=Q(1)
    for v in values:r*=v
    return r
def reduce(state,a,b):
    rows=state["rows"];n=len(rows);square=a==b;extra=1 if square else 2
    pad=lambda e:e+(0,)*extra
    a1=pad(a);b1=pad(b);m=tuple(x+y for x,y in zip(a1,b1))
    u=mon(n+extra,{n:1});v=u if square else mon(n+extra,{n+1:1})
    out=[{pad(e):c for e,c in p.items()} for p in rows]
    def plus(e,f):return tuple(x+y for x,y in zip(e,f))
    replacement={plus(u,a1):Q(2),plus(u,u):Q(-1)} if square else {plus(u,b1):Q(1),plus(v,a1):Q(1),plus(u,v):Q(-1)}
    for i,p in enumerate(out):
        c=p.pop(m,Q(0))
        if c:out[i]=add(p,replacement,c)
    carriers=[{u:Q(1),a1:Q(-1)}] if square else [{u:Q(1),a1:Q(-1)},{v:Q(1),b1:Q(-1)}]
    # Elementary target row operations remove standalone carrier monomials.
    clears=[]
    for j,factor in enumerate([a1] if square else [a1,b1]):
        for i in range(n):
            c=out[i].get(factor,Q(0))
            if c:out[i]=add(out[i],carriers[j],c);clears.append([i,j,str(c)])
    out+=carriers
    points=[p+[evaluate_monomial(a,p)] if square else p+[evaluate_monomial(a,p),evaluate_monomial(b,p)] for p in state["points"]]
    move=dict(a=list(a),b=list(b),clears=clears,square=square)
    return dict(rows=out,moves=state["moves"]+[move],points=points)
def signature(state):
    return tuple(tuple(sorted(p.items())) for p in state["rows"])
def metric(state):
    rows=state["rows"];degrees=[sum(e) for p in rows for e in p]
    high={e for p in rows for e in p if sum(e)>3}
    return (max(degrees),sum((sum(e)-3)**2 for e in high),len(high),
            sum(sum(e)>3 for p in rows for e in p),sum(sum(e)>=2 for p in rows for e in p))
def serialize(state):
    return dict(dimension=len(state["rows"]),degree=max(sum(e) for p in state["rows"] for e in p),
      nonlinear_counts=[sum(sum(e)>=2 for e in p) for p in state["rows"]],
      coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in sorted(p.items())] for p in state["rows"]],
      points=[[str(q) for q in p] for p in state["points"]],moves=state["moves"])
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--beam",type=int,default=40);ap.add_argument("--max-dimension",type=int,default=19);args=ap.parse_args()
    t=time.monotonic();queues={3:[initial()]};log=[];processed=0
    for n in range(3,args.max_dimension+1):
        states=queues.pop(n,[])
        unique={}
        for state in states:
            sig=signature(state)
            if sig not in unique:unique[sig]=state
        states=sorted(unique.values(),key=metric)[:args.beam]
        if not states:continue
        completed=[s for s in states if metric(s)[0]<=3]
        row=dict(dimension=n,unique_candidates=len(unique),retained=len(states),best_metric=metric(states[0]),elapsed_seconds=round(time.monotonic()-t,3))
        log.append(row);print(json.dumps(row),flush=True)
        receipt=dict(status="running",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parameters=vars(args),log=log,processed=processed,boundary="Bounded beam discovery; not an optimum certificate.")
        if completed:
            best=min(completed,key=lambda s:sum(sum(e)>=2 for p in s["rows"] for e in p))
            receipt.update(status="complete",witness=serialize(best))
            (HERE/"cubic-witness.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
            print("CUBIC",n,"nonlinear",sum(receipt["witness"]["nonlinear_counts"]),flush=True);return
        for state in states:
            processed+=1
            top=metric(state)[0]
            exponents=sorted({e for p in state["rows"] for e in p if sum(e)==top})
            for m in exponents:
                for a in product(*(range(q+1) for q in m)):
                    b=tuple(q-r for q,r in zip(m,a))
                    if min(sum(a),sum(b))<2 or a>b:continue
                    new=reduce(state,a,b);nn=len(new["rows"])
                    if nn<=args.max_dimension:queues.setdefault(nn,[]).append(new)
        (HERE/"search-checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    receipt["status"]="bounded-search-exhausted"
    (HERE/"search-checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

