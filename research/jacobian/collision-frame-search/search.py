#!/usr/bin/env python3
"""Exact scalar searches in the affine collision-frame stabilizer.

Each line H -> S(t)^(-1) H(S(t)X) is analyzed at all algebraic coefficient
zero strata. The optional bounded beam is discovery only, not orbit coverage.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import comb
from pathlib import Path
import argparse,hashlib,json,time
import sympy as sp

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/"sparse-nine"
DIRECTIONS=((0,1),(0,2),(1,2),(2,1))
T=sp.Symbol("t")

def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+".research/.scratch")
    temp.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
    temp.replace(path)
def load_base():
    data=json.loads((PREV/"normalized-49-witness.json").read_text())
    return tuple({tuple(v["exponent"]):Q(v["coefficient"]) for v in row}
                 for row in data["coefficients"])
def dump_map(outputs):
    return [[dict(exponent=list(e),coefficient=str(c)) for e,c in sorted(row.items())]
            for row in outputs]
def load_map(data):
    return tuple({tuple(v["exponent"]):Q(v["coefficient"]) for v in row} for row in data)
def identifier(outputs):
    return hashlib.sha256(json.dumps(dump_map(outputs),sort_keys=True).encode()).hexdigest()
def support(outputs):
    return sum(sum(sum(e)>=2 for e in row) for row in outputs)
def height(outputs):
    values=[abs(c.numerator).bit_length()+c.denominator.bit_length()
            for row in outputs for c in row.values()]
    return max(values),sum(values)
def canonical(outputs):
    # Quotient only by explicit rational diagonal conjugation fixing e1.
    entries=[]
    for i,row in enumerate(outputs):
        for e,c in sorted(row.items(),key=lambda v:(sum(v[0]),v[0])):
            if sum(e)>=2:
                entries.append(((e[1]-int(i==1),e[2]-int(i==2)),c))
    a=b=Q(1)
    for (p,c),(q,d) in combinations(entries,2):
        determinant=p[0]*q[1]-p[1]*q[0]
        if abs(determinant)==1:
            a=(1/c)**(q[1]//determinant)*(1/d)**(-p[1]//determinant)
            b=(1/c)**(-q[0]//determinant)*(1/d)**(p[0]//determinant)
            break
    normalized=tuple({e:c*a**(e[1]-int(i==1))*b**(e[2]-int(i==2))
                      for e,c in row.items()} for i,row in enumerate(outputs))
    assert support(normalized)==support(outputs)
    return normalized,(a,b)

def line_coefficients(outputs,i,j):
    source=[]
    for row in outputs:
        expanded=defaultdict(lambda:defaultdict(Q))
        for e,c in row.items():
            for r in range(e[i]+1):
                f=list(e);f[i]-=r;f[j]+=r
                expanded[tuple(f)][r]+=c*comb(e[i],r)
        source.append({e:dict(p) for e,p in expanded.items()})
    for e,p in source[j].items():
        target=source[i].setdefault(e,{})
        for degree,c in p.items():target[degree+1]=target.get(degree+1,Q(0))-c
    cleaned=[]
    for row in source:
        out={}
        for e,p in row.items():
            p={k:c for k,c in p.items() if c}
            if p:out[e]=p
        cleaned.append(out)
    for row_index,row in enumerate(cleaned):
        linear={e:p for e,p in row.items() if sum(e)<=1}
        unit=tuple(int(k==row_index) for k in range(3))
        assert linear=={unit:{0:Q(1)}}
    return cleaned

def evaluate(rows,t):
    return tuple({e:c for e,p in row.items()
                  if (c:=sum(v*t**d for d,v in p.items()))}
                 for row in rows)

def analyze(outputs,i,j,certificate=False):
    rows=line_coefficients(outputs,i,j)
    coefficients=[]
    metadata=[]
    for row_index,row in enumerate(rows):
        for e,p in sorted(row.items()):
            if sum(e)<2:continue
            polynomial=sp.Poly.from_dict({(d,):sp.Rational(c.numerator,c.denominator)
                                          for d,c in p.items()},T,domain=sp.QQ)
            coefficients.append(polynomial)
            metadata.append(dict(row=row_index,exponent=list(e),coefficient=str(polynomial.as_expr())))
    factors={}
    for p in coefficients:
        for f,multiplicity in sp.factor_list(p)[1]:
            f=f.monic();factors[str(f.as_expr())]=f
    strata=[]
    for key,f in factors.items():
        zero_indices=[n for n,p in enumerate(coefficients) if p.rem(f).is_zero]
        entry=dict(factor=key,degree=f.degree(),support=len(coefficients)-len(zero_indices))
        if certificate:entry["zero_indices"]=zero_indices
        if f.degree()==1:
            root=Q(str(-f.nth(0)/f.nth(1)))
            entry["rational_parameter"]=str(root)
        strata.append(entry)
    strata.sort(key=lambda r:(r["support"],r["factor"]))
    minimum=min([len(coefficients)]+[r["support"] for r in strata])
    candidates=[]
    for row in strata:
        if "rational_parameter" in row:
            parameter=Q(row["rational_parameter"])
            candidate=evaluate(rows,parameter)
            assert support(candidate)==row["support"]
            candidates.append((candidate,parameter))
    result=dict(direction=[i,j],base_support=support(outputs),generic_support=len(coefficients),
                minimum_over_algebraic_parameters=minimum,strata=strata)
    if certificate:result["coefficients"]=metadata
    return result,candidates

def matrix_identity():
    return [[Q(int(i==j)) for j in range(3)] for i in range(3)]
def mmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def step_matrix(i,j,t,diagonal):
    a,b=diagonal
    S=matrix_identity();S[i][j]=t
    D=[[Q(1),Q(0),Q(0)],[Q(0),a,Q(0)],[Q(0),Q(0),b]]
    return mmul(S,D)
def dump_matrix(M):return [[str(c) for c in row] for row in M]
def load_matrix(M):return [[Q(c) for c in row] for row in M]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth",type=int,default=1)
    parser.add_argument("--beam",type=int,default=48)
    parser.add_argument("--slack",type=int,default=40)
    a=parser.parse_args()
    base=load_base()
    state_path=HERE/"beam-checkpoint.json"
    initial_id=identifier(base)
    if a.depth==1:
        certificates=[]
        for i,j in DIRECTIONS:
            result,candidates=analyze(base,i,j,True)
            certificates.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ("strata","coefficients")}),flush=True)
        save(HERE/"elementary-lines.json",dict(status="complete",base_sha256=initial_id,
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    lines=certificates))
        return
    if state_path.exists():
        state=json.loads(state_path.read_text())
        assert state["source_sha256"]==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        assert state["base_sha256"]==initial_id and state["parameters"]==dict(depth=a.depth,beam=a.beam,slack=a.slack)
        queue=state["queue"];next_layer=state["next_layer"];seen=set(state["seen"])
        best=state["best"];layer=state["layer"];processed=state["processed"]
    else:
        initial=dict(id=initial_id,map=dump_map(base),matrix=dump_matrix(matrix_identity()),word=[],support=49)
        queue=[initial];next_layer=[];seen={initial_id};best=initial;layer=0;processed=0
        state=dict(schema="collision-frame-beam.v1",base_sha256=initial_id,
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   parameters=dict(depth=a.depth,beam=a.beam,slack=a.slack))
    def checkpoint(status):
        state.update(status=status,queue=queue,next_layer=next_layer,seen=sorted(seen),
                     best=best,layer=layer,processed=processed)
        save(state_path,state)
    while layer<a.depth:
        while queue:
            parent=queue.pop(0)
            outputs=load_map(parent["map"])
            for i,j in DIRECTIONS:
                result,candidates=analyze(outputs,i,j)
                for candidate,t in candidates:
                    if support(candidate)>best["support"]+a.slack:continue
                    candidate,diagonal=canonical(candidate)
                    key=identifier(candidate)
                    if key in seen:continue
                    seen.add(key)
                    M=mmul(load_matrix(parent["matrix"]),step_matrix(i,j,t,diagonal))
                    child=dict(id=key,map=dump_map(candidate),matrix=dump_matrix(M),
                               word=parent["word"]+[dict(direction=[i,j],parameter=str(t),
                                                        diagonal=list(map(str,diagonal)))],
                               support=support(candidate))
                    next_layer.append(child)
                    if child["support"]<best["support"]:
                        best=child
                        save(HERE/"best-witness.json",best)
                        print(json.dumps(dict(event="improvement",support=best["support"],
                                              word=best["word"])),flush=True)
            processed+=1
            checkpoint("running")
            print(json.dumps(dict(layer=layer,processed=processed,pending=len(queue),
                                  new=len(next_layer),best=best["support"])),flush=True)
        layer+=1
        next_layer.sort(key=lambda r:(r["support"],height(load_map(r["map"])),r["id"]))
        queue=next_layer[:a.beam];next_layer=[]
        checkpoint("complete" if layer==a.depth else "running")
    print(json.dumps(dict(status="complete",best_support=best["support"],
                          processed=processed,distinct=len(seen),parameters=state["parameters"])),flush=True)
if __name__=="__main__":main()

