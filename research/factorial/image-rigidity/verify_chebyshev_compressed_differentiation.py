#!/usr/bin/env python3
"""Exact rational tests of the uniform compression proof, not its premise."""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
if not __debug__: raise SystemExit("Optimized Python is not supported.")
def matmul(A,B):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def coeff_product(nodes):
    p=[F(1)]
    for x in nodes:
        q=[F(0)]*(len(p)+1)
        for k,a in enumerate(p): q[k]-=x*a; q[k+1]+=a
        p=q
    return p
def deriv(p): return [k*p[k] for k in range(1,len(p))]
def value(p,x):
    v=F(0)
    for a in reversed(p): v=v*x+a
    return v
def det_I_coeff(K):
    d=len(K); powers=K; tr=[]
    for k in range(1,d+1):
        tr.append(sum((powers[i][i] for i in range(d)),F(0)))
        if k<d: powers=matmul(powers,K)
    a=[F(1)]
    for k in range(1,d+1):
        a.append(-sum((tr[i-1]*a[k-i] for i in range(1,k+1)),F(0))/k)
    return a
def expected(nodes,h):
    p=deriv(coeff_product(nodes)); c=nodes[h]; den=value(p,c); out=[]
    while p:
        out.append(value(p,c)/den); p=deriv(p)
    return out
def check(nodes,h,corrupt=False):
    c=nodes[h]; inds=[i for i in range(len(nodes)) if i!=h]
    S=[sum((1/(x-y) for j,y in enumerate(nodes) if i!=j),F(0)) for i,x in enumerate(nodes)]
    K=[[S[i]+1/(nodes[i]-c) if i==j else 1/(nodes[i]-nodes[j]) for j in inds] for i in inds]
    if corrupt: K[0][0]+=1
    assert det_I_coeff(K)==expected(nodes,h)
    # General-node version of J K+Kt J: diagonal 2(x-c)S+2, offdiagonal 1.
    for i,ni in enumerate(inds):
        for j,nj in enumerate(inds):
            x,y=nodes[ni],nodes[nj]
            got=(x-c)*K[i][j]+(y-c)*K[j][i]
            want=2*(x-c)*S[ni]+2 if i==j else F(1)
            assert got==want
    return len(inds)**2
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--output",required=True); a=ap.parse_args()
    out=Path(a.output)
    if out.exists(): raise SystemExit("Refusing overwrite.")
    cases=entries=0
    for N in range(2,12):
        nodes=[F(i*i+3*i-2, i+1) for i in range(N)]
        assert len(set(nodes))==N
        for h in range(N):
            entries+=check(nodes,h); cases+=1
    # A second irregular grid tests signs and nonsymmetric node layouts.
    for N in (3,5,8):
        nodes=[F((-1)**i*(i+2),i+1) for i in range(N)]
        for h in range(N): entries+=check(nodes,h); cases+=1
    neg=0
    for N in (2,3,5,8):
        try: check([F(i) for i in range(N)],N//2,True)
        except AssertionError: neg+=1
        else: raise AssertionError("Corruption accepted.")
    receipt={"status":"PASS","scope":"exact rational controls for arbitrary-node compression and energy identities",
             "node_grid_and_selected_node_cases":cases,"energy_entry_checks":entries,
             "diagonal_corruptions_rejected":neg,
             "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with out.open("x") as f: json.dump(receipt,f,sort_keys=True,indent=2); f.write("\n")
    print(json.dumps(receipt,sort_keys=True))
if __name__=="__main__": main()
