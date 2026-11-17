#!/usr/bin/env python3
"""Independent Fraction replay of the complete finite all-pole certificate.

No imports from the producer, discovery, legacy algebra or SymPy. A closed
monomial bracket formula and separate Gaussian elimination verify the packet.
"""
from collections import defaultdict,Counter
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json
import time

HERE=Path(__file__).resolve().parent

def add(out,m,c):
    if c:
        out[m]=out.get(m,F(0))+c
        if not out[m]: del out[m]
def read(rows): return {tuple(m):F(c) for m,c in rows}
def dz(p,n):
    for _ in range(n): p={(i,j,k-1):c*k for (i,j,k),c in p.items() if k}
    return p
def delta(p,n):
    for _ in range(n):
        out={}
        for (i,j,k),c in p.items():
            add(out,(i+1,j,k),c*(3*i-6*j))
            if j: add(out,(i,j-1,k),c*2*j)
        p=out
    return p
def multiply(p,q):
    out={}
    for m,c in p.items():
        for n,d in q.items(): add(out,tuple(i+j for i,j in zip(m,n)),c*d)
    return out
def pi(p,q,n):
    out={}
    for j in range(n+1):
        term=multiply(delta(dz(p,n-j),j),dz(delta(q,n-j),j))
        for m,c in term.items(): add(out,m,(-1)**j*comb(n,j)*c)
    return out
def bracket_terms(m,n):
    i,j,k=m; ii,jj,kk=n
    z=k+kk-1
    if z<0: return ()
    return (((i+ii+1,j+jj,z),3*(k*(ii-2*jj)-(i-2*j)*kk)),
            ((i+ii,j+jj-1,z),2*(k*jj-j*kk)))
def bracket(p,q):
    out={}
    for m,c in p.items():
        for n,d in q.items():
            for w,v in bracket_terms(m,n): add(out,w,c*d*v)
    return out
def period(p,lam): return sum(c*lam.get(m,F(0)) for m,c in p.items())
def monomial_period(m,n,lam): return sum(v*lam.get(w,F(0)) for w,v in bracket_terms(m,n))

def main():
    start=time.monotonic()
    packet=json.loads((HERE/"combined_linear_packet.json").read_text())
    receipt=json.loads((HERE/"combined_laurent.json").read_text())
    for name,h in receipt["hashes"].items(): assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,name
    for name,h in receipt["source_hashes"].items(): assert hashlib.sha256((HERE.parents[2]/name).read_bytes()).hexdigest()==h,name
    S,T=read(packet["S"]),read(packet["T"])
    original=json.loads((HERE/"inputs.json").read_text())
    original_receipt=json.loads((HERE/"receipt.json").read_text())
    assert hashlib.sha256((HERE/"inputs.json").read_bytes()).hexdigest()==original_receipt["inputs_sha256"]
    assert S=={tuple(m[:3]):F(c) for m,c in original["S"]}
    assert T=={tuple(m[:3]):F(c) for m,c in original["T"]}
    assert bracket(S,T)=={(0,0,0):F(1)}
    lam=read(packet["functional"])
    sm=list(map(tuple,packet["s_monomials"])); tm=list(map(tuple,packet["t_monomials"]))
    expected=lambda d,z:[(i,j,k) for k in range(z+1) for i in range(-15,d-4*k+1) for j in range(d-4*k-i+1)]
    assert sm==expected(28,3) and tm==expected(23,2)
    split=len(sm); total=split+len(tm)
    columns=[]
    for side,mons in ((0,sm),(1,tm)):
        for m in mons:
            col=bracket({m:F(1)},T) if side==0 else bracket(S,{m:F(1)})
            columns.append({w:c for w,c in col.items() if w[0]>=-3})
    rhs={m:-c/24 for m,c in pi(S,T,3).items() if m[0]>=-3}
    base={i:F(c) for i,c in packet["base"]}
    kernel=[{i:F(c) for i,c in row} for row in packet["kernel"]]
    def apply(vec):
        out={}
        for i,c in vec.items():
            for m,v in columns[i].items(): add(out,m,c*v)
        return out
    assert apply(base)==rhs
    for vec in kernel: assert not apply(vec)
    incidence=Counter(i for vec in kernel for i in vec)
    free=[]
    for vec in kernel:
        candidates=[i for i,c in vec.items() if c==1 and incidence[i]==1]
        assert candidates
        free.append(min(candidates))
    assert len(set(free))==1191
    # Complementary columns must have full column rank. Use a separate
    # sparse forward elimination, keeping the highest monomial as pivot.
    pivots={}
    for i in range(total):
        if i in set(free): continue
        row=dict(columns[i])
        while row:
            p=max(row)
            if p not in pivots:
                c=row[p]
                pivots[p]={m:v/c for m,v in row.items()}
                break
            c=row[p]
            for m,v in pivots[p].items(): add(row,m,-c*v)
        else: raise AssertionError("complementary pivot columns dependent")
    assert len(pivots)==3719 and total==4910
    print("PASS independent affine completeness",total,len(pivots),len(kernel),flush=True)
    # Ambient bilinear form, computed from direct monomial bracket terms.
    B=defaultdict(list)
    for i,m in enumerate(sm):
        for j,n in enumerate(tm):
            if 13<=m[0]+n[0]<=20 and 0<=m[2]+n[2]-1<=4 and m[1]+n[1]<=7:
                c=monomial_period(m,n,lam)
                if c: B[i].append((j,c))
    def bilinear(u,v): return sum(c*d*v.get(j,F(0)) for i,c in u.items() for j,d in B.get(i,()))
    L=[]
    for side,mons in ((0,sm),(1,tm)):
        for m in mons:
            if not(0<=m[0]<=20 and m[1]<=9): L.append(F(0)); continue
            p={m:F(1)}
            L.append(period(pi(p,T,3) if side==0 else pi(S,p,3),lam)/24)
    def splitvec(vec): return ({i:c for i,c in vec.items() if i<split},{i-split:c for i,c in vec.items() if i>=split})
    bp=splitvec(base); kp=list(map(splitvec,kernel))
    value=bilinear(*bp)+sum(L[i]*c for i,c in base.items())+period(pi(S,T,5),lam)/1920
    assert value==1
    for vec,(u,v) in zip(kernel,kp):
        assert bilinear(bp[0],v)+bilinear(u,bp[1])+sum(L[i]*c for i,c in vec.items())==0
    tinc=defaultdict(list)
    for i,(_,v) in enumerate(kp):
        for j,c in v.items(): tinc[j].append((i,c))
    rows=[]
    for u,_ in kp:
        middle={}
        for i,c in u.items():
            for j,v in B.get(i,()): add(middle,j,c*v)
        row={}
        for j,c in middle.items():
            for k,v in tinc[j]: add(row,k,c*v)
        rows.append(row)
    for i,row in enumerate(rows):
        for j,c in row.items(): assert c+rows[j].get(i,F(0))==0
    count=0
    for side,d,zcap in ((0,22,1),(1,17,0)):
        for k in range(zcap+1):
            for i in range(d-4*k+1):
                for j in range(d-4*k-i+1):
                    m=(i,j,k)
                    col=bracket({m:F(1)},T) if side==0 else bracket(S,{m:F(1)})
                    assert period(col,lam)==0
                    count+=1
    assert count==637
    out={"status":"PASS","arithmetic":"Python Fraction; separate monomial bracket and forward elimination", "projected_rank":len(pivots),"kernel_dimension":len(kernel),"current_columns":count,"linear_variations":len(kernel),"quadratic_coefficient_pairs":len(kernel)*(len(kernel)+1)//2,"period":str(value),
         "elapsed_seconds":time.monotonic()-start,"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"packet_sha256":hashlib.sha256((HERE/"combined_linear_packet.json").read_bytes()).hexdigest(),
         "boundary":"Complete independent replay of the finite cutoff certificate. The infinite-pole cutoff and cubic-shear transport are written proofs, not formal verification."}
    (HERE/"combined_independent.json").write_text(json.dumps(out,indent=2)+"\n")
    print("PASS independent full strong cocycle",value,"seconds",out["elapsed_seconds"],flush=True)

if __name__=="__main__":main()
