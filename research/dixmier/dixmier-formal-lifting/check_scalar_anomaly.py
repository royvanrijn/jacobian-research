#!/usr/bin/env python3
"""Exact controls for QCL1. Universal degree and all-rank claims use written proofs."""
from pathlib import Path
from itertools import product
from math import factorial
import hashlib
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
x,p,y,q,h=sp.symbols("x p y q h")
coords=(x,p,y,q)
J=sp.diag(sp.Matrix([[0,-1],[1,0]]),sp.Matrix([[0,-1],[1,0]]))

def contraction(f,g,k):
    total=0
    # Counts of p_left/x_right, x_left/p_right in each canonical pair.
    for counts in product(range(k+1),repeat=4):
        if sum(counts)!=k:
            continue
        ap,ax,bq,by=counts
        left=sp.diff(f,p,ap,x,ax,q,bq,y,by)
        right=sp.diff(g,x,ap,p,ax,y,bq,q,by)
        total+=(-1)**(ax+by)*left*right/sp.prod(factorial(a) for a in counts)
    return sp.expand(factorial(k)*total)

def poisson(f,g):
    return contraction(f,g,1)

def normalized_moyal(f,g):
    return sp.expand(poisson(f,g)+h**2*contraction(f,g,3)/24)

def zero(f):
    assert sp.cancel(f)==0,sp.cancel(f)

def gram_schmidt(B):
    n=B.rows
    residual=[sp.eye(n).row(i) for i in range(n)]
    chosen=[]
    def beta(a,b):
        return sp.cancel((a*B*b.T)[0])
    for i in range(0,n,2):
        u=residual[i]
        v=residual[i+1].applyfunc(lambda z:sp.cancel(-z/beta(u,residual[i+1])))
        chosen.extend([u,v])
        for j in range(i+2,n):
            w=residual[j]
            residual[j]=(w+beta(w,v)*u-beta(w,u)*v).applyfunc(sp.cancel)
    return sp.Matrix.vstack(*chosen)

def main():
    f=x**3+p*y**2+x*p*q
    g=p**3+x**2*q+p*y*q
    cubic_correction=sp.expand(normalized_moyal(f,g)-poisson(f,g))
    assert sp.Poly(cubic_correction,*coords).total_degree()==0
    assert sp.expand(cubic_correction/h**2).has(h)==False
    assert cubic_correction!=0

    u1,u2=x+q,y+p
    canonical=[u1,p+u1**3+u1*u2**2,u2,q+u1**2*u2+u2**3]
    relation_count=0
    for i in range(4):
        for j in range(i+1,4):
            zero(normalized_moyal(canonical[i],canonical[j])-J[i,j])
            relation_count+=1

    C=sp.Matrix([[0,2,3,-1],[-2,0,1,4],[-3,-1,0,5],[1,-4,-5,0]])
    B=J+h**2*C
    A=gram_schmidt(B)
    for entry in A*B*A.T-J:
        zero(entry)
    assert A.subs(h,0)==sp.eye(4)
    assert all(sp.denom(sp.cancel(entry)).subs(h,0)!=0 for entry in A)

    H=sp.Matrix([0,x**3+x*y**2,0,x**2*y+y**3])
    variables=sp.Matrix(coords)
    N=H.jacobian(coords)
    assert N*N==sp.zeros(4)
    assert (N*variables-3*H).applyfunc(sp.expand)==sp.zeros(4,1)
    assert (N*H).applyfunc(sp.expand)==sp.zeros(4,1)
    inverse_sub={z:z-H[i] for i,z in enumerate(coords)}
    for i,z in enumerate(coords):
        zero((z+H[i]).subs(inverse_sub,simultaneous=True)-z)

    z1,z2,z3,r1,r2,r3=sp.symbols("z1 z2 z3 r1 r2 r3")
    base=sp.Matrix([z1+z2**3,z2+z3**3,z3])
    basevars=(z1,z2,z3)
    Nbase=base.jacobian(basevars)-sp.eye(3)
    assert Nbase**3==sp.zeros(3) and Nbase**2!=sp.zeros(3)
    inverse_jacobian=sp.eye(3)-Nbase+Nbase**2
    assert base.jacobian(basevars)*inverse_jacobian==sp.eye(3)
    momentum=inverse_jacobian.T*sp.Matrix([r1,r2,r3])
    degree=max(sp.Poly(a,z1,z2,z3,r1,r2,r3).total_degree() for a in momentum)
    assert degree==5

    receipt={
        "claim":"QCL1","status":"PASS",
        "nonzero_scalar_cubic_Moyal_correction":str(cubic_correction),
        "mixed_cubic_full_rank_two_relations":relation_count,
        "exact_rational_skew_matrix_normalization":True,
        "normalization_regular_and_identity_at_zero":True,
        "homogeneous_square_zero_and_inverse_controls":True,
        "triangular_index_three_cotangent_degree":int(degree),
        "coverage_boundary":"Exact symbolic identities only. The all-rank scalar-anomaly criterion, cubic bounded lifting, all-degree homogeneous symplectic invertibility and the nilpotency-index degree formula use written proofs. No noninvertible cubic symplectic map or bounded Long lift is constructed.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("SCALAR_ANOMALY_AND_CUBIC_LIFTS.md","check_scalar_anomaly.py")}
    }
    (HERE/"scalar_anomaly_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
