#!/usr/bin/env python3
"""Exact rational tests of the uniform compression proof, not its premise."""
import argparse, hashlib, json
from math import comb, factorial
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
def bernstein_coefficients(p,a,b):
    degree=len(p)-1
    shifted=[sum((p[k]*comb(k,j)*a**(k-j)*(b-a)**j
                  for k in range(j,degree+1)),F(0)) for j in range(degree+1)]
    return [sum((shifted[j]*F(comb(k,j),comb(degree,j))
                 for j in range(k+1)),F(0)) for k in range(degree+1)]
def strict_cauchy_energy_counterexample():
    # Q has four negative simple roots. P=xQ/12000, P'=h.
    Q=list(map(F,(12000,24000,12400,2200,121)))
    Qp=deriv(Q); Qpp=deriv(Qp)
    R=[3*x for x in Qp]
    for j,x in enumerate(Qpp): R[j+1]+=x
    intervals=[(F(-19,2),F(-47,5)),(F(-11,2),F(-27,5)),
               (F(-13,5),F(-5,2)),(F(-4,5),F(-3,4))]
    rows=[]
    for i,(a,b) in enumerate(intervals):
        expected=-1 if i%2==0 else 1
        assert a<b<0 and value(Q,a)*value(Q,b)<0
        qp=bernstein_coefficients(Qp,a,b); r=bernstein_coefficients(R,a,b)
        assert all(expected*x>0 for x in qp+r)
        rows.append({"interval":[str(a),str(b)],
                     "Q_endpoint_values":[str(value(Q,a)),str(value(Q,b))],
                     "Q_prime_Bernstein":[str(x) for x in qp],
                     "R_Bernstein":[str(x) for x in r],"common_sign":expected})
    h=[F(1),F(4),F(31,10),F(11,15),F(121,2400)]
    A=[[F(-4),F(2),F(0),F(0)],
       [F(-31,10),F(0),F(3),F(0)],
       [F(-11,15),F(0),F(0),F(4)],
       [F(-121,2400),F(0),F(0),F(0)]]
    assert det_I_coeff(A)==[h[j]*factorial(j) for j in range(5)]
    B=matmul(A,A)
    for i in range(4):
        for j in range(4): B[i][j]+=2*A[i][j]+(F(11,10) if i==j else 0)
    assert any(x for row in B for x in row)
    assert matmul(B,B)==[[F(0)]*4 for _ in range(4)]
    return {"polynomial_Q_ascending":[str(x) for x in Q],
            "compressed_matrix":[[str(x) for x in row] for row in A],
            "det_I_minus_bA":"(1+2*b+11*b^2/10)^2",
            "minimal_polynomial":"(lambda^2+2*lambda+11/10)^2",
            "q_of_A_nonzero_and_square_zero":True,
            "positive_energy_ratio":"R(x)/Q_prime(x), R=3Q_prime+xQ_double_prime",
            "exact_root_and_positive_diagonal_certificates":rows}
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
    counter=strict_cauchy_energy_counterexample()
    receipt={"status":"PASS","scope":"exact rational controls for compression and energy identities, plus a strict-energy Cauchy counterexample to simplicity",
             "node_grid_and_selected_node_cases":cases,"energy_entry_checks":entries,
             "diagonal_corruptions_rejected":neg,
             "strict_cauchy_energy_counterexample":counter,
             "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with out.open("x") as f: json.dump(receipt,f,sort_keys=True,indent=2); f.write("\n")
    print(json.dumps(receipt,sort_keys=True))
if __name__=="__main__": main()
