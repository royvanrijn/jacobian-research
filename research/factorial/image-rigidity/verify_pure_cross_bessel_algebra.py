#!/usr/bin/env python3
"""Exact supplemental controls; the uniform arguments are written proofs."""
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
if not __debug__:
    raise SystemExit("Run ordinary Python: optimized execution is rejected.")
def add(a,b):
    c=dict(a)
    for k,v in b.items():
        c[k]=c.get(k,F(0))+v
        if not c[k]: del c[k]
    return c
def scale(a,v): return {k:v*x for k,x in a.items() if v*x}
def mul(a,b):
    c={}
    for u,x in a.items():
        for v,y in b.items():
            w=tuple(i+j for i,j in zip(u,v))
            c[w]=c.get(w,F(0))+x*y
    return {k:v for k,v in c.items() if v}
def deriv(a,i):
    c={}
    for k,v in a.items():
        if k[i]:
            j=list(k); j[i]-=1
            c[tuple(j)]=v*k[i]
    return c
def shift(a,i,n=1):
    return {tuple(x+n if j==i else x for j,x in enumerate(k)):v for k,v in a.items()}
one={(0,0):F(1)}
def cnst(x): return scale(one,F(x))
def cheb(N):
    U=[one,{(1,0):F(2)}]; T=[one,{(1,0):F(1)}]
    for j in range(2,N+1):
        U.append(add(scale(shift(U[-1],0),2),scale(U[-2],-1)))
        T.append(add(scale(shift(T[-1],0),2),scale(T[-2],-1)))
    return U,T
def G(H):
    out={}; j=0
    while H:
        out=add(out,shift(H,1,j)); H=deriv(H,0); j+=1
    return out
def eq(a,b):
    if a!=b: raise AssertionError("polynomial identity failed")
def zero(a): eq(a,{})
def canonical_hash(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(",",":")).encode()).hexdigest()
def ode_lhs(W,n):
    Wb=deriv(W,1); Wbb=deriv(Wb,1)
    A=add(one,{(2,0):F(-1)})
    lhs=add(add(scale(shift(Wbb,1,4),-1),
                scale(shift(Wb,1,3),-3)),scale(shift(shift(Wb,0),1,2),-2))
    lhs=add(lhs,scale(shift(W,1,2),n*(n+2)))
    lhs=add(lhs,scale(shift(shift(W,0),1),-1))
    return add(lhs,mul(A,W))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--output",required=True); args=ap.parse_args()
    output=Path(args.output)
    if output.exists(): raise SystemExit("Refusing to overwrite an existing receipt.")
    U,T=cheb(27); count=0; transformed=0
    A=add(one,{(2,0):F(-1)})
    for n in range(1,26):
        W=G(U[n]); Wb=deriv(W,1); Wbb=deriv(Wb,1)
        lhs=ode_lhs(W,n)
        rhs=add(mul(A,U[n]),scale(shift(T[n+1],1),-(n+1)))
        eq(lhs,rhs); count+=1
        lhs=add(add(mul(add({(1,0):F(1)}, {(0,1):F(-n)}),W),
                    shift(Wb,1,2)),scale(T[n+1],-1))
        eq(lhs,G(U[n-1])); count+=1
        rhs=add(add(scale(shift(W,1),2*(n+1)),G(U[n-1])),scale(T[n+1],2))
        eq(G(U[n+1]),rhs); count+=1
        # R(c,t)=t^n W(c,1/t), now variable 1 denotes t.
        R={(i,n-j):v for (i,j),v in W.items()}
        Rt=deriv(R,1); Rtt=deriv(Rt,1)
        l=shift(Rtt,1)
        l=add(l,scale(shift(shift(Rt,0),1),-2))
        l=add(l,scale(Rt,-(2*n+1)))
        l=add(l,shift(mul(add({(2,0):F(1)},cnst(-1)),R),1))
        l=add(l,scale(shift(R,0),2*n+1))
        r=add(scale(shift(mul(A,U[n]),1,n+1),-1),
              scale(shift(T[n+1],1,n),n+1))
        eq(l,r); transformed+=1
        # V=e^-ct R: direct algebraic gauge after factoring the exponential.
        DcR=add(deriv(R,1),scale(shift(R,0),-1))
        Dc2R=add(deriv(DcR,1),scale(shift(DcR,0),-1))
        vexpr=add(shift(Dc2R,1),scale(DcR,-(2*n+1)))
        vexpr=add(vexpr,scale(shift(R,1),-1))
        eq(vexpr,l); transformed+=1
    # Typed broad-shortcut counterexamples.
    n=2; W=G(U[n])
    eq(W,add(add({(2,0):F(4)},cnst(-1)),add({(1,1):F(8)},{(0,2):F(8)})))
    p=[F(1),F(4),F(31,10),F(11,15),F(121,2400)]
    transform=[p[j]*math.factorial(j) for j in range(5)]
    q=[F(1),F(2),F(11,10)]
    square=[sum((q[i]*q[j-i] for i in range(3) if 0<=j-i<3),F(0)) for j in range(5)]
    if transform!=square: raise AssertionError("typed transform mismatch")
    vals=[]
    for x in (0,-1,-2,-5,-10):
        vals.append(int(2400*sum((a*x**j for j,a in enumerate(p)),F(0))))
    if vals != [2400,-1399,816,-3975,100400]: raise AssertionError("sign witness mismatch")
    # Negative controls deliberately corrupt indispensable coefficients.
    negative=0
    for index,n in enumerate((2,3,7,13,25)):
        W=G(U[n])
        broken=add(W,{(n-index,index):F(1)})
        uncorrupted_rhs=add(mul(A,U[n]),scale(shift(T[n+1],1),-(n+1)))
        try: eq(ode_lhs(broken,n),uncorrupted_rhs)
        except AssertionError: negative+=1
        else: raise AssertionError("corrupted native coefficient passed the actual ODE")
    out={"status":"PASS","scope":"supplemental exact transcription checks, not a uniform proof",
         "universal_polynomial_identity_checks":count,"reciprocal_and_gauge_checks":transformed,
         "typed_counterexample_sign_values":vals,"negative_controls_rejected":negative,
         "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with output.open("x",encoding="utf-8") as f:
        json.dump(out,f,indent=2,sort_keys=True); f.write("\n")
    print(json.dumps(out,sort_keys=True))
if __name__=="__main__": main()
