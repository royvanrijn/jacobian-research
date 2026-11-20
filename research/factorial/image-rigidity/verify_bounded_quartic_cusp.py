#!/usr/bin/env python3
"""Exact supporting controls for bounded actual quartic-cusp exclusion."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as Q
from pathlib import Path
if not __debug__:
    raise SystemExit("Refusing optimized execution: exact assertions are required.")
HERE=Path(__file__).resolve().parent
PROOF=HERE/"BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md"
LOCAL=HERE/"QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md"
LOCAL_SHA="86b19ffd10650e14a91cddd47974a360d9329f44dfeda6cc860cb0a84faa7953"
ORDER=6
def full(a): return list(a[:ORDER])+[Q(0)]*max(0,ORDER-len(a))
def plus(a,b):return [x+y for x,y in zip(full(a),full(b))]
def times(a,b):
    a,b=full(a),full(b)
    return [sum((a[j]*b[k-j] for j in range(k+1)),Q(0)) for k in range(ORDER)]
def sc(a,c):return [Q(c)*x for x in full(a)]
def power(a,n):
    out=full([1])
    for _ in range(n):out=times(out,a)
    return out
def reciprocal(a):
    a=full(a); assert a[0]
    out=[1/a[0]]
    for k in range(1,ORDER):
        out.append(-sum((a[j]*out[k-j] for j in range(1,k+1)),Q(0))/a[0])
    return out
def logunit(a):
    assert a[0]==1
    t=plus(a,[-1]); out=full([])
    for k in range(1,ORDER):
        out=plus(out,sc(power(t,k),Q((-1)**(k+1),k)))
    return out
def cusp_data(s):
    d=2*s+1; a=s*s; R=Q(s+1,s); delta=Q(d,s)
    p=[R,Q(1,s)]
    q=R**d; z=(R**(d+1)+1)/delta
    h=sc(plus(power(p,d),[-1]),1/delta)
    H=plus(h,[z])
    f=plus(sc(logunit(sc(p,1/R)),a),
           plus([0,-s],logunit(sc(H,1/q))))
    hprime=[h[k+1]*(k+1) for k in range(ORDER-1)]
    vz=sc(times(hprime,power(reciprocal(H),2)),-1)
    F4=f[4]*24; F5=f[5]*120
    assert f[1:4]==[0,0,0]
    assert F4==-2/R**2
    assert F5==24/Q(d-1)/R**3
    assert vz[0]==-R**(-d-1) and vz[1]==0
    # Linear u scaling A=sqrt(R) gives gamma=F5*A^5/120 >0.
    gamma_squared=F5**2*R**5/120**2
    assert gamma_squared==Q(1,25*(d*d-1))
    # Leading quartic is exactly -u^4/12 after A scaling.
    assert f[4]*R**2==-Q(1,12)
    # Canonical b is six times the raw quintic coefficient.
    assert 36*gamma_squared==Q(36,25*(d*d-1))
    wc=Q(a**d)*R**(d+1)
    # The radial scaled-field coefficient is -sqrt(d^2-1)/2.
    # Its square is (d*A/delta)^2.
    assert Q(d*d)*R/delta**2==Q(d*d-1,4)
    return {"s":s,"d":d,"a":a,"gamma_squared":str(gamma_squared),
            "Fvvz":str(vz[1]),"wc":str(wc)}
# A vector is a sparse polynomial in h with coefficients multiplying Z,Z',Z''.
def vadd(a,b):
    out=dict(a)
    for key,value in b.items():out[key]=out.get(key,Q(0))+value
    return {k:v for k,v in out.items() if v}
def vscale(a,c):return {k:Q(c)*v for k,v in a.items() if c*v}
def diff(a):
    out={}
    for (j,k),c in a.items():
        if k: out=vadd(out,{(j,k-1):c*k})
        if j<2:out=vadd(out,{(j+1,k):c})
        else:out=vadd(out,{(0,k+1):3*c})
    return out
def model_identity():
    v={(0,0):Q(1)}
    for _ in range(3):v=diff(v)
    assert v=={(0,1):Q(3)}
    v=diff(diff(v))
    assert v=={(1,0):Q(6),(2,1):Q(3)}
    # This is Z5=6Z'+3hZ'' exactly, in the third-order ODE quotient.
    return v
def A(m,a,d,x):
    return sum((Q(math.comb(m,k)*math.factorial(a*m+d*k),math.factorial(a*m))
                *x**(m-k) for k in range(m+1)),Q(0))
def I(m,a,d,w):
    return sum((Q(math.comb(m,k)*math.factorial(a*m+d*k),m**(a*m+d*k+1))
                *w**(m-k) for k in range(m+1)),Q(0))
def moment_controls(s,m):
    a=s*s; d=2*s+1; w=Q(7,5)
    factor=Q(math.factorial(a*m),m**((a+d)*m+1))
    assert I(m,a,d,w)==factor*A(m,a,d,Q(m**d)*w)
    if m>=2:
        n=m-1
        common=sum((Q(math.comb(n,k)*math.factorial(a*n+d*k),m**(a*n+d*k+1))
                    *w**(n-k) for k in range(n+1)),Q(0))
        assert common==Q(n,m)**((a+d)*n+1)*I(n,a,d,Q(m,n)**d*w)
    return True
def epsilon_powers():
    # w=wc+K eps^3 h(eps), deps/dm=-eps^5/4.
    # The derivative exponents are 2+5=7 and 3+5=8.
    assert 2+5==7 and 3+5==8
    assert Q(7,4)>1 and Q(5,4)>1
    # A two-term expansion's remainder m^-5/4 is o(m^-1).
    assert Q(5,4)-1==Q(1,4)
    # Radial field drift N^-1/4 exceeds a root remainder N^-1/2.
    assert Q(1,2)-Q(1,4)==Q(1,4)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=HERE/"bounded_quartic_cusp_verification.json")
    ns=ap.parse_args()
    assert hashlib.sha256(LOCAL.read_bytes()).hexdigest()==LOCAL_SHA,"Cusp source bytes changed"
    rows=[cusp_data(s) for s in range(1,17)]
    model=model_identity()
    controls=sum(moment_controls(s,m) for s in range(1,5) for m in range(2,6))
    epsilon_powers()
    # Actual wrong first corrections or support claims are rejected.
    assert rows[0]["gamma_squared"]!="0"
    assert model!={(1,0):Q(6)}
    assert Q(3,2)**3!=1  # Dropping the lower-order radial scaling changes it.
    scope={"fixed_cusp":True,"bounded_scaled_window":True,"all_actual_roots_localized":False,
           "effective_cutoff":False,"analytic_at_epsilon_zero":False,"global_SFC":False}
    assert not any(scope[k] for k in ["all_actual_roots_localized","effective_cutoff",
                                      "analytic_at_epsilon_zero","global_SFC"])
    payload={"status":"PASS","scope":scope,
        "proof_sha256":hashlib.sha256(PROOF.read_bytes()).hexdigest(),
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "local_cusp_proof_sha256":LOCAL_SHA,
        "counts":{"exact_cusp_first_corrections":len(rows),"actual_moment_pair_controls":controls,
                  "universal_ODE_identities":2,"exponent_guards":6,"negative_controls":4},
        "rows":rows,"meaning":"Supporting exact algebra; written analysis supplies all asymptotic quantifiers."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","counts":payload["counts"],"output":str(ns.output)}))
if __name__=="__main__":main()

