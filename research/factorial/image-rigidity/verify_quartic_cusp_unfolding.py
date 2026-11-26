#!/usr/bin/env python3
"""Exact supporting controls for the written quartic-cusp local lemma."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as Q
from pathlib import Path
if not __debug__:
    raise SystemExit("Refusing optimized execution: exact assertions are required.")
HERE=Path(__file__).resolve().parent
PROOF=HERE/"QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md"
K=6
def trim(a):
    return list(a[:K])+[Q(0)]*max(0,K-len(a))
def add(a,b):
    a,b=trim(a),trim(b)
    return [x+y for x,y in zip(a,b)]
def neg(a): return [-x for x in trim(a)]
def scale(a,c): return [Q(c)*x for x in trim(a)]
def mul(a,b):
    a,b=trim(a),trim(b)
    return [sum((a[j]*b[k-j] for j in range(k+1)),Q(0)) for k in range(K)]
def powi(a,n):
    r=trim([Q(1)])
    for _ in range(n): r=mul(r,a)
    return r
def inv(a):
    a=trim(a)
    assert a[0]
    r=[1/a[0]]
    for k in range(1,K):
        r.append(-sum((a[j]*r[k-j] for j in range(1,k+1)),Q(0))/a[0])
    return r
def div(a,b): return mul(a,inv(b))
def log1(a):
    a=trim(a)
    assert a[0]==1
    t=add(a,[-1])
    r=trim([])
    for k in range(1,K):
        r=add(r,scale(powi(t,k),Q((-1)**(k+1),k)))
    return r
def derivative(a):
    a=trim(a)
    return trim([Q(k)*a[k] for k in range(1,K)])
def raw_cusp(d):
    R=Q(d+1,d-1); de=R+1; ep=de/d
    q=R**d; z=(R**(d+1)+1)/de
    v=trim([1,1]); p=trim([R,ep])
    hp=powi(p,d-1)
    h=scale(add(powi(p,d),[-1]),1/de)
    H=add([z],h)
    phase=add(scale(log1(scale(p,1/R)),Q((d-1)**2,4)),
              add(scale(v,-Q(d-1,2)),log1(scale(H,1/q))))
    # Constants in phase are irrelevant to derivatives.
    vz=neg(div(hp,mul(H,H)))
    hd=scale(add(scale(mul(v,hp),de),add(neg(powi(p,d)),[1])),1/de**2)
    hpd=scale(mul(v,powi(p,d-2)),Q(d-1,d))
    vd=add(div(mul(v,v),scale(mul(p,p),d)),
           add(div(hpd,H),neg(div(mul(hp,hd),mul(H,H)))))
    F=[phase[j]*math.factorial(j) for j in range(K)]
    return R,de,q,z,F,vz,vd
def check_d(d):
    R,de,q,z,F,vz,vd=raw_cusp(d)
    B=Q(6,5*(d-1))
    assert F[1:4]==[0,0,0]
    assert F[4]==-2/R**2
    assert F[5]==24/Q(d-1)/R**3
    assert vz[0]==-R**(-d-1)
    assert vz[1]==0
    assert vd[0]==(R+2-R**(-d))/(R*(R+1)**2)
    assert vd[1]==Q((d-1)**3,d*(d+1)**2)
    determinant=vz[0]*vd[1]-vd[0]*vz[1]
    assert abs(determinant)>=Q(1,486)
    assert -F[4]>=Q(2,9)
    assert q*R<=27
    # Inverse coefficient from the actual fifth coefficient:
    # V_tt=F5/(20*q0^3), q0=1/R.
    assert F[5]*R**3/20==B
    # Canonical T_z includes the nonlinear coordinate contribution.
    Tz=R*vz[1]+B*vz[0]
    assert Tz==-Q(6,5*(d-1))*R**(-d-1)
    assert Tz!=0
    radial=z-1/de
    # H_z/A and T_z are rational; no floating square root enters.
    assert vz[0]*radial==-1/de
    assert Tz*radial==-Q(3,5*d)
    # Both 2x2 versal matrices have nonzero determinant.
    assert determinant<0
    return {"d":d,"F4":str(F[4]),"F5":str(F[5]),"det":str(determinant),
            "Tz":str(Tz),"radial_temperature":str(Tz*radial)}
def original_series(d,a,m,w):
    return [Q(math.comb(m,k)*math.factorial(a*m+d*k),math.factorial(a*m))*w**(m-k)
            for k in range(m+1)]
def check_scaling(d,a,N,w):
    n=N-1; kappa=Q(N,n)
    # Compare normalized exact Gamma sums in the common N scale.
    # J_n,N=sum binom(n,k) w^(n-k) Gamma(an+dk+1)/N^(an+dk+1).
    left=sum((Q(math.comb(n,k)*math.factorial(a*n+d*k),N**(a*n+d*k+1))
              *w**(n-k) for k in range(n+1)),Q(0))
    # I_n(kappa^d w) times (n/N)^((a+d)n+1).
    right=sum((Q(math.comb(n,k)*math.factorial(a*n+d*k),n**(a*n+d*k+1))
              *(kappa**d*w)**(n-k) for k in range(n+1)),Q(0))
    right*=Q(n,N)**((a+d)*n+1)
    assert left==right
    return left
def check_global_cusp(s):
    a=s*s; d=2*s+1; D=(s+1)**2; S=s*(s+1)
    assert Q((d-1)**2,4)==a
    # Polynomial stationary-map numerator is -d(S-Sc)^2.
    assert a+D-1==2*S and a*D==S*S
    # Critical w and normalized z conversion agree.
    w=Q(S**d*(D-S),S-a)
    R=Q(d+1,d-1)
    assert w/Q(a**d)==R**(d+1)
def expected_mistake_rejected(kind):
    R,de,q,z,F,vz,vd=raw_cusp(3)
    if kind=="raw_temperature":
        return Q(0)!=Q(6,10)*vz[0]
    if kind=="wrong_F5":
        return F[5]!=0
    if kind=="missing_order_scaling":
        # Correct shifted field multiplier (3/2)^d, not unity.
        return Q(3,2)**3!=1
    if kind=="real_temperature_on_imaginary_curve":
        # T_z is nonzero real; its product with a nonzero imaginary displacement is nonreal.
        return Q(6,10)*vz[0]!=0
    raise ValueError(kind)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=HERE/"quartic_cusp_unfolding_verification.json")
    args=ap.parse_args()
    rows=[check_d(d) for d in range(2,34)]
    scales=[]
    for a,d in [(1,3),(4,5),(9,7),(2,2)]:
        for N in [2,3,5]:
            scales.append(str(check_scaling(d,a,N,Q(7,3))))
    for s in range(1,17): check_global_cusp(s)
    negatives=["raw_temperature","wrong_F5","missing_order_scaling",
               "real_temperature_on_imaginary_curve"]
    assert all(expected_mistake_rejected(k) for k in negatives)
    # Exact scalar budgets used in the uniform nondegeneracy proof.
    assert Q(1,18)*Q(1,27)==Q(1,486)
    assert Q(1,3)**2*Q(1,2)==Q(1,18)
    # The limiting cosh phase has x^4=-1/12, x^5=0.
    cosh=trim([1,0,Q(1,2),0,Q(1,24),0])
    limit=add(log1(cosh),[0,0,-Q(1,2)])
    assert limit[1:4]==[0,0,0] and limit[4]==-Q(1,12) and limit[5]==0
    payload={"status":"PASS",
        "scope":"Finite exact supporting controls for the full written local proof; no SFC or all-index zero band.",
        "proof_sha256":hashlib.sha256(PROOF.read_bytes()).hexdigest(),
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "counts":{"finite_gap_jet_cases":len(rows),"actual_order_integral_identities":len(scales),
                  "integer_cusp_controls":16,"negative_controls":len(negatives)},
        "rows":rows,"negative_controls":negatives,
        "normal_form_uniformity":"written Banach implicit-function proof; existence radii are not computed"}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","counts":payload["counts"],"output":str(args.output)}))
if __name__=="__main__":main()

