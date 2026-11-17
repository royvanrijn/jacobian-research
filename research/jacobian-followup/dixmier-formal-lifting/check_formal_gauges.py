#!/usr/bin/env python3
"""Exact controls for QFG1; the group classification and all-order escape are written proofs."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import hashlib
import importlib.util
import itertools
import json
import sympy as s

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zero(expr,label):
    if s.cancel(s.expand(expr))!=0:
        raise AssertionError(label+": "+str(expr))


def main():
    x,y,h=s.symbols("x y h")
    f=s.Function("f")(x,y)
    theta=lambda g:x**2*s.diff(g,x)-(1+2*x*y)*s.diff(g,y)
    a=1+x*f
    XX=x*a
    YY=(y-f)/a**2
    J=(1+2*x*f+theta(f))/a**2
    zero(XX*(1+XX*YY)-x*(1+x*y),"arbitrary formal base preserves r")
    zero(s.det(s.Matrix([[s.diff(XX,x),s.diff(XX,y)],
                         [s.diff(YY,x),s.diff(YY,y)]]))-J,"arbitrary base Jacobian")
    zero(theta(XX)-J*XX**2,"arbitrary base theta-X")
    zero(theta(YY)+J*(1+2*XX*YY),"arbitrary base theta-Y")

    # Exact cotangent+gradient example, with rational formal units.
    XX=x*(1+h*x)
    YY=(y-h)/(1+h*x)**2
    jac=s.Matrix([[s.diff(XX,x),s.diff(XX,y)],
                  [s.diff(YY,x),s.diff(YY,y)]])
    inverse=jac.inv()
    G=h*x**2*y
    pp=[inverse[0,0],inverse[1,0],
        inverse[0,0]*s.diff(G,x)+inverse[1,0]*s.diff(G,y)]
    zz=[inverse[0,1],inverse[1,1],
        inverse[0,1]*s.diff(G,x)+inverse[1,1]*s.diff(G,y)]
    gens=[[0,0,XX],[0,0,YY],pp,zz]
    def comm_first(A,B):
        ax,ay,ac=A
        bx,by,bc=B
        return [h*(ax*s.diff(bx,x)+ay*s.diff(bx,y)-bx*s.diff(ax,x)-by*s.diff(ax,y)),
                h*(ax*s.diff(by,x)+ay*s.diff(by,y)-bx*s.diff(ay,x)-by*s.diff(ay,y)),
                h*(ax*s.diff(bc,x)+ay*s.diff(bc,y)-bx*s.diff(ac,x)-by*s.diff(ac,y))]
    for i,j in itertools.combinations(range(4),2):
        val=comm_first(gens[i],gens[j])
        expected=-h if (i,j) in [(0,2),(1,3)] else 0
        zero(val[0],"full first-order relation momentum x")
        zero(val[1],"full first-order relation momentum y")
        zero(val[2]-expected,"full first-order scalar relation")
    lam=1/s.det(jac)
    b_image=[s.cancel(XX**2*pp[k]-(1+2*XX*YY)*zz[k]) for k in range(3)]
    for actual,expected in zip(b_image,[lam*x**2,-lam*(1+2*x*y),lam*theta(G)]):
        zero(actual-expected,"full centralizer normal form example")
    zero(s.det(jac)*(h**3/s.det(jac))-h**3,"normalized nonadmission translation")

    # Reuse the retained exact PBW calculus, after binding its entire control receipt.
    retained=json.loads((HERE/"painleve_long_controls.json").read_text())
    if retained["status"]!="PASS":
        raise AssertionError("QPL1 receipt status")
    for name,expected in retained["hashes"].items():
        if sha(HERE/name)!=expected:
            raise AssertionError("QPL1 input drift: "+name)
    spec=importlib.util.spec_from_file_location("pi_controls",HERE/"check_painleve_long.py")
    pc=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pc)
    mon,add,mul,power,scale=pc.mon,pc.add,pc.mul,pc.power,pc.scale
    BB=mon(c=1)
    R=mon(a=1)
    X=mon(b=1)
    hh=mon(d=1)
    a0=mon(a=1,b=-1)
    w=add(power(a0,3),scale(power(a0,2),-5),scale(a0,9),mon(k=-2),
          mon(a=1,b=1,c=1,k=-1),mon(a=1,b=2,d=1,k=-F(3,2)))
    I=scale(add(mul(power(X,3),w),mul(w,power(X,3))),F(1,2))
    H=power(I,2)
    global_H=pc.globalize(H)
    if max(e[2] for e in H)!=2:
        raise AssertionError("Hamiltonian momentum order two")
    CK=add(mon(b=-2,c=1,k=-F(1,2)),mon(a=1,b=-3,c=1,k=F(1,6)),
        mon(a=-1,b=-1,c=1,k=-F(1,15)),mon(b=-4,k=-F(3,4)),
        mon(a=1,b=-5,k=F(1,5)),mon(a=-1,b=-3,k=F(3,5)),
        mon(a=-2,b=-2,k=F(1,5)),mon(a=-1,d=1,k=F(1,10)))
    dB=add(mon(b=-1,c=1,k=3),mon(a=1,b=-2,c=1,k=-F(3,2)))
    def delta(p):
        out={}
        for (a,b,c,d),coef in p.items():
            posdelta=add(mon(a=a-1,b=b,d=d,k=a),
                         mon(a=a,b=b-1,d=d,k=F(3*b,2)),
                         mon(a=a+1,b=b-2,d=d,k=-F(b,2)))
            delta_power=add(*[mul(mul(power(BB,j),dB),power(BB,c-1-j))
                              for j in range(c)])
            term=add(mul(posdelta,power(BB,c)),mul(mon(a=a,b=b,d=d),delta_power))
            out=add(out,scale(term,coef))
        return out
    pc.equal(add(mul(hh,delta(I)),pc.comm(CK,I)),{},"exact K-I commutation")
    pc.equal(add(mul(hh,delta(H)),pc.comm(CK,H)),{},"exact K-H commutation")
    pc.equal(pc.comm(R,H),{},"exact R-H commutation")
    pc.equal(pc.comm(I,X),mon(a=1,b=6,d=1,k=-3),"exact invariant position field")
    def deriv(p):
        out=pc.comm(H,p)
        if not all(e[3]>=1 for e in out):
            raise AssertionError("normalized commutator divisibility")
        return {(a,b,c,d-1):coef for (a,b,c,d),coef in out.items()}
    pos_iter=X
    pos_samples=[]
    for n in range(5):
        leading={e:c for e,c in pos_iter.items() if e[2]==n}
        expected=mon(a=2*n,b=9*n+1,c=n,k=6**n*prod(5*j+1 for j in range(n)))
        pc.equal(leading,expected,"all-coefficient leading position iterate "+str(n))
        if max(e[2] for e in pos_iter)!=n:
            raise AssertionError("position iterate order")
        pos_samples.append({"n":n,"B_order":n,"leading_coefficient":str(next(iter(expected.values())))})
        if n<4:
            pos_iter=deriv(pos_iter)

    seed=json.loads((HERE/"painleve_long_seed.json").read_text())
    S0={tuple(e):F(c) for e,c in seed["local"]["S"] if e[3]==0}
    op_iter=S0
    full_samples=[]
    for n in range(4):
        leading={e:c for e,c in op_iter.items() if e[2]==n+2}
        expected=mon(a=2*n+1,b=9*n+2,c=n+2,
                     k=-F(6**n,3)*prod(5*j-6 for j in range(n)))
        pc.equal(leading,expected,"leading full-frame iterate "+str(n))
        full_samples.append({"n":n,"S_B_order":n+2,"coefficient":str(next(iter(expected.values())))})
        if n<3:
            op_iter=deriv(op_iter)

    locked=["FORMAL_GAUGES_AND_MOMENTUM_ESCAPE.md","check_formal_gauges.py",
            "PAINLEVE_GLOBAL_LONG_LIFT.md","painleve_long_controls.json",
            "BOUNDED_GAUGE_RIGIDITY.md","RELATIVE_PERIOD_TORSION.md",
            "QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md"]
    out={"claim":"QFG1","status":"PASS",
        "arbitrary_formal_base_factorization_and_Jacobian":"PASS",
        "full_cotangent_gradient_example_relations":6,
        "centralizer_and_admission_normalization_example":"PASS",
        "global_invariant_H_polynomial_terms":len(global_H),
        "exact_R_and_K_preservation":"PASS",
        "position_flow_leading_samples":pos_samples,
        "transported_full_frame_leading_samples":full_samples,
        "contact_scope":"every N>=1, by the written h-valuation and leading-order argument",
        "hashes":{name:sha(HERE/name) for name in locked},
        "coverage_boundary":"Exact generic Jacobian identities, one full formal-unit example and selected invariant-flow commutators only. The full classification, admission criterion, unique-comparison obstruction and all-order unbounded momentum use written proofs. Exact R and K do not themselves force finite-momentum gauges. No bounded Long candidate or all-representative nonlifting is certified.",
        "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"formal_gauge_controls.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"claim":"QFG1","status":"PASS",
        "full_example_relations":6,"exact_R_K_invariant_flow":"PASS",
        "position_iterate_orders":[x["B_order"] for x in pos_samples],
        "transported_S_orders":[x["S_B_order"] for x in full_samples]},indent=2))


if __name__=="__main__":
    main()
