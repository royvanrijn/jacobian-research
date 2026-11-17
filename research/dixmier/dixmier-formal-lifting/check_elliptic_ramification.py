#!/usr/bin/env python3
"""Exact identities for QHR1; the all-degree geometric implications are written."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent

def main():
    r,ss,t,d,lam,X,Y,m=s.symbols("r ss t d lam X Y m")
    Delta=4*ss**2-4*r*ss**3+128*t-144*r*ss*t-432*r*r*t*t
    vx=lambda g:s.expand(s.diff(g,ss)*s.diff(Delta,t)-s.diff(g,t)*s.diff(Delta,ss))
    x=3*r*ss-4
    y=108*r*r*t+18*r*ss-16
    A=(4*r*ss*t-32*t-ss**2)/(20*r*r)
    assert s.expand(y*y+4*x**3+27*r*r*Delta)==0
    assert s.cancel(s.diff(Delta,r)-vx(A)-s.Rational(8,5)*Delta/r)==0
    assert s.diff(Delta,t)==-8*y
    assert s.det(s.Matrix([[s.diff(x,ss),s.diff(x,t)],
                           [s.diff(y,ss),s.diff(y,t)]]))==324*r**3
    # Smooth generic cubic, singular special level lambda=0.
    cubic=-4*X**3-27*r*r*lam
    assert s.discriminant(cubic,X)==-314928*lam**2*r**4
    assert s.factor(cubic.subs(lam,0))==-4*X**3
    # Differential reduction: y*dX = (2/5)d(X*y) -(3/5)c*dX/y.
    c=s.symbols("c")
    curve=-4*X**3-c
    assert s.expand(curve-s.Rational(2,5)*(curve+X*s.diff(curve,X)/2)+s.Rational(3,5)*c)==0

    # The three-dimensional Keller core has determinant one. Its cotangent
    # lift is ordinarily equivalent after stabilization, but need not be canonical.
    bx,by,bb=s.symbols("bx by bb")
    bw=1+bx*by
    core=[2*bx-3*bx*bx*by-bx**3*bb,
          -(bw**3*bb+by*by*bw*(4+3*bx*by))/2,
          by+3*bx*bw*bw*bb+3*bx*by*by*(4+3*bx*by)]
    J=s.Matrix(core).jacobian([bx,by,bb])
    assert s.expand(J.det())==1
    for value in J*J.adjugate()-s.eye(3):
        assert s.expand(value)==0

    # Source ordinary-coordinate inverse, without treating it as canonical.
    xx,q,p,z,beta,d0,yy=s.symbols("xx q p z beta d0 yy")
    aa=1-3*xx*q
    cc=(1+3*xx*q)/2
    M=beta+9*q*q
    pback=3*q*q*M+2*aa*d0
    zback=cc*M-3*xx*xx*d0
    assert s.expand(3*xx*xx*pback+2*aa*zback-9*q*q-beta)==0
    assert s.expand(cc*pback-3*q*q*zback-d0)==0
    assert s.expand(-9*xx*xx*q*q-2*aa*cc)==-1

    v,a=s.symbols("v a")
    newY=a-v
    newBeta=5*v*v-3*v*a-r*v**3
    uv=(1/v)*newY
    w=1+uv
    Rcore=2/v-3*newY/v**2-newBeta/v**3
    Score=newY+3*w*w*newBeta/v+3*newY*newY*(4+3*uv)/v
    Tcore=-(w**3*newBeta+newY*newY*w*(4+3*uv))/2
    S=2*v+4*a-3*r*a*a
    T=(r*a**3-a*a-v*a)/2
    assert s.cancel(Rcore-r)==0
    assert s.cancel(Score-S)==0
    assert s.cancel(Tcore-T)==0
    assert s.expand(Delta.subs({ss:S,t:T})-4*v*v*((3*r*a-2)**2-8*r*v))==0
    assert s.cancel((s.diff(S,a)*s.diff(T,v)-s.diff(S,v)*s.diff(T,a))/v)==1

    # Arbitrary conjugate power, with the d coefficient kept symbolic.
    # On d^m the full bracket is m*Delta_r*d^(m-1).
    assert s.simplify(s.diff(d**m,d)-m*d**(m-1))==0
    f=s.Function("f")(ss,t,d)
    full=s.diff(f,ss)*s.diff(Delta,t)-s.diff(f,t)*s.diff(Delta,ss)+s.diff(f,d)*s.diff(Delta,r)
    assert s.simplify(full.subs(f,d**3).doit()-3*d*d*s.diff(Delta,r))==0
    # Spectators commute with Delta, while r and spectator positions attain the bound.
    u,b=s.symbols("u b")
    def bracket(g,h):
        return s.expand(s.diff(g,d)*s.diff(h,r)-s.diff(g,r)*s.diff(h,d)+
                        s.diff(g,ss)*s.diff(h,t)-s.diff(g,t)*s.diff(h,ss)+
                        s.diff(g,b)*s.diff(h,u)-s.diff(g,u)*s.diff(h,b))
    assert bracket(u,Delta)==bracket(b,Delta)==0
    assert bracket(r,u)==0 and bracket(b,u)==1
    assert bracket(d,r)==1 and bracket(ss,t)==1
    # Sharpness: a Hamiltonian not preserving Delta can be locally nilpotent.
    assert bracket(ss,Delta)!=0
    assert bracket(ss,bracket(ss,bracket(ss,Delta)))==0

    old=json.loads((HERE/"position_avoidance_controls.json").read_text())
    for name,digest in old["hashes"].items():
        path=HERE/name if name!="verify_long_pbw.py" else HERE.parent/"dixmier-audit"/name
        assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name
    paths=[HERE/"ELLIPTIC_RAMIFICATION_AND_STABLE_POLARIZATION.md",Path(__file__),
           HERE/"RAMIFICATION_AND_POSITION_AVOIDANCE.md",HERE/"position_avoidance_controls.json",
           HERE.parent/"dixmier-audit"/"verify_long_pbw.py"]
    result={
        "claim":"QHR1","status":"PASS",
        "Weierstrass_identity":True,
        "elliptic_primitive_identity":"Delta_r=V(A)+8*Delta/(5*r)",
        "primitive":"(4*r*s*t-32*t-s^2)/(20*r^2)",
        "generic_cubic_discriminant":"-314928*lambda^2*r^4",
        "nonzero_holomorphic_class":"8*lambda/(5*r) times omega; written proof",
        "ordinary_source_chart_inverse":True,"Keller_core_cotangent_fiber_inverse":True,"actual_Long_core_chart":True,
        "discriminant_chart_factorization":True,
        "arbitrary_conjugate_power_derivative":True,
        "spectator_Poisson_controls":True,"nonpreserving_LND_control":True,
        "retained_input_hashes":"PASS",
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "coverage_boundary":"Exact algebraic identities only. The nonproperness component, no algebraic primitive, all-degree Poisson centralizer, Hamiltonian LND classification, spectator bound and exclusion of all stable polynomial canonical positional forms are written proofs. No finite computation establishes those universal implications.",
        "independent_replay":False,"formal_verification":False,"external_review":False
    }
    (HERE/"elliptic_ramification_controls.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: Weierstrass and elliptic primitive identities, source chart, discriminant, conjugate derivative, spectator and sharpness controls")
    print("Universal geometric and stable-polarization implications use the written proof.")

if __name__=="__main__":
    main()
