#!/usr/bin/env python3
"""Exact finite controls for QEI1; universal filtered arguments are written."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s

HERE=Path(__file__).resolve().parent

def main():
    spec=importlib.util.spec_from_file_location("qei_ore",HERE/"check_curve_and_module.py")
    o=importlib.util.module_from_spec(spec);spec.loader.exec_module(o)
    r,v,a,h,x=s.symbols("r v a h x")
    S=2*v+4*a-3*r*a*a
    T=(r*a**3-a*a-v*a)/2
    lam=4*v*v*((3*r*a-2)**2-8*r*v)
    b=16*r*r*lam
    xx=8*r*v; yy=xx*(3*r*a-2)
    xp=16-12*r*S;yp=432*r*r*T+72*r*S-64
    assert s.cancel(xp-(xx**3+4*b)/xx**2)==0
    assert s.cancel(yp+yy*(xx**3-8*b)/xx**3)==0
    assert s.expand(yy**2-xx**3-b)==0
    assert s.expand(yp**2-xp**3+27*b)==0
    z,B=s.symbols("z B")
    assert s.diff(z+4*B/z**2,z)==1-8*B/z**3

    X={0:3*r*x-4}
    Y={1:-108*r*r*h,0:18*r*x-16}
    X2=o.mul(X,X,x);X3=o.mul(X2,X,x)
    energy=o.add(o.mul(Y,Y,x),o.scale(X3,4))
    Delta=o.scale(energy,-1/(27*r*r))
    A={1:(-4*r*h*x+32*h)/(20*r*r),0:(-2*r*h-x*x)/(20*r*r)}
    Dr=lambda Q:{i:s.diff(c,r) for i,c in Q.items()}
    covariant=lambda Q:o.add(Dr(Q),o.scale(o.comm(A,Q,x),-1/h))
    assert o.equal(o.comm(X,Y,x),{0:324*h*r**3})
    assert o.equal(covariant(Delta),o.scale(Delta,s.Rational(8,5)/r))
    assert o.equal(covariant(X),o.scale(X,s.Rational(6,5)/r))
    assert o.equal(covariant(Y),o.scale(Y,s.Rational(9,5)/r))
    zeta=s.symbols("zeta")
    X0={i:q.subs(r,zeta**5)*zeta**(-6) for i,q in X.items()}
    Y0={i:q.subs(r,zeta**5)*zeta**(-9) for i,q in Y.items()}
    A0={i:q.subs(r,zeta**5) for i,q in A.items()}
    cover=lambda Q:o.add({i:s.diff(q,zeta) for i,q in Q.items()},
                         o.scale(o.comm(A0,Q,x),-5*zeta**4/h))
    assert not cover(X0) and not cover(Y0)
    assert o.equal(o.comm(X0,Y0,x),{0:324*h})
    energy0=o.add(o.mul(Y0,Y0,x),o.scale(o.mul(o.mul(X0,X0,x),X0,x),4))
    assert o.equal({i:q.subs(r,zeta**5) for i,q in Delta.items()},
                   o.scale(energy0,-zeta**8/27))

    # Universal oscillator identities in an independent scalar presentation.
    c=s.symbols("c",nonzero=True)
    U={0:x};V={1:-c}
    L=o.add(o.mul(V,V,x),{0:4*x**3})
    C1=o.comm(L,U,x);C2=o.comm(L,V,x)
    assert o.equal(C1,o.scale(V,-2*c))
    assert o.equal(C2,{0:12*c*x*x})
    assert o.equal(o.comm(L,C1,x),{0:-24*c*c*x*x})

    aa,aa0,bb,bb1,bb0=s.symbols("aa aa0 bb bb1 bb0")
    Us={0:aa*x+aa0}
    Vs={1:-bb*c,0:bb1*x+bb0}
    candidate=o.add(o.mul(Vs,Vs,x),o.scale(o.mul(o.mul(Us,Us,x),Us,x),4))
    assert s.expand(o.comm(Us,Vs,x)[0]-aa*bb*c)==0
    # After the leading constraints force aa=bb=1, the cross terms force zero shifts.
    reduced=o.add({i:q.subs({aa:1,bb:1}) for i,q in candidate.items()},o.scale(L,-1))
    assert s.expand(reduced[1]+2*c*(bb1*x+bb0))==0
    reduced0=s.expand(reduced[0].subs({bb1:0,bb0:0}))
    assert s.Poly(reduced0,x).coeff_monomial(x*x)==12*aa0

    iterates=[]
    Q=U
    for n in range(7):
        # c is a scalar; the leading weighted degree grows at each commutator.
        wt=max(2*s.degree(q,x)+3*j for j,q in Q.items())
        assert wt==2+n
        iterates.append({"n":n,"weighted_degree":int(wt)})
        Q=o.scale(o.comm(L,Q,x),1/c)

    for receipt in ["position_avoidance_controls.json","elliptic_ramification_controls.json"]:
        old=json.loads((HERE/receipt).read_text())
        for name,digest in old["hashes"].items():
            p=HERE/name if name!="verify_long_pbw.py" else HERE.parent/"dixmier-audit"/name
            assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,(receipt,name)
    paths=[HERE/"ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md",Path(__file__),
           HERE/"RAMIFICATION_AND_POSITION_AVOIDANCE.md",HERE/"position_avoidance_controls.json",
           HERE/"ELLIPTIC_RAMIFICATION_AND_STABLE_POLARIZATION.md",HERE/"elliptic_ramification_controls.json",
           HERE/"check_curve_and_module.py"]
    result={
        "claim":"QEI1","status":"PASS",
        "Long_generic_fiber_is_negative_Velu_three_isogeny":True,
        "elliptic_curve_parameters":{"input":"b=16*r^2*lambda","output":"-27*b"},
        "exact_quantum_normalizer_eigenweights":["6/5","9/5","8/5"],
        "exact_quantum_elliptic_identity":True,"ramified_constant_commutator":True,
        "universal_oscillator_commutator_controls":True,"finite_stabilizer_coefficient_controls":True,
        "formal_fixed_energy_iterate_samples":iterates,
        "retained_QPA1_QHR1_hash_chain":"PASS",
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        "coverage_boundary":"Finite exact identities only. The isogeny interpretation, polynomial centralizer and commutator quotient, full source normalizer, all-degree finite bounds, stabilizer rigidity, forced Long energy-order jump and formal unboundedness use written proofs. No bounded full quantum lift or all-order nonlifting theorem is certified.",
        "independent_replay":False,"formal_verification":False,"external_review":False
    }
    (HERE/"quantum_elliptic_energy_controls.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: explicit three-isogeny, exact quantum normalizer, oscillator bounds identities, stabilizer controls and seven formal-flow coefficients")
    print("Universal finite-degree and energy-order statements use the written proof.")

if __name__=="__main__":
    main()
